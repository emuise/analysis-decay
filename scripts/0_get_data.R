library(tidyverse)
library(wdpar)
library(terra)
library(tidyterra)
library(sf)
# remotes::install_github("https://github.com/emuise/budR")
library(budR) # has my keys in it
source(here::here("scripts", "0a_make_folders.R"))

terraOptions(memfrac = 0.90,
             tempdir = "F:\\scratch2")

# bc boundary
bcb_loc <- here::here(shapefile_loc, "bcb.shp")

if (!file.exists(bcb_loc)) {
  # clean and validate bc bounds
  bcb <- bcmaps::bc_bound_hres() %>%
    st_make_valid() %>%
    st_combine() %>%
    st_union() %>%
    st_make_valid()
  
  write_sf(bcb, bcb_loc)
}

bcb <- vect(bcb_loc)

# forests
forests_loc <- here::here(raster_loc, "forests.dat")

if (!dir.exists(dirname(forests_loc))) {
  dir.create(dirname(forests_loc), recursive = T)
}

if (!file.exists(forests_loc)) {
  vlce_template <- rast("F:/mosaiced/VLCE2.0/LC_Class_terr_bc.dat")
  
  forest_rcl <- keys$vlce %>%
    mutate(class_name = ifelse(forest == "Forest", class_val, NA)) %>%
    select(class_val, class_name) %>%
    as.matrix(ncol = 2)
  
  forests <- classify(vlce_template, forest_rcl) %>%
    trim() %>%
    extend(bcb)
  
  writeRaster(
    forests,
    forests_loc,
    filetype = "envi",
    datatype = "INT1U",
    overwrite = T
  )
}

forests <- rast(forests_loc)



bcb_wgs84 <- bcb %>%
  project("epsg:4326")

# some are in wgs84 and thus need to have a different cropping method
# trying cropping to bcb_wgs84 ext plus a boundary, also used to determine
# which provinces/states have large cities for the distance raster

bcb_wgs84_buff <- buffer(bcb_wgs84, 25000)




# bc protected areas
pa_loc <- here::here(shapefile_loc, "bc_pa_filt.shp")

if (!file.exists(pa_loc)) {
  # get canadian PAs from wdpar
  cad_pa <-
    wdpa_fetch("CAN",
               wait = T,
               download_dir = here::here("data"))
  
  # get bc terrestrial PA
  bc_pa <- cad_pa %>%
    filter(SUB_LOC == "CA-BC") %>%
    st_transform(3005)
  
  bc_terr_pa <- bc_pa %>%
    st_make_valid() %>%
    st_intersection(bcb %>%
                      st_as_sf())
  
  # filter pa based off bolton et al. (2018)
  # less than 100 ha; in IUCN class Ia Ib II and IV
  bc_pa_filt <- bc_terr_pa %>%
    mutate(area = st_area(.) %>%
             as.numeric() %>%
             round(digits = 2)) %>%
    filter(IUCN_CAT %in% c("Ia", "Ib", "II", "IV")) %>%
    group_by(WDPAID, IUCN_CAT) %>%
    summarize(
      area = sum(area),
      minyear = min(STATUS_YR),
      maxyear = max(STATUS_YR)
    ) %>%
    filter(area > 1e6)
  
  write_sf(bc_pa_filt, pa_loc)
  rm(cad_pa, bc_pa, bc_terr_pa, bc_pa_filt)
}

bc_pa_filt <- vect(pa_loc)



# bc bec zones dissolved
bec_loc <- here::here(shapefile_loc, "bec_terr_agg.shp")


if (!file.exists(bec_loc)) {
  bec <- bcmaps::bec() %>%
    vect()
  
  bec_agg <- aggregate(bec,
                       by = "ZONE",
                       dissolve = T)
  
  bec_terr <- intersect(bec_agg, bcb)
  writeVector(bec_terr, bec_loc)
  rm(bec, bec_agg)
}

bec_terr <- vect(bec_loc)

# iter expand function
iter_expand <- function(raster) {
  # project and crop to bcb twice
  # crops to bcb twice because of projection issues
  if (crs(raster) == crs("epsg:4326")) {
    raster <- raster %>%
      crop(bcb_wgs84_buff) %>%
      project(
        "epsg:3005",
        threads = T,
        gdal = T,
        by_util = T
      )
  }
  
  bcb_rast <- bcb %>%
    rasterize(raster, touches = T)
  
  projed <- raster %>%
    crop(bcb_rast, mask = T)
  
  # count how many pixels there should be
  old <- bcb_rast %>%
    freq() %>%
    pull(count)
  
  print("value to match")
  print(old)
  print("----")
  
  # set a base to expand from as the OG raster
  expanded <- projed
  
  
  # loop the focal analysis until the number of pixels in bcb is the same
  # focal, then mask to bcb, then classify and count number of cells
  # if they are the same, break the loop and return the masked output
  for (i in 1:50) {
    expanded <-
      focal(
        expanded,
        w = 5,
        fun = "mean",
        na.rm = F,
        na.policy = "only"
      )
    
    # mask to bcb, classify, then count number of pixels
    masked <- mask(expanded, bcb_rast)
    
    classed <-
      classify(masked, cbind(minmax(projed)[1], minmax(projed)[2], 1), include.lowest = T)
    
    new <- freq(classed) %>%
      pull(count)
    
    print(new)
    
    # if num pixels same, break the loop, return the final masked raster
    if (unique(old == new)) {
      print(paste0("done after ", i, " loops"))
      break
    }
    
    
  }
  masked
}



# travel time to city; Nelson et al., 2019; data is for 2015 (perfect!)
travel_loc <-
  here::here(raster_loc, "travel_time_to_cities_bc.dat")

if (!file.exists(travel_loc)) {
  travel_city <-
    rast("https://figshare.com/ndownloader/files/14189825")
  
  travel_rcl <- subst(travel_city, 65535, NA)
  
  masked <- iter_expand(travel_rcl)
  
  travel_bc <-
    resample(masked, forests, method = "cubicspline", threads = T)
  
  names(travel_bc) <- "travel_time_to_cities"
  
  writeRaster(travel_bc,
              travel_loc,
              filetype = "envi",
              overwrite = T)
}

travel_bc <- rast(travel_loc)


# population density and count!
# density (2015) from https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11
# count (2015) from https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-rev11
# need usgs account to download or i would've done programmatically

# i do not understand why we need to match based on both, but this is following
# Duncanson et al., (2023)

pops_locs <-
  here::here(raster_loc, c("pop_count.dat", "pop_density.dat"))


if (!all(file.exists(pops_locs))) {
  pops <-
    list.files("F:/global", full.names = T, pattern = ".tif$") %>%
    map(rast)
  
  pops_bc <- pops %>%
    map(iter_expand, .progress = T) %>%
    rast() %>%
    resample(forests, method = "cubicspline", threads = T)
  
  names(pops_bc) <- c("pop_count", "pop_density")
  
  writeRaster(pops_bc[[1]],
              pops_locs[1],
              overwrite = T,
              filetype = "envi")
  
  writeRaster(pops_bc[[2]],
              pops_locs[2],
              overwrite = T,
              filetype = "envi")
  
}


bc_pops <- rast(pops_locs)

distance_loc <- here::here("data", "rasters", "city_distance.dat")

if (!file.exists(distance_loc)) {
  download_loc <- here::here(scratch, "distance_city.zip")
  
  if (!file.exists(download_loc)) {
    distance_url <-
      "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_SMOD_GLOBE_R2023A/GHS_SMOD_E2015_GLOBE_R2023A_54009_1000/V1-0/GHS_SMOD_E2015_GLOBE_R2023A_54009_1000_V1_0.zip"
    download.file(distance_url, download_loc)
  }
  
  unzipped <- unzip(download_loc, exdir = scratch)
  
  ghs <- unzipped %>% str_subset(".tif$") %>%
    rast()
  
  rcl <- rbind(cbind(0, 23, NA), cbind(24, 30, 1))
  # keep only urban centre grid cells
  
  ghs_rcl <- classify(ghs, rcl)
  
  ghs_84 <-
    ghs_rcl %>% project(
      "epsg:4326",
      method = "near",
      threads = T,
      gdal = T,
      by_util = T
    )
  
  canus <-
    geodata::gadm(c("Canada", "US"), path = scratch, level = 1)
  
  isos <- bcb_wgs84_buff %>% intersect(canus) %>% pull(ISO_1)
  
  canus_filt <- canus %>%
    filter(ISO_1 %in% isos)
  
  bc_ghs <- ghs_84 %>%
    crop(canus_filt, mask = T) %>%
    trim() %>%
    project(
      "epsg:3005",
      threads = T,
      method = "near",
      gdal = T,
      by_util = T
    )
  
  ghs_30 <- bc_ghs %>%
    resample(forests, method = "cubicspline", threads = T)
  
  distance_city <- distance(bc_ghs) %>%
    resample(forests, method = "cubicspline", threads = T) %>%
    crop(bcb, mask = T)
  
  names(distance_city) = "distance_city"
  
  writeRaster(distance_city,
              distance_loc,
              filetype = "envi",
              overwrite = T)
}

distance_city <- rast(distance_loc)


# climate data; from chris mulverhill derived from climate NA
# mean annual precipitation, mean annual temperature
# mean coldest month temperature
# mean warmest month temperature
# all climate normals (1990-2020)

clim_locs <-
  list.files(path = raster_loc,
             pattern = "^clim.*\\.dat$",
             full.names = T)

if (length(clim_locs) != 4) {
  clim_rasts <- list.files(path = "F:/mosaiced/climatena_1k",
                           pattern = ".dat$",
                           full.names = T) %>%
    map(rast)
  
  names <- map(clim_rasts, sources) %>%
    map(basename) %>%
    map(str_replace, pattern = "BC", replacement = "clim") %>%
    unlist() %>%
    here::here(raster_loc, .)
  
  clim_rasts <- clim_rasts %>%
    map(iter_expand)
  
  clim_rasts <- rast(clim_rasts)
  
  names(clim_rasts) <- names %>%
    basename() %>%
    tools::file_path_sans_ext()
  
  clim_res <- clim_rasts %>%
    resample(forests, method = "cubicspline", threads = T)
  
  map2(
    writeRaster,
    .x = clim_res %>% as.list(),
    .y = names,
    filetype = "envi",
    overwrite = T,
    .progress = "Save"
  )
}

clim_locs <-
  list.files(path = raster_loc,
             pattern = "^clim.*\\.dat$",
             full.names = T)

clim_rasts <- rast(clim_locs)
names(clim_rasts) <- sources(clim_rasts) %>%
  basename() %>%
  tools::file_path_sans_ext()


check_extreme_clim = F
if (check_extreme_clim) {
  # tile the climate things so i can look at the data
  grid_ncol <- 5
  grid_nrow <- 5
  
  # make tiles to grab samples from
  grid <-
    rast(ncol = grid_ncol,
         nrow = grid_ncol,
         ext = ext(clim_rasts))
  
  dir.create(scratch, showWarnings = F)
  
  tiles <-
    makeTiles(
      clim_rasts,
      grid,
      filename = here::here(scratch, "til_.tif"),
      overwrite = T
    )
  
  rename_tiles <- function(tile) {
    names(tile) <- names(clim_rasts)
    tile
  }
  
  til_dfs <- tiles %>%
    map(rast, .progress = "rasts") %>%
    map(rename_tiles) %>%
    map(as.data.frame, .progress = "dfs") %>%
    map(tibble)
  
  outlier_MAP <- til_dfs %>%
    map_dfr(filter, clim_MAP < 0, .progress = "filter")
  # 6 values; safe to just clear them by reclassifying them to NA or 0 since they're
  # impossible anyway
  
  outlier_MAT <- til_dfs %>%
    map_dfr(filter, clim_MAT < -10, .progress = "filter")
  # -12.5 seems to be the extreme end of MAT
  
  outlier_MCMT <- til_dfs %>%
    map_dfr(filter, clim_MCMT < -20, .progress = "filter")
  # -23 seems to be the exteme end of MCMT
  
  outlier_MWMT <- til_dfs %>%
    map_dfr(filter, clim_MWMT < 10, .progress = "filter")
  
  keep_na <- function(tibble) {
    tibble %>%
      filter_all(any_vars(is.na(.)))
  }
  
  outlier_MAT %>%
    filter(clim_MAT > -20) %>%
    ggplot(aes(x = clim_MAT)) +
    geom_density()
  # -12.5 is the barrier
  
  outlier_MCMT %>%
    filter(clim_MCMT > -25) %>%
    ggplot(aes(x = clim_MCMT)) +
    geom_density()
  # -23 is the barrier
  
  outlier_MWMT %>%
    filter(clim_MWMT > -25) %>%
    ggplot(aes(x = clim_MWMT)) +
    geom_density()
  # -5 is the barrier; checked by mapping it in arcgis.
  # these ~-5 values are in mt fairweather
  # which does not have any nearby british columbian area.
  # these values seem reasonable, having checked the average climate for the
  # base of mt fairweather
  
  # checking for extreme high precip
  outlier_MAP <- til_dfs %>%
    map_dfr(filter, clim_MAP > 5000, .progress = "filter")
  
  outlier_MAP %>%
    filter() %>%
    ggplot(aes(x = clim_MAP)) +
    geom_density()
  
  
  outlier_NA <- til_dfs %>%
    map_dfr(keep_na, .progress = T)
  
  # NA values in temps are always consistent, and have real values in the
  # MAP
}

minmax <- minmax(clim_rasts)
# passes the sniff test after filtering for only forested pixels

mins <- minmax %>%
  as.data.frame() %>%
  tibble() %>%
  mutate(minmax = c("min", "max")) %>%
  filter(minmax == "min")

if (mins$clim_MAP < 0) {
  clim_rasts[[1]] <-
    classify(clim_rasts[[1]], cbind(-Inf, 0, NA)) # can't have minimum precipitation below zero
}
if (mins$clim_MAT < -12.5) {
  clim_rasts[[2]] <-
    classify(clim_rasts[[2]], cbind(-Inf,-12.5, NA)) # manual check of MAT
}
if (mins$clim_MCMT < -23) {
  clim_rasts[[3]] <-
    classify(clim_rasts[[3]], cbind(-Inf,-23, NA)) # manual check of MCMT
}
if (mins$clim_MWMT < -5) {
  clim_rasts[[4]] <-
    classify(clim_rasts[[4]], cbind(-Inf,-5, NA)) # manual check of MWMT
}
# end climate data sanity checks

footprint_loc <-
  here::here(raster_loc, "bc_albers_footprint.dat")

# human footprint cropping projecting etc
if (!file.exists(footprint_loc)) {
  footprint <-
    rast("Z:/_CanadaLayers/Rasters/canada_human_footprint/cum_threat2020.02.18.tif")
  
  bcb_lamb <- bcb %>%
    project(footprint)
  
  crop_footprint <- footprint %>%
    crop(bcb_lamb) %>%
    project("epsg:3005",
            threads = T,
            gdal = T,
            by_util = T)
  
  footprint_albers <- iter_expand(crop_footprint) %>%
    resample(forests, method = "cubicspline", threads = T)
  
  names(footprint_albers) = "cad_footprint"
  
  writeRaster(footprint_albers,
              footprint_loc,
              filetype = "envi",
              overwrite = T)
}

footprint <- rast(footprint_loc)

# making zone all layers for covariates
if (list.files(raster_cov, pattern = ".dat$") %>% length() != 16) {
  topo <-
    list.files("F:/mosaiced/topo",
               full.names = T,
               pattern = ".dat$") %>%
    map(rast) %>%
    map(crop, forests, .progress = "crop") %>%
    rast()
  
  names(topo) <- list.files("F:/mosaiced/topo", pattern = ".dat$")
  
  r_gen <-
    list.files(here::here(raster_loc),
               full.names = T,
               pattern = ".dat$") %>%
    map(rast) %>%
    rast()
  
  pa_loc <- here::here(shapefile_loc, "bc_pa_filt.shp")
  
  bc_pa_filt <- vect(pa_loc)
  
  pa_rast <- rasterize(bc_pa_filt, forests)
  
  names(pa_rast) <- "pa"
  
  names(r_gen) <-
    list.files(here::here(raster_loc), pattern = ".dat$")
  
  #all <- c(rast(arc_tifs), rast(topo), r_gen)
  
  all <- c(topo, r_gen, pa_rast)
  
  names(all) <- str_replace(names(all), ".dat", "")
  
  all_masked <- mask(all, forests)
  
  map2(.x = as.list(all_masked)[1:2], .y = names(all_masked)[1:2], \(x, y) {
    savename <- here::here(raster_loc, glue::glue("{y}.dat"))
    print(savename)
    if (!file.exists(savename))
      writeRaster(x, savename, overwrite = T, filetype = "envi")
  })
  
  bec_loc <- here::here(shapefile_loc, "bec_terr_agg.shp")
  
  bec_vect <- vect(bec_loc)
  
  zones <- bec_vect %>% pull(ZONE)
  
  
  for (zone in zones) {
    print(zone)
    zone_save <-
      here::here(raster_cov, paste0(zone, ".dat"))
    
    if (!file.exists(zone_save)) {
      zone_vect <- bec_vect %>%
        filter(ZONE == zone)
      
      zone_rast <- rasterize(zone_vect, forests, field = "ZONE") %>%
        trim()
      
      zone_all <- crop(all_masked, zone_rast, mask = T) %>%
        trim()
      
      #zone_csv <- as.data.frame(zone_all)
      
      #write_csv(zone_csv, here::here("data", "csv", paste0(zone, ".csv")))
      
      writeRaster(zone_all,
                  zone_save,
                  filetype = "envi",
                  overwrite = T)
    }
    
  }
}

rm(all, all_masked, bc_pa_filt, pa_rast, r_gen, topo)


#### structure and dhi

struct_locs <- here::here("F://", "mosaiced", "structure")

struct_varnames <- list.files(struct_locs)

struct_rasts <- map(struct_varnames, \(x) {
  sname <- here::here(mosaic_mask_loc, glue::glue("{x}.dat"))
  print(sname)
  if (!file.exists(sname)) {
    r <- rast(here::here(struct_locs, x, glue::glue("{x}_2015.dat")))
    masked <- r %>%
      crop(y = forests, mask = T) %>%
      writeRaster(., sname, filetype = "envi", overwrite = T)
  }
  rast(sname)
}, .progress = "Structure Masking")

dhi_loc <- here::here("F://", "mosaiced", "DHI_nomask")

dhi_files <-
  list.files(dhi_loc, pattern = ".tif$", full.names = T)

dhi_rasts <- map(dhi_files, \(file) {
  x <- basename(file) %>%
    tools::file_path_sans_ext() %>%
    str_split("-") %>%
    unlist() %>%
    str_subset("DHI")
  
  sname <- here::here(mosaic_mask_loc, glue::glue("{x}.dat"))
  print(sname)
  if (!file.exists(sname)) {
    r <- rast(here::here(file))
    masked <- r %>%
      crop(y = forests, mask = T) %>%
      writeRaster(., sname, filetype = "envi", overwrite = T)
  }
  rast(sname)
}, .progress = "DHI Masking")

all_rasts <- c(struct_rasts, dhi_rasts)

all_names <- all_rasts %>%
  map(sources) %>%
  unlist() %>%
  basename() %>%
  tools::file_path_sans_ext()

names(all_rasts) <- all_names

bec_loc <- here::here(shapefile_loc, "bec_terr_agg.shp")

bec_vect <- vect(bec_loc)

zones <- bec_vect %>% pull(ZONE)

for (zone in zones) {
  print(zone)
  zone_loc <-
    here::here(raster_vars, zone)
  
  dir.create(zone_loc, showWarnings = F)
  
  zone_vect <- bec_vect %>%
    filter(ZONE == zone)
  
  zone_save <-
    here::here(raster_loc,
               "zone_masks",
               glue::glue("{zone}_mask.tif"))
  
  dir.create(dirname(zone_save), showWarnings = F)
  
  if (!file.exists(zone_save)) {
    print(glue::glue("Making {zone} mask and saving"))
    zone_rast <- rasterize(zone_vect, forests, field = "ZONE") %>%
      trim()
    
    writeRaster(zone_rast, zone_save)
  }
  zone_rast <- rast(zone_save)
  
  save_names <-
    here::here(zone_loc,
               glue::glue("{all_names}.dat"))
  
  # zone_all <- all_rasts %>%
  #   map(crop,
  #       y = zone_rast,
  #       mask = T,
  #       .progress = "Zone Cropping/masking") %>%
  #   map(trim, .progress = "Trim")
  #
  # map2(
  #   .x = zone_all,
  #   .y = save_names,
  #   .f = writeRaster,
  #   filetype = "envi",
  #   overwrite = T,
  #   .progress = "Save"
  # )
  
  map2(
    .x = all_rasts,
    .y = save_names,
    .f = \(x, y) {
      print(y)
      if (!file.exists(y)) {
        crop(x, zone_rast, mask = T) %>%
          trim() %>%
          writeRaster(
            x = .,
            filename = y,
            filetype = "envi",
            overwrite = T
          )
      }
    },
    .progress = "Crop, Trim, Save"
  )
  
}
