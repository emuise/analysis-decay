library(terra)
library(tidyverse)
library(tidyterra)

shapefile_loc <- here::here("data", "shapefiles")
raster_loc <- here::here("data", "rasters")

forests <- rast(here::here(raster_loc, "forests.dat"))

bcb <- vect(here::here("data", "shapefiles", "bcb.shp"))

topo <-
  list.files("F:/mosaiced/topo",
             full.names = T,
             pattern = ".dat$") %>%
  map(rast) %>%
  map(crop, forests) %>%
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

bec_loc <- here::here(shapefile_loc, "bec_terr_agg.shp")

bec_vect <- vect(bec_loc)

zones <- bec_vect %>% pull(ZONE)

dir.create(here::here(raster_loc, "zone_all"))
dir.create(here::here("data", "csv"))

zone <- "CWH"

for (zone in zones) {
  print(zone)
  zone_save <-
    here::here(raster_loc, "zone_all", paste0(zone, ".dat"))
  
  if (!file.exists(zone_save)) {
    zone_vect <- bec_vect %>%
      filter(ZONE == zone)
    
    zone_rast <- rasterize(zone_vect, forests, field = "ZONE") %>%
      trim()
    
    zone_all <- crop(all_masked, zone_rast, mask = T)
    
    #zone_csv <- as.data.frame(zone_all)
    
    #write_csv(zone_csv, here::here("data", "csv", paste0(zone, ".csv")))
    
    writeRaster(zone_all,
                zone_save,
                filetype = "envi",
                overwrite = T)
  }
  
}

rm(all, all_masked, bc_pa_filt, pa_rast, r_gen, topo)

dat_info <- list.files(here::here(raster_loc, "zone_all"), pattern = ".dat$", full.names = T) %>%
  map_dfr(file.info) %>%
  arrange(size) %>%
  rownames_to_column(var = "name") %>%
  tibble()

ordered_rasts <- dat_info %>% 
  pull(name)

for(rast in ordered_rasts) {
  r <- rast(rast)
  
  df <- as.data.frame(r, xy = T)
}