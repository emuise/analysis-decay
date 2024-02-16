# bc intact forest landscapes (IFL)
ifl_loc <- here::here(shapefile_loc, "ifl_2020.shp")

if (!file.exists(ifl_loc)) {
  ifl_url <- "https://intactforests.org/shp/IFL_2020.zip"
  ifl_dest <- here::here("data", "ifl.zip")
  
  download.file(ifl_url, ifl_dest)
  
  unzip(ifl_dest, exdir = shapefile_loc)
  rm(ifl_url, ifl_dest)
}

bc_ifl_loc <- here::here(shapefile_loc, "bc_ifl_2020.shp")

if (!file.exists(bc_ifl_loc)) {
  ifl_2020 <- ifl_loc %>%
    vect() %>% # projection in terra much faster on this dataset?
    terra::project("epsg:3005")
  
  bc_ifl <- ifl_2020 %>%
    intersect(bcb)
  
  writeVector(bc_ifl, bc_ifl_loc)
}

bc_ifl <- vect(bc_ifl_loc)

pa_ifl <- intersect(bc_ifl, bc_pa_filt)