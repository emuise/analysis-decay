library(tidyverse)
library(terra)

all_rasts_loc <- here::here("data", "rasters", "structure", "all_rasts_bec.dat")

if (!file.exists(all_rasts_loc)) {
  # bec_rast <- bec() %>%
  #   rasterize(y = vlce, field = "ZONE") %>%
  #   crop(bcb_rast) %>%
  #   mask(bcb_rast)
  # 
  # names(bec_rast) <- "bec"
  
  struct_locs <- here::here("F://", "mosaiced", "structure")
  struct_rasts <- list.files(
    struct_locs,
    recursive = T,
    pattern = ".dat$",
    full.names = T
  ) %>%
    map(rast) %>%
    rast()
  
  struct_names <- str_split(sources(struct_rasts), pattern = "/") %>%
    lapply("[[", 4) %>% # get fourth index from names
    unlist()
  
  names(struct_rasts) <- struct_names
  
  struct_rasts <- struct_rasts %>%
    crop(forests, mask = T)
  
  dhi_rasts <-
    list.files(here::here("F://", "mosaiced", "DHI"),
               pattern = "2010s.tif$",
               full.names = T) %>%
    map(rast) %>%
    rast()
  
  dhi_names <- str_split(sources(dhi_rasts), pattern = "/") %>%
    lapply("[[", 4) %>% # get fourth index from names
    unlist() %>%
    str_split(pattern = "[[:punct:]]") %>%
    lapply("[[", 2) %>% # get second index from split names
    unlist()
  
  names(dhi_rasts) <- dhi_names
  
  dhi_rasts <- dhi_rasts %>%
    crop(bcb_rast, mask = T)
  
  # dhi_mm <- minmax(dhi_rasts)
  #
  # dhi_rasts_mm <- (dhi_rasts - dhi_mm[1,]) / (dhi_mm[2,] - dhi_mm[1,])
  
  all_rasts <- c(vlce, struct_rasts, dhi_rasts) %>%
    mask(bcb_rast)
  
  writeRaster(
    all_rasts,
    all_rasts_loc,
    overwrite = T,
    filetype = "envi",
    datatype = "FLT4S"
  )
}