# propensity scoring using confounders
library(terra)
library(sf)
library(tidyverse)

rasters <- list.files(raster_loc, full.names = T, pattern = ".dat$") %>% 
  map(rast) %>%
  map(crop, forests, .progress = "crop") %>%
  map(extend, forests, .progress = "extend")
