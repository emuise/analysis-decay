library(tidyverse)

shapefile_loc <- here::here("data", "shapefiles")
raster_loc <- here::here("data", "rasters")
scratch <- here::here(raster_loc, "scratch")
dir.create(shapefile_loc, recursive = T, showWarnings = F)
dir.create(raster_loc, recursive = T, showWarnings = F)
dir.create(scratch, recursive = T, showWarnings = F)

zone_all_loc <- here::here(raster_loc, "zone_all")
csv_loc <- here::here("data", "csv")

raster_cov <- here::here(zone_all_loc, "covariates")
csv_cov <- here::here(csv_loc, "covariates")

dir.create(raster_cov, showWarnings = F)
dir.create(csv_cov, showWarnings = F)

raster_vars <- here::here(zone_all_loc, "variables")
csv_vars <- here::here(csv_loc, "variables")

dir.create(raster_vars, showWarnings = F)
dir.create(csv_vars, showWarnings = F)

csv_all <- here::here(csv_loc, "all")
dir.create(csv_all, showWarnings = F)


mosaic_mask_loc <- here::here(raster_loc, "bc_mosaic_masked")
dir.create(mosaic_mask_loc, showWarnings = F)

split_df_loc <- here::here(csv_loc, "split_dfs")
dir.create(split_df_loc)
