import xarray as xr
import os
import numpy as np
import rioxarray as rio

clim_dir = "E:/Sync/Masters/analysis_03_decay/data/rasters"

clim = [os.path.join(clim_dir, file) for file in os.listdir(clim_dir) if file.endswith(".dat") and file.startswith("clim")]

cov_paths = clim + ["E:/Sync/Masters/analysis_03_decay/data/rasters/DEM.dat", "E:/Sync/Masters/analysis_03_decay/data/rasters/slope.dat"]

for file in cov_paths:
    rast = xr.open_dataset(file, engine = "rasterio").to_array()
    # Clip the raster with the forest raster
    rast.quantile(np.arange(0, 1.1, 0.1))