# Script to explore tiff files

library(raster)
library(rgdal)
library(ggplot2)
library(dplyr)

# tif_file <- 'data/workfiles/soilgrids_ocs_uncertainty_central-iowa.tif'
tif_file <- 'data/get_soilgrids_data/Iowa_ocs_0-30_Q0.05_2.tif'
GDALinfo(tif_file)
tif <- raster(tif_file)
# hist(tif)
plot(tif)
# tif2 <- reclassify(tif, cbind(5000, 50000, NA))
# # hist(tif2)
# plot(tif2)

ocs_file <- 'data/workfiles/soilgrid_ocs_mean_iowa_crs-aea.tif'
GDALinfo(ocs_file)
ocs <- raster(ocs_file)
ocs
plot(ocs)
hist(ocs)

ocs_file <- 'data/soil/soilgrid_ocs_mean.tif'
GDALinfo(ocs_file)
ocs <- raster(ocs_file)
ocs
plot(ocs)

crop_file <- 'data/crop_mask/CDL_2020_19/CDL_2020_19.tif'
# GDALinfo(crop_file)
crop <- raster(crop_file)
crop

# It's hard to see what the area coverage is just from looking at the description.
# Lets plot to have a better idea
plot(ocs)
plot(crop)

# Subset to Iowa
# We will use the crop raster data as a mask to subset the ocs data

# Subset to part of Iowa

