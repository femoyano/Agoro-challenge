# Script to explore tiff files

library(raster)
library(rgdal)
library(ggplot2)
library(dplyr)

# tif_file <- 'data/climate/1988_2017/kc_1988_2017.tif'

ocs_file <- 'data/soil/soilgrid_ocs_mean.tif'
GDALinfo(ocs_file)
ocs <- raster(ocs_file)
ocs

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

