library(terra)
library(rgdal)
library(ggplot2)
library(dplyr)

# Subset cropped areas
ocs_file <- 'data/workfiles/ocs_mean_iowa_onlycrop.tif'
# GDALinfo(ocs_file)
ocs <- rast(ocs_file)
# ocs
# plot(ocs)
# hist(ocs)

agg <- 4
ocs2 <- aggregate(ocs, fact=agg, fun="mean", na.rm=TRUE, cores=2)
# Fifth, we write out to file
writeRaster(ocs2, paste0('data/workfiles/ocs_mean_iowa_onlycrop_agg', agg,'.tif'))
