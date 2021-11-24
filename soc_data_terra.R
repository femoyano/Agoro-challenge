# Script to subset soc data to cropped areas in Iowa and to a small part of Iowa (for testing)

library(terra)
library(rgdal)
library(ggplot2)
library(dplyr)

# Subset cropped areas ---
crop_file <- 'data/crop_mask/CDL_2020_19/CDL_2020_19.tif'
# GDALinfo(crop_file)
crop <- rast(crop_file)
cr_vals <- values(crop)
cr_uni <- unique(cr_vals)
rm(cr_vals)

# Codes for land use obtained from https://www.nass.usda.gov/Research_and_Science/Cropland/metadata/metadata_ia20.htm
notcrop <- c(0, 7, 8, 9, 15, 16, 17, 18, 19, 20, 40, 62, 63, 64, 65,
             73, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
             91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
             104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114,
             115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 
             126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136,
             137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 
             148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158,
             159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 
             170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
             181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
             192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202,
             203, 251, 252, 253, 255)
iscrop <- c(1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 21, 22, 23, 24, 25,
            26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
            41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 
            55, 56, 57, 58, 59, 60, 61, 66, 67, 68, 69, 70, 71, 72, 
            74, 75, 76, 77, 204, 205, 206, 207, 208, 209, 210, 211, 
            212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 
            223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233,
            234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 
            245, 246, 247, 248, 249, 250, 254)

cr_yes <- cr_uni[cr_uni %in% iscrop]
cr_no <- cr_uni[cr_uni %in% notcrop]

# Set all non-crop values to NA
crop[crop %in% cr_yes] <- 1
crop[crop %in% cr_no] <- 0

# plot(crop)
# writeRaster(crop, 'cropland-iowa.tif')


# Now we use the new crop raster to mask the ocs raster -----------

# First we read in the ocs mean and error data (Subset to Iowa was done in QGIS)
ocs_file <- 'data/workfiles/soilgrid_ocs_mean_iowa_crs-aea.tif'
ocs <- rast(ocs_file)
# GDALinfo(ocs_file)
# plot(ocs)
# hist(ocs)

ocserr_file <- 'data/workfiles/Iowa_ocs_0-30_Q0.05.tif'
ocserr <- rast(ocserr_file)
# GDALinfo(ocserr_file)
# hist(ocserr)
# plot(ocserr)


# Second we need to change the resolution to match ocs data and subset the ocs data using the crop mask
crop2 <- resample(crop, ocs, method='near')
ocs2 <- mask(x = ocs, mask = crop2, maskvalue=0)

crop2 <- resample(crop, ocserr, method='near')
ocserr2 <- mask(x = ocserr, mask = crop2, maskvalue=0)

# Lets convert 0.05 percentile to standard deviation assuming a normal distribution of the error:
ocserr2 <- resample(ocserr2, ocs2, method='near')
ocssd2 <- (ocserr2-ocs2)/(-1.645)  #  -1.645 value of Z at 005 percentile

# We can create an aggregated (low res) version for testing purposes
ocs3 <- aggregate(ocs2, fact=10, fun="mean", na.rm=TRUE, cores=2)
plot(ocs3)

# Fifth, we write out to file
writeRaster(ocs2, 'data/workfiles/ocs_mean_iowa_onlycrop.tif')
writeRaster(ocs3, 'data/workfiles/ocs_mean_iowa_onlycrop_agg.tif')
