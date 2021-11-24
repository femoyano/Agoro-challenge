# Get ocs uncertainty

tif_file <- 'data/get_soilgrids_data/Iowa_ocs_0-30_Q0.05.tif'
GDALinfo(tif_file)
tif <- raster(tif_file)
hist(tif)
plot(tif)