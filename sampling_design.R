# This script optimizes a SOC sampling design taking into account stratification and target confidence interval 

library(raster)
library(rgdal)

# Definitions and settings
conf_lev   <- 0.95  # confidence level required
max_ci     <- 0.15  # maximum allowed confidence interval as proportion of the mean
n_mps      <- 10    # minimum number of sampling points allowed per strata
min_strata <- 2     # minimum number of strata

# Required functions
Get_kmeans <- function(d, c) {
  # Function that returns a clustering based on k_means
  set.seed(42)
  kmeans(x = na.exclude(ocs_df$d), centers = c, iter.max = 500, nstart = 5) #, algorithm = "Lloyd")
}

Optim_sampling <- function(d, conf_lev, max_ci, n_mps, min_strata) {
  # This function obtains the minimum amount of total samples necessary to attain the value of max_ci
  
  
}

# Read in the data of SOC predictions over the target area
ocs_file <- 'data/workfiles/ocs_mean_iowa_onlycrop_agg.tif'
ocs_r <- raster(ocs_file)

# Get the values as a matrix
ocs_df <- data_frame(ocs_m = getValues(ocs_r), ocs_clus = NA)

ocs_df$ocs_clus[which(!is.na(ocs_df$ocs_m))] <- as.integer(km_ocs$cluster)

ocskm_r <-raster(ocs_r)
values(ocskm_r) <- ocs_df$ocs_clus
# mycolor <- c("#fef65b","#ff0000", "#daa520","#0000ff","#0000ff","#00ff00","#cbbeb5",
             # "#c3ff5b", "#ff7373", "#00ff00", "#808080")
mycolor <- rainbow(5)
plot(ocskm_r, col = mycolor)

