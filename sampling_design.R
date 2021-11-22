# This script optimizes a SOC sampling design taking into account stratification and target confidence interval 

library(raster)
library(rgdal)
library(tidyverse)

# Definitions and settings
conf_lev   <- 0.95  # confidence level required
max_ci     <- 0.15  # maximum allowed confidence interval as proportion of the mean
n_mps      <- 10    # minimum number of sampling points allowed per strata
min_strata <- 2     # minimum number of strata
max_strata <- 100   # maximum number of strata
max_n      <- 1000  # maximum number of total samples
sample_inc <- 5     # sample number increment per iteration. For large grids (e.g. Iowa) larger value can save processing time.

# Required functions
get_km_clust <- function(d, c) {
  # Function returns a clustering based on k_means
  # d: data to be clustered
  # c: number of centriods (clusters)
  set.seed(42)
  kmeans(x = na.exclude(d), centers = c, iter.max = 500, nstart = 5) #, algorithm = "Lloyd")
}

get_stratstats <- function(df) {
  df_out <- df %>%
    group_by(stratum) %>%
    summarize(stratum = mean(stratum), S = sd(ocs_m), V = sd(ocs_m)^2, N = length(ocs_m)) %>%
    drop_na(S)
}

get_sampvariance <- function(H, N_H, N, S_H, n_H) {
  # Function returns the sampling variance for a stratified design
  # H: scalar for number of strata
  # N_H: vector size H of size (number of grid points) of stratum h
  # N: scalar for total grid points (sum of all strata) 
  # S_H: vector size H of standard deviation in each stratum
  # n_H: vector size H of sample size allocated to stratum h
  sampvar <- 0
  for(i in 1:H) {
    v_h <- (N_H[i]/N)^2 * S_H[i]^2 / n_H[i]
    sampvar <- sampvar+v_h
  }
  return(sampvar)
}

get_confint <- function(conf_lev, sampvar, n_tot) {
  # Two-sided confidence interval calculation
  z <- qnorm(1-(1-conf_lev)/2)  # use of the Z value assumes a large enough number of samples (> 30)
  se <- sqrt(sampvar) / sqrt(n_tot) # standard error
  ci <- z * se * 2  # Multiplied by 2 to get the full interval length
}

get_optsampalloc <- function(H, N_H, S_H, n_tot) {
  # Optimal sample distribution (sizes) for all strata h n’_h = n * (N_h * S_h) / (SUM(h=1 -> H) {N_h * S_h})
  # H: scalar for number of strata
  # N_H: vector size H of size (total number of grid points) of stratum h
  # S_H: vector size H of standard deviations for each stratum

  sum_sdn <- 0  # sum of standard deviations times number of samples
  for(i in 1:H) {
    x <- N_H[i] * S_H[i]
    sum_sdn <- sum_sdn + x
  }
  opt_H <- rep(NA, length = H)  # vector to hold outputs
  for(i in 1:H) {
    opt_H[i] <- n_tot * (N_H[i] * S_H[i]) / sum_sdn
  }
  return(as.integer(opt_H))
}


optim_sampling <- function(ocs_pr, conf_lev, max_ci, n_mps, min_strata, sample_inc) {
  # This function obtains the minimum amount of total samples necessary to attain the value of max_ci
  # It also returns the sampling density for each stratum 
  # Returns a list containing: 
  #   1. dataframe with grid values and optimized strata values
  #   2. dataframe with strata numbers and optimized samples per strata
  
  # Working variables
  n_opt       <- NA  # Optimal number of total samples
  nstrata_opt <- NA  # Optimal number of strata
  nstrata     <- NA  # Number of strata
  n_tot       <- NA  # Number of total samples
  rci         <- NA  # Confidence interval relative to the mean
  ocs_avg     <- mean(ocs_pr$ocs_m, na.rm = TRUE)
  
  iter1 <- 1 
  
  while(TRUE) {
    
    if(iter1==1) {nstrata <- min_strata} else {nstrata <- nstrata+1}
    
    # Create strata
    clust <- get_km_clust(d = ocs_pr$ocs_m, c = nstrata)
    ocs_df$stratum[which(!is.na(ocs_df$ocs_m))] <- as.integer(clust$cluster)
    
    # Get statistics for each stratum
    str_stats <- get_stratstats(ocs_df)
    
    iter2 <- 1
    
    while(TRUE) {
      
      # Calculate the total sampling size
      if(iter2==1) {n_tot <- n_mps * nstrata} else {n_tot <- n_tot + sample_inc}
      
      # Get the optimal allocation
      opt_alloc <- get_optsampalloc(H = nstrata, N_H = str_stats$N, S_H = str_stats$S, n_tot = n_tot)
      
      if(any(opt_alloc < n_mps)) next  # if we don't get the minimum n per stratum, skip to next iteration
      
      
      
      if(n >= max_n) {
        print(paste("Reached maximum number of allowed samples at stratum", nstrata))
        break
        }
        
      iter2 <- iter2+1
    }
    
    if(strata == max_strata) {
      print(paste("Reached maximum number of allowed strata:", strata))
      break
    }
    iter1 <- iter1+1
  }
  
  strat_dat <- data.frame(stratum = , n = )
  return(list(ocs_pr, strat_dat))
}


# Data Analysis ============

# Read in the data of SOC predictions over the target area
ocs_file <- 'data/workfiles/ocs_mean_iowa_onlycrop_agg.tif'
ocs_r <- raster(ocs_file)

# Get the values as a matrix
ocs_df <- data.frame(ocs_m = getValues(ocs_r), stratum = NA)  # add a variable to hold stratum values

# Call the sampling optimization function
samp_opt <- optim_sampling(ocs_pr = ocs_df, conf_lev = conf_lev, max_ci = max_ci,
                           n_mps = n_mps, min_strata = min_strata, sample_inc = sample_inc)


# Transform back into raster
ocskm_r <-raster(ocs_r)
values(ocskm_r) <- ocs_df$ocs_clus
# mycolor <- c("#fef65b","#ff0000", "#daa520","#0000ff","#0000ff","#00ff00","#cbbeb5",
             # "#c3ff5b", "#ff7373", "#00ff00", "#808080")
mycolor <- rainbow(5)
plot(ocskm_r, col = mycolor)

