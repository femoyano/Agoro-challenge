# This script optimizes a SOC sampling design taking into account stratification and target confidence interval 

library(raster)
library(rgdal)
library(tidyverse)

# Definitions and settings
conf_lev   <- 0.95  # confidence level required
max_cif    <- 0.15  # maximum allowed confidence interval as fraction of the mean
n_mps      <- 3     # minimum number of sampling points allowed per strata
min_strata <- 2     # minimum number of strata
max_strata <- 100   # maximum number of strata
max_n      <- 1000  # maximum number of total samples
sample_inc <- 1     # sample number increment per iteration. For large grids (e.g. Iowa) larger value can save processing time.

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

get_sampvar <- function(H, N_H, N, S_H, n_H) {
  # Function returns the sampling variance for a stratified design
  # H: scalar for number (of strata
  # N_H: vector size H of size (number of grid points) of stratum h
  # N: scalar for total grid points (sum of all strata) 
  # S_H: vector size H of standard deviation in each stratum
  # n_H: vector size H of sample size allocated to stratum h
  sampvar <- 0
  for(i in 1:H) {
    v_h <- (N_H[i]/N)^2 * (S_H[i])^2 / n_H[i]
    sampvar <- sampvar+v_h
  }
  return(sampvar)
}

get_sampvar2 <- function(H, N_H, N, S_H, n_H) {
  # Same as get_sampvar but with a finite population correction term (not significant in this case)
  sampvar <- 0
  for(i in 1:H) {
    v_h <- (N_H[i]/N)^2 * ((N_H[i] - n_H[i]) / N_H[i]) * ((S_H[i])^2 / n_H[i])
    sampvar <- sampvar+v_h
  }
  return(sampvar)
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
  opt_n <- rep(NA, length = H)  # vector to hold outputs
  for(i in 1:H) {
    opt_n[i] <- n_tot * (N_H[i] * S_H[i]) / sum_sdn
  }
  return(as.integer(round(opt_n)))
}

get_confint <- function(conf_lev, sampvar, n_tot) {
  # Two-sided confidence interval calculation
  z <- qnorm(1-(1-conf_lev)/2)  # use of the Z value assumes a large enough number of samples (> 30)
  se <- sqrt(sampvar) / sqrt(n_tot) # standard error
  ci <- z * se * 2  # Multiplied by 2 to get the full interval length
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
  cif         <- NA  # Confidence interval as fraction of the mean
  ocs_avg     <- mean(ocs_pr$ocs_m, na.rm = TRUE)
  N           <- length(!is.na(ocs_pr$ocs_m))
  
  iter1 <- 1 
  
  while(TRUE) {
    
    if(iter1==1) {nstrata <- min_strata} else {nstrata <- nstrata+1}
    if(nstrata > max_strata) {
      cat(paste("Reached maximum number of allowed strata:", strata))
      break
    }
    
    # Create strata
    clust <- get_km_clust(d = ocs_pr$ocs_m, c = nstrata)
    ocs_pr$stratum[which(!is.na(ocs_pr$ocs_m))] <- as.integer(clust$cluster)
    
    # Get statistics for each stratum
    str_stats <- get_stratstats(ocs_pr)
    
    iter2 <- 1
    target_reached <- FALSE
    
    while(TRUE) {
      
      # Calculate the total sampling size
      if(iter2==1) {n_tot <- n_mps * nstrata} else {n_tot <- n_tot + sample_inc}
      if(n_tot > max_n) {
        cat(paste0("Reached maximum number of allowed samples at stratum: ", nstrata))
        break
      }
      
      # Get the optimal allocation
      opt_n <- get_optsampalloc(H = nstrata, N_H = str_stats$N, S_H = str_stats$S, n_tot = n_tot)
      if(any(opt_n < n_mps)) { # if we don't get the minimum n per stratum, skip to next iteration
        iter2 <- iter2+1
        next
        }  
      str_stats$opt_n <- opt_n
      n_tot <- sum(str_stats$opt_n)  #redefine in case of rounding differences
      
      # Get the overall sampling variance as a function of stratified sampling
      samp_var <- get_sampvar2(H = nstrata, N_H = str_stats$N, N = N, S_H = str_stats$S, n_H = str_stats$opt_n)
      
      # Get confidence interval
      ci <- get_confint(conf_lev = conf_lev, sampvar = samp_var, n_tot = n_tot)
      cif <- ci / ocs_avg
      
      # Check if cif is same or lower than the maximum allowed confidence interval
      if(cif <= max_cif) {
        target_reached <- TRUE
        break
      }
      
      iter2 <- iter2+1
    }
    
    if(target_reached) {
      n_opt <- n_tot
      nstrata_opt <- nstrata
      
      cat(" Confidence interval equal or less than the target (" ,  max_cif * 100, ") was found.\n")
      cat(" CI = ", ci, "\n",
          "CIrel = ", cif * 100, "%\n",
          "Mean SOC stocks = ", ocs_avg, "\n",
          "Total strata = ", nstrata, "\n",
          "Total samples = ", n_tot, "\n",
          "Samples per strata = ", str_stats$opt_n, "\n"
          )
      return(list(ocs_out = ocs_pr, str_stats = str_stats))
    }
    iter1 <- iter1+1
  }

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

