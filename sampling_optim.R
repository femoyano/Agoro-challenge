# This script optimizes a SOC sampling design taking into account stratification and target confidence interval 

library(terra)
library(rgdal)
library(tidyverse)

# Settings ======
test_id     <- "test2"
ocs_file    <- 'data/workfiles/soilgrids_ocs_mean_iowa.tif'
ocserr_file <- 'data/workfiles/soilgrids_ocserr_sd_iowa.tif'
use_prederr <- TRUE

conf_lev   <- 0.95  # confidence level required
max_cif    <- 0.15  # maximum allowed confidence interval as fraction of the mean
n_mps      <- 3     # minimum number of sampling points allowed per strata
min_H      <- 2     # minimum number of strata
max_H      <- 10    # number of strata to test
max_n      <- 1000  # maximum number of total samples
sample_inc <- 1     # sample number increment per iteration. In case processing time becomes critical.


# Function definitions ========
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
    summarize(stratum = mean(stratum), stratum_sd = sd(X), stratum_var = sd(X)^2,
              stratum_size = length(X), stratum_errsd = mean(X_sd)) %>%
    drop_na(stratum_sd)
}

get_sampvar <- function(H, N_H, N, S_H, n_H) {
  # Function returns the sampling variance for a stratified design
  # H: scalar for number of strata
  # N_H: vector of stratum sizes (number of grid points)
  # N: scalar for total grid points (sum of all strata) 
  # S_H: vector of standard deviation of ocs by stratum
  # n_H: vector of sample size allocated to each stratum
  sampvar <- 0
  for(i in 1:H) {
    v_h <- (N_H[i]/N)^2 * (S_H[i])^2 / n_H[i]
    sampvar <- sampvar+v_h
  }
  return(sampvar)
}

get_sampvar2 <- function(H, N_H, N, S_H, n_H) {
  # Same as get_sampvar but with a finite population correction term (only important for small grids)
  sampvar <- 0
  for(i in 1:H) {
    v_h <- (N_H[i]/N)^2 * ((N_H[i] - n_H[i]) / N_H[i]) * ((S_H[i])^2 / n_H[i])
    sampvar <- sampvar+v_h
  }
  return(sampvar)
}

get_sampvar3 <- function(H, N_H, N, S_H, E_H, n_H) {
  # Function returns the sampling variance for a stratified design
  # H: scalar for number of strata
  # N_H: vector of stratum sizes (number of grid points)
  # N: scalar for total grid points (sum of all strata) 
  # S_H: vector of standard deviation of ocs by stratum
  # E_H: vector of mean standard deviation in ocs prediction error by stratum
  # n_H: vector of sample size allocated to each stratum
  sampvar <- 0
  for(i in 1:H) {
    v_h <- (N_H[i]/N)^2 * ((S_H[i])^2 + E_H[i]^2) / n_H[i]
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
  t <- qt(1-(1-conf_lev)/2, n_tot-1)  # the t-score from students distribution
  se <- sqrt(sampvar) / sqrt(n_tot) # standard error
  ci <- t * se * 2  # Multiplied by 2 to get the full interval length
}

optim_sampling <- function(df_in, conf_lev, max_ci, n_mps, min_H, sample_inc) {
  # This function obtains the minimum amount of total samples necessary to attain the value of max_ci
  # It also returns the sampling density for each stratum 
  # Returns a list containing: 
  #   1. dataframe with grid values and optimized strata values
  #   2. dataframe with strata numbers and optimized samples per strata
  
  # Working variables
  H_all   <- tibble(nstrata=min_H:max_H, n_opt = NA, svar = NA, cif = NA, ci = NA)  # table holding output values
  nH      <- NA  # Number of strata
  n_tot   <- NA  # Number of total samples
  ocs_avg <- mean(df_in$X, na.rm = TRUE)  # overall mean organic carbon stock
  N_tot   <- sum(!is.na(df_in$X))  # N total number of grid points
  H_levs  <- min_H:max_H  # stratification levels to test
  
  for(i in 1:length(H_levs)) {
    
    nH <- H_levs[i]
    cat("Optimizing at stratification:", nH, "\n")
    
    # if(nH == 20) browser()  # debug line
    
    # Create strata
    clust <- get_km_clust(d = df_in$X_raw, c = nH)
    df_in$stratum[which(!is.na(df_in$X_raw))] <- as.integer(clust$cluster)
    
    # Get statistics for current stratification
    H_stats <- get_stratstats(df_in)
    
    iter1 <- 1
    # target_reached <- FALSE  # not used
    
    while(TRUE) {
      
      # Calculate the total sampling size
      if(iter1==1) {n_tot <- n_mps * nH} else {n_tot <- n_tot + sample_inc}
      
      # print(n_tot)  # for debug
      
      # Stop if maximum allowed samples reached
      if(n_tot > max_n) {
        cat(paste0("Reached maximum number of allowed samples at stratification: ", nH, "\n"))
        break
      }
      
      # if(n_tot == 64 & nH == 6) browser()  # debug line
      
      # Get the optimal allocation
      opt_n <- get_optsampalloc(H = nH, N_H = H_stats$stratum_size, S_H = H_stats$stratum_sd, n_tot = n_tot)
      opt_n[opt_n < n_mps] <- n_mps  # if we don't get the minimum n per stratum, add where necessary
      H_stats$opt_n <- opt_n
      n_tot2 <- sum(H_stats$opt_n)  # recalculate to fix any difference after rounding
      
      # Get the overall sampling variance as a function of stratified sampling
      if(use_prederr) {
        samp_var <- get_sampvar3(H = nH, N_H = H_stats$stratum_size, N = N_tot, S_H = H_stats$stratum_sd,
                                 E_H = H_stats$stratum_errsd, n_H = H_stats$opt_n)
      } else {
        samp_var <- get_sampvar(H = nH, N_H = H_stats$stratum_size, N = N_tot, S_H = H_stats$stratum_sd,
                                n_H = H_stats$opt_n)  
      }
      
      # Get confidence interval
      ci <- get_confint(conf_lev = conf_lev, sampvar = samp_var, n_tot = n_tot2)
      cif <- ci / ocs_avg
      
      # Check if cif is same or lower than the maximum allowed confidence interval
      if(cif <= max_cif) {
        # target_reached <- TRUE
        H_all$n_opt[i] <- n_tot2
        H_all$svar[i] <- samp_var
        H_all$cif[i] <- cif * 100
        H_all$ci[i] <- ci
        break
      }
      
      iter1 <- iter1+1
      
    }
    
    # if(target_reached) {} 
    
  }
  
  # First select strata with lowest n_opt
  H_sel <- H_all[H_all$n_opt == min(H_all$n_opt, na.rm=TRUE),] 
  # If more than one selected, take the case with more strata (err on the safe side)
  H_best  <- H_sel[H_sel$nstrata == max(H_sel$nstrata, na.rm=TRUE),]
  
  # Recalculate clusters and sampling for the selected case ----- (deduplicate code?)
  
  # Create strata
  clust <- get_km_clust(d = df_in$X, c = H_best$nstrata[1])
  df_in$stratum[which(!is.na(df_in$X))] <- as.integer(clust$cluster)
  
  # Get statistics for best stratification
  Hb_stats <- get_stratstats(df_in)
  
  # Get the optimal allocation
  
  opt_n <- get_optsampalloc(H = H_best$nstrata[1], N_H = Hb_stats$stratum_size, S_H = Hb_stats$stratum_sd, n_tot = H_best$n_opt[1])
  opt_n[opt_n < n_mps] <- n_mps  # if we don't get the minimum n per stratum, add where necessary
  Hb_stats$opt_n <- opt_n
  
  cat(" Confidence interval equal or less than the target (" ,  max_cif * 100, "% ) was found.\n")
  cat(" CI = ", H_best$ci, "\n",
      "CIrel = ", H_best$cif * 100, "%\n",
      "Mean SOC stocks = ", ocs_avg, "\n",
      "Total strata = ", H_best$nstrata, "\n",
      "Total samples = ", H_best$n_opt, "\n",
      "Samples per strata = ", Hb_stats$opt_n, "\n"
  )
  
  return(list(Strata_out = df_in, Strata_stats = Hb_stats, Strata_all = H_all))
  
}


# Load and execute ============

# Read in the data of SOC predictions over the target area
r_ocs_m <- rast(ocs_file)
r_ocs_err <- rast(ocserr_file)

# Get the values as a matrix
ocs_raw <- data.frame(
  X_raw = as.vector(terra::values(r_ocs_m)),
  X_sd_raw = as.vector(terra::values(r_ocs_err)),
  stratum = NA )  # add a variable to hold stratum values

# Check the data
hist(ocs_raw$X_raw)
hist(ocs_raw$X_sd_raw)

# Log transform to get normal distributions
ocs_tran <- ocs_raw
small <- 0.01 # add to avoid 0 values for log transform
ocs_tran <- ocs_tran %>%
  mutate(X = log(X_raw+small), X_sd =  log(X_sd_raw+small))

# Call the sampling optimization function
optim_out <- optim_sampling(df_in = ocs_tran, conf_lev = conf_lev, max_ci = max_ci,
                           n_mps = n_mps, min_H = min_H, sample_inc = sample_inc)

print(optim_out[[2]])
print(optim_out[[3]])

# Create raster of strata
ocs_H <- terra::setValues(r_ocs_m, optim_out[[1]][[2]])

# Save results ======
info_file <- paste0("data/workfiles/", test_id, "_strata_optim_info.csv")
write_csv(optim_out[[2]], info_file)
strata_file <- paste0("data/workfiles/", test_id, "_ocs_strata.tif")
terra::writeRaster(x = ocs_H, strata_file, overwrite = TRUE)
