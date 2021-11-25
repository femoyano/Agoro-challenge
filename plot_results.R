# Plot results

library(raster)
library(rgdal)
library(tidyverse)
library(viridis)  # better colors for everyone
library(ggthemes) # theme_map()

test_id    <- "test1"

info_file <- paste0("data/workfiles/", test_id, "_strata_optim_info.csv")
strata_info <- read_csv(info_file)
strata_file <- paste0("data/workfiles/", test_id, "_ocs_strata.tif")
strata_r <- raster(strata_file)

# # mycolors <- c("#fef65b","#ff0000", "#daa520","#0000ff","#0000ff","#00ff00","#cbbeb5",
#              # "#c3ff5b", "#ff7373", "#00ff00", "#808080")
mycolors <- rainbow(nrow(strata_info))
# mycolors <- terrain.colors(nrow(strata_info))

plot(strata_r, col = mycolors, legend = FALSE)
legend(title = "Strata", "topright", legend = strata_info$stratum, fill = mycolors)
