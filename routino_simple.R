########################################################################################
# Project           : AHAH Index 
# Program name      : routino_simple.R
# Author            : Konstantinos Daras - konstantinos.daras@gmail.com
# Date created      : 3 Oct 2016
# Purpose           : Calculates network distance between each postcode and nearest POI.
# 
# Revision History  --------------------------------------------------------------------
# Date         Author      Ref    Revision (Date in YYYYMMDD format) 
# 3 Oct 2016   KD          1      Original code. 
########################################################################################

library(data.table)
library(FNN)
library(splitstackshape)

# Local paths
setwd("/media/kd/Data/Dropbox/Repos/iahc_index/off_l_net/")

# NSPL Lookup table
#---------------------------
# pcd: Postcode trimmed, 
# pcd2: Postcode - 2 spaces,
# pcds: Postcode - 1 space,
# lat: Latitute,
# long: Longitute,
# oseast1m: Easting coord,
# osnrth1m: Northing coord
# --------------------------
nspl_file <- "nspl_pcd.csv"

# POIs file
#---------------------------
# postcode: Postcode trimmed,
# poi_id: POI identification code (Character),
# pcd2: Postcode - 2 spaces,
# pcds: Postcode - 1 space,
# lat: Latitute,
# long: Longitute,
# oseast1m: Easting coord,
# osnrth1m: Northing coord
# --------------------------
pois_file <- "pois_pcd.csv"

# Output file
out_file <- "output.csv"

# Set Routino parameters
transport <- "motorcar"
# routino_profile <- "/media/kd/Data/Dropbox/Repos/netdist_pc/routino-profiles.xml"
# transport <- "foot"
routino_profile <- "/media/kd/Data/Dropbox/Repos/iahc_index/netdist_pc/routino-profile-norules.xml"
rout_path <- "/media/kd/Data/routino_data"   # OSM database path

# Load Postcode lookup table
nspl_data <- read.csv(nspl_file, sep=",")
nspl_data$pcd <- as.character(nspl_data$pcd)

#Load POIs Practice table
pois_pcd <- read.csv(pois_file, sep = ",")
pois_pcd$postcode <- as.character(pois_pcd$postcode)

# Create unique pcd codes with lat/long cols
ll_clean <- pois_pcd[!duplicated(pois_pcd$postcode),]
ll_clean$postcode <- as.character(ll_clean$postcode)

# Prepare easting/northing cols for knn 
ll_mtx <- ll_clean[,c(7,8)]
pois_mtx <- nspl_data[,c(6,7)]

# Run knn for k=10 neighbours
knn_res <- knnx.index(ll_mtx, pois_mtx, k=10)

# Loop for each postcode 
ptm <- proc.time()
for (r in 1:nrow(nspl_data)){

  long_from <- nspl_data[r,5]
  lat_from <- nspl_data[r,4]
  
  near_pcs <- data.frame("postcode" = character(), "long"= numeric(), "lat"= numeric())
  near_pcs$postcode <-as.character(near_pcs$postcode)

  for (n in 1:ncol(knn_res)){
    row_loc <- knn_res[r,n]
    near_pcs[n,1] <- ll_clean[row_loc,1]
    near_pcs[n,2] <- ll_clean[row_loc,6]
    near_pcs[n,3] <- ll_clean[row_loc,5]
  }

  for (neig in 1:nrow(near_pcs)){
    long_to <- near_pcs[neig,2]
    lat_to <- near_pcs[neig,3]
    
    router <- paste0("routino-router --transport=",
                     transport, 
                     " --prefix=gb --quickest --profiles=",
                     routino_profile,
                     " --lon1=", 
                     long_from, 
                     " --lat1=", 
                     lat_from, 
                     " --lon2=",
                     long_to, 
                     " --lat2=",
                     lat_to, 
                     " --output-text-all --output-stdout --quiet --dir=",
                     rout_path)

    sysmsg <- system(router, intern=TRUE)  # Send the routing command
    sysmsg <- as.data.frame(sysmsg)
    if (nrow(sysmsg)>0){
      sysmsg <- sysmsg[-c(1,2,3,4,5,6,7),] 
      sysmsg <- as.data.frame(sysmsg)
      sysmsg <- cSplit(sysmsg, "sysmsg", "\t")
      sysmsg <- as.data.frame(sysmsg)
      
      # Get the distance/time of the route
      km <- max(as.numeric(sysmsg$sysmsg_07))
      t <- max(as.numeric(sysmsg$sysmsg_08))
      
      near_pcs[neig,"dist"] <- km
      near_pcs[neig,"time"] <- t
    }
    else {
      near_pcs <- near_pcs[0,]
    }
  }
  
  if (nrow(near_pcs)>0) {
    # Get the minimum distance
    dst <- min(as.numeric(as.character(near_pcs$dist)))
    pc <- near_pcs[near_pcs$dist==dst,1][1]
    km <- dst
    t <- near_pcs[near_pcs$dist==dst,5][1]
    
    nspl_data[r,"poi_pc"] <- pc
    nspl_data[r,"dist"] <- km
    nspl_data[r,"time"] <- t
    
    print(paste("Postcode:",nspl_data[r,1], "|", km, "km", "|", t, "mins"))
  }
  else {
    print(paste("Postcode:",nspl_data[r,1], "| NO ROUTE AVAILABLE"))
  }
}
proc.time() - ptm

write.csv(nspl_data, out_file, quote = FALSE, row.names = FALSE)
