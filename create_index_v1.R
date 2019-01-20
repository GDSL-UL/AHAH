########################################################################################
# Project           : AHAH Index 
# Program name      : create_index_v1.R
# Author            : Konstantinos Daras - konstantinos.daras@gmail.com
# Date created      : 3 Oct 2016
# Purpose           : Calculates the domains and the overall AHAH index.
# 
# Revision History  --------------------------------------------------------------------
# Date         Author      Ref    Revision (Date in YYYYMMDD format) 
# 3 Oct 2016   KD          1      Original code. 
########################################################################################

library(nFactors)

# Quantile calcs
idr <- function(x){((x-quantile(x, 0.5, na.rm=TRUE)))/
    ((quantile(x, 0.95, na.rm=TRUE))-(quantile(x, 0.05, na.rm=TRUE)))}

# Exponential transformation  ~>   X=-23 ln(1-R(1-exp(-100/23)))
exp_trans <- function(x,y){-23*log(1-(x/nrow(y))*(1-exp(-100/23)), base = exp(1))}
exp_default <- function(x,y){(x-0.5)/nrow(y)} # Use of Rankit rank-based normalisation


# Dataset 
data <- read.csv("/ahahinputs.csv", sep = ",")


# Ranking
data$gpp_dist <- rank(data$gpp_dist,ties.method= "min")
data$ed_dist <- rank(data$ed_dist,ties.method= "min")
data$dent_dist <- rank(data$dent_dist,ties.method= "min")
data$pharm_dist <- rank(data$pharm_dist,ties.method= "min")
data$gamb_dist <- rank(data$gamb_dist,ties.method= "min")
data$gamb_dist <- rank(-data$gamb_dist) # Invert ranking
data$ffood_dist <- rank(data$ffood_dist,ties.method= "min")
data$ffood_dist <- rank(-data$ffood_dist) # Invert ranking
data$pubs2_dist <- rank(data$pubs2_dist,ties.method= "min")
data$pubs2_dist <- rank(-data$pubs2_dist) # Invert ranking
data$off2_dist <- rank(data$off2_dist,ties.method= "min")
data$off2_dist <- rank(-data$off2_dist) # Invert ranking
data$tobac_dist<- rank(data$tobac_dist,ties.method= "min")
data$tobac_dist <- rank(-data$tobac_dist) # Invert ranking
data$leis_dist <- rank(data$leis_dist,ties.method= "min")
data$green900 <- rank(data$green900,ties.method= "min")
data$green900 <- rank(-data$green900) # Invert ranking
data$no2 <- rank(data$no2,ties.method= "min")
data$pm10 <- rank(data$pm10,ties.method= "min")
data$so2 <- rank(data$so2,ties.method= "min")

# Exponential transformation
data$gpp_dist <- exp_default(data$gpp_dist, data)
data$gpp_dist <- qnorm(data$gpp_dist, mean = 0, sd = 1)
data$ed_dist <- exp_default(data$ed_dist, data)
data$ed_dist <- qnorm(data$ed_dist, mean = 0, sd = 1)
data$dent_dist <- exp_default(data$dent_dist, data)
data$dent_dist <- qnorm(data$dent_dist, mean = 0, sd = 1)
data$pharm_dist <- exp_default(data$pharm_dist, data)
data$pharm_dist <- qnorm(data$pharm_dist, mean = 0, sd = 1)
data$gamb_dist <- exp_default(data$gamb_dist, data)
data$gamb_dist <- qnorm(data$gamb_dist, mean = 0, sd = 1)
data$ffood_dist <- exp_default(data$ffood_dist, data)
data$ffood_dist <- qnorm(data$ffood_dist, mean = 0, sd = 1)
data$pubs2_dist <- exp_default(data$pubs2_dist, data)
data$pubs2_dist <- qnorm(data$pubs2_dist, mean = 0, sd = 1)
data$leis_dist <- exp_default(data$leis_dist, data)
data$leis_dist <- qnorm(data$leis_dist, mean = 0, sd = 1)
data$green900 <- exp_default(data$green900, data)
data$green900 <- qnorm(data$green900, mean = 0, sd = 1)
data$off2_dist <- exp_default(data$off2_dist, data)
data$off2_dist <- qnorm(data$off2_dist, mean = 0, sd = 1)
data$tobac_dist <- exp_default(data$tobac_dist, data)
data$tobac_dist <- qnorm(data$tobac_dist, mean = 0, sd = 1)
data$no2 <- exp_default(data$no2, data)
data$no2 <- qnorm(data$no2, mean = 0, sd = 1)
data$pm10 <- exp_default(data$pm10, data)
data$pm10 <- qnorm(data$pm10, mean = 0, sd = 1)
data$so2 <- exp_default(data$so2, data)
data$so2 <- qnorm(data$so2, mean = 0, sd = 1)

# Write data 
saveRDS(data,"norm_data_all_variables.rds")

# Load data
data <- readRDS("norm_data_all_variables.rds")

# Calculate weights

# Domain scores
data$r_domain <- (0.20 * data$gamb_dist +
                  0.20 * data$ffood_dist +
                  0.20 * data$pubs2_dist +
                  0.20 * data$off2_dist +
                  0.20 * data$tobac_dist)

data$h_domain <- (0.20 * data$gpp_dist +
                  0.20 * data$ed_dist +
                  0.20 * data$dent_dist +
                  0.20 * data$pharm_dist +
                  0.20 * data$leis_dist)

data$e_domain <- (0.25 * data$green900 + 
                  0.25 * data$no2 + 
                  0.25 * data$pm10 + 
                  0.25 * data$so2)

# Domain ranks
data$r_rank <- rank(data$r_domain,ties.method= "min")
#data$r_rank <- rank(-data$r_rank) # Inverse ranking
data$h_rank <- rank(data$h_domain,ties.method= "min")
#data$h_rank <- rank(-data$h_rank) # Inverse ranking
data$e_rank <- rank(data$e_domain,ties.method= "min")
#data$e_rank <- rank(-data$e_rank) # Inverse ranking


# Exp domains
data$r_exp <- exp_trans(data$r_rank,data)
data$h_exp <- exp_trans(data$h_rank,data)
data$e_exp <- exp_trans(data$e_rank,data)

# AHAH score
data$ahah <- ((1/3) * data$r_exp + 
              (1/3) * data$h_exp +
              (1/3) * data$e_exp)


write.csv(data,"ahahdomainsindex.csv",quote = FALSE, row.names = FALSE)


