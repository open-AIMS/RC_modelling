# Run the spatio-temporal FRK model and save model outputs into a list 
# @author Julie Vercelloni
# August 2025 

rm(list=ls())
#setwd("c:/Users/jvercell/OneDrive - Australian Institute of Marine Science/AIMS/01_Research projects/ReefCloud/SP_models/FRK_dev/for_andrew/scripts")

setwd(paste0(here(), "/scripts"))

# Load R package 
source("../R/packages.R")
source("../R/functions.R")

# data.grp.tier.ready <- read.csv("data/data.grp.tier.csv") %>%
#  filter(GROUP_DESC == "HARD CORAL") %>%
#  dplyr::select(P_CODE, Tier2, Tier3, Tier4, Tier5, SITE_NO, LATITUDE, LONGITUDE, fYEAR, COUNT, TOTAL, max_cyc, max_cyc.lag1, max_cyc.lag2,
#                    max_dhw, max_dhw.lag1, max_dhw.lag2)

# write.csv(data.grp.tier.ready, file = "data/data.grp.tier.ready.csv", row.names = F)

# Read the data
data.grp.tier.ready <- read.csv("../data/data.grp.tier.ready.csv")

HexPred_reefid2 <- st_read("../data/hexpred.shp") |>
 rename(max_cyc_lag1 = mx_cy_1,
        max_cyc_lag2 = mx_cy_2,
        max_dhw_lag1 = mx_dh_1,
        max_dhw_lag2 = mx_dh_2)

# Model run 
obj_frk <- frk_prep(data.grp.tier.ready, HexPred_reefid2) 

# Fit FRK model 
M <- FRK(f = COUNT ~ 1 + (1 | reefid) + max_cyc + max_cyc_lag1 + max_cyc_lag2 +
                   max_dhw + max_dhw_lag1 + max_dhw_lag2, 
         data = list(obj_frk$STObj), 
         BAUs = obj_frk$ST_BAUs, 
         basis = obj_frk$basis, 
         response = "binomial", 
         link = "logit", 
         K_type = "precision", 
         method = "TMB", 
         est_error = FALSE)

saveRDS(list(data.grp.tier = data.grp.tier.ready,
             M = M,
             obj_frk = obj_frk, 
             HexPred_reefid2 = HexPred_reefid2),
        file = "../data/FRK_fit.RData")
