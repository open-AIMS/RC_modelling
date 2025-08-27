# Run the explorative analysis of the use of different basis functions
# @author Julie Vercelloni
# August 2025 

#setwd("c:/Users/jvercell/OneDrive - Australian Institute of Marine Science/AIMS/01_Research projects/ReefCloud/SP_models/FRK_dev/for_andrew/scripts")
setwd(paste0(here(), "/scripts"))

rm(list=ls())
# Load R package 
source("../R/packages.R")
source("../R/functions.R")

# Read model outputs 
mod_ <- readRDS("../data/FRK_fit.RData")
.x <- mod_
obj_frk = .x$obj_frk
hexpred = .x$HexPred_reefid2

# Find locations of temporal knots 
orign_knots <- obj_frk$basis@Basis2@df    
p2_1 <- show_basis(obj_frk$basis@Basis2)

########## 1. Select 5 knots 
# Spatial remains similar 
G_spatial <- auto_basis(manifold = plane(), 
                        data= obj_frk$ST_BAUs,
                        #nres = 2L, # for development 
                        nres = 3L, # for final run 
                        prune= 0)

knots <- sample(orign_knots$loc1,5)

G_temporal <- local_basis(manifold=real_line(),      # functions on real line
                          loc = matrix(knots),   # location parameter
                          scale = rep(2,length(knots)),          # scale parameter
                          type = "Gaussian")
basis <- TensorP(G_spatial,G_temporal)
p2_2 <- show_basis(basis@Basis2)

start_time <- Sys.time()
M_5K <- FRK(f = COUNT ~ 1 + (1 | reefid) + max_cyc + max_cyc_lag1 + max_cyc_lag2 + max_dhw + max_dhw_lag1 + max_dhw_lag2, 
            data = list(obj_frk$STObj), 
            BAUs = obj_frk$ST_BAUs, 
            basis = basis, 
            response = "binomial", 
            link = "logit", 
            K_type = "precision", 
            method = "TMB", 
            est_error = FALSE)
end_time <- Sys.time()
time_5K <- end_time - start_time

pred_5K <- predict(M_5K, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_5K <- as.data.frame(pred_5K$MC$mu_samples) %>% 
  mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = obj_frk$ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  mutate(knots = "Model 5K") %>%
tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc, knots),
                    names_to = "draw", 
                    values_to = "pred"
)

saveRDS(list(M = M_5K,
             post_dist_df = post_dist_df_5K,
             time = time_5K,
             basis = basis),
        file = "../data/Basis_function_5K.RData")

########## 2. Select 3 knots 
knots <- sample(orign_knots$loc1,3)

G_temporal <- local_basis(manifold=real_line(),      # functions on real line
                          loc = matrix(knots),   # location parameter
                          scale = rep(2,length(knots)),          # scale parameter
                          type = "Gaussian")
basis <- TensorP(G_spatial,G_temporal)
p2_3 <- show_basis(basis@Basis2)

start_time <- Sys.time()
M_3K <- FRK(f = COUNT ~ 1 + (1 | reefid) + max_cyc + max_cyc_lag1 + max_cyc_lag2 + max_dhw + max_dhw_lag1 + max_dhw_lag2, 
            data = list(obj_frk$STObj), 
            BAUs = obj_frk$ST_BAUs, 
            basis = basis, 
            response = "binomial", 
            link = "logit", 
            K_type = "precision", 
            method = "TMB", 
            est_error = FALSE)
end_time <- Sys.time()
time_3K <- end_time - start_time

pred_3K <- predict(M_3K, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_3K <- as.data.frame(pred_3K$MC$mu_samples) %>% 
  mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = obj_frk$ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  mutate(knots = "Model 3K") %>%
tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc, knots),
                    names_to = "draw", 
                    values_to = "pred"
)

saveRDS(list(M = M_3K,
             post_dist_df = post_dist_df_3K,
             time = time_3K,
             basis = basis),
        file = "../data/Basis_function_3K.RData")

########## 3. Select 1 knot 
knots <- sample(orign_knots$loc1,1)

G_temporal <- local_basis(manifold=real_line(),      # functions on real line
                          loc = matrix(knots),   # location parameter
                          scale = rep(2,length(knots)),          # scale parameter
                          type = "Gaussian")
basis <- TensorP(G_spatial,G_temporal)

p2_4 <- show_basis(basis@Basis2)

start_time <- Sys.time()
M_1K <- FRK(f = COUNT ~ 1 + (1 | reefid) + max_cyc + max_cyc_lag1 + max_cyc_lag2 + max_dhw + max_dhw_lag1 + max_dhw_lag2, 
            data = list(obj_frk$STObj), 
            BAUs = obj_frk$ST_BAUs, 
            basis = basis, 
            response = "binomial", 
            link = "logit", 
            K_type = "precision", 
            method = "TMB", 
            est_error = FALSE)
end_time <- Sys.time()
time_1K <- end_time - start_time

pred_1K <- predict(M_1K, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_1K <- as.data.frame(pred_1K$MC$mu_samples) %>% 
  mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = obj_frk$ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  mutate(knots = "Model 1K") %>%
  tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc, knots),
                    names_to = "draw", 
                    values_to = "pred"
)

saveRDS(list(M = M_1K,
             post_dist_df = post_dist_df_1K,
             time = time_1K,
             basis = basis),
        file = "../data/Basis_function_1K.RData")