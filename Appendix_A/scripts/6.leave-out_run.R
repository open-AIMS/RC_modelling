# Run the leave-out data analysis
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

X <- mod_$data.grp.tier %>%
  mutate(id = row_number()) %>%
  mutate(Transect = paste(REEF,TRANSECT_NO, sep="_"),
         Reef = str_replace(sub(" Site.*", "", REEF), " ", "_"))

HexPred_reefid2 <- mod_$HexPred_reefid2 

# Prepare the lists 

test.validation <- c("1. rm(20% obs)", "2. rm(20% reef)", "3. rm(20% site)", "4. rm(20% transect)", "5. rm(5YRS)")

# Make training and testing datasets for each validation test 

for (i in 1:length(test.validation)){
  
  test.name <- test.validation[i]
  if (test.name == test.validation[1]) {
    train <- X %>% sample_frac(.80)
    test  <- anti_join(X, train, by = 'id')
  } else if (test.name == test.validation[2]){
    reef_pick <-  sample(unique(X$Reef), ceiling(length(unique(X$Reef))*.8), replace = F) 
    train <-  X %>% filter(Reef %in% reef_pick)
    test <- anti_join(X, train, by = 'id')
  } else if (test.name == test.validation[3]){
    site_pick <-  sample(unique(X$REEF), ceiling(length(unique(X$REEF))*.8)) 
    train <-  X %>% filter(REEF %in% site_pick)
    test <-  anti_join(X, train, by = 'id')
  } else if(test.name == test.validation[4]){
    transect_pick <-  sample(unique(X$Transect), ceiling(length(unique(X$Transect))*.8)) 
    train <-  X %>% filter(Transect %in% transect_pick)
    test <-  anti_join(X, train, by = 'id')  
  }else{
    tal_transect <- X %>% group_by(Transect) %>% tally() %>% mutate(year_keep = n - 5) %>% filter(!year_keep < 1) %>% arrange(Transect)
    
    train <- X %>% inner_join(tal_transect) %>%   arrange(Transect) %>%
      tidyr::nest(.by = Transect)  %>% 
      mutate(n = tal_transect$year_keep) %>%
      mutate(samp = map2(data, n, sample_n)) %>%
      dplyr::select(-c(data,n)) %>%
      tidyr::unnest(samp) %>%
      data.frame()
    
    test <-  anti_join(X, train, by = 'id')   
  }

# Model run 
obj_frk <- frk_prep(train, HexPred_reefid2) 

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

# Predict on test data
pred <- predict(M, type = c("mean"))
  
# Extracting posterior distributions of predictive locations 
  
post_dist_df <- as.data.frame(pred$MC$mu_samples) %>% 
    mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
    mutate(Tier5 = obj_frk$ST_BAUs@data$Tier5) %>%
    mutate(id_loc = row_number()) %>%
    tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc),
                        names_to = "draw", 
                        values_to = "pred"
    )

# Mean test dataset 

test_mean <- test %>% 
 group_by(Tier5, fYEAR) %>%
 summarize(mean_cov = mean(COUNT / TOTAL))

test_mean$Tier5 <- as.character(test_mean$Tier5)
test_mean$fYEAR <- as.numeric(as.character(test_mean$fYEAR))
  
#test_pred <- inner_join(test_mean, post_dist_df, by = c("Tier5","fYEAR"))
  
test_pred_sum <- post_dist_df %>% group_by(fYEAR,Tier5) %>% 
    ggdist::median_hdci(pred)%>%
    inner_join(test_mean) %>%
    mutate(p_se = .upper - .lower) %>%
    mutate(Test = test.name)
  
# Measures
  
list.indicators <- c(test.name, round(coverage95(test_pred_sum$mean_cov, test_pred_sum$.lower, test_pred_sum$.upper),2),
                           round(IS95(test_pred_sum$mean_cov, test_pred_sum$.lower, test_pred_sum$.upper),2),
                           round(RMSPE(test_pred_sum$mean_cov, test_pred_sum$pred),2),
                           round(crps(test_pred_sum$mean_cov, data.frame(test_pred_sum$pred, test_pred_sum$p_se))$CRPS,2))

saveRDS(list(train_data = train,
             test_data = test,
             M = M,
             post_dist_df = post_dist_df,
             test_pred_sum = test_pred_sum, 
             list.indicators = list.indicators),
        file = paste0("../data/",
                      "leave_out",
                      "_", test.name, ".RData"))
}
