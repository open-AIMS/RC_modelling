#' @title Extract model predictions 
#' @author Julie Vercelloni

## Import data and average by tiers 
dat <- read.csv(paste0(title_of_run,"/data/reef_data_NAs_with_cov_",surveys,".csv")) %>%
  filter(!is.na(COUNT)) %>%
  filter(fGROUP == "HCC") %>%
  arrange(REEF_NAME, SITE_NO, TRANSECT_NO, fYEAR) %>%
  group_by(fYEAR, tier) %>%
  summarize(COUNT_t = sum(COUNT), TOTAL_t = sum(TOTAL)) %>%
  dplyr::select(tier, fYEAR, COUNT_t, TOTAL_t)

# Import predictive layer

hexpred <- st_read(paste0(title_of_run,"/data/hexpred_cov.shp"))

hexpred_unique <- hexpred %>%
                group_by(tier) %>%
                filter(row_number() == 1) %>%
                dplyr::select(tier) 
      
# Read model_outputs 

model_list <- list.files(paste0(title_of_run,"/model_outputs/full_run"), recursive = TRUE, pattern = "Rdata")

##################################
########################################### PREDICTIONS 
##################################

model_frk <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "FRK"))

for (l in 1:nrow(model_frk)){
  if(nrow(model_frk)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_frk[l,]))

# Making spatial prediction and saving data table 
obj_frk <- mod$obj_frk
pred_sum_sf <- predictions_FRK(model.out = mod$fit) 
head(pred_sum_sf)

# Making plots and saving them 

# Mean
p_pred <- plot_predictions(pred_sum_sf)

ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/pred",str_remove(model_frk[l,], ".Rdata"),".png"),
       width=6, height=6)

# Uncertainty 
p_unc <- plot_predictions_unc(pred_sum_sf)

ggsave(p_unc, filename = paste0(title_of_run,"/report/extra/pred_unc",str_remove(model_frk[l,], ".Rdata"),".png"),
       width=6, height=6)

# Save predictions
save(pred_sum_sf, file = paste0(title_of_run,"/model_outputs/predictions/",str_remove(model_frk[l,], ".Rdata"),".Rdata"),
          row.names = F)

##################################
########################################### TRENDS AT FINE SPATIAL SCALE 
##################################

### Coral cover trajectories of tier with data 

pred_FRK_data <- pred_sum_sf %>% filter(tier %in% unique(dat$tier))

pred_FRK_data$tier <- as.character(pred_FRK_data$tier)
dat$tier <- as.character(dat$tier)
pred_FRK_data$fYEAR <- as.character(pred_FRK_data$fYEAR)
dat$fYEAR <- as.character(dat$fYEAR)

pred_with_data <- pred_FRK_data %>% data.frame() %>% full_join(dat) 

# Sample 9 tiers randomly for vizulisations 
unique_tier <- sample(unique(dat$tier), 9)

p_pred <- plot_traj(pred_with_data, unique_tier)$p_traj
ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/trends_data_",str_remove(model_frk[l,], ".Rdata"),".png"),
       width=6, height=6)

### Coral cover trajectories of tier without data 
pred_FRK_nodata <- pred_sum_sf %>% filter(!tier %in% unique(dat$tier))

# Sample 9 tiers randomly for vizulisations 
unique_tier <- sample(unique(pred_FRK_nodata$tier), 9)

p_pred <- plot_traj(pred_FRK_nodata, unique_tier)$p_traj_no
ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/trends_nodata_",str_remove(model_frk[l,], ".Rdata"),".png"),
       width=6, height=6)

}


###########################################
########################################### TRENDS AT BROAD SPATIAL SCALE
###########################################

# Import GRMF layer
GRMF <- st_read(paste0(title_of_run,"/data/true_sp_field.shp")) 

GRMF_all <- st_join(hexpred_unique, GRMF) %>%
    data.frame %>%
    filter(!is.na(Year)) %>%
    ungroup() %>%
    group_by(Year) %>% 
    median_hdci(HCC) %>%
    mutate(fYEAR = Year + min(as.numeric(dat$fYEAR)) - 1) %>%
    mutate(true = plogis(HCC)) #%>%
    #st_as_sf() 

##################################
############ FRK model(s)  
##################################

model_frk <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "FRK"))

for (l in 1:nrow(model_frk)){
  if(nrow(model_frk)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_frk[l,]))

# Making spatial prediction and saving data table 
obj_frk <- mod$obj_frk
pred_sum_sf <- predictions_FRK_broad(model.out = mod$fit) %>% data.frame()
head(pred_sum_sf)

# Making plots and saving them 
pred_sum_sf$fYEAR <- as.factor(pred_sum_sf$fYEAR)

p_pred <- plot_traj_broad(pred_sum_sf, GRMF_all)
ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/pred_broad",str_remove(model_frk[l,], ".Rdata"),".png"),
       width=6, height=4)

# Save predictions 
write.csv(pred_sum_sf, file = paste0(title_of_run,"/model_outputs/predictions/",str_remove(model_frk[l,], ".Rdata"),"_broad.csv"),
          row.names = F)

}

