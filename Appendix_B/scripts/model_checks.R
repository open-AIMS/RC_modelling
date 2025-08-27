#' @title Extract model predictive performance measures
#' @author Julie Vercelloni

# Extract computing time for each model 
model_list_run <- list.files(paste0(title_of_run,"/model_outputs/full_run"), recursive = TRUE, pattern = "Rdata") %>% 
  data.frame() 

computing_time_list <- list()

for (k in 1:nrow(model_list_run)){
load(paste0(title_of_run,"/model_outputs/full_run/", model_list_run[k,]))
computing_time_list[[k]] <- cbind(round(as.numeric(mod$computing_time, units = "mins"),3), str_remove(model_list_run[k,], ".Rdata"))
}

computing_time <- data.frame(do.call(rbind,computing_time_list))
colnames(computing_time) <- c("time", "Model")

##################################
############ Leave_out tier analysis using true count  
################################## SPATIAL 

#Import data
dat <- read.csv(paste0(title_of_run,"/data/reef_data_NAs_with_cov_",surveys,".csv")) %>%
  mutate(tier_cat = case_when(is.na(COUNT) ~ "tier_test", TRUE ~ "tier_train")) 

# Group_by tier and year 
dat_grp <- dat %>%
  filter(!is.na(LONGITUDE)) %>% # needed because no values of TRUE_COUNT in every tier 
  #filter(fGROUP == "HCC") %>% 
  group_by(tier, fYEAR) %>%
  mutate(TRUE_COUNT_tier = mean(TRUE_COUNT),
          TOTAL_tier = mean(TOTAL)) %>%
  filter(row_number() == 1) %>%
  dplyr::select(tier, fYEAR, TRUE_COUNT_tier, TOTAL_tier, tier_cat) %>%
  mutate(COVER_tier = TRUE_COUNT_tier / TOTAL_tier) %>%
  data.frame() 

dat_grp$fYEAR <- as.character(dat_grp$fYEAR) 
dat_grp$tier <- as.character(dat_grp$tier) 

# Read model_outputs from model fit 

model_list_pred <- list.files(paste0(title_of_run,"/model_outputs/predictions"), recursive = TRUE, pattern = "Rdata") %>% 
  data.frame() 

# List to save outputs 
pred_dat_true.list <- list()
list.indicators_tier <- list()

for (j in 1:nrow(model_list_pred)){ 
  
load(paste0(title_of_run,"/model_outputs/predictions/", model_list_pred[j,]))

pred_dat <- pred_sum_sf %>%
  st_drop_geometry() %>%
  data.frame() 

pred_dat$fYEAR <- as.character(pred_dat$fYEAR) 
pred_dat$tier <- as.character(pred_dat$tier) 

 # Join with true values 

pred_dat_true.list[[j]] <- pred_dat %>%
 inner_join(dat_grp) %>%
 mutate(Model = str_remove(model_list_pred[j,], ".Rdata")) %>%
 mutate(Diff = COVER_tier - pred)

 # Measures

test_pred_sum <- pred_dat_true.list[[j]] %>%
 mutate(true = TRUE_COUNT_tier / TOTAL_tier) %>%
 filter(tier_cat == "tier_test")

list.indicators_tier[[j]] <- c(round(coverage95(test_pred_sum$true, test_pred_sum$.lower, test_pred_sum$.upper),2),
                              round(IS95(test_pred_sum$true, test_pred_sum$.lower, test_pred_sum$.upper),2),
                              round(RMSPE(test_pred_sum$true, test_pred_sum$pred),2),
                              round(crps(test_pred_sum$true, data.frame(test_pred_sum$pred, test_pred_sum$Unc))$CRPS,2), 
                              str_remove(model_list_pred[j,], ".Rdata"), title_of_run)


}

pred_dat_true <- data.frame(do.call(rbind,pred_dat_true.list)) %>%
                 mutate(Title_of_run = title_of_run)

# Vizualisation 
pal_model <- lacroix_palette("Pamplemousse", n = length(unique(pred_dat_true$Model)), type = "discrete")

p_check <- ggplot(pred_dat_true %>% filter(tier_cat == "tier_test"), aes(x = COVER_tier, y = pred)) + 
 geom_point(alpha=0.4) + 
 facet_wrap(~Model, nrow = 1) +
 geom_abline(linetype = "dashed") +
 labs(title = unique(pred_dat_true$Title_of_run)) +
 xlab("True values") +
 ylab("Predicted values") +
 theme_pubr() +
 scale_color_manual(values = pal_model)

ggsave(p_check, filename = paste0(title_of_run,"/report/extra/check_true.png"),
       width=8, height=6)


# Get performances measures 
indicators_table_tier <- data.frame(do.call(rbind,list.indicators_tier)) 
colnames(indicators_table_tier) <- c("Cvg", "IS", "RMSPE", "CRPS",
                                     "Model", "Title_of_run")

indicators_table_tier  <- indicators_table_tier %>%  left_join(computing_time) %>%
                          mutate(across(c(Cvg, IS, RMSPE, CRPS, time), as.numeric))

write.csv(indicators_table_tier, file = paste0(title_of_run,"/model_outputs/leave_out/table_performances_tier_true.csv"), row.names = F)

##################################
############ Leave_out tier analysis using true count  
################################## TEMPORAL TRENDS 

# Coral cover trajectories of tier with data 
dat_plot <- dat %>%
filter(fGROUP == "HCC") %>%
filter(!is.na(COUNT)) 

dat_plot$fYEAR <- as.character(dat_plot$fYEAR)

# Pick 4 tier randomly 
tier.observed <- dat_plot  %>%
  dplyr::select(tier) %>%
  distinct() %>%
  sample_n(size = 4) %>%
  pull(tier)

hexpred_unique <- st_read(paste0(title_of_run,"/data/hexpred_cov.shp")) %>%
group_by(tier) %>%
filter(row_number() == 1) %>%
dplyr::select(tier) 

hexpred_unique$tier <- as.character(hexpred_unique$tier) 

pred_dat_true_sf <- pred_dat_true %>%
left_join(hexpred_unique, by = "tier") %>%
st_as_sf() 

dat_plot <- dat_plot %>%
  filter(tier %in% tier.observed) %>%
  filter(fDEPTH == 10)

# average across tiers 
dat_plot_tier <- dat_plot %>% 
  group_by(tier, fYEAR) %>%
  summarize(COUNT_t = sum(COUNT, na.rm = TRUE), TOTAL_t= sum(TOTAL, na.rm = TRUE)) %>%
  mutate(COVER_t = (COUNT_t / TOTAL_t) *100)

dat_pred <- pred_dat_true_sf %>%
  filter(tier %in% tier.observed) 

p_data <- ggplot() + 
  geom_ribbon(data = dat_pred, aes(x=fYEAR,ymin=.lower*100, ymax=.upper*100, group=1),alpha=.2, fill ="#0072B2FF")  +
  geom_point(data = dat_plot_tier, 
                   aes(x = fYEAR, y = COVER_t), fill = "grey66", col = "black", size = 2.5, alpha = .6, shape = 20) + 
 geom_line(data = dat_pred, aes(x=fYEAR, y=pred*100, group=1),linewidth=.4) + 
  facet_wrap(~tier, ncol=2) +
  ylab("Coral cover (%)") + xlab("Year") +  theme_pubr() +
  theme(axis.text.x = element_text(size=11, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=11),axis.title.y=element_text(size=13),
        axis.title.x=element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white')) +
  scale_x_discrete(breaks= nth_element(unique(dat_pred$fYEAR),1,4))

#ggsave(p_traj, filename = paste0(title_of_run,"/report/extra/pred_traj_tier_data", unique(dat_pred[[i]]$Model),".png"),
 #      width=6, height=8)


# Coral cover trajectories of tier without data 
# Pick 4 tiers randomly 
dat_plot_no <- dat %>%
filter(fGROUP == "HCC") %>%
filter(is.na(COUNT)) 

dat_plot_no$fYEAR <- as.character(dat_plot_no$fYEAR)

tier.without <- dat_plot_no %>%
  dplyr::select(tier) %>%
  distinct() %>%
  sample_n(size = 4) %>%
  pull(tier)

dat_plot_no <- dat_plot_no %>%
  filter(tier %in% tier.without) 

# average across tiers 
dat_plot_no_tier <- dat_plot_no %>% 
  group_by(tier, fYEAR) %>%
  summarize(COUNT_t = sum(TRUE_COUNT, na.rm = TRUE), TOTAL_t= sum(TOTAL, na.rm = TRUE)) %>%
  mutate(COVER_t = (COUNT_t / TOTAL_t) *100)

dat_pred_without <- pred_dat_true_sf %>%
  filter(tier %in% tier.without) 

p_nodata <- ggplot() + 
  geom_ribbon(data = dat_pred_without, aes(x=fYEAR,ymin=.lower*100, ymax=.upper*100, group=1),alpha=.2, fill ="#D68A8A")  +
  geom_point(data = dat_plot_no_tier, 
                   aes(x = fYEAR, y = COVER_t), fill = "red", col = "black", size = 2.5, alpha = .6, shape = 20) + 
 geom_line(data = dat_pred_without, aes(x=fYEAR, y=pred*100, group=1),linewidth=.4) + 
  facet_wrap(~tier, ncol=2) +
  ylab("Coral cover (%)") + xlab("Year") +  theme_pubr() +
  theme(axis.text.x = element_text(size=11, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=11), axis.title.y = element_text(size=13),
        axis.title.x = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white')) +
  scale_x_discrete(breaks= nth_element(unique(dat_pred_without$fYEAR),1,4))

#ggsave(p_traj, filename = paste0(title_of_run,"/report/extra/pred_traj_tier_no_data", unique(dat_pred_without[[i]]$Model),".png"),
#       width=6, height=8)

p_traj <- p_data / p_nodata +  
  plot_layout(guides = "collect", axis_titles = "collect") +  
  plot_annotation(tag_levels = "a", tag_suffix = ')') &
  theme(plot.tag.position  = c(.035, 1.020)) 

p_traj 

ggsave(p_traj, filename = paste0(title_of_run,"/report/extra/trends_", unique(pred_dat_true_sf$Model),".png"),
       width=6, height=6)

