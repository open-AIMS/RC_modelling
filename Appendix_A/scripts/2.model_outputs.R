# Extract model outputs including model goodness-of-fit, residuals and coral cover trajectories at multiple spatial scales
# Save figures 
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

################################
################################
################################ R²

mean_obs <- mod_$data.grp.tier %>% 
  group_by(fYEAR, Tier5) %>%
  summarize(mean_cover = mean(COUNT / TOTAL)*100) %>%
  mutate(fYEAR = as.factor(fYEAR),
         Tier5 = as.factor(Tier5))

# Extract predictions 

pred <- predict(mod_$M, type = c("mean"))
  
# Extracting posterior distributions of predictive locations 
  
post_dist_df <- as.data.frame(pred$MC$mu_samples) %>% 
    mutate(fYEAR = mod_$obj_frk$ST_BAUs@data$fYEAR) %>%
    mutate(Tier5 = mod_$obj_frk$ST_BAUs@data$Tier5) %>%
    mutate(id_loc = row_number()) %>%
    tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc),
                        names_to = "draw", 
                        values_to = "pred")

mean_pred <- post_dist_df %>% group_by(fYEAR,Tier5) %>% 
  ggdist::median_hdci(pred)%>%
  inner_join(mod_$HexPred_reefid2 %>% group_by(Tier5) %>% slice(1) %>% dplyr::select(geometry,Tier5)) %>% 
  st_as_sf(sf_column_name = "geometry") %>%
  mutate(pred = pred*100)

mean_all <- mean_obs %>% inner_join(mean_pred %>% mutate(Tier5 = as.factor(Tier5),
                                                        fYEAR = as.factor(fYEAR)) %>% dplyr::select(fYEAR, Tier5, pred))

# Compute R²
r2 <- cor(mean_all$mean_cover, mean_all$pred)^2
r2_label <- paste0("R² = ", round(r2, 2))

model_fit <- ggpubr::ggscatter(mean_all, x = "mean_cover", y = "pred",
          add = NULL,                   
          shape = 16,                
          alpha = 0.4,                 
          size = 2.5) +                 
  geom_abline(slope = 1, intercept = 0, 
              color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = Inf, y = -Inf, 
           label = r2_label, hjust = 1.1, vjust = -1.5, size = 5) +
  labs(
    x = "Observed Mean Coral Cover (%)",
    y = "Predicted Coral Cover (%)"
  ) 

ggsave(model_fit, filename = "../figures/model_fit.png",
       width=10, height=8)


################################
################################
################################ Model residuals 
res <- make_dharma_res(mod_)

png("../figures/dharma_residuals.png", width = 800, height = 600)
plot(res)
dev.off()

################################
################################
################################ Coral trends

mod_$HexPred_reefid2$Tier5 <- as.factor(mod_$HexPred_reefid2$Tier5)
  
pred_sum_sf <- post_dist_df %>% group_by(fYEAR, Tier5) %>% 
    median_hdci(pred) %>%
    inner_join(mod_$HexPred_reefid2 %>% group_by(Tier5) %>% 
    summarize() %>% dplyr::select(geometry, Tier5)) %>% 
    st_as_sf(sf_column_name = "geometry") %>%
    mutate(Unc = .upper - .lower) %>%
    mutate(tier_fYEAR = paste0(Tier5,fYEAR)) %>%
    dplyr::select(fYEAR, Tier5, pred, .lower, .upper, Unc, tier_fYEAR)

################################ Fine-scale trends
 unique_hex <- mod_$HexPred_reefid2 %>%
     group_by(Tier5) %>%
     filter(row_number()==1) %>%
     dplyr::select(Tier5, geometry)

 hexpred_nodata <- unique_hex %>% 
     filter(!Tier5 %in%  unique(mean_obs$Tier5))

# hexpred_data <- unique_hex %>% 
#     filter(Tier5 %in%  unique(mean_obs$Tier5))

### Locations of data-tiers
# sf_use_s2(FALSE) 
# GBR_shp_cut <- st_crop(GBR_shp, st_bbox(tier4_cut))

# map_data_tier <- ggplot() + 
#   geom_sf(data = GBR_shp_cut, fill="transparent") +
#   geom_sf(data = tier4_cut, fill="transparent") + 
#   geom_sf(data = hexpred_nodata , fill= "gray88", col = "gray70")    +    
#   geom_sf(data = hexpred_data , fill="#72b0d3")    +          
#   theme_pubr() +
#   xlab("") + ylab("")

# ggsave(plot =  map_data_tier, width=6, height=8, file = "figures/viz_data_tier.png")

### Coral cover trajectories of tier5 with data 

pred_FRK_data <- pred_sum_sf %>% filter(Tier5 %in% unique(mean_obs$Tier5))

pred_FRK_data$Tier5 <- as.character(pred_FRK_data$Tier5)
pred_FRK_data$fYEAR <- as.character(pred_FRK_data$fYEAR)

pred_with_data <- pred_FRK_data %>% data.frame() %>% full_join(mean_obs) 

# Remove tiers with less than three replicated years 

pred_with_data_tal <- pred_with_data %>%
 filter(! is.na(mean_cover)) %>%
 group_by(fYEAR, Tier5) %>%
 count() %>%
 group_by(Tier5) %>%
 count() %>%
 filter(n < 4)

pred_data_transect <- pred_with_data %>% 
                   filter(!Tier5 %in% pred_with_data_tal$Tier5) 

p_data_all <- ggplot() + 
   geom_ribbon(data = pred_with_data %>% filter(!Tier5 %in% pred_with_data_tal$Tier5), 
                      aes(x=as.numeric(fYEAR),ymin=.lower*100, ymax=.upper*100, group=1),alpha=.3, fill ="#72b0d3") +
  geom_point(data = pred_data_transect, 
                   aes(x = as.numeric(fYEAR), y = mean_cover), fill = "grey66", col = "black", size = 2.5, alpha = .6, shape = 20) + 
 geom_line(data = pred_with_data %>% filter(!Tier5 %in% pred_with_data_tal$Tier5), 
                    aes(x=as.numeric(fYEAR), y=pred*100, group=1),size=.8) + 
  facet_wrap(~Tier5, ncol=4) +
  ylab("Coral cover (%)") + xlab("Year") +  theme_pubr() +
  theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=10),axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white')) +
    scale_y_continuous(
      breaks = seq(0, 80, by = 40)
    )

ggsave(p_data_all, filename = paste0("../figures/trend_data_all_.png"),
       width=12, height=14)


### Coral cover trajectories of tier5 without data 

# Get 9 cells for common viz 
pred_without_data <- pred_sum_sf %>% filter(Tier5 %in% unique(hexpred_nodata$Tier5))
tier5_nodata_viz <- sample(unique(pred_without_data$Tier5), 25)

nodata_tier5 <- pred_without_data %>%
  filter(Tier5 %in% tier5_nodata_viz)

p_nodata <- ggplot(nodata_tier5) + 
  geom_ribbon(aes(x=as.numeric(fYEAR),ymin=.lower*100, ymax=.upper*100, group=1),alpha=.9, fill= "#D68A8A") +
  geom_line(aes(x=as.numeric(fYEAR), y=pred*100, group=1), linewidth=.6) +
  facet_wrap(~Tier5) + 
  ylab("Coral cover (%)") + xlab("Year") +  theme_pubr() +
  theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=10),axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white')) +
    scale_y_continuous(
      breaks = seq(0, 80, by = 40)
    )
ggsave(p_nodata, filename = paste0("../figures/trend_data_no_data.png"),
       width=12, height=14)

################################ Regional trends

# Weight predictions by reef areas 
reef_areas <- mod_$HexPred_reefid2 %>% 
  group_by(Tier5) %>%
  filter(row_number() == 1) %>%
  mutate(reef_ar = reef_ar / 1000000) %>% # areas in square kilometer
  dplyr::select(Tier5, reef_ar) %>% 
  st_drop_geometry() 

post_dist_df <- left_join(post_dist_df,  reef_areas )

tot_area <- sum(reef_areas$reef_ar) 
 
######## Plot temporal trends at different spatial scales 

## The whole region, sum over all tiers

post_dist_region <- post_dist_df %>% 
  mutate(weighted_pred = pred * reef_ar) %>%
  group_by(fYEAR, draw) %>%
  summarize(cover = sum(weighted_pred, na.rm = TRUE),
            cover_prop = cover / tot_area) 

pred_region <-  post_dist_region %>% group_by(fYEAR) %>% 
  ggdist::median_hdci(cover_prop) %>%
  dplyr::select(fYEAR:.upper) %>%
  data.frame() 

p1 <-  ggplot() +
  geom_ribbon(data = pred_region %>% data.frame(), aes(x = fYEAR, ymin=.lower*100, ymax=.upper*100, group=1), alpha=.2, fill="#00FFC2")+
  geom_line(data = pred_region %>% data.frame(), aes(x=fYEAR, y=cover_prop*100,group=1), col="black", linewidth=1.1)+
  xlab("Year") +ylab("Coral cover (%)") + 
  ylim(0, 40) + 
  labs(subtitle = "All-tiers") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    legend.position = "none",
    plot.subtitle = element_text(size = 16)
  )

## Plot regions keeping data-tier only 

# Filter data tier
data_tier <- unique(mod_$data.grp.tier$Tier5)

# Update sum reef areas 
reef_areas_datatier <- reef_areas %>% 
  filter(Tier5 %in% data_tier)

tot_area_datatier <- sum(reef_areas_datatier$reef_ar) 

post_dist_region_data <- post_dist_df %>% 
  filter(Tier5 %in% data_tier) %>%
  mutate(weighted_pred = pred * reef_ar) %>%
  group_by(fYEAR, draw) %>%
  summarize(cover = sum(weighted_pred, na.rm = TRUE),
            cover_prop = cover / tot_area_datatier ) 

pred_region_data <-  post_dist_region_data %>% group_by(fYEAR) %>% 
  ggdist::median_hdci(cover_prop) %>%
  dplyr::select(fYEAR:.upper) %>%
  data.frame() 

p2 <- ggplot() +
  geom_ribbon(data = pred_region_data %>% data.frame(), aes(x = fYEAR, ymin=.lower*100, ymax=.upper*100, group=1), alpha=.2, fill="#72b0d3")+
  geom_line(data = pred_region_data %>% data.frame(), aes(x=fYEAR, y=cover_prop*100,group=1), col="black", linewidth=1.1)+
  xlab("Year") +ylab("Coral cover (%)") + 
  ylim(0, 40) + 
  labs(subtitle = "Data-tiers (7.1%)") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    legend.position = "none",
    plot.subtitle = element_text(size = 16)
  )


pcut <- p1 + p2 +  
  plot_layout(guides = "collect", axis_titles = "collect") +  
  plot_annotation(tag_levels = "a", tag_suffix = ')')&
  theme(plot.tag = element_text(size = 18))
pcut 

ggsave("../figures/WetTropics_region_trends_cut.png", plot = pcut, width = 10, height = 6, dpi = 300)
