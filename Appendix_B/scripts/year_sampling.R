#' @title Select sampling years and replot data viz
#' @author Julie Vercelloni

# Import full data
dat_full <- read.csv(paste0(title_of_run, "/data/reef_data_NAs_with_cov_", surveys, ".csv"))

# Focus on HCC group
dat_hcc <- dat_full %>%
  filter(fGROUP == "HCC") %>%
  arrange(REEF_NAME, SITE_NO, TRANSECT_NO, fYEAR)

sample_years_with_condition <- function(years, min_n = 5, max_n = 30, max_gap = 5, max_tries = 100) {
  years <- sort(unique(years))
  n_years <- length(years)
  
  if (n_years < min_n) return(NULL)
  
  tries <- 0
  while (tries < max_tries) {
    n_pick <- sample(min_n:min(max_n, n_years), 1)
    candidate <- sort(sample(years, n_pick))
    
    if (all(diff(candidate) <= max_gap)) {
      return(candidate)
    }
    
    tries <- tries + 1
  }
  
  return(NULL)  # If no valid combo found after max_tries
}

# Sample valid years per reef
reef_years <- dat_hcc %>%
  group_by(REEF_NAME) %>%
  summarise(selected_years = list(sample_years_with_condition(fYEAR)), .groups = "drop") %>%
  unnest(cols = c(selected_years)) %>%
  mutate(reef_year = paste(REEF_NAME, selected_years, sep=""))

dat_hcc_masked <- dat_hcc %>%
  mutate(reef_year = paste(REEF_NAME, fYEAR, sep="")) %>%
  mutate(COUNT = case_when(
    reef_year %in% reef_years$reef_year ~ COUNT,
    TRUE ~ NA_real_
  )) %>%
  dplyr::select(!reef_year)

head(dat_hcc_masked)

# Combine with non-HCC rows (unchanged)
dat_other <- dat_full %>%
  filter(fGROUP != "HCC")

# Final data
dat_final <- bind_rows(dat_hcc_masked, dat_other) %>%
  arrange(REEF_NAME, SITE_NO, TRANSECT_NO, fYEAR) %>%
  mutate(COVER = COUNT / TOTAL)

write_csv(dat_final, file=paste0(title_of_run,"/data/reef_data_NAs_with_cov_", surveys, ".csv"))

## Data viz 

p_vis_data_fixed <- ggplot(dat_final  %>% filter(!is.na(COUNT)) %>% filter(fGROUP == "HCC")) + 
  geom_line(aes(x = fYEAR, y = (COUNT/TOTAL)*100, group = interaction(as.factor(TRANSECT_NO), as.factor(SITE_NO), REEF_NAME),
  col = as.factor(SITE_NO)), 
   show.legend = FALSE) + 
  facet_wrap(~REEF_NAME, ncol=4) + theme_bw() +
  labs(x = "Year", y = "Coral cover") +
  ylab("Coral cover") + xlab("Year")+theme_bw()+
  theme(axis.text.x = element_text(size=10, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=10),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10, margin = margin())) + 
  ggtitle("")

ggsave(filename = paste0(title_of_run,"/report/extra/trend_data_",surveys,".png"),
       plot = p_vis_data_fixed, width=13, height=12)  


reef_data.synthetic_fixed_ready_site <- dat_final %>%
 filter(!is.na(COUNT)) %>% filter(fGROUP == "HCC") %>% 
 group_by(REEF_NAME, SITE_NO, fYEAR, fDEPTH, fGROUP) %>%
 summarize(COUNT_sum = sum(COUNT),
           TOTAL_sum = sum(TOTAL)) %>%
 mutate(COVER_site = COUNT_sum / TOTAL_sum) %>%
  mutate(fYEAR = as.numeric(as.character(fYEAR))) 

if (length(unique(reef_data.synthetic_fixed_ready_site$fYEAR)) >10){
p_heat <- ggplot(reef_data.synthetic_fixed_ready_site) +
  geom_tile(aes(x = fYEAR, y = as.factor(SITE_NO),
                fill = COVER_site * 100)) +
  scale_fill_viridis(
    name = "Coral cover (%)", 
    option = "plasma", 
    begin = 0, 
    end = ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE)), 
    limits = c(0, ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE) * 100)), 
    na.value = "grey90"
  ) +
  scale_x_continuous(
    breaks = seq(min(reef_data.synthetic_fixed_ready_site$fYEAR, na.rm = TRUE),
                 max(reef_data.synthetic_fixed_ready_site$fYEAR, na.rm = TRUE),
                 by = 4)
  ) +
  facet_wrap(~REEF_NAME, ncol = 4) +
  labs(x = "Year", y = "Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
         axis.text.y = element_text(size=10),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10, margin = margin()),
        legend.position = "bottom")
}else{
p_heat <- ggplot(reef_data.synthetic_fixed_ready_site) +
  geom_tile(aes(x = fYEAR, y = as.factor(SITE_NO),
                fill = COVER_site * 100)) +
  scale_fill_viridis(
    name = "Coral cover (%)", 
    option = "plasma", 
    begin = 0, 
    end = ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE)), 
    limits = c(0, ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE) * 100)), 
    na.value = "grey90"
  ) +
  facet_wrap(~REEF_NAME, ncol = 4) +
  labs(x = "Year", y = "Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
         axis.text.y = element_text(size=10),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10, margin = margin()),
        legend.position = "bottom")
}

ggsave(filename = paste0(title_of_run,"/report/extra/tile_data_",surveys,".png"),
       plot = p_heat, width=13, height=12)  
