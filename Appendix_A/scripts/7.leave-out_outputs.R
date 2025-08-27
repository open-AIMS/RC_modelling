# Extract outputs of the leave-out data analysis and visualize with a spider plot 
# @author Julie Vercelloni
# August 2025 

#setwd("c:/Users/jvercell/OneDrive - Australian Institute of Marine Science/AIMS/01_Research projects/ReefCloud/SP_models/FRK_dev/for_andrew/scripts")
setwd(paste0(here(), "/scripts"))

rm(list=ls())

# Load R package 
source("../R/packages.R")
source("../R/functions.R")

mod_out <-
  list.files("../data", pattern = "\\.RData$", full.names = TRUE) %>%
  tibble(file = .) %>% 
  filter(str_detect(file, "leave_out")) %>% 
  pull(file) %>%  
  map(readRDS)  

# Get performances measures 
indicators_table_tier <- mod_out %>% 
  map(~ tibble(indicator = .x$list.indicators)) %>% 
  bind_cols()

indicators_table_tier <- as_tibble(t(indicators_table_tier))
colnames(indicators_table_tier) <- c("Test","Cvg", "IS", "RMSPE", "CRPS", "AIC")

indicators_table_tier  <- indicators_table_tier %>%  
                          mutate(across(c(Cvg, IS, RMSPE, CRPS), as.numeric)) %>%
                          data.frame() 

# Viz model performances
scale <- 4.5 # scale is an adjustment for a 4k screen

# Define the earthy and contrasted color palette
colors <- c(
  "#67AFCB", 
  "#A23B72", 
  "#F3C13A", 
  "#3B7D4E",
  "#E78C8C"
)

tab <- indicators_table_tier %>% dplyr::select(Test, RMSPE, `Cvg`, `CRPS`, `IS`)

p_perf <- ggspider(tab, 
         axis_name_offset = 0.15,
         background_color = "gray98", 
         fill_opacity = 0.15, 
         polygon = FALSE) +
  labs(col = "") +
  scale_color_manual(values = colors) +  
  theme(
    legend.position = "top", 
    plot.title = element_text(size = 5 * scale, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 4 * scale, hjust = 0.5),
    plot.caption = element_text(size = 3 * scale),
    legend.text = element_text(size = 2 * scale),
    legend.title = element_text(size = 2.5 * scale)
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))


ggsave(plot = p_perf , width=5.5, height=6, file = "../figures/viz_pperf.png")

