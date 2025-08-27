# Extract model outputs including covariates effects and conditional plots
# @author Julie Vercelloni
# August 2025 

rm(list=ls())
setwd(paste0(here::here(), "/Appendix_A/scripts"))

rm(list=ls())
# Load R package 
source("../R/packages.R")
source("../R/functions.R")

# Read model outputs 
mod_ <- readRDS("../data/FRK_fit.RData")

# Conditional plots
mod_M <- mod_$M
mean_cov <- colMeans(mod_M@X)
BAUs <- mod_$obj_frk$ST_BAUs
obj_frk <- mod_$obj_frk
hexpred <- mod_$HexPred_reefid2

################################ Method 1: regional level 

# Extract covariate effects 
coef_ <- cov_extract(
        model.out = mod_$M,
        obj_frk = mod_$obj_frk) %>% 
        filter(Type == "fixed") %>%
        mutate(param = ifelse(param == "Intercept", "(Intercept)", param))

# Plot coef_ 
coef_plot <- coef_ %>%
 filter(!param == "(Intercept)") %>%
 mutate(param = case_when(param == "max_cyc" ~ "Cyclone Exposure",
                          param == "max_cyc_lag1" ~ "Cyclone Exposure (lag 1)",
                          param == "max_cyc_lag2" ~ "Cyclone Exposure (lag 2)",
                          param == "max_dhw" ~ "Heat Stress",
                          param == "max_dhw_lag1" ~ "Heat Stress (lag 1)",
                          param == "max_dhw_lag2" ~ "Heat Stress (lag 2)",
                          TRUE ~ "NA"))

p_coef <- ggplot(data = coef_plot, aes(x = `50%`, y = param)) + 
  geom_point(size = 3, shape = 21, fill = "black", color = "black") + 
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.1) + 
  theme_pubr() + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  xlab("Effect Size") + ylab("") +
  scale_x_continuous(limits = c(-0.055, 0.015), 
                     breaks = seq(-0.055, 0.015, by = 0.02))

ggsave("../figures/coef_.png", plot = p_coef, width = 8, height = 6, dpi = 300)

coef_scaled <- coef_ %>%
  mutate(
    mean_val = mean_cov[param],  
    effect_50 = `50%` * mean_val,
    effect_low = `2.5%` * mean_val,
    effect_high = `97.5%` * mean_val
  ) %>%
  dplyr::select(param, effect_50, effect_low, effect_high)

# DHW
var <- "max_dhw"
range_vals <- seq(0, max(hexpred$max_dhw, na.rm = TRUE)) 
p1 <- plot_cond(range_vals, var, xlabel = "Heat stress (DHW)", ylabel = "Coral Cover (%)")

# DHW (lag1)
var <- "max_dhw_lag1"
range_vals <- seq(0, max(hexpred$max_dhw, na.rm = TRUE)) 
p2 <- plot_cond(range_vals, var, xlabel = "Heat stress (DHW, lag 1)", ylabel = "Coral Cover (%)")

# Cyclone exposure
var <- "max_cyc"
range_vals <- seq(0, max(hexpred$max_cyc, na.rm = TRUE)) 
p3 <- plot_cond(range_vals, var, xlabel = "Cyclone exposure (in hours)", ylabel = "Coral Cover (%)")

# Cyclone exposure (lag1)
var <- "max_cyc_lag1"
range_vals <- seq(0, max(hexpred$max_cyc, na.rm = TRUE)) 
p4 <- plot_cond(range_vals, var, xlabel = "Cyclone exposure (in hours, lag 1)", ylabel = "Coral Cover (%)")

p_lay <- "
AB
CD
"
p_cond <- (p3 + p4 + p1 + p2) +  
  plot_layout(design = p_lay, guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ')') &
  theme(plot.tag.position  = c(.2, .96)) 

for (i in c(2, 4)) p_cond[[i]] <- p_cond[[i]] + theme(plot.tag.position  = c(.15, .96))
p_cond

ggsave("../figures/condtional_effects.png", plot = p_cond, width = 8, height = 6, dpi = 300)

################################ Method 2: BAU level 

################## Same data 
sf::sf_use_s2(FALSE)

olddata <- mod_M@BAUs@data 
p_old <- predict_newdata(mod_M, olddata, hexpred, obj_frk, title_n = "Old data")
p_pred <- predict_year(mod_M, olddata, hexpred, obj_frk, title_n = "")$p
p_unc <- predict_year(mod_M, olddata, hexpred, obj_frk, title_n = "")$p_unc

ggsave("../figures/pred_new.png", plot = p_pred, width = 10, height = 14, dpi = 300)
ggsave("../figures/pred_unc.png", plot = p_unc, width = 10, height = 14, dpi = 300)

################## Cyclone exposure 
## Low Cyclone
N <- nrow(olddata)

new_BAU_data <- olddata %>%
  mutate(max_cyc      = rnorm(N, 2, 0.001),
         max_cyc_lag1 = rnorm(N, 0, 0.001),
         max_cyc_lag2 = rnorm(N, 0, 0.001),
         max_dhw      = rnorm(N, 0, 0.001),
         max_dhw_lag1 = rnorm(N, 0, 0.001),
         max_dhw_lag2 = rnorm(N, 0, 0.001)) 

p_cyc_low <- predict_newdata(mod_M, new_BAU_data, hexpred, obj_frk, title_n = "a) Low Cyclone Exposure")

new_BAU_data <- olddata %>%
  mutate(max_cyc      = rnorm(N, 44, 0.001),
         max_cyc_lag1 = rnorm(N, 0, 0.001),
         max_cyc_lag2 = rnorm(N, 0, 0.001),
         max_dhw      = rnorm(N, 0, 0.001),
         max_dhw_lag1 = rnorm(N, 0, 0.001),
         max_dhw_lag2 = rnorm(N, 0, 0.001)) 

p_cyc_high <- predict_newdata(mod_M, new_BAU_data, hexpred, obj_frk, title_n = "b) High Cyclone Exposure")

################## Heat stress
## Low 
new_BAU_data <- olddata %>%
  mutate(max_cyc      = rnorm(N, 0, 0.001),
         max_cyc_lag1 = rnorm(N, 0, 0.001),
         max_cyc_lag2 = rnorm(N, 0, 0.001),
         max_dhw      = rnorm(N, 2, 0.001),
         max_dhw_lag1 = rnorm(N, 0, 0.001),
         max_dhw_lag2 = rnorm(N, 0, 0.001)) 

p_dhw_low <- predict_newdata(mod_M, new_BAU_data, hexpred, obj_frk, title_n = "c) Low Heat Stress Event")

## High 
new_BAU_data <- olddata %>%
  mutate(max_cyc      = rnorm(N, 0, 0.001),
         max_cyc_lag1 = rnorm(N, 0, 0.001),
         max_cyc_lag2 = rnorm(N, 0, 0.001),
         max_dhw      = rnorm(N, 12, 0.001),
         max_dhw_lag1 = rnorm(N, 0, 0.001),
         max_dhw_lag2 = rnorm(N, 0, 0.001)) 

p_dhw_high <- predict_newdata(mod_M, new_BAU_data, hexpred, obj_frk, title_n = "d) High Heat Stress Event")

## All together 
p_cond_tier5_cyc <-  (p_cyc_low + p_cyc_high) +
  plot_layout(nrow = 1, guides = "collect") & 
  theme(legend.position = 'none') 
p_cond_tier5_cyc

ggsave("../figures/condtional_effects_tier5_cyc.png", plot = p_cond_tier5_cyc, width = 8, height = 6, dpi = 300)

p_cond_tier5_dhw <- (p_dhw_low + p_dhw_high) +
  plot_layout(nrow = 1, guides = "collect") & 
  theme(legend.position = 'none') 
 p_cond_tier5_dhw

ggsave("../figures/conditional_effects_tier5_dhw.png", plot = p_cond_tier5_dhw, width = 8, height = 6, dpi = 300)
ggsave("../figures/conditional_effects_tier5_data.png", plot = p_old, width = 8, height = 6, dpi = 300)
