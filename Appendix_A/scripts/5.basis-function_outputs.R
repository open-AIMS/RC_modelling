# Extract outputs of the exploration of different basis functions
# @author Julie Vercelloni
# August 2025

#setwd("c:/Users/jvercell/OneDrive - Australian Institute of Marine Science/AIMS/01_Research projects/ReefCloud/SP_models/FRK_dev/for_andrew/scripts")
setwd(paste0(here(), "/scripts"))

rm(list=ls())

# Load R package 
source("../R/packages.R")
source("../R/functions.R")

mod_ <- readRDS("../data/FRK_fit.RData")
pred <- predict(mod_$M, type = c("mean"))
hexpred = mod_$HexPred_reefid2

# Extracting posterior distributions of predictive locations 
  
post_dist_df_10K <- as.data.frame(pred$MC$mu_samples) %>% 
    mutate(fYEAR = mod_$obj_frk$ST_BAUs@data$fYEAR) %>%
    mutate(Tier5 = mod_$obj_frk$ST_BAUs@data$Tier5) %>%
    mutate(id_loc = row_number()) %>%
    mutate(knots = "Model full") %>%
    tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc, knots),
                        names_to = "draw", 
                        values_to = "pred")

# Read model outputs 
mod_bf <-
  list.files("../data", pattern = "\\.RData$", full.names = TRUE) %>%
  tibble(file = .) %>% 
  filter(str_detect(file, "Basis_function")) %>% 
  pull(file) %>%  
  map(readRDS)  

############# Extract all together 
post_dist_df <- bind_rows(post_dist_df_10K, mod_bf[[1]]$post_dist_df, mod_bf[[2]]$post_dist_df, mod_bf[[3]]$post_dist_df)

# Broad scale level 
reef_areas <- hexpred %>% 
  group_by(Tier5) %>%
  filter(row_number() == 1) %>%
  mutate(reef_ar = reef_ar / 1000000) %>% # areas in square kilometer
  dplyr::select(Tier5, reef_ar) %>% 
  st_drop_geometry() 

post_dist_df <- left_join(post_dist_df,  reef_areas )

tot_area <- sum(reef_areas$reef_ar) 
 
post_dist_df <- post_dist_df %>% 
  mutate(weighted_pred = pred * reef_ar) %>%
  group_by(fYEAR, knots, draw) %>%
  summarize(cover = sum(weighted_pred, na.rm = TRUE),
            cover_prop = cover / tot_area) 

pred_tier4 <-  post_dist_df %>% group_by(fYEAR, knots) %>% 
  ggdist::median_hdci(cover_prop) %>%
  dplyr::select(fYEAR:.upper) %>%
  data.frame() 

pred_tier4$knots <- as.factor(pred_tier4$knots)

pred_tier4 <- pred_tier4 %>%
 mutate(knots_plot = case_when(knots == "Model full" ~ "Default",
                               knots == "Model 5K" ~ "5",
                               knots == "Model 3K" ~ "3",
                               knots == "Model 1K" ~ "1",
                               TRUE ~ "NA") )

pred_tier4$knots_plot <- factor(pred_tier4$knots_plot, levels = c("Default", "5",
                                                                "3", "1"))
plot_traj_broad <- ggplot() +
geom_ribbon(data = pred_tier4 %>% data.frame(), aes(x = fYEAR, ymin=.lower*100, ymax=.upper*100, group=1), alpha=.2, fill="#00FFC2")+
  geom_line(data = pred_tier4 %>% data.frame(), aes(x=fYEAR, y=cover_prop*100,group=1), col="black", linewidth=1.1)+
  xlab("Year") +ylab("Coral cover")+ theme_pubr() +
  ylim(0, 40) + 
  facet_wrap(~knots_plot) +
  theme(axis.text.x = element_text(size=13),legend.position = "none",
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15)) 

ggsave(plot_traj_broad, width=8, height=6, file = "../figures/basisfunction_tier4.png")

# Viz basis function locations 

p2_1 <- show_basis(mod_$obj_frk$basis@Basis2) + xlab("Year") + ylab("Weight") +
  theme(axis.text.x = element_text(size=13),legend.position = "none",
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15)) 
p2_2 <- show_basis(mod_bf[[3]]$basis@Basis2) + xlab("Year") + ylab("Weight") +
  theme(axis.text.x = element_text(size=13),legend.position = "none",
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15)) 
p2_3 <- show_basis(mod_bf[[2]]$basis@Basis2) + xlab("Year") + ylab("Weight") +
  theme(axis.text.x = element_text(size=13),legend.position = "none",
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15)) 
p2_4 <- show_basis(mod_bf[[1]]$basis@Basis2) + xlab("Year") + ylab("Weight") +
  theme(axis.text.x = element_text(size=13),legend.position = "none",
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15)) 

p_lay <- "
AB
CD
"
p_basis <- (p2_1 + p2_2 + p2_3 + p2_4) +  
  plot_layout(design = p_lay, axis_titles = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ')') &
  theme(plot.tag.position  = c(.075, 1.037)) 

for (i in c(1, 3)) p_basis [[i]] <- p_basis [[i]] + theme(plot.tag.position  = c(0.125, 1.037))
p_basis

ggsave(plot =  p_basis , width=6, height=5, file = "../figures/viz_basisfunction.png")

#### AIC table
model_name <- c("Default", "5",
                "3", "1")
AIC_values <- c(AIC(mod_$M), AIC(mod_bf[[3]]$M), AIC(mod_bf[[2]]$M), AIC(mod_bf[[1]]$M))

AIC_table <- data.frame(cbind(model_name, round(AIC_values,2))) 
colnames(AIC_table) <- c("Model", "AIC")

AIC_table <- write.csv(AIC_table, file = "../figures/basisfunction_aic.csv", row.names = F)

#### Effect of disturbances
# Full table 
coef_table_all <- coef_uncertainty(mod_$M, percentiles = c(2.5, 50, 97.5), nsim = 400, 
  random_effects = TRUE)%>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value")%>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[1])

coef_table_5K <- coef_uncertainty(mod_bf[[3]]$M, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = FALSE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[2])

coef_table_3K <- coef_uncertainty(mod_bf[[2]]$M, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = TRUE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[3])

coef_table_1K<- coef_uncertainty(mod_bf[[1]]$M, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = FALSE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[4])

# All table together 

coef_table <- rbind(coef_table_all[-1,], coef_table_5K[-1,],
                    coef_table_3K[-1,], coef_table_1K[-1,]) %>% 
  mutate(param = case_when(param == "max_dhw" ~ "Heat stress",
                           param == "max_dhw_lag1" ~ "Heat stress (lag 1)",
                           param == "max_dhw_lag2" ~ "Heat stress (lag 2)",
                           param == "max_cyc" ~ "Cyclone exposure",
                           param == "max_cyc_lag1" ~ "Cyclone exposure (lag 1)",
                           param == "max_cyc_lag2" ~ "Cyclone exposure (lag 2)",
                          TRUE ~ "Intercept"))

coef_table$Name <- as.factor(coef_table$Name)
coef_table$Name_plot <- factor(coef_table$Name, levels = c("Default", "5",
                                                        "3", "1", "none"))
# Make viz
p_coef <- ggplot()+ geom_point(data = coef_table, aes(y=param, x=`50%`)) +
  geom_errorbar(data = coef_table, aes(y = param, xmin = `2.5%`, xmax = `97.5%`), width=.1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + facet_wrap(~Name_plot) +
  theme_pubr() + labs(x = "Effect size", y = "")

ggsave(plot =  p_coef , width=6, height=5, file = "../figures/basisfunction_coef.png")

