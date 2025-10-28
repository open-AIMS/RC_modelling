#' @title Extract disturbance effects
#' @author Julie Vercelloni

# Read model_outputs 

model_list <- list.files(paste0(title_of_run,"/model_outputs/full_run"), recursive = TRUE, pattern = "Rdata")

##################################
############ FRK model(s)  
##################################

model_frk <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "FRK"))

for (l in 1:nrow(model_frk)){
  if(nrow(model_frk)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_frk[l,]))
mod_ <- mod

# Extract covariate effects 
coef_ <- cov_extract(
        model.out = mod_$fit,
        obj_frk = mod_$obj_frk) %>% 
        filter(Type == "fixed") %>%
        mutate(param = ifelse(param == "Intercept", "(Intercept)", param))

# Plot coef_ 
coef_plot <- coef_ %>%
 filter(!param == "(Intercept)") %>%
 mutate(param = case_when(param == "CYC" ~ "Cyclone Exposure",
                          param == "DHW" ~ "Heat Stress",
                          TRUE ~ "NA"))

p_coef <- ggplot(data = coef_plot, aes(x = `50%`, y = param)) + 
  geom_point(size = 3, shape = 21, fill = "black", color = "black") + 
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.1) + 
  theme_pubr() + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  xlab("Effect Size") + ylab("") +
  scale_x_continuous(limits = c(-0.50, 0.005), 
                     breaks = seq(-0.50, 0.01, by = 0.1))

ggsave(p_coef, filename = paste0(title_of_run,"/report/extra/p_coef",str_remove(model_frk[l,], ".Rdata"),".png"),
       width=5, height=4)

}