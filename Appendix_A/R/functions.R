frk_prep <- function(data.grp.tier, HexPred_reefid2) {

      # Convert fYEAR to Date (start of year)
      data.grp.tier$Year <- as.Date(paste0(as.character(data.grp.tier$fYEAR), "-01-01"))
      data.grp.tier$k_Z <- data.grp.tier$TOTAL  # number of trials
      
      lon_idx <- which(names(data.grp.tier) == "LONGITUDE")
      lat_idx <- which(names(data.grp.tier) == "LATITUDE")
      
      STObj <- stConstruct(x = data.grp.tier,
                           space = c(lon_idx, lat_idx),
                           time = "Year",
                           interval = TRUE)
      
      # Convert HexPred_reefid2 to sp object for BAU construction
      HexPred_sp <- as_Spatial(HexPred_reefid2)
      nHEX <- nrow(subset(HexPred_sp, fYEAR == min(HexPred_sp@data$fYEAR)))
      nYEAR <- length(unique(HexPred_sp@data$fYEAR))
      
      HexPred_sp@data$n_spat <- rep(1:nHEX, each = nYEAR)
      BAUs_spat <- subset(HexPred_sp, fYEAR == min(HexPred_sp@data$fYEAR))
      coordnames(BAUs_spat) <- c("LONGITUDE", "LATITUDE")
      
      ST_BAUs <- FRK::auto_BAUs(manifold = STplane(),
                           data = STObj,
                           spatial_BAUs = BAUs_spat,
                           tunit = "years")
      
      ST_BAUs <- ST_BAUs[, 1:nYEAR, 1:2]  # Remove extra year
      ST_BAUs$fYEAR <- as.character(ST_BAUs$t + (min(as.numeric(HexPred_sp@data$fYEAR)) - 1))
      ST_BAUs$n_spat <- rep(1:nHEX, nYEAR)
      
      ST_BAUs@data$fYEAR <- as.integer(ST_BAUs@data$fYEAR)
      ST_BAUs@data <- dplyr::left_join(ST_BAUs@data, HexPred_sp@data, by = c("fYEAR", "n_spat"))
      
      ST_BAUs$fs <- 1
      ST_BAUs@sp@proj4string <- CRS()
      
      ST_BAUs@data$reefid <- as.factor(as.character(ST_BAUs@data$reefid))
      ST_BAUs@data$yearid <- as.factor(ST_BAUs@data$fYEAR)
      
      # Remove overlapping covariates from STObj
      overlapping_fields <- intersect(names(STObj@data), names(ST_BAUs@data))
      STObj@data[, overlapping_fields] <- NULL
      
      # Create basis functions
      basis <- FRK::auto_basis(STplane(),
                          ST_BAUs,
                          tunit = "years",
                          #nres = 2L, #for dev
                          nres = 3L,
                          regular = TRUE)

  obj_frk <- list("ST_BAUs" = ST_BAUs, "STObj" = STObj, "basis" = basis)
  return(obj_frk)
}

make_dharma_res <- function(mod_) {
  response <- mod_$data.grp.tier$COUNT
  preds <- simulate(mod_$M, nsim = 1000) 
  DHARMa::createDHARMa(
    simulatedResponse = preds,
    observedResponse = response,
    integerResponse = TRUE
  )
}

cov_extract <- function(model.out, obj_frk){

coef_table_all <- coef_uncertainty(model.out, percentiles = c(2.5, 50, 97.5), nsim = 400, 
  random_effects = TRUE)%>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value")%>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value)

return(coef_table_all)
}

build_table <- function(vary_var, effect_col, range_vals = NULL) {
  fixed_vals <- coef_scaled %>% 
                filter(! param == vary_var) %>%
                dplyr::select(param, !!sym(effect_col)) %>%
                deframe()

  map_dfr(range_vals, function(x) {
    tibble(
      !!!as.list(fixed_vals),
      !!vary_var       := x
    ) %>%
      relocate(!!sym(vary_var), .after = `(Intercept)`)
  })
}

predict_cover <- function(effect_col) {
  coef_value <- coef_scaled %>%
    filter(param == var) %>%
    pull(all_of(effect_col))

  table <- build_table(
    vary_var   = var,
    effect_col = effect_col,
    range_vals = range_vals
  ) %>%
    mutate(new_var = .data[[var]] * coef_value) %>%
    data.frame() %>%
    dplyr::select(!var)
  plogis(rowSums(table))
}

plot_cond <- function(range_vals, var, xlabel, ylabel){
cond_table <- tibble(
  !!sym(var) := range_vals,
  value_coral_50  = predict_cover("effect_50"),
  value_coral_low = predict_cover("effect_low"),
  value_coral_high = predict_cover("effect_high")
)

p <- ggplot(cond_table, aes(x = !!sym(var))) +
  geom_ribbon(aes(ymin = value_coral_low * 100, ymax = value_coral_high * 100), 
              fill = "gray66", alpha = 0.3) +
  geom_smooth(aes(y = value_coral_50 * 100), color = "black", size = 1, se = FALSE) +
  labs(
    x = xlabel,
    y = ylabel
  ) +
  theme_pubr() +
  scale_x_continuous(
    limits = c(0, max(range_vals)),
    breaks = seq(0, max(range_vals), by = 4)
  ) +
  ylim(0, 40)

return(p)
}

predict_newdata <- function(mod_M, new_BAU_data, hexpred, obj_frk, title_n){
  
  newpred <- predict(mod_M, new_BAU_data = new_BAU_data)
  
  post_dist_df <- as.data.frame(newpred$MC$mu_samples) %>% 
    mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
    mutate(Tier5 = obj_frk$ST_BAUs@data$Tier5) %>%
    mutate(id_loc = row_number()) %>%
    tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc),
                        names_to = "draw", 
                        values_to = "pred")    
  
  # Spatial maps 
  hexpred$Tier5 <- as.factor(hexpred$Tier5)
  
  pred_sum_sf <- post_dist_df %>% group_by(Tier5) %>% 
    median_hdci(pred) %>%
    inner_join(hexpred %>% group_by(Tier5) %>% 
                 summarize() %>% dplyr::select(geometry, Tier5)) %>% 
    st_as_sf(sf_column_name = "geometry") %>%
    mutate(Unc = .upper - .lower) %>%
    dplyr::select(Tier5, pred, .lower, .upper, Unc) %>%
    mutate(Name = "Dummy")
  
  pal_pred <- LaCroixColoR::lacroix_palette("Pamplemousse", n = 100, type = "continuous")

  p <- ggplot() + 
    geom_sf(data = pred_sum_sf, aes(fill = pred*100), col = "transparent") +
    coord_sf(expand = F) +
    scale_fill_gradientn(
      colours = pal_pred,    
      name = "Predicted coral cover (%)",
      limits = c(0,50), 
    ) + 
    theme_pubr() +
    xlab("") +
    ylab("") +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      strip.background = element_blank(),
      strip.text = element_text(size = 12)
    ) +
    ggtitle(title_n) +
    theme(plot.title = element_text(size = 13, hjust = 0.5),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
           axis.text.y = element_text(size = 10),
           axis.title.x = element_text(size = 12),
           axis.title.y = element_text(size = 12))  +
  scale_y_continuous(
    breaks = c(-16, -17, -18),
    labels = c("16", "17", "18")
  ) +
  scale_x_continuous(
    breaks = c(145.4, 145.8, 146.2, 146.6),
    labels = c("145.4","145.8", "146.2", "146.6"))

p_insert <- ggplot(pred_sum_sf, aes(x = pred * 100, y = "", fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_gradientn(
    colours = pal_pred,
    name = "Predicted coral cover (%)",
    limits = c(0, 50)
  ) +
  labs(x = NULL, y = NULL) +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

p_all <- p +
  inset_element(p_insert, 0.6, 0.6, 1, 1) 

  return(p_all)
}

predict_year <- function(mod_M, new_BAU_data, hexpred, obj_frk, title_n){
  
  newpred <- predict(mod_M, new_BAU_data = new_BAU_data)
  
  post_dist_df <- as.data.frame(newpred$MC$mu_samples) %>% 
    mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
    mutate(Tier5 = obj_frk$ST_BAUs@data$Tier5) %>%
    mutate(id_loc = row_number()) %>%
    tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc),
                        names_to = "draw", 
                        values_to = "pred")    
  
  # Spatial maps 
  hexpred$Tier5 <- as.factor(hexpred$Tier5)
  
  pred_sum_sf <- post_dist_df %>% group_by(fYEAR, Tier5) %>% 
    median_hdci(pred) %>%
    inner_join(hexpred %>% group_by(Tier5) %>% 
               summarize() %>% dplyr::select(geometry, Tier5)) %>% 
    st_as_sf(sf_column_name = "geometry") %>%
    mutate(Unc = .upper - .lower) %>%
    mutate(tier_fYEAR = paste0(Tier5,fYEAR)) %>%
    dplyr::select(fYEAR, Tier5, pred, .lower, .upper, Unc, tier_fYEAR)

  pal_pred <- LaCroixColoR::lacroix_palette("Pamplemousse", n = 100, type = "continuous")

  p <- ggplot() + 
    geom_sf(data = pred_sum_sf, aes(fill = pred*100), col = "transparent") +
    coord_sf(expand = F) +
    scale_fill_gradientn(
      colours = pal_pred,    
      name = "Coral Cover (%)",
      limits = c(0,50), 
    ) + 
    facet_wrap(~ fYEAR) +
    theme_pubr() +
    xlab("Longitude") +
    ylab("Latitude") +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      strip.background = element_blank(),
      strip.text = element_text(size = 12)
    ) +
    ggtitle(title_n) +
    theme(plot.title = element_text(size = 13, hjust = 0.5),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
           axis.text.y = element_text(size = 10),
           axis.title.x = element_text(size = 12),
           axis.title.y = element_text(size = 12)) +
  scale_y_continuous(
    breaks = c(-16, -17, -18),
    labels = c("16", "17", "18")
  ) +
  scale_x_continuous(
    breaks = c(145.4, 145.8, 146.2, 146.6),
    labels = c("145.4","145.8", "146.2", "146.6"))

pal_unc<- wes_palette("Zissou1", 100, type = "continuous")

  p_unc <- ggplot() + 
    geom_sf(data = pred_sum_sf, aes(fill = Unc*100), col = "transparent") +
    coord_sf(expand = F) +
    scale_fill_gradientn(
      colours = pal_unc,    
      name = "Uncertainty range (%)",
      limits = c(0,80), 
    ) +    
    facet_wrap(~ fYEAR) +
    theme_pubr() +
    xlab("Longitude") +
    ylab("Latitude") +
    theme(
     # axis.text = element_blank(),
     # axis.ticks = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      strip.background = element_blank(),
      strip.text = element_text(size = 12)
    ) +
    ggtitle(title_n) +
    theme(plot.title = element_text(size = 13, hjust = 0.5),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
           axis.text.y = element_text(size = 10),
           axis.title.x = element_text(size = 12),
           axis.title.y = element_text(size = 12)) +
  scale_y_continuous(
    breaks = c(-16, -17, -18),
    labels = c("16", "17", "18")
  ) +
  scale_x_continuous(
    breaks = c(145.4, 145.8, 146.2, 146.6),
    labels = c("145.4","145.8", "146.2", "146.6"))

 return(list(p = p, p_unc = p_unc))
}

##################################
############ predictive indicators  
##################################

## 95% coverage
coverage95 <- function(z, lower, upper) {
  abs(0.95 - (sum((z < upper) & (z > lower)) / length(z)))
}

## 95% interval score
IS95 <- function(true, lower, upper) {
  alpha = 0.05
  pred95l <- lower 
  pred95u <- upper
  ISs <- (pred95u - pred95l) + 2/alpha * (pred95l - true) * (true < pred95l) +
    2/alpha * (true - pred95u) * (true > pred95u)
  mean(ISs)
}

## Root-mean-squared prediction error
RMSPE <- function(z,pred) {
  Y <- (z - pred)^2
  sqrt(mean(Y))
}

crps <- function(obs, pred, ...)
  ## Tilmann Gneiting's crps code, assumes pred is either a vector of length
  ## 2 (mu, sig) or a matrix of mu and sig if each forcast is different
{
  if(is.null( dim(pred)) & length(pred)==2){mu <- pred[1];
  sigma <- pred[2]} else {
    mu<- as.numeric( pred[,1] ); sigma <- as.numeric( pred[,2]) }
  
  z <- (obs-mu)/sigma ## center and scale
  
  crps<- sigma * (z*(2*pnorm(z,0,1)-1) + 2*dnorm(z,0,1) - 1/sqrt(pi))
  ign <-  0.5*log(2*pi*sigma^2) + (obs - mu)^2/(2*sigma^2)
  pit <- pnorm(obs, mu,sigma )
  
  return(list(crps = crps, CRPS = mean(crps), ign = ign, IGN = mean(ign), pit = pit) )
  
}

print_attr_plots <- function(appendix_dir, name_plot) {


  # List files matching the pattern in the specified path
  file_mod <- list.files(path = appendix_dir, pattern = name_plot, full.names = TRUE)
  
  # Generate plots for each file
  p_mod <- lapply(file_mod, function(file_mod) {
    ggplot2::ggplot() +
      ggplot2::annotation_custom(grid::rasterGrob(png::readPNG(file_mod), interpolate = TRUE)) +
      ggplot2::theme_void() 
  })
  
  # Combine the plots using patchwork
  combined_mod <- patchwork::wrap_plots(p_mod, ncol = 1)
  return(combined_mod)
}

