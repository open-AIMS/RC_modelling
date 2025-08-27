##################################
############ Create folders structure of the pipeline 
##################################

make_folders <- function(title_of_run){
  
  wd <- getwd() 
  
  # lEVEL 1 - name of run 
  pathway = paste0(wd,"/", title_of_run, "/")
  dir.create(pathway)
  ifelse(dir.exists(title_of_run)!=TRUE,print("directory already exists - is it a new simulation?"), FALSE)
  
  # LEVEL 2 - subdirectories 
  
  create_subdir = function (x) {
    for (i in x){
      if(!dir.exists(paste0(pathway, i))){dir.create(paste0(pathway, i))}
    }}
  
  
  create_subdir(c("data","model_outputs", "report"))
  
  # LEVEL 3 - subsubdirectories 
  
  # within model_outputs 
  pathway_2 = paste0(wd,"/", title_of_run, "/model_outputs/")
  
  create_subsubdir = function (x) {
    for (i in x){
      if(!dir.exists(paste0(pathway_2, i))){dir.create(paste0(pathway_2, i))}
    }}
  
  create_subsubdir(c("full_run","leave_out", "predictions"))
  
  # within report 
  pathway_3 = paste0(wd,"/", title_of_run, "/report/")

  create_subsubdir2 = function (x) {
  for (i in x){
    if(!dir.exists(paste0(pathway_3, i))){dir.create(paste0(pathway_3, i))}
  }}

  create_subsubdir2("extra")

}

##################################
############ Generate and assign settings 
##################################

generateSettings <- function(nreefs, nyears){

## Config of the spatio-temporal model
config_sp <- list(
  seed = 1,
  crs = 4326,
  model = "Exp",
  psill = 1,
  range = 15,
  nugget = 0,
  alpha = 2,
  kappa = 1,
  variance = 1,
  patch_threshold = 1.75,
  reef_width = 0.01,
  years = 1:nyears,
  dhw_weight = 0.6,
  cyc_weight = 0.39,
  other_weight = 0.01,
  hcc_cover_range = c(0.1, 0.7),
  hcc_growth = 0.3,
  sc_cover_range = c(0.01, 0.1),
  sc_growth =  0.3
)
assign("config_sp", config_sp, envir = .GlobalEnv)

## Config of sampling design for large scale details
config_lrge <- list(n_locs = nreefs, n_sites = 2, seed = 123)
assign("config_lrge", config_lrge, envir = .GlobalEnv)

## Config for sampling details for fine scale details 
config_fine <- list(
  years =  1:nyears,
  Number_of_transects_per_site = 5,
  Depths = 1,
#  Number_of_frames_per_transect = 100,
#  Points_per_frame = 5,
  ## Note, the following are on the link scale
  hcc_site_sigma = 0.5, # variability in Sites within Locations
  hcc_transect_sigma = 0.2, # variability in Transects within Sites
  hcc_sigma = 0.1, # random noise

  sc_site_sigma = 0.05, # variability in Sites within Locations
  sc_transect_sigma = 0.02, # variability in Transects within Sites
  sc_sigma = 0.01, # random noise

  ma_site_sigma = 0.5, # variability in Sites within Locations
  ma_transect_sigma = 0.2, # variability in Transects within Sites
  ma_sigma = 0.1 # random noise
)
assign("config_fine", config_fine, envir = .GlobalEnv)

## Generate point-based data 
config_pt <- list(
#  Depths = 2,
#  Depth_effect_multiplier = 2,
#  Number_of_transects_per_site = 5,
  Number_of_frames_per_transect = 100,
  Points_per_frame = 50
)
assign("config_pt", config_pt, envir = .GlobalEnv)

config_list <- list(config_sp, config_lrge, config_fine, config_pt)
save(config_list, file = paste0(title_of_run,"/lists_of_parameters.RData"))

## Type of monitoring surveys
surveys <- "fixed"
assign("surveys", surveys, envir = .GlobalEnv)

## Choose the grid size in degree
grid_size <- .1
assign("grid_size", grid_size, envir = .GlobalEnv)
}

#################################
##################################
############ Information about the depth(s) 
##################################

depth_info <- function(x) {
if (x > 1){
  Depth_info = seq(3, 10, length=x)
}else{
  Depth_info = 10
}
return(Depth_info)
}

##################################
############ Select nth element from a vector 
##################################

nth_element <- function(vector, starting_position, n) { 
  vector[seq(starting_position, length(vector), n)] 
  }

##################################
############ Extract covariates 
##################################

extract_cov <- function(predictive_layer, cov_name) {

predictive_layer <- predictive_layer %>%
   st_transform(crs = 4326)

intersectlist <- st_intersects(predictive_layer, cov_name %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(4326)))
  
datlist <- list()
for (i in 1:nrow(pred_layer)){
  intersect <-  intersectlist[[i]]
  datlist[[i]] <-  cov_name[intersect,] %>%
       group_by(Year) %>%
       summarize(mean_value = mean(Value)) %>%
       mutate(tier = i + 999) %>%
       mutate(fYEAR = Year + min(dat$fYEAR) - 1)
  }
  
dat_all <- do.call(rbind, datlist) %>%
    merge(pred_layer) %>%
    st_as_sf() %>%
    dplyr::select(tier, fYEAR, mean_value,geometry)

return(dat_all)
}

##################################
############ FRK model 
##################################

frk_prep <- function(dat){
  
  ## Construct STIDF object from data
  dat$Year <- as.Date(paste0(as.character(dat$fYEAR),"-01-01")) 
  dat$k_Z <- dat$TOTAL                                         
  lon_idx <- which(names(dat) == "LONGITUDE")                  
  lat_idx <- which(names(dat) == "LATITUDE")
  STObj <- stConstruct(x = dat,                               
                       space = c(lon_idx, lat_idx), 
                       time = "Year",                      
                       interval = TRUE)     
  
  ## Predictive layer
  HexPred_sp <- as_Spatial(hexpred)                                   
  nHEX <- nrow(subset(HexPred_sp, fYEAR == unique(dat$fYEAR)[1]))       
  nYEAR <- length(unique(HexPred_sp@data$fYEAR))            
  
  HexPred_sp@data$n_spat <- rep(1:nHEX, each = nYEAR)   
  BAUs_spat <- subset(HexPred_sp, fYEAR == unique(dat$fYEAR)[1])        
  coordnames(BAUs_spat) <- c("LONGITUDE", "LATITUDE")
  
  nrow(BAUs_spat@data)
  
  ## Construct spatio-temporal BAUs (will not contain covariate information for now)
  ST_BAUs <- auto_BAUs(manifold = STplane(),
                       data = STObj,
                       spatial_BAUs = BAUs_spat,
                       tunit = "years")
  
  ST_BAUs <- ST_BAUs[, 1:nYEAR, 1:2]                 
  ST_BAUs$fYEAR <- as.character(ST_BAUs$t + min(dat$fYEAR)[1]-1)    
  ST_BAUs$n_spat <- rep(1:nHEX, nYEAR)              
  
  HexPred_sp@data$fYEAR <- as.character(HexPred_sp@data$fYEAR) 
  HexPred_sp@data$tier <- as.factor(HexPred_sp@data$tier) 
  
  ST_BAUs@data <- left_join(ST_BAUs@data, HexPred_sp@data , by = c("fYEAR","n_spat")) 
  
  ST_BAUs$fs <- 1                  
  ST_BAUs@sp@proj4string  <- CRS()  
  
  head(ST_BAUs@data)
  head(HexPred_sp@data)
  
  ## Covariates must only be in BAUs, so remove covariates associated with data
  overlapping_fields <- intersect(names(STObj@data), names(ST_BAUs@data))
  STObj@data[,overlapping_fields] <- NULL
  
  ## Create basis functions
  basis <- auto_basis(STplane(),
                      STObj,
                      tunit = "years",
                      #nres = 2L, # for development
                      nres = 3L, # for final run
                      regular = TRUE)

  
  obj_frk <- list("ST_BAUs" = ST_BAUs, "STObj" = STObj, "basis" = basis)
  return(obj_frk)
  
}

predictions_FRK <- function(model.out){
  
  pred <- predict(model.out, type = c("mean"))
  
  # Extracting posterior distributions of predictive locations 
  
  post_dist_df <- as.data.frame(pred$MC$mu_samples) %>% 
    mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
    mutate(tier = obj_frk$ST_BAUs@data$tier) %>%
    mutate(id_loc = row_number()) %>%
    tidyr::pivot_longer(!c(fYEAR,tier,id_loc),
                        names_to = "draw", 
                        values_to = "pred"
    )
  
  # Summary predictions at tier level
  hexpred$tier <- as.factor(hexpred$tier)
  
  pred_sum_sf <- post_dist_df %>% group_by(fYEAR,tier) %>% 
    median_hdci(pred)%>%
    inner_join(hexpred %>% group_by(tier) %>% 
    summarize() %>% dplyr::select(geometry,tier)) %>% 
    st_as_sf(sf_column_name = "geometry") %>%
    mutate(Unc = .upper - .lower) %>%
    mutate(tier_fYEAR = paste0(tier,fYEAR)) %>%
    dplyr::select(fYEAR, tier, pred, .lower, .upper, Unc, tier_fYEAR)

  return(pred_sum_sf)
}

predictions_FRK_broad <- function(model.out){
  
  pred <- predict(model.out, type = c("mean"))
  
  # Extracting posterior distributions of predictive locations 
  
  post_dist_df <- as.data.frame(pred$MC$mu_samples) %>% 
    mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
    mutate(tier = obj_frk$ST_BAUs@data$tier) %>%
    mutate(id_loc = row_number()) %>%
    tidyr::pivot_longer(!c(fYEAR,tier,id_loc),
                        names_to = "draw", 
                        values_to = "pred"
    )
  
  # Summary predictions at tier level
  hexpred$tier <- as.factor(hexpred$tier)
  
  pred_sum_sf <- post_dist_df %>% group_by(fYEAR) %>% 
    median_hdci(pred)

  return(pred_sum_sf)
}

##################################
############ plotting 
##################################

plot_predictions <- function(dat){

 if (length(unique(dat$fYEAR)) <= 16) {
  ncols <- 4
} else {
  ncols <- 6
}

pal_pred <- lacroix_palette("Pamplemousse", n = 100, type = "continuous")
  
p_pred <- ggplot() + 
  geom_sf(data = pred_sum_sf, aes(fill = pred*100), col = "transparent") +
  facet_wrap(~fYEAR, ncol = ncols) +   scale_fill_gradientn(
    colours = pal_pred,    
    name = "Predicted coral cover (%)" 
  ) + 
  theme_pubr() +
  xlab("longitude") +
  ylab("latitude") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )

return(p_pred)
}

plot_predictions_unc <- function(dat){

 if (length(unique(dat$fYEAR)) <= 16) {
  ncols <- 4
} else {
  ncols <- 6
}

pal_unc<- wes_palette("Zissou1", 100, type = "continuous")

  p_unc <- ggplot() + 
    geom_sf(data = pred_sum_sf, aes(fill = Unc*100), col = "transparent") +
    facet_wrap(~fYEAR, ncol = ncols) +  scale_fill_gradientn(
      colours = pal_unc,
      name = "Uncertainty range (%)" 
    )  + 
    theme_pubr() + 
    xlab("longitude") +
    ylab("latitude") +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      strip.background = element_blank(),
      strip.text = element_text(size = 12)
    )

return(p_unc)
}

plot_traj <- function(dat, dat_pred){
p_traj <- ggplot() + 
  geom_line(data = dat, aes(x = fYEAR, y = (COUNT/TOTAL)*100, group = interaction(as.factor(TRANSECT_NO), REEF_NAME)), 
            show.legend = FALSE, linewidth=.1, col="grey30") + 
  geom_ribbon(data = dat_pred, aes(x=fYEAR,ymin=.lower*100, ymax=.upper*100, group=1),alpha=.2, fill ="#72b0d3") +
  geom_line(data = dat_pred, aes(x=fYEAR, y=pred*100, group=1),size=.4) +
  facet_wrap(~tier, ncol=3) +
  ylab("Cover") + xlab("Year")+
  theme_pubr() + 
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),
        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white')) +
  scale_x_discrete(breaks= nth_element(unique(dat$fYEAR),1,4))

return(p_traj)
}

plot_traj_broad <- function(dat, GRMF_all){
plot_traj_broad <- ggplot() +
  geom_ribbon(data = pred_sum_sf %>% data.frame(), aes(x = fYEAR, ymin=.lower, ymax=.upper, group=1), alpha=.2, fill="#00FFC2")+
  geom_line(data = pred_sum_sf %>% data.frame(), aes(x=fYEAR, y=pred,group=1), col="black", linewidth=1.1)+
  geom_point(data = pred_sum_sf %>% data.frame(), aes(x=fYEAR, y=pred), col="black", size=2.1)+
  geom_line(data = GRMF_all, aes(x=as.factor(fYEAR), y=true, group=1),linewidth = .4, col = "blue", linetype = "dashed") +
  xlab("Year") +ylab("Coral cover")+
  theme_pubr() + 
  ylim(0,1) + 
  theme(axis.text.x = element_text(size=13),legend.position = "none",
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15))+
  scale_x_discrete(breaks= nth_element(unique(dat$fYEAR),1,4))

return(plot_traj_broad)
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
RMSPE <- function(z, pred) {
  Y <- (z - pred)^2
  sqrt(mean(Y))
}

crps <- function(obs, pred, ...){
  if(is.null( dim(pred)) & length(pred)==2){mu <- pred[1];
  sigma <- pred[2]} else {
    mu<- as.numeric( pred[,1] ); sigma <- as.numeric( pred[,2]) }
  
  z <- (obs-mu)/sigma ## center and scale
  
  crps<- sigma * (z*(2*pnorm(z,0,1)-1) + 2*dnorm(z,0,1) - 1/sqrt(pi))
  ign <-  0.5*log(2*pi*sigma^2) + (obs - mu)^2/(2*sigma^2)
  pit <- pnorm(obs, mu,sigma )
  
  return(list(crps = crps, CRPS = mean(crps, na.rm = TRUE), ign = ign, IGN = mean(ign), pit = pit) )
  
}

#################################
############# automated report
#################################

## Grab png files associated with name_plot and add them together 

print_attr_plots_level1 <- function(list_folders, name_plot) {


  # List files matching the pattern in the specified path
  file_mod <- list.files(pattern = name_plot, full.names = TRUE)
  
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

print_attr_plots_level2 <- function(list_folders, name_plot) {

  # List files matching the pattern in the specified path
  file <- list.files(paste0(list_folders, FIGURES_PATH),
    pattern = name_plot,
    full.names = TRUE
  )

  # Generate plots for each file
  p <- lapply(file, function(file) {
    ggplot2::ggplot() +
      ggplot2::annotation_custom(grid::rasterGrob(png::readPNG(file), interpolate = TRUE)) +
      ggplot2::theme_void()
  })
  
  # Combine the plots using patchwork
  combined_mod <- patchwork::wrap_plots(p, ncol = 2)
  return(combined_mod)
}


print_attr_plots_level3 <- function(list_folders, name_plot, subtitle = FALSE, plot_labels = FALSE, label = "", n_cols) {

  # Debug: print the list_folders argument
  cat("####", list_folders, "\n\n")
  
  # Define file paths
  file <- list.files(paste0(list_folders, FIGURES_PATH),
    pattern = name_plot,
    full.names = TRUE
  )
  
  # Generate plots for each file
  p <- lapply(file, function(file) {
    plot <- ggplot2::ggplot() +
      ggplot2::annotation_custom(grid::rasterGrob(png::readPNG(file), interpolate = TRUE)) +
      ggplot2::theme_void() #+ 
      #ggplot2::theme(plot.subtitle = element_text(size = 18, hjust = 0.5) )
    
    # Add subtitle if subtitle = TRUE
    if (subtitle) {
      plot <- plot + ggplot2::labs(
        subtitle = sub(
          pattern = paste0("^.*",str_replace(name_plot, "\\[.*", ""), "(.*)\\.png$"), # Use name_plot dynamically
          replacement = "\\1",
          x = file
        )
      )
    }
   return(plot)
  })
  
  # Combine the plots using patchwork
  combined_mod <- patchwork::wrap_plots(p, ncol = n_cols)
  
  # Print plots with labels if plot_labels = TRUE
  if (plot_labels) {
    div_label <- paste0("#fig-pred-data", list_folders)
    cat("\n::: {", div_label, "}\n", sep = "")
    
    print(p) # Print combined plot
    cat(paste0("\n\n", label, " ", list_folders, "\n"), sep = "")
    cat("\n:::\n")
  }
  
  # Return combined plot
  return(combined_mod)
}

# Adapted from https://github.com/Ringomed/ggvanced/blob/main/R/ggspider.R
ggspider_2 <- function(p_data,
                     ci_data = NULL,
                     polygon = TRUE,
                     scaled = FALSE,
                     draw_axis = TRUE,
                     n_labels = 5,
                     zero_centered = FALSE,
                     subset = NULL,
                     reorder_axis = NULL,
                     background_color = "gray99",
                     area_fill = TRUE,
                     fill_opacity = 0.05,
                     central_distance = 0.2,
                     axis_name_offset = 0.2,
                     digit_rounding = 2,
                     axis_label_font_size = NULL,
                     axis_label_font_face = NULL,
                     axis_name_font_size = NULL,
                     axis_name_font_face = NULL
){

  legend_title <- names(p_data)[[1]]
  p_data <- p_data %>% dplyr::rename(group = 1) %>% dplyr::mutate(group = factor(group))

  if(zero_centered == TRUE){
    zero_tibble <- tibble::as_tibble(as.list(setNames(rep(0, ncol(p_data)), names(p_data)))) %>% dplyr::mutate(group = "zero_centered")
    p_data <- p_data %>% dplyr::bind_rows(zero_tibble)
  }
  
  if(!is.null(reorder_axis)){
    p_data <- p_data %>% dplyr::mutate(dplyr::across(dplyr::all_of(reorder_axis), ~ -.))

    if(!is.null(ci_data)){
      ci_data <- ci_data %>%
        dplyr::mutate(c = min,
               min = dplyr::case_when(parameter %in% reorder_axis ~ -max, TRUE ~ min),
               max = dplyr::case_when(parameter %in% reorder_axis ~ -c, TRUE ~ max))
    }
  }


  if(!is.null(ci_data)){
    ci_data <- ci_data %>% dplyr::rename(group = 1, parameter = 2)
  }

  circle_coords <- function(r, n_axis = ifelse(polygon == TRUE, ncol(p_data) - 1, 100)){
    fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
    x <- r*cos(fi)
    y <- r*sin(fi)

    tibble::tibble(x, y, r)
  }

  step_1 <- purrr::map_df(seq(0, 1, 0.25) + central_distance, circle_coords) %>%
      ggplot2::ggplot(ggplot2::aes(x, y)) +
      ggplot2::geom_polygon(data = circle_coords(1 + central_distance), alpha = 1, fill = background_color, lty = 2) +
      ggplot2::geom_path(ggplot2::aes(group = r), lty = 2, alpha = 0.5) +
      ggplot2::theme_void()


  axis_coords <- function(n_axis){
    fi <- seq(0, (1 - 1/n_axis)*2*pi, (1/n_axis)*2*pi) + pi/2
    x1 <- central_distance*cos(fi)
    y1 <- central_distance*sin(fi)
    x2 <- (1 + central_distance)*cos(fi)
    y2 <- (1 + central_distance)*sin(fi)

    tibble::tibble(x = c(x1, x2), y = c(y1, y2), id = rep(1:n_axis, 2))
  }

  if(is.null(ci_data)){
    text_data <- p_data %>%
      dplyr::select(-group) %>%
      purrr::map_df(~ min(.) + (max(.) - min(.)) * seq(0, 1, 1/(n_labels - 1))) %>%
      dplyr::mutate(r = seq(0, 1, 1/(n_labels - 1))) %>%
      tidyr::pivot_longer(-r, names_to = "parameter", values_to = "value")
  }else{
    text_data <-
      p_data %>%
      dplyr::select(-group) %>%
      purrr::map_df(~ min(.) + (max(.) - min(.)) * seq(0, 1, 1/(n_labels - 1))) %>%
      dplyr::mutate(r = seq(0, 1, 1/(n_labels - 1))) %>%
      tidyr::pivot_longer(-r, names_to = "parameter", values_to = "value") %>%
      dplyr::select(r, parameter) %>%
      dplyr::left_join(
        ci_data %>%
        dplyr::group_by(parameter) %>%
        dplyr::reframe(min = min(min), max = max(max)) %>%
        dplyr::group_by(parameter) %>%
        dplyr::mutate(data = list(tibble::tibble(r = seq(0, 1, 1/(n_labels - 1)), value = min + (max-min)*r))) %>%
        tidyr::unnest("data") %>%
        dplyr::select(r, parameter, value) %>%
        dplyr::arrange(r), by = c("parameter", "r")
      )
  }

  text_coords <- function(r, n_axis = ncol(p_data) - 1){
    fi <- seq(0, (1 - 1/n_axis)*2*pi, (1/n_axis)*2*pi) + pi/2 + 0.01*2*pi/r
    x <- r*cos(fi)
    y <- r*sin(fi)

    tibble::tibble(x, y, r = r - central_distance)
  }

  labels_data <- purrr::map_df(seq(0, 1, 1/(n_labels - 1)) + central_distance, text_coords) %>%
    dplyr::bind_cols(text_data %>% dplyr::select(-r))


  if(!is.null(reorder_axis)){
    labels_data <- labels_data %>% dplyr::mutate(value = case_when(parameter %in% reorder_axis ~ -value,
                                                                   TRUE ~ value))
  }


  rescaled_coords <- function(r, n_axis){
    fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
    tibble::tibble(r, fi) %>% dplyr::mutate(x = r*cos(fi), y = r*sin(fi)) %>% dplyr::select(-fi)
  }

  if(!is.null(ci_data)){
    pdata_long_rescaled <- ci_data %>%
      dplyr::select(group, parameter, y = min) %>%
      mutate(measure = "min") %>%
      bind_rows(ci_data %>%
                   dplyr::select(group, parameter, y = max) %>%
                   mutate(measure = "max")) %>%
      bind_rows(p_data %>%
                  tidyr::pivot_longer(-1, names_to = "parameter", values_to = "mean") %>%
                  dplyr::select(group, parameter, y = mean) %>%
                  mutate(measure = "mean")) %>%
      group_by(parameter) %>%
      dplyr::mutate(y = scales::rescale(y))

    rescaled_min <- pdata_long_rescaled %>% filter(measure == "min") %>%
      dplyr::select(-measure) %>%
      pivot_wider(names_from = "parameter", values_from = "y") %>%
      dplyr::mutate(copy = dplyr::pull(., 2)) %>% #da se moze geom_path spojiti opet na pocetnu tocku
      tidyr::pivot_longer(-group, names_to = "parameter", values_to = "value") %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(coords = rescaled_coords(value + central_distance, ncol(p_data) - 1)) %>%
      tidyr::unnest(cols = c(coords)) %>%
      dplyr::select(group, parameter, xmin = x, ymin = y)

    rescaled_max <- pdata_long_rescaled %>% filter(measure == "max") %>%
      dplyr::select(-measure) %>%
      pivot_wider(names_from = "parameter", values_from = "y") %>%
      dplyr::mutate(copy = dplyr::pull(., 2)) %>% #da se moze geom_path spojiti opet na pocetnu tocku
      tidyr::pivot_longer(-group, names_to = "parameter", values_to = "value") %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(coords = rescaled_coords(value + central_distance, ncol(p_data) - 1)) %>%
      tidyr::unnest(cols = c(coords)) %>%
      dplyr::select(group, parameter, xmax = x, ymax = y)

    rescaled_ci <- rescaled_min %>% dplyr::left_join(rescaled_max, by = dplyr::join_by(group, parameter))

    rescaled_ci_alt <- rescaled_min %>% rename(x = xmin, y = ymin) %>% dplyr::bind_rows(rescaled_max %>% rename(x = xmax, y = ymax))
  }else{
    pdata_long_rescaled <- p_data %>%
                  tidyr::pivot_longer(-1, names_to = "parameter", values_to = "mean") %>%
                  dplyr::select(group, parameter, y = mean) %>%
                  mutate(measure = "mean") %>%
      group_by(parameter) %>%
      dplyr::mutate(y = scales::rescale(y))
  }

  rescaled_data <- pdata_long_rescaled %>% filter(measure == "mean") %>%
    dplyr::select(-measure) %>%
    pivot_wider(names_from = "parameter", values_from = "y") %>%
    dplyr::mutate(copy = dplyr::pull(., 2)) %>% #da se moze geom_path spojiti opet na pocetnu tocku
    tidyr::pivot_longer(-group, names_to = "parameter", values_to = "value") %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(coords = rescaled_coords(value + central_distance, ncol(p_data) - 1)) %>%
    tidyr::unnest(cols = c(coords))

  rescaled_data <- if(is.null(subset) & zero_centered == FALSE){
    rescaled_data
  } else if(is.null(subset) & zero_centered == TRUE) {
    rescaled_data %>% dplyr::filter(group != "zero_centered")
  } else {
      rescaled_data %>% dplyr::filter(group %in% subset)
  }

  step_1 +
    {if(draw_axis == TRUE) ggplot2::geom_line(data = axis_coords(ncol(p_data) - 1), ggplot2::aes(x, y, group = id), alpha = 0.3)} +
    {if(!is.null(ci_data)) ggplot2::geom_segment(data = rescaled_ci,
                          ggplot2::aes(x = xmin, y = ymin, xend = xmax, yend = ymax, col = group, lwd = 1),
                          alpha = 0.5, lineend = "square", show.legend = FALSE,
                          arrow = arrow(ends = "both", angle = 90, length = unit(.1,"cm")))} +
    scale_linewidth(range = c(0, 4)) +
    ggplot2::geom_point(data = rescaled_data, ggplot2::aes(x, y, group = group, col = group), size = 2, stroke = 2) +
    ggplot2::geom_path(data = rescaled_data, ggplot2::aes(x, y, group = group, col = group), size = 1) +
    {if(area_fill == TRUE) ggplot2::geom_polygon(data = rescaled_data, ggplot2::aes(x, y, group = group, col = group, fill = group), size = 1, alpha = fill_opacity, show.legend = FALSE)} +
    {if(scaled == TRUE){
      ggplot2::geom_text(data = labels_data %>% dplyr::filter(parameter == labels_data$parameter[[1]]), ggplot2::aes(x, y, label = r), alpha = 0.65,
                         family = theme_get()$text[["family"]],
                         size = ifelse(is.null(axis_label_font_size), theme_get()$text[["size"]]/2.75, axis_label_font_size),
                         fontface = ifelse(is.null(axis_label_font_face), "plain", axis_label_font_face))
    }else{
        ggplot2::geom_text(data = labels_data, ggplot2::aes(x, y, label = round(value, digit_rounding)), alpha = 0.65,
                           family = theme_get()$text[["family"]],
                           size = ifelse(is.null(axis_label_font_size), theme_get()$text[["size"]]/2.75, axis_label_font_size),
                           fontface = ifelse(is.null(axis_label_font_face), "plain", axis_label_font_face))
      }
    } +
    ggplot2::geom_text(data = text_coords(1 + central_distance + axis_name_offset), ggplot2::aes(x, y), label = labels_data$parameter[1:(ncol(p_data)-1)],
                       family = theme_get()$text[["family"]],
                       size = ifelse(is.null(axis_name_font_size), theme_get()$text[["size"]]/2.75, axis_name_font_size),
                       fontface = ifelse(is.null(axis_name_font_face), "plain", axis_name_font_face)) +
    ggplot2::labs(col = legend_title) +
    ggplot2::theme(legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 12),
                   legend.title = ggplot2::element_text(size = 12))
}
