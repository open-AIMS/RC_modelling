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

plot_traj <- function(pred_df, unique_tier){
p_traj <- ggplot() + 
   geom_ribbon(data = pred_df %>% filter(tier %in% unique_tier), 
                      aes(x=as.numeric(fYEAR),ymin=.lower*100, ymax=.upper*100, group=1),alpha=.3, fill ="#72b0d3") +
  geom_point(data = pred_df %>% filter(tier %in% unique_tier), 
                   aes(x = as.numeric(fYEAR), y = (COUNT_t / TOTAL_t)*100), fill = "grey66", col = "black", size = 2.5, alpha = .6, shape = 20) + 
 geom_line(data = pred_df %>% filter(tier %in% unique_tier), 
                    aes(x=as.numeric(fYEAR), y=pred*100, group=1),size=.8) + 
  facet_wrap(~tier, ncol=3) +
  ylab("Coral cover (%)") + xlab("Year") +  theme_pubr() +
  theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=10),axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white')) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 20)
    )

p_traj_no <- ggplot(pred_df %>% filter(tier %in% unique_tier)) + 
  geom_ribbon(aes(x=as.numeric(fYEAR),ymin=.lower*100, ymax=.upper*100, group=1),alpha=.9, fill= "#D68A8A") +
  geom_line(aes(x=as.numeric(fYEAR), y=pred*100, group=1), linewidth=.6) +
  facet_wrap(~tier) + 
  ylab("Coral cover (%)") + xlab("Year") +  theme_pubr() +
  theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=10),axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white')) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 20)
    )

return(list(p_traj = p_traj, p_traj_no = p_traj_no))
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
############# extract disturbance effects
#################################


cov_extract <- function(model.out, obj_frk){

coef_table_all <- FRK::coef_uncertainty(model.out, percentiles = c(2.5, 50, 97.5), nsim = 400, 
  random_effects = TRUE)%>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value")%>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value)

return(coef_table_all)
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

