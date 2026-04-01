#############################################
######### Run predictive model
#############################################

mid_quant_75 <- function(y){
  y <- as.numeric(as.character(y))
  y <- y[!is.na(y)]
  
  # Case 1: all zeros → return 0
  if (length(y) == 0 || all(y == 0)) {
    return(0)
  }
  
  # Case 2: otherwise compute mid-quantile
  ymid <- Qtools::midquantile(y, probs = 3/4)
  return(ymid$y)
}

select_covariates <- function(x) {

  variables_name_full <- names(x)
  variables_name_full <- grep("^max", variables_name_full, value = TRUE)
  
  filtered_data <- x |>
    dplyr::select(all_of(variables_name_full)) |>
    sf::st_drop_geometry() |>
    dplyr::summarise(across(everything(), ~ mid_quant_75(.x))) |>
    tidyr::pivot_longer(everything(), names_to = "column", values_to = "q75_value") |>
    dplyr::filter(q75_value != 0) |>
    dplyr::pull(column)


  return(filtered_data)
}


select_lowest_lag <- function(var1, var2, group1, group2) {
    if(group1 == group2 && group1 %in% c("max_cyc", "max_dhw")) {
      lag1 <- ifelse(str_detect(var1, "lag"), as.numeric(str_extract(var1, "\\d+")), 0)
      lag2 <- ifelse(str_detect(var2, "lag"), as.numeric(str_extract(var2, "\\d+")), 0)
      if (lag1 <= lag2) return(var1) else return(var2)
    } else {
      # Keep both if not in same group
      return(c(var1, var2))
    }
  }
  
filter_non_collinear <- function(df, vars, threshold = 0.7) {
  
  vars_valid <- vars[
  sapply(df[vars], function(x) {
    x <- x[complete.cases(df[vars])]
    sd(x, na.rm = TRUE) > 0
  })
]
  m_coll <- cor(df[vars_valid], use = "pairwise.complete.obs")
  
  corr_long <- as.data.frame(as.table(m_coll)) %>%
    filter(Var1 != Var2) %>%
    filter(abs(Freq) > threshold) %>%
    arrange(desc(abs(Freq))) %>%
    mutate(group = case_when(
      str_detect(Var1, "cyc") & str_detect(Var2, "cyc") ~ "max_cyc",
      str_detect(Var1, "dhw") & str_detect(Var2, "dhw") ~ "max_dhw",
      TRUE ~ "other"
    ))

  if(nrow(corr_long) == 0) return( vars_valid)
  
  corr_long <- corr_long %>%
    rowwise() %>%
    mutate(keep_var = list(select_lowest_lag(as.character(Var1),
                                             as.character(Var2),
                                             group,
                                             group))) %>%
    ungroup() %>%
    unnest(cols = c(keep_var))
  
  corr_long <- corr_long %>%
    mutate(keep_var = as.character(keep_var),
           Var1 = as.character(Var1),
           Var2 = as.character(Var2))
  
  non_flagged_vars <- setdiff(vars, unique(c(corr_long$Var1, corr_long$Var2)))
  
  final_vars <- unique(c(corr_long$keep_var, non_flagged_vars))
  
  return(final_vars)
}

make_reefid <- function(tier.sf.joined, HexPred_sf, reef_layer.sf) {
  sf::sf_use_s2(TRUE) |> suppressMessages()

  covs.hexpred_tier_sf <- HexPred_sf |>
    dplyr::left_join(tier.sf.joined, by = c("Tier5" = "Tier5")) |>
    dplyr::filter(fYEAR == min(fYEAR)) |>
    droplevels() |>
    sf::st_as_sf() |>
    sf::st_cast("POLYGON") |>
    sf::st_transform(crs = sf::st_crs(reef_layer.sf)) |>
    suppressMessages() |>
    suppressWarnings()

  # Check CRS units
  testthat::expect_equal(sf::st_crs(covs.hexpred_tier_sf)$units, "m")
  testthat::expect_equal(sf::st_crs(reef_layer.sf)$units, "m")

  Reef_layer_tier5_84 <- reef_layer.sf |>
    sf::st_crop(covs.hexpred_tier_sf) |>
    sf::st_cast("POLYGON") |>
    dplyr::mutate(reefid = dplyr::row_number()) |>
    dplyr::select(-GRIDCODE) |>
    sf::st_transform(crs = 4326) |>
    sf::st_buffer(dist = 450) |> # careful with units and version of sf
    suppressMessages() |>
    suppressWarnings()

  covs.hexpred_tier_sf_84 <- covs.hexpred_tier_sf |>
    sf::st_transform(crs = 4326)

  # Join the shapefiles
  sf::sf_use_s2(FALSE) |> suppressMessages()

  covs.hexpred_tier_sf_v2_prep <- covs.hexpred_tier_sf_84 |>
    sf::st_join(Reef_layer_tier5_84) |>
    dplyr::select(Tier5, reefid, geometry) |>
    suppressMessages() |>
    suppressWarnings()

  return(covs.hexpred_tier_sf_v2_prep)
}

rm_obs_outside <- function(data.grp.tier, HexPred_reefid2) {

  data.grp.tier.sf <- data.grp.tier |>
    sf::st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

  within_check <- sf::st_within(data.grp.tier.sf, HexPred_reefid2)

  inside_indices <- which(lengths(within_check) > 0)

  data.grp.tier.filtered <- data.grp.tier.sf[inside_indices, ] %>%
    dplyr::mutate(
      LONGITUDE = sf::st_coordinates(.)[, 1],
      LATITUDE  = sf::st_coordinates(.)[, 2]
    ) %>%
    sf::st_drop_geometry()
  return(data.grp.tier.filtered)
}

frk_prep_3L <- function(data.grp.tier, HexPred_reefid2) {

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

   ## Prep FRK model inputs - 3 resolutions 
frk_prep_2L <- function(data.grp.tier, HexPred_reefid2) {

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
                          nres = 2L, #for dev
                          #nres = 3L,
                          regular = TRUE)
  obj_frk <- list("ST_BAUs" = ST_BAUs, "STObj" = STObj, "basis" = basis)
  return(obj_frk)
}

#############################################
######### Plot coral cover observations
#############################################

data_viz <- function(dat, title_n){

  dat_hov <- dat %>% 
    arrange(desc(LATITUDE), LONGITUDE) %>%
    group_by(fYEAR, Tier5) %>% 
    summarize(mean_cover = mean(COUNT/TOTAL) * 100, .groups = "drop") %>%
    arrange(Tier5, fYEAR)

  pal_hov <- viridisLite::plasma(100, begin = 0, end = .70)

  hov_plot <- ggplot(dat_hov,
                     aes(x = as.numeric(as.character(fYEAR)),
                         y = Tier5,
                         fill = mean_cover)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_gradientn(
      colours = pal_hov,
      name = "Coral cover (%)",
      limits = c(0, 70)
    ) +
    theme_minimal() +
    labs(x = "Year", y = "Data-tiers") +
    ggtitle(sub("_.*", "", title_n)) +
    theme(
      plot.title = element_text(size = 13, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 13),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 13),
      axis.title = element_text(size = 15),
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(0.5, "cm")
    ) +
    scale_x_continuous(breaks = seq(2015, 2023, by = 2))

  return(hov_plot)
}


#############################################
######### Plot disturbances
#############################################

summarise_disturbance <- function(dat, var, source_label) {
  
  n_tier5 <- n_distinct(dat$Tier5)
  
  dat %>%
    group_by(fYEAR) %>%
    summarize(
      mean_val  = mean(.data[[var]], na.rm = TRUE),
      sd_val    = sd(.data[[var]], na.rm = TRUE),
      n_nonzero = sum(.data[[var]] != 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      se_val = sd_val / sqrt(n_tier5),
      prop   = n_nonzero / n_tier5,
      Source = source_label
    )
}

# Generic plot function
make_disturbance_plot <- function(mod, var, dist_label, fill_name, palette_name, direction = 1) {
  
  dat_coral <- mod$data.grp.tier
  hexpred   <- mod$HexPred_reefid2
  
  full_cov_plots <- hexpred %>%
    inner_join(
      tiers.lookup %>% dplyr::select(-reef_area, -tier_id),
      by = "Tier5"
    ) 
  
  # All tiers
  dat_all <- summarise_disturbance(full_cov_plots, var, "all-tier")
  
  # Data tiers only
  dat_tier_year <- dat_coral %>%
    mutate(tier_year = paste(Tier5, fYEAR, sep = "_")) %>%
    pull(tier_year)
  
  dat_data <- full_cov_plots %>%
    mutate(tier_year = paste(Tier5, fYEAR, sep = "_")) %>%
    filter(tier_year %in% dat_tier_year) %>%
    summarise_disturbance(var, "data-tier")
  
  # Combine and prepare for plotting
  dat_plot <- bind_rows(dat_all, dat_data) %>%
    filter(prop != 0) %>%
    dplyr::select(fYEAR, mean_val, prop, Source) %>%
    group_by(fYEAR, Source) %>%
    group_split() %>%
    map_dfr(~ add_row(
      .x,
      fYEAR = first(.x$fYEAR),
      mean_val = NA,
      prop = 1 - sum(.x$prop),
      Source = first(.x$Source)
    )) %>%
    filter(as.numeric(as.character(fYEAR)) >= 2010)
  
  # Plot
  ggplot(dat_plot, aes(x = Source, y = prop, fill = mean_val)) +
    geom_col(width = 0.8, colour = "white") +
    facet_wrap(~fYEAR, nrow = 3) +
    scico::scale_fill_scico(
      name = fill_name,
      palette = palette_name,
      na.value = "gray98",
      direction = direction
    ) +
    scale_y_continuous(
      labels = scales::percent_format(),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = NULL,
      y = paste0("Proportion of tiers with ", dist_label, " > 0")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(size = 11),
      legend.position = "bottom",
      panel.spacing = unit(0.9, "lines")
    )
}

#############################################
######### Plot islands
#############################################

plot_islands_custom_buffer <- function(add_islands) {
  
  # # Compute bounding boxes with individual buffers
  island_bboxes <- add_islands %>%
    rowwise() %>%
    mutate(
      xmin = lon - buffer_deg,
      xmax = lon + buffer_deg,
      ymin = lat - buffer_deg,
      ymax = lat + buffer_deg
    ) %>%
    select(NAME, xmin, xmax, ymin, ymax, buffer_deg)
  
  island_bboxes_sf <- add_islands %>%
  rowwise() %>%
  mutate(
    xmin = lon - buffer_deg,
    xmax = lon + buffer_deg,
    ymin = lat - buffer_deg,
    ymax = lat + buffer_deg,
    geometry = st_as_sfc(st_bbox(c(
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax
    ), crs = 4326))
  ) %>%
  st_as_sf() %>%
  select(NAME)

  unique_hex_island <- st_join(
  unique_hex,
  island_bboxes_sf,
  join = st_intersects,
  left = FALSE   # keep only hexes inside islands
   )

  # Generate plots
  plots <- lapply(1:nrow(island_bboxes), function(i) {
    
    island <- island_bboxes[i, ]
    
    ggplot() +
      geom_sf(data = asm_crop,
              fill = "grey95",
              colour = "black",
              linewidth = 0.6) +
      geom_sf(data = tier.sf.joined,
              fill = NA,
              colour = "black",
              alpha = 0.2) +
      geom_sf(data = unique_hex_island,
              aes(fill = NAME),
              alpha = 0.8,
              linewidth = 0.05) +
      geom_point(data = tier4_south_showed,
                 aes(x = lon, y = lat),
                 shape = 21, fill = "red", colour = "black", size = 1.5, alpha = 0.6) +
      coord_sf(xlim = c(island$xmin, island$xmax),
               ylim = c(island$ymin, island$ymax)) +
      ggpubr::theme_pubr() +
      labs(
        title = island$NAME,
        x = "Longitude",
        y = "Latitude"
      ) +
      annotation_scale(
        location = "bl",
        width_hint = 0.3
      ) +
      theme(
        axis.text  = element_text(size = 10, angle = 45, hjust = 1),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 16),
        legend.position = "none"
      ) +
      scale_fill_manual(values = c(
  "Tutuila"       = "#5b2a86",  
  "Ofu-Olosega"   = "#2c7fb8",  
  "Ta‘ū"          = "#41b6c4",  
  "Rose Atoll"    = "#7fcdbb",  
  "Swains Island" = "#f768a1"
))
  })
  
  names(plots) <- island_bboxes$NAME
  return(list(plots = plots, island_bboxes_sf = island_bboxes_sf, unique_hex_island = unique_hex_island))
}

#############################################
######### Spatial weighted aggregation
#############################################

make_contrasts <- function(pred_tierIndex, tier_col) {
    cellmeans_wide_list <- pred_tierIndex |>
      group_split(!!sym(tier_col)) |>
      map(~ .x |> 
        pivot_wider(
          names_from = fYEAR,
          values_from = cover_prop
        )
     )
  
  # Now process contrasts on each list
  predictions_list <- map(cellmeans_wide_list, ~ process_contrasts(.x, tier_col = tier_col))
  return(predictions_list)
}

process_contrasts <- function(cellmeans_wide, tier_col) {
  predictions_i <- cellmeans_wide |>
    mutate(iter = seq_len(n())) |>
    pivot_longer(cols = contains("20"), names_to = "year") |>
    mutate(year = as.integer(year)) |>
    arrange(iter, year) |>
    group_by(iter) |>
    mutate(diff = value / lag(value, n = 1),
           diff_id = factor(paste0("diff_", seq_len(n())))) |>
    ungroup()

  max_iter <- max(predictions_i$iter)

  plot_data <- predictions_i |>
    drop_na() |>
    mutate(cat_up = ifelse(diff > 1, 1, 0),
           cat_down = ifelse(diff < 1, 1, 0)) |>
    group_by(diff_id, year) |>
    summarise(prob_up = sum(cat_up) / max_iter,
              prob_down = sum(cat_down) / max_iter,
              .groups = "drop")

  fold_change <- predictions_i |>
    drop_na() |>
    group_by(year, diff_id) |>
    summarise(fold_change = mean(diff), .groups = "drop") |>
    left_join(plot_data, by = join_by(year, diff_id))

  direction_arrow <- fold_change |>
    group_by(year) |>
    mutate(arrow = case_when(
      fold_change > 1 & prob_up >= 0.9 ~ "Up",
      fold_change < 1 & prob_down >= 0.9 ~ "Down",
      TRUE ~ "Flat"
    )) |>
    select(year, fold_change, prob_up, prob_down, arrow)

  predictions_i |>
    filter(!is.na(value)) |>
    group_by(year, !!sym(tier_col), model_name) |>
    ggdist::median_hdci(value) |>
    select(-.width, -.interval) |>
    left_join(direction_arrow, by = "year")
}

get_sum_area <- function(post_dist_df_all, tier_col, group = NULL) {
  if (is.null(group)) {
    sum_df <- post_dist_df_all  |> 
      dplyr::group_by(Tier5) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::group_by(!!sym(tier_col)) |>
      dplyr::summarise(sum_area = sum(reef_area), .groups = "drop")
  } else {
    sum_df <- post_dist_df_all |>
      dplyr::filter(tier_type == group) |> 
      dplyr::group_by(Tier5) |>
      dplyr::slice_head(n = 1) |>
      dplyr::group_by(!!sym(tier_col)) |>
      dplyr::summarise(sum_area = sum(reef_area), .groups = "drop")
  }

  return(sum_df)
}

extract_info_region <- function(post_dist_df_all, tier_col) {
## Extract sum_area 

sum_area <- get_sum_area(post_dist_df_all, tier_col) %>%
  arrange(!!sym(tier_col))

## Extract year range 

year_range <- post_dist_df_all %>% group_by(!!sym(tier_col)) %>%
  mutate(
    fYEAR_numeric = as.numeric(as.character(fYEAR))) %>%
  summarize(
    year_range = paste0(min(fYEAR_numeric, na.rm = TRUE), "–", max(fYEAR_numeric, na.rm = TRUE))
  ) %>%
  dplyr::select(!!sym(tier_col), year_range) %>%
  arrange(!!sym(tier_col))


# Extract % of data tier 
data_perc <- post_dist_df_all %>% 
  group_by(!!sym(tier_col)) %>%
  dplyr::count(tier_type) %>%
  dplyr::mutate(prop = (n / sum(n))*100) %>%
  dplyr::select(!!sym(tier_col), tier_type, prop) %>%
  pivot_wider(names_from = tier_type, values_from = prop, values_fill = 0) %>%
  arrange(!!sym(tier_col))

# Extract % of model types 
model_perc <- post_dist_df_all %>% 
  group_by(!!sym(tier_col)) %>%
  dplyr::count(model_name) %>%
  dplyr::mutate(prop = (n / sum(n))*100) %>%
  dplyr::select(!!sym(tier_col), model_name, prop) %>%
  pivot_wider(names_from = model_name, values_from = prop, values_fill = 0) %>%
  arrange(!!sym(tier_col))


# List of tables to join
tables <- list(sum_area, year_range, data_perc, model_perc)

# Iteratively left join all by `tier_col`
all_info <- reduce(tables, left_join, by = tier_col)  %>%
            dplyr::rename(
              Size.area = sum_area, Year.range = year_range, data.tier = data,
              new.tier = new, FRK.prop = FRK
            ) 
return(all_info)
}

spatial_agg_prep <- function() {

        files <- list.files(
        path = "../data",
        pattern = "FRK", 
        full.names = TRUE
      )
      files <- files[!grepl('TIER', files, perl = TRUE)]

      # ---- Load and prepare all model outputs ----
      data.list <- vector('list', length(files))
      post_dist_df_list <- list()
      data_tier_list <- list()

      for (i in seq_along(data.list)) {
        GROUP <- "HARD CORAL"
        tier <- stringr::str_extract(files[i], "(?<=FRK_).*?(?=\\.RData)")
        obj <- readRDS(files[i]) 

        # Temporary renaming for compatibility
        if ("data.sub" %in% names(obj)) {
          names(obj)[names(obj) == "data.sub"] <- "data.grp.tier"
        }

        post_dist_df_list[[i]] <- obj$post_dist_df %>%
        mutate(depth = tier)

        data_tier_list[[i]] <- unique(obj$data.grp.tier$Tier5) 
      }

      # Add tier type variable
      for (i in seq_along(post_dist_df_list)) {
        post_dist_df_list[[i]] <- post_dist_df_list[[i]] |>
          dplyr::mutate(
            tier_type = ifelse(
              as.character(Tier5) %in% data_tier_list[[i]],
              "data", "new"
            )
          )
      }

      rm(data.list)

      # Keep only valid elements
      post_dist_df_list <- post_dist_df_list |> purrr::keep(~ "model_name" %in% names(.x))

      # Format data types and structure
      post_dist_df_list <- map(post_dist_df_list, ~ .x |>
        dplyr::mutate(
          fYEAR = as.factor(fYEAR),
          Tier5 = as.factor(Tier5),
          id_loc = as.integer(id_loc),
          draw = as.character(draw),
          pred = as.numeric(pred),
          model_name = as.character(model_name),
          tier_type = as.character(tier_type)
        ) |>
        dplyr::select(fYEAR, Tier5, id_loc, draw, pred, model_name, tier_type, depth)
      )
    
      # Bind and weight all tiers
      post_dist_df_all <- dplyr::bind_rows(post_dist_df_list) %>%
        dplyr::left_join(tiers.lookup) %>%
        dplyr::mutate(
          reef_area = reef_area / 1000000,
          weighted_pred = pred * reef_area
        )

    # Bind all Tier5 rows
      post_dist_df_tier5 <- dplyr::bind_rows(post_dist_df_list) %>%
        dplyr::left_join(tiers.lookup) %>%
        dplyr::mutate(
          reef_area = reef_area / 1000000
        )

      rm(post_dist_df_list)

      # Drop NA rows (likely wrongly saved)
      post_dist_df_tier5 <- post_dist_df_tier5 %>% dplyr::filter(if_all(everything(), ~ !is.na(.)))
      post_dist_df_all <- post_dist_df_all %>% dplyr::filter(if_all(everything(), ~ !is.na(.)))

      return(list(post_dist_df_tier5 = post_dist_df_tier5, post_dist_df_all = post_dist_df_all))
    }

spatial_agg <- function(post_dist_df_tier5, post_dist_df_all, tiers.lookup, tier.sf, tier_col) {

        if (tier_col == "Tier5") {
        
        # Local scale
      
        pred_tierIndex <- post_dist_df_tier5 %>%
            dplyr::rename(cover_prop = pred) %>%
            dplyr::select(fYEAR, !!sym(tier_col), draw, model_name, cover_prop)

        predictions <- make_contrasts(pred_tierIndex, tier_col) 

        pred_tierIndex <- dplyr::bind_rows(predictions) %>%
            dplyr::rename(
              Median = value, Lower = .lower, Upper = .upper,
              Fold.Change = fold_change, P.up = prob_up,
              P.down = prob_down, Change = arrow, Model.name = model_name,
              Year = year
            ) %>%
            dplyr::select(
              !!sym(tier_col), Year, Median, Lower, Upper, Fold.Change, P.up, P.down, Change, Model.name
            )


        } else if (tier_col == "Tier4") {
          
          # ecoregion scale
          sum_area <- get_sum_area(post_dist_df_all, tier_col)

          pred_tierIndex <- post_dist_df_all %>%
            dplyr::group_by(fYEAR, draw, !!sym(tier_col), model_name) %>%
            dplyr::summarise(
              cover = sum(weighted_pred, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            dplyr::left_join(sum_area, by = tier_col) %>%
            dplyr::mutate(
              cover_prop = cover / sum_area
            ) %>%
            dplyr::select(fYEAR, !!sym(tier_col), draw, model_name, cover_prop)

        predictions <- make_contrasts(pred_tierIndex, tier_col)

        pred_tierIndex <- dplyr::bind_rows(predictions) %>%
            dplyr::mutate(
              Median = value,
              Lower  = .lower,
              Upper  = .upper
             ) %>%
           dplyr::rename(
             Fold.Change = fold_change,
             P.up        = prob_up,
             P.down      = prob_down,
            Change      = arrow,
            Model.name  = model_name,
            Year        = year
      ) %>%
      dplyr::select(
      !!sym(tier_col), Year, Median, Lower, Upper,
      Fold.Change, P.up, P.down, Change, Model.name
     )
      } else {
     stop("tier_col must be either 'Tier5' or 'Tier4'")
    
     }

      # Extra information of the regions

         info_region <- extract_info_region(post_dist_df_all, tier_col)

        # ---- Save results ----
        readr::write_csv(
          pred_tierIndex,
          file = paste0("../data/output_", tier_col, ".csv"),
          quote = "none"
        )

             # ---- Save results ----
        readr::write_csv(
          info_region,
          file = paste0("../data/info_", tier_col, ".csv"),
          quote = "none"
        )
    
     return(list(info_region = info_region, pred_table = pred_tierIndex)) 

}


spatial_agg_island <- function(post_dist_df_all, tiers.lookup, tier.sf, unique_hex_island) {

        tier_col <- "NAME"
        
        # Add Island name 
        post_dist_df_all_island <- post_dist_df_all %>%
        left_join(unique_hex_island)

          
        # get reef area of each islands
        sum_area <- get_sum_area(post_dist_df_all_island, tier_col)

        pred_tierIndex <- post_dist_df_all_island %>%
            dplyr::group_by(fYEAR, draw, !!sym(tier_col), model_name) %>%
            dplyr::summarise(
              cover = sum(weighted_pred, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            dplyr::left_join(sum_area, by = tier_col) %>%
            dplyr::mutate(
              cover_prop = cover / sum_area
            ) %>%
            dplyr::select(fYEAR, !!sym(tier_col), draw, model_name, cover_prop)

        predictions <- make_contrasts(pred_tierIndex, tier_col)

        pred_tierIndex <- dplyr::bind_rows(predictions) %>%
            dplyr::mutate(
              Median = value,
              Lower  = .lower,
              Upper  = .upper
             ) %>%
           dplyr::rename(
             Fold.Change = fold_change,
             P.up        = prob_up,
             P.down      = prob_down,
            Change      = arrow,
            Model.name  = model_name,
            Year        = year
      ) %>%
      dplyr::select(
      !!sym(tier_col), Year, Median, Lower, Upper,
      Fold.Change, P.up, P.down, Change, Model.name
     )

      # Extra information of the regions

         info_region <- extract_info_region(post_dist_df_all_island, tier_col)

        # ---- Save results ----
        readr::write_csv(
          pred_tierIndex,
          file = paste0("../data/output_", tier_col, ".csv"),
          quote = "none"
        )

             # ---- Save results ----
        readr::write_csv(
          info_region,
          file = paste0("../data/info_", tier_col, ".csv"),
          quote = "none"
        )
    
     return(list(info_region = info_region, pred_table = pred_tierIndex)) 

}

#############################################
######### Attribution of disturbances 
#############################################

# Function to get scaled effect sizes
build_table <- function(vary_var, effect_col, range_vals = NULL, coef_scaled) {
  
  fixed_vals <- coef_scaled %>% 
                filter(param != vary_var) %>%
                dplyr::select(param, !!sym(effect_col)) %>%
                deframe()
  
  map_dfr(range_vals, function(x) {
    tibble(
      !!!as.list(fixed_vals),
      !!vary_var := x
    ) #%>%
     # relocate(!!sym(vary_var), .after = `(Intercept)`)
  })
}

cov_extract <- function(model.out) {
  coef_uncertainty(
    model.out,
    percentiles = c(2.5, 50, 97.5),
    nsim = 400,
    random_effects = TRUE
  ) %>%
    data.frame() %>%
    tibble::rownames_to_column("rowname") %>%
    tidyr::pivot_longer(
      cols = -rowname,
      names_to = "param",
      values_to = "value"
    ) %>%
    mutate(
      Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")
    ) %>%
    tidyr::pivot_wider(
      names_from = rowname,
      values_from = value
    )
}

get_coef_data <- function(mod) {

  mean_cov <- Matrix::colMeans(mod$M@X)

  coef_ <- cov_extract(mod$M) %>%
    filter(Type == "fixed") %>%
    mutate(param = ifelse(param == "Intercept", "(Intercept)", param))
  
  coef_ %>%
    mutate(
      mean_val   = mean_cov[param],
      effect_50  = `50%`  * mean_val,
      effect_low = `2.5%` * mean_val,
      effect_high= `97.5%`* mean_val
    ) %>%
    dplyr::select(param, effect_50, effect_low, effect_high)
}

# Function to relabel parameters for plotting
label_params <- function(df) {
  df %>%
    filter(param != "(Intercept)") %>%
    mutate(
      param = dplyr::case_when(
        param == "max_cyc"       ~ "Cyclone Exposure",
        param == "max_cyc.lag1"  ~ "Cyclone Exposure (lag 1)",
        param == "max_cyc.lag2"  ~ "Cyclone Exposure (lag 2)",
        param == "max_dhw"       ~ "Heat Stress",
        param == "max_dhw.lag1"  ~ "Heat Stress (lag 1)",
        param == "max_dhw.lag2"  ~ "Heat Stress (lag 2)",
        TRUE ~ param
      )
    )
}


compute_decline <- function(cond_table){

  # name of the first column (predictor)
  xvar <- names(cond_table)[1]

  # row with min and max disturbance
  row_min <- cond_table %>% dplyr::slice_min(.data[[xvar]], n = 1)
  row_max <- cond_table %>% dplyr::slice_max(.data[[xvar]], n = 1)

  tibble::tibble(
    mean_decline = (row_max$value_coral_50 - row_min$value_coral_50), #/ row_max$value_coral_50,
    low_decline  = (row_max$value_coral_low - row_min$value_coral_low), #/ row_max$value_coral_low,
    high_decline = (row_max$value_coral_high - row_min$value_coral_high)# / row_max$value_coral_high
  )
}

predict_cover <- function(effect_col, label, coef_scaled, range_vals) {
  
  coef_value <- coef_scaled %>%
    filter(param == label) %>%
    pull(all_of(effect_col))
  
  table <- build_table(
    vary_var   = label,
    effect_col = effect_col,
    range_vals = range_vals,
    coef_scaled = coef_scaled   
  ) %>%
    mutate(new_var = .data[[label]] * coef_value) %>%
    dplyr::select(!label)
  
  plogis(rowSums(table))
}

plot_cond <- function(range_vals, label, coef_scaled, xlabel, ylabel){
  
  cond_table <- tibble(
    !!sym(label) := range_vals,
    value_coral_50  = predict_cover("effect_50", label, coef_scaled, range_vals),
    value_coral_low = predict_cover("effect_low", label, coef_scaled, range_vals),
    value_coral_high = predict_cover("effect_high", label, coef_scaled, range_vals)
  )

  xlabel_final <- if (str_detect(xlabel, "Cyclone")) {
    paste0(xlabel, ", in hours")
     } else {
    paste0(xlabel, ", in DHW")
   }
  
  p <- ggplot(cond_table, aes(x = !!sym(label))) +
    geom_ribbon(aes(ymin = value_coral_low * 100, ymax = value_coral_high * 100), 
                fill = "gray66", alpha = 0.3) +
    geom_smooth(aes(y = value_coral_50 * 100), color = "black", size = 1, se = FALSE) +
    labs(
      x = xlabel_final,
      y = ylabel
    ) +
    theme_pubr() +
    scale_x_continuous(
      limits = c(0, max(range_vals)),
      breaks = seq(0, max(range_vals), by = 5)
    ) +
    ylim(0, round(max(cond_table$value_coral_high*100),1) + 1)
  
  return(list(p = p, cond_table = cond_table))
}

# ---- plot_conditional_time ----
plot_conditional_time <- function(mod_obj, mod_name, save_dir = "../figures") {
  
  # Extract coefficients
  coef_scaled <- get_coef_data(mod_obj) %>%
    label_params() 

  # Extract prediction table
  hexpred <- mod_obj$HexPred_reefid2
  
  # List of variables to plot with labels
  plot_vars <- list(
    max_dhw = "Heat Stress",
   # max_dhw.lag1 = "Heat Stress (lag 1)",
  #  max_dhw.lag2 = "Heat Stress (lag 2)",
    max_cyc = "Cyclone Exposure"
  #  max_cyc.lag1 = "Cyclone Exposure (lag 1)",
  #  max_cyc.lag2 = "Cyclone Exposure (lag 2)"
  )
  
  # Loop over variables
  plots <- purrr::imap(plot_vars, function(label, var_name) {

    range_vals <- seq(0, max(hexpred[[var_name]], na.rm = TRUE), length.out = 100)
    
    p <- plot_cond(range_vals, label, coef_scaled, xlabel = label, ylabel = "Coral Cover (%)")$p
    cond_table <- plot_cond(range_vals, label, coef_scaled, xlabel = label, ylabel = "Coral Cover (%)")$cond_table
    # Save plot
    file_name <- file.path("../figures/", paste0(mod_name, "_", var_name, ".png"))
    ggsave(file_name, p, width = 6, height = 4)
    
    return(list(p =p, cond_table = cond_table))
  })  
  return(plots)
}

predict_newdata <- function(mod_M, new_BAU_data, hexpred, obj_frk, title_n, max_cover){
  
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
  
  pal_pred <- LaCroixColoR::lacroix_palette("Pamplemousse", n = 30, type = "continuous")

  p <- ggplot() + 
    geom_sf(data = pred_sum_sf, aes(fill = pred*100), col = "transparent") +
    coord_sf(expand = F) +
    scale_fill_gradientn(
      colours = pal_pred,    
      name = "Predicted coral cover (%)",
      limits = c(0,max_cover), 
    ) + 
    theme_pubr() +
    xlab("Longitude") +
    ylab("Latitude") +
    theme(
      legend.position = "bottom",
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
      annotation_north_arrow(
            location = "tr",         
            which_north = "true",
            style = north_arrow_fancy_orienteering,
            pad_x = unit(0.1, "cm"),
            pad_y = unit(.02, "cm")
     ) + annotation_scale(
           location = "bl",        
           width_hint = 0.3
  )

# p_insert <- ggplot(pred_sum_sf, aes(x = pred * 100, y = "", fill = ..x..)) +
#   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
#   scale_fill_gradientn(
#     colours = pal_pred,
#     name = "Predicted coral cover (%)",
#     limits = c(0, 70)
#   ) +
#   labs(x = NULL, y = NULL) +
#   theme_pubr() +
#   theme(
#     legend.position = "none",
#     axis.title = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.line = element_blank()
#   )

# p_all <- p +
#   inset_element(p_insert, 0.6, 0.6, 1, 1) 

  return(list(p, pred_sum_sf)) # p_insert,  
}

predict_year_island <- function(mod_M, new_BAU_data, hexpred, obj_frk, 
                                unique_hex_island, title_n, max_cover){ 
  
  newpred <- predict(mod_M, new_BAU_data = new_BAU_data)
  
  post_dist_df <- as.data.frame(newpred$MC$mu_samples) %>% 
    mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
    mutate(Tier5 = obj_frk$ST_BAUs@data$Tier5) %>%
    mutate(id_loc = row_number()) %>%
    tidyr::pivot_longer(!c(fYEAR, Tier5, id_loc),
                        names_to = "draw", 
                        values_to = "pred")    
  
  # Spatial maps 
  hexpred$Tier5 <- as.factor(hexpred$Tier5)
  
  pred_sum_sf <- post_dist_df %>% 
    group_by(fYEAR, Tier5) %>% 
    median_hdci(pred) %>%
    inner_join(hexpred %>% 
                 group_by(Tier5) %>% 
                 summarize() %>% 
                 dplyr::select(geometry, Tier5)) %>% 
    st_as_sf(sf_column_name = "geometry") %>%
    mutate(Unc = .upper - .lower,
           tier_fYEAR = paste0(Tier5, fYEAR)) %>%
    dplyr::select(fYEAR, Tier5, pred, .lower, .upper, Unc, tier_fYEAR) %>%
  left_join(
    unique_hex_island %>% st_drop_geometry(),
    by = "Tier5"
  )
  
  # palettes
  pal_pred <- LaCroixColoR::lacroix_palette("Pamplemousse", n = 100, type = "continuous")
  pal_unc  <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
  
  # ---- split by island ----
  plots <- pred_sum_sf %>%
    split(.$NAME) %>%
    purrr::imap(function(df_island, island_name) {
      
      # prediction plot
      p <- ggplot() + 
        geom_sf(data = df_island, aes(fill = pred * 100), col = "transparent") +
        coord_sf(expand = FALSE) +
        scale_fill_gradientn(
          colours = pal_pred,
          name = "Coral Cover (%)",
          limits = c(0, max_cover)
        ) +
        facet_wrap(~ fYEAR) +
        theme_pubr() +
        labs(
          title = paste0(title_n, " - ", island_name),
          x = "Longitude",
          y = "Latitude"
        ) +
        theme(
          legend.position = "top",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 13, hjust = 0.5),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12)
        )
      
      # uncertainty plot
      p_unc <- ggplot() + 
        geom_sf(data = df_island, aes(fill = Unc * 100), col = "transparent") +
        coord_sf(expand = FALSE) +
        scale_fill_gradientn(
          colours = pal_unc,
          name = "Uncertainty range (%)",
          limits = c(0, 100)
        ) +
        facet_wrap(~ fYEAR) +
        theme_pubr() +
        labs(
          title = paste0(title_n, " - ", island_name),
          x = "Longitude",
          y = "Latitude"
        ) +
        theme(
          legend.position = "top",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 13, hjust = 0.5),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12)
        )
      
      list(p = p, p_unc = p_unc)
    })
                                }

plot_diff <- function(pred_high, pred_low,  title_n){
diff_pred <- pred_high %>%
  st_drop_geometry() %>%                # remove geometry before join
  select(Tier5, pred_high = pred) %>%
  left_join(
    pred_low %>%
      st_drop_geometry() %>%
      select(Tier5, pred_low = pred),
    by = "Tier5"
  ) %>%
  mutate(abs_diff = pred_high- pred_low,
         rel_diff = (pred_high - pred_low) / ((pred_high + pred_low)/2))

diff_sf <- pred_high %>%
  select(Tier5, geometry) %>%
  left_join(diff_pred, by = "Tier5")

#pal_pred <- LaCroixColoR::lacroix_palette("Pamplemousse", n = 100, type = "continuous")

p <- ggplot() + 
  geom_sf(data = diff_sf, 
          aes(fill = abs_diff * 100), 
          col = "transparent") +
  
  coord_sf(expand = FALSE) +
  
  scico::scale_fill_scico(
  palette = "vik",
  midpoint = -5,
  name = "Absolute change (%)"
  ) +
  
  theme_pubr() +
  xlab("Longitude") +
  ylab("Latitude") +

  ggtitle(title_n) +

  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  ) +
  
  annotation_north_arrow(
    location = "tr",         
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    pad_x = unit(0.1, "cm"),
    pad_y = unit(.02, "cm")
  ) +
  
  annotation_scale(
    location = "bl",        
    width_hint = 0.3
  )
return(p)
}

#############################################
######### Model validation 
#############################################

compute_model_fit <- function(mod, model_name) {

mean_obs <- mod$data.grp.tier  %>%
 group_by(Tier5, fYEAR) %>%
 summarize(COUNT_sum = sum(COUNT),
           TOTAL_sum = sum(TOTAL)) %>%
 ungroup() %>% 
    group_by(fYEAR, Tier5) %>%
    summarize(
      mean_cover = mean(COUNT_sum / TOTAL_sum) * 100,
      .groups = "drop"
    ) %>%
    mutate(
      fYEAR = as.factor(fYEAR),
      Tier5 = as.factor(Tier5)
    )

  mean_pred <- mod$pred_sum_sf %>%
    mutate(
      pred  = pred * 100,
      fYEAR = as.factor(fYEAR),
      Tier5 = as.factor(Tier5)
    ) %>%
    dplyr::select(fYEAR, Tier5, pred)

  mean_all <- inner_join(mean_obs, mean_pred, by = c("fYEAR", "Tier5"))

  r2 <- cor(mean_all$mean_cover, mean_all$pred)^2
  r2_label <- paste0("R² = ", round(r2, 2))

  plot <- ggplot(data = mean_all) + 
    geom_point(
    aes(
    x = mean_cover,
    y = pred),
    shape = 16,
    alpha = 0.4,
    size = 2.5
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "red",
      linetype = "dashed",
      size = 1
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = r2_label,
      hjust = 1.1,
      vjust = -1.5,
      size = 5
    ) +
    labs(
      x = "Observed Mean Coral Cover (%)",
      y = "Predicted Coral Cover (%)",
      title = str_extract(model_name, "^[^_]+")
    )  + 
  theme_pubr()

  list(
    data = mean_all,
    r2   = r2,
    plot = plot
  )
}


make_FRK_dharma_res <- function(mod) {
  # Observed response
  response_df <- mod$data.grp.tier %>% 
   # group_by(fYEAR, Tier5) %>% 
   # summarize(COUNT = sum(COUNT), k = sum(TOTAL))  %>% 
    mutate(fYEAR = as.factor(fYEAR), Tier5 = as.factor(Tier5))

  # Simulated response 
  pred <- predict(mod$M, type = c("mean"))
  pred_sum <-  mod$pred_sum_sf %>% 
    mutate(fYEAR = as.factor(mod$obj_frk$ST_BAUs@data$fYEAR)) %>%
    mutate(Tier5 = as.factor(mod$obj_frk$ST_BAUs@data$Tier5)) %>%
    semi_join(response_df, by = c("fYEAR", "Tier5")) 

  simulation_df <- as.data.frame(pred$MC$mu_samples) %>% 
    mutate(fYEAR = as.factor(mod$obj_frk$ST_BAUs@data$fYEAR)) %>%
    mutate(Tier5 = as.factor(mod$obj_frk$ST_BAUs@data$Tier5)) %>%
    semi_join(response_df, by = c("fYEAR", "Tier5")) 
  
  # Make sure the order is preserved 
  response_df <- simulation_df %>%
    dplyr::select(fYEAR, Tier5) %>%
    dplyr::left_join(response_df, by = c("fYEAR", "Tier5"))

  observedResponse <- response_df$COUNT
  simulation_df <- as.matrix(simulation_df[, 1:400]) 
 
  # Flatten sampled probabilities, number of trials; simulate from binomial distribution; then reshape back to correct dimensions
  probs <- as.vector(simulation_df)
  trials <- rep(response_df$TOTAL, times = ncol(simulation_df))
  sims <- rbinom(length(probs), size = trials, prob = probs)
  simulatedResponse <- matrix(sims, nrow = nrow(simulation_df), ncol = ncol(simulation_df))
  
  DHARMa::createDHARMa(
    simulatedResponse = simulatedResponse,
    observedResponse = observedResponse,
    fittedPredictedResponse = pred_sum$pred,
    integerResponse = TRUE
  )

}

###############################
#### ggplot2-style DHARMa plots
# from https://github.com/open-AIMS/stats_helpers/blob/main/R/dharma.R
###############################

check_dots <- DHARMa:::checkDots
dispersion_fct <- DHARMa:::testDispersion
ensure_dharma <- DHARMa:::ensureDHARMa
ensure_predictor <- DHARMa:::ensurePredictor
get_p_val <- DHARMa:::getP
outliers_fct <- DHARMa:::testOutliers
uniformity_fct <- DHARMa:::testUniformity
test_categorical <- DHARMa:::testCategorical
test_quantiles <- DHARMa:::testQuantiles

plot_qq_unif <- function(sim_listrix, test_uniformity = FALSE,
                         test_outliers = FALSE, test_dispersion = TRUE, ...) {
  sim_listrix <- ensure_dharma(sim_listrix, convert = "Model")
  res_ <- data.frame(y = sim_listrix$scaledResiduals) %>%
    dplyr::mutate(x = seq_len(dplyr::n()) / (dplyr::n() + 1))
  res_ <- qqplot(res_$x, res_$y, plot.it = FALSE) %>%
    as.data.frame
  p_ <- ggplot(data = res_) +
    geom_point(mapping = aes(x = x, y = y), shape = 16,
               colour = "skyblue", alpha = 0.8) +
    geom_abline(slope = 1, linetype = 2) +
    labs(x = "Expected", y = "Observed") +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(plot.title = element_text(size = 9))
  if (test_uniformity) {
    tmp_ks <- uniformity_fct(sim_listrix, plot = FALSE)
    k_lab_a <- paste("KS test: p =", round(tmp_ks$p.value, digits = 5))
    k_lab_b <- paste("Deviation:", ifelse(tmp_ks$p.value < 0.05, "Significant",
                                          "N.S."))
    p_ <- p_ +
      annotate("text", x = 0, y = 0.98, hjust = 0, vjust = 0.5, label = k_lab_a,
               colour = ifelse(tmp_ks$p.value < 0.05, "tomato", "black"),
               size = 3) +
      annotate("text", x = 0, y = 0.90, hjust = 0, vjust = 0.5, label = k_lab_b,
               colour = ifelse(tmp_ks$p.value < 0.05, "tomato", "black"),
               size = 3)
  }
  if (test_outliers) {
    tmp_ot <- outliers_fct(sim_listrix, plot = FALSE)
    pval_ot <- ceiling(tmp_ot$p.value * 10) / 10
    ot_lab_a <- paste("Outlier test: p =", pval_ot)
    ot_lab_b <- paste("Deviation:", ifelse(pval_ot < 0.05, "Significant",
                                           "N.S."))
    p_ <- p_ +
      annotate("text", x = 0, y = 0.78, hjust = 0, vjust = 0.5, size = 3,
               label = ot_lab_a, colour = ifelse(pval_ot < 0.05,
                                                 "tomato", "black")) +
      annotate("text", x = 0, y = 0.70, hjust = 0, vjust = 0.5, size = 3,
               label = ot_lab_b, colour = ifelse(pval_ot < 0.05,
                                                 "tomato", "black"))
    }
  if (test_dispersion) {
    tmp_ds <- dispersion_fct(sim_listrix, alternative = "two.sided",
                             plot = FALSE)
    pval_disp <- ceiling(tmp_ds$p.value * 10) / 10
    ds_lab_a <- paste("Dispersion test: p =", pval_disp)
    ds_lab_b <- paste("Deviation:", ifelse(pval_disp < 0.05, "significant",
                                           "N.S."))
    p_ <- p_ +
      annotate("text", x = 1, y = 0.12, hjust = 1, vjust = 0.5, size = 3,
               label = ds_lab_a, colour = ifelse(pval_disp < 0.05,
                                                 "tomato", "black")) +
      annotate("text", x = 1, y = 0.04, hjust = 1, vjust = 0.5, size = 3,
               label = ds_lab_b, colour = ifelse(pval_disp < 0.05,
                                                 "tomato", "black"))
  }
  p_
}

gg_dharma <- function(x, title = title_n, ...) {
  a_ <- plot_qq_unif(x)
  a_ + patchwork::plot_annotation(title = title)
}

#############################################
######### Import figures
#############################################

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

