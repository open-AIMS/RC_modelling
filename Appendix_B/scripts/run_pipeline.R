#' @title Run the modelling pipeline
#' @author Julie Vercelloni

# ##################################
# ############ 1. Generate synthetic data
# ##################################

source(paste0("scripts/make_data_",surveys,"_synthos_custom.R"))

# ##################################
# ############ 2. Create grid
# ##################################

source("scripts/make_grid.R")

# ##################################
# ############ 3. Create predictive layer
# ##################################

source("scripts/make_predictive_layer.R")

# ##################################
# ############ 4. Random-year sampling
# ##################################

if (sparse) {
  source("scripts/year_sampling.R")
} else {
  message("Sparse sampling not selected. Skipping year sampling script.")
}

# ##################################
# ############ 5. FRK model
# ##################################

source("scripts/FRK_model.R")

# ##################################
# ############ 6. Model predictions (incl. regional level)
# ##################################

source("scripts/model_predictions.R")

# ##################################
# ############ 7.Model predictive performances 
# ##################################

source("scripts/model_checks.R")