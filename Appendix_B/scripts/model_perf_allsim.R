#' @title Generate model predictive performance viz across all simulations
#' @author Julie Vercelloni

setwd("c:/Users/jvercell/OneDrive - Australian Institute of Marine Science/AIMS/01_Research projects/ReefCloud/SP_models/FRK_dev/git_content/Appendix_B")

rm(list=ls())

# # Loading packages 
library(dplyr)
library(readr)
library(ggvanced)

source("R/functions.R")

# Grep names of simulations 
path <- "./"
all_dirs <- list.dirs(path, full.names = FALSE, recursive = FALSE)
exclude_patterns <- c("R", "scripts", "extra", "resources", ".git")
included_dirs <- all_dirs[!sapply(all_dirs, function(dir) any(grepl(paste(exclude_patterns, collapse = "|"), dir)))] 

path_fig.list <- list()
for (i in 1:length(included_dirs)){
path_fig.list[[i]] <- paste0(path, included_dirs[i], "/model_outputs/leave_out/table_performances_tier_true.csv")
}
file_paths <- do.call(rbind, path_fig.list)

combined_table <- file_paths %>%
  map_dfr(~ read_csv(.x) ) %>%
  filter(!grepl("^\\d", Title_of_run)) %>%
  data.frame() %>%
  rename(CvgErr = Cvg)
  
combined_table$Title_of_run <- factor(combined_table$Title_of_run, levels = c("High-resolution", "Medium-resolution", 
                                                                               "Low-resolution", "Temporally-sparse" ))
# Viz model performances
scale <- 4.5 # scale is an adjustment for a 4k screen

# Define the earthy and contrasted color palette
grand_budapest_colors <- c(
  "#67AFCB", 
  "#A23B72", 
  "#F3C13A", 
  "#3B7D4E"
)

p_perf <- ggspider(combined_table %>% dplyr::select(Title_of_run, RMSPE, `CvgErr`, `CRPS`, `IS`), # , time
         axis_name_offset = 0.20,
         background_color = "gray98", 
         fill_opacity = 0.15, 
         polygon = FALSE) +
  labs(col = "") +
  scale_color_manual(values = grand_budapest_colors) +  
  theme(
    legend.position = "top", 
    plot.title = element_text(size = 5 * scale, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 4 * scale, hjust = 0.5),
    plot.caption = element_text(size = 3 * scale),
    legend.text = element_text(size = 2 * scale),
    legend.title = element_text(size = 2.5 * scale)
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
p_perf
ggsave(p_perf, filename = "model_perf_all.png",
       width=6.5, height=7)

