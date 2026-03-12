
source(here('R/MultiSeroFunctions.R'))
source(here('R/Functions.R'))

# Read files
# Import multisero fits 
fit_full_model <- readRDS(here('Results/full_model_fits.rds'))
onnv_results <- readRDS(here('Results/ONNV_INLAResults.rds'))


# Read preprocessed data
preprocessed_data_full_model <- readRDS(here('Results/preprocessed_data_full_model.rds'))
meta_data_with_labels <- read.csv(here('Results/final_meta_data_with_labels.csv'))
# shapefile with district geometries 
cam_shapefile_districts_merged <- readRDS(here('Data/cam_shapefile_districts_merged.rds'))

fit_full_model <- readRDS(here('Results/full_model_fits.rds'))

# extract chains and parameters
chains_full <- fit_full_model$draws(format='df')
chains_df_full <- as.data.frame(chains_full)


age_prev_model_fits <- plot_age_seroprevalence_model_fits(
  result = onnv_results_pop_grid, 
  data = onnv_results_pop_grid$data_filtered,
  model_data = meta_data_with_labels,
  chains_df = chains_df_full,
  infM = preprocessed_data_full_model$data$infM,
  pathogen_col = "a" # a == ONNV 
)
print(age_prev_model_fits)


age_gender_prev_model_fits <- plot_age_seroprevalence_model_fits_by_gender(
  result = onnv_results_pop_grid, 
  data = onnv_results_pop_grid$data_filtered,
  model_data = meta_data_with_labels,
  chains_df = chains_df_full,
  infM = preprocessed_data_full_model$data$infM,
  pathogen_col = "a" # a == ONNV 
)

age_gender_prev_model_fits <- plot_age_seroprevalence_model_fits_by_gender_bars_points(
  result = onnv_results_pop_grid, 
  data = onnv_results_pop_grid$data_filtered,
  model_data = meta_data_with_labels,
  chains_df = chains_df_full,
  infM = preprocessed_data_full_model$data$infM,
  pathogen_col = "a" # a == ONNV 
)


ggsave(here("Results/Fig3.png"), 
       plot = age_gender_prev_model_fits[[1]],    # swap for your actual plot object name
       width = 18, 
       height = 12, 
       units = "in", 
       dpi = 300,
       bg = "white")
