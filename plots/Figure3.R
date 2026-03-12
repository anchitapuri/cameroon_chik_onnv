
source(here('R/MultiSeroFunctions.R'))
source(here('R/Functions.R'))

# Read files
# Import multisero fits 
fit <- readRDS(here('Results/full_model_fits.rds'))
# Import INLA model fits 
onnv_results <- readRDS(here('Results/ONNV_INLAResults.rds'))


# Read preprocessed data
preprocessed_data <- readRDS(here('Results/preprocessed_data_full_model.rds'))
meta_data_with_labels <- read.csv(here('Results/final_meta_data_with_labels.csv'))
# shapefile with district geometries 
cam_shapefile_districts_merged <- readRDS(here('Data/cam_shapefile_districts_merged.rds'))


# Extract chains 
chains <- fit$draws(format='df')
chains_df <- as.data.frame(chains)


age_prev_model_fits <- plot_age_sex_seroprevalence_model_fits(
  year_intro = onnv_results_pop_grid$year,
  result = onnv_results_pop_grid, 
  data = onnv_results_pop_grid$data_filtered,
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = "a" # a == ONNV 
)
print(age_prev_model_fits)
