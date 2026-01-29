# Fit geostatistical models to predict ONNV FOI and prevelance across Cameroon
# using the stochastic partial differential equation (SPDE) approach and the R-INLA package.

library(here)
library(emdbook)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(matrixStats)
library(data.table)
library(dplyr)
library(scales)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(INLA)

# --- Source functions
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/Spatial Analysis/Functions.R'))

# Get Cameroon boundary
cameroon <- ne_countries(country = "Cameroon", returnclass = "sf")

# Population and mosquito rasts
anopheles_funestus <- rast('2010_Anopheles_funestus_CMR.tiff')
anopheles_gambiae <- rast('2010_Anopheles_gambiae_ss_CMR.tiff')
cam_pop <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_ppp_2020_UNadj.tif")


# Cameroon population by age
cameroon_age_2025 <- read.csv('CameroonAge2025.csv')
cameroon_age_2025 <- cameroon_age_2025 %>%
  mutate(total = M + F)

# ----- Read preprocessed data with coords 
meta_data_with_coords <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/meta_data_with_coords.rds')
nrow(meta_data_with_coords)

meta_data_with_labels <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/final_meta_data_with_labels.csv')
nrow(meta_data_with_labels)

meta_data_with_labels$Easting <- meta_data_with_coords$Easting
meta_data_with_labels$Northing <- meta_data_with_coords$Northing

model_data <- meta_data_with_labels

# --- Run INLA model for ONNV (with historic year of intro 1900)
onnv_results <- run_inla(
   year_intro = 1900,
  data = model_data,
  cameroon = cameroon,
  positive_col = "ONNV_pos")



# --- Index of prediction and estimation stacks 
index_pred_onnv <- inla.stack.index(best_model$stk.full, "pred")$data
length(index_pred_onnv)
index_est_onnv <- inla.stack.index(best_model$stk.full, "est")$data
length(index_est_onnv)


# --- Extract and plot FOI
foi_onnv <- extract_and_plot_foi(best_model, best_model$coop, pathogen_name = "ONNV")


# --- Prob of seropositive proportion 
age_mid <- c(
  2.5, 7.5, 12.5, 17.5,
  22.5, 27.5, 32.5, 37.5, 42.5, 47.5,
  52.5, 57.5, 62.5, 67.5, 72.5, 77.5,
  82.5, 87.5, 92.5, 97.5, 100
)
cameroon_age_2025$total
w_age <- cameroon_age_2025$total / sum(cameroon_age_2025$total)

sero_onnv <- plot_predicted_seroprevalence(
  foi_result = foi_onnv,
  age_mid = age_mid,
  age_weights = w_age,
  pathogen_name = "ONNV"
)

# --- Annual Infections 
age_weights <- cameroon_age_2025$total / sum(cameroon_age_2025$total)

# Population per pixel (from gridded data)
total_pop_per_pixel <- extract(cam_pop, foi_onnv$coop)

# Distribute pixel population across age groups
pop_data <- outer(total_pop_per_pixel, age_weights)





                         
# --- Model fits 
plot_age_seroprevalence_model_fits(best_model$year,best_model, model_data, "ONNV_pos")


                         
# --- Save best model results 
saveRDS(best_model, 'ONNV_INLAResults.rds', compress = "gzip")



# --- Mosquito and population proportion vs proportion positive 