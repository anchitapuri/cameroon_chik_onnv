# Fit geostatistical models to predict ONNV FOI and prevelance across Cameroon
# using the stochastic partial differential equation (SPDE) approach and the R-INLA package.


# --- Source functions
source(here('Functions.R'))

# Get Cameroon boundary
cameroon <- ne_countries(country = "Cameroon", returnclass = "sf")

# Population and mosquito rasts
anopheles_funestus <- rast('2010_Anopheles_funestus_CMR.tiff')
anopheles_gambiae <- rast('2010_Anopheles_gambiae_ss_CMR.tiff')
cam_pop <- rast("cmr_ppp_2020_UNadj.tif")


# Cameroon population by age
cameroon_age_2025 <- read.csv('CameroonAge2025.csv')
cameroon_age_2025 <- cameroon_age_2025 %>%
  mutate(total = M + F)

# ----- 
cameroon_data <- readRDS('sf_meta_data_with_coords_pw.rds')
nrow(cameroon_data)


# --- Compare models with difference covariates ----
model_comparison <- compare_models(
  year_intro = 1900,
  data = model_data,
  cameroon,
  anopheles_funestus,
  anopheles_gambiae, 
  positive_col = "ONNV_pos"
)


# Plot DIC and WAIC to compare models 
plot_model_comparison(model_comparison$comparison)


# Run best model (with best covariate/s)
best_model <- run_inla_model_comparision(
  year_intro = 1900,
  data = model_data,
  cameroon = cameroon,
  anopheles_funestus = anopheles_funestus,  
  anopheles_gambiae = anopheles_gambiae,
  positive_col = "ONNV_pos",
  covariates = "baseline"
)



# --- index of prediction and estimation stacks ---
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


                         
# --- MODEL FITS ---- 
plot_age_seroprevalence_model_fits(best_model$year,best_model, model_data, "ONNV_pos")


                         
# --- Save best model results ----
saveRDS(best_model, 'ONNV_INLAResults.rds', compress = "gzip")
