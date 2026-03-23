# Fit geostatistical models to predict ONNV FOI and prevelance across Cameroon
# using the stochastic partial differential equation (SPDE) approach and the R-INLA package.
library(ggraph)
library(igraph)
library(lhs)
library(matrixStats)
library(mvtnorm)
library(matrixcalc)
library(here)
library(mixR)
library(raster)
library(INLA)
library(PBSmapping)  # For convUL function
library(sp)
library(gstat)
library(sf)
library(rnaturalearth)
library(viridis)
library(patchwork)  # for combining plots
library(spatstat)
library(dplyr)
library(tidyr)
library(rworldmap)
library(exactextractr)
library(rworldxtra)
library(centr)
library(terra)
library(colorspace)
library(ggspatial)

# --- Source functions
source(here('R/Functions.R'))

# Population and mosquito rasts
anopheles_funestus <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/2010_Anopheles_funestus_CMR.tiff')
anopheles_gambiae <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/2010_Anopheles_gambiae_ss_CMR.tiff')
cam_pop <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_ppp_2020_UNadj.tif")
sum(values(cam_pop), na.rm = TRUE)

# Cameroon population by age
cameroon_age_2025 <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/CameroonAge2025.csv')
cameroon_age_2025 <- cameroon_age_2025 %>%
  mutate(total = M + F)

cameroon_age_2025$total
w_age <- cameroon_age_2025$total / sum(cameroon_age_2025$total)
sum(cameroon_age_2025$total)

age_groups <- data.frame(
  age_string = cameroon_age_2025$Age,
  stringsAsFactors = FALSE
)
# Extract lower and upper bounds from strings like "0-4", "5-9", etc.
age_groups$age_lower <- as.numeric(sub("-.*", "", age_groups$age_string))
age_groups$age_upper <- as.numeric(sub(".*-", "", age_groups$age_string))

# Handle the last age group if it's something like "80+" 
if (grepl("\\+", age_groups$age_string[nrow(age_groups)])) {
  age_groups$age_lower[nrow(age_groups)] <- as.numeric(sub("\\+", "", age_groups$age_string[nrow(age_groups)]))
  age_groups$age_upper[nrow(age_groups)] <- 120  # or whatever max age you want
}

# ----- Read preprocessed data with coords 
meta_data_with_coords <- readRDS(here('Results/meta_data_with_coords.rds'))
colnames(meta_data_with_coords)
nrow(meta_data_with_coords)

# ----- Read labels data
meta_data_with_labels <- read.csv(here('Results/meta_data_with_labels.csv'))
colnames(meta_data_with_labels)
nrow(meta_data_with_labels)

meta_data_onnv_samples <- read.csv(here('Results/meta_data_onnv_samples_with_labels.csv'))
nrow(meta_data_onnv_samples)


# Convert coords to Easting and Northing
sp_vill <- SpatialPoints(cbind(meta_data_with_labels$Longitude, meta_data_with_labels$Latitude))
points_to_extract <- terra::vect(sp_vill)
length(points_to_extract)

data_points <- meta_data_with_coords %>%
  st_drop_geometry() %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

data_utm <- st_transform(data_points, crs = 32633)
coords_utm <- st_coordinates(data_utm) / 1000  # Convert to km
colnames(coords_utm) <- c("Easting", "Northing")

# Add Easting and Northing to dataframe
meta_data_with_coords$Easting <- coords_utm[, "Easting"]
meta_data_with_coords$Northing <- coords_utm[, "Northing"]

meta_data_with_labels$Easting <- meta_data_with_coords$Easting
meta_data_with_labels$Northing <- meta_data_with_coords$Northing

meta_data_onnv_samples$Easting <- meta_data_with_coords$Easting
meta_data_onnv_samples$Northing <- meta_data_with_coords$Northing

range(meta_data_with_labels$AgeInYears, na.rm = TRUE)


model_data <- meta_data_with_labels
model_data_onnv_samples <- meta_data_onnv_samples

sum(is.na(model_data$AgeInYears))
sum(model_data$AgeInYears == 0, na.rm = TRUE)
sum(is.na(model_data$Easting) & !is.na(model_data$Northing))
sum(is.na(model_data$ONNV_pos))

# --- Run INLA model for ONNV (with historic year of intro 1900)
onnv_results_pop_grid <- run_inla(
  year_intro = 1900,
  data = model_data,
  cam_pop = cam_pop,
  positive_col = "ONNV_pos")


# --- Run INLA model for ONNV (with historic year of intro 1900)
onnv_results_pop_grid_onnv_samples <- run_inla(
  year_intro = 1900,
  data = model_data_onnv_samples,
  cam_pop = cam_pop,
  positive_col = "ONNV_pos")

# --- Save prediction results 
saveRDS(onnv_results_pop_grid, here('Results/ONNV_INLAResults.rds'))
saveRDS(onnv_results_pop_grid_onnv_samples, here('Results/ONNV_INLAResults_ONNV_samples.rds'))

# Read saved results
onnv_results_pop_grid <- readRDS(here('Results/ONNV_INLAResults.rds'))



# --- SPATIAL PREDICTIONS: # Overall cameroon estimates ---- 

# overall FOI
foi_summary <- onnv_results_pop_grid$output$summary.fixed
est_cameroonwide_foi <- list(
  mean = exp(foi_summary$mean),
  ciL  = exp(foi_summary$`0.025quant`),
  ciU  = exp(foi_summary$`0.975quant`)
)
est_cameroonwide_foi


foi_summary_onnv_samples <- onnv_results_pop_grid_onnv_samples$output$summary.fixed
est_cameroonwide_foi_onnv_samples <- list(
  mean = exp(foi_summary_onnv_samples$mean),
  ciL  = exp(foi_summary_onnv_samples$`0.025quant`),
  ciU  = exp(foi_summary_onnv_samples$`0.975quant`)
)
est_cameroonwide_foi_onnv_samples

# -- cameroon wide summary 
compute_foi_metrics <- function(foi_val, age_groups, w_age, cam_pop, total_cameroon_pop) {
  
  avg_susceptible_by_age <- numeric(nrow(age_groups))
  
  for (j in 1:nrow(age_groups)) {
    a_lower <- age_groups$age_lower[j]
    a_upper <- age_groups$age_upper[j]
    age_width <- a_upper - a_lower
    
    if (foi_val > 1e-10) {
      avg_susceptible_by_age[j] <- (1 / (foi_val * age_width)) *
        (exp(-foi_val * a_lower) -
           exp(-foi_val * a_upper))
    } else {
      avg_susceptible_by_age[j] <- 1
    }
  }
  
  avg_susceptible  <- sum(avg_susceptible_by_age * w_age)
  avg_seroprev     <- 1 - avg_susceptible
  infections       <- total_cameroon_pop * foi_val * avg_susceptible
  
  list(
    avg_susceptible = avg_susceptible,
    avg_seroprev    = avg_seroprev,
    infections      = infections
  )
}

# --- Run for mean, ciL, ciU
total_cameroon_pop <- sum(values(cam_pop), na.rm = TRUE)

metrics_mean <- compute_foi_metrics(est_cameroonwide_foi$mean, age_groups, w_age, cam_pop, total_cameroon_pop)
metrics_ciL  <- compute_foi_metrics(est_cameroonwide_foi$ciL,  age_groups, w_age, cam_pop, total_cameroon_pop)
metrics_ciU  <- compute_foi_metrics(est_cameroonwide_foi$ciU,  age_groups, w_age, cam_pop, total_cameroon_pop)

# ---summary 
cameroon_summary <- data.frame(
  metric    = c("FOI", "Avg susceptible", "Avg seroprev", "Expected infections"),
  mean      = c(est_cameroonwide_foi$mean,  metrics_mean$avg_susceptible, metrics_mean$avg_seroprev, metrics_mean$infections),
  ciL       = c(est_cameroonwide_foi$ciL,   metrics_ciL$avg_susceptible,  metrics_ciL$avg_seroprev,  metrics_ciL$infections),
  ciU       = c(est_cameroonwide_foi$ciU,   metrics_ciU$avg_susceptible,  metrics_ciU$avg_seroprev,  metrics_ciU$infections)
)

cameroon_summary_fmt <- cameroon_summary
cameroon_summary_fmt[, c("mean", "ciL", "ciU")] <- lapply(
  cameroon_summary[, c("mean", "ciL", "ciU")],
  function(x) formatC(x, format = "f", digits = 3)
)

print(cameroon_summary_fmt)

# number of individuals with a history of ONNV infection (ie seropositive)
seropos_mean <- total_cameroon_pop * metrics_mean$avg_seroprev
seropos_ciL  <- total_cameroon_pop * metrics_ciL$avg_seroprev
seropos_ciU  <- total_cameroon_pop * metrics_ciU$avg_seroprev

cat(sprintf(
  "We estimated that an average of %s (95%% CI: %s-%s) 
  individuals get infected each year and that in 2020, 
  %s (95%% CI: %s-%s) individuals had a history of ONNV infection, 
  representing %.1f%% (95%% CI: %.1f%%-%.1f%%) of the population.",
  formatC(metrics_mean$infections, format = "f", digits = 0),
  formatC(metrics_ciL$infections,  format = "f", digits = 0),
  formatC(metrics_ciU$infections,  format = "f", digits = 0),
  formatC(seropos_mean, format = "f", digits = 0),
  formatC(seropos_ciL,  format = "f", digits = 0),
  formatC(seropos_ciU,  format = "f", digits = 0),
  metrics_mean$avg_seroprev * 100,
  metrics_ciL$avg_seroprev  * 100,
  metrics_ciU$avg_seroprev  * 100
))



# --- Cameroon Wide prediction, aggregated by region 
foi_onnv <- predicted_foi(onnv_results_pop_grid, onnv_results_pop_grid$coop, pathogen_name = "ONNV")
sero_onnv <- predicted_seroprevalence( foi_result = foi_onnv, model = onnv_results_pop_grid,
  age_groups = age_groups,
  age_weights = w_age,
  crs = 32633,
  pathogen_name = "ONNV"
)
infections_onnv <- predicted_annual_infections(
  foi_result = foi_onnv,
  model = onnv_results_pop_grid,
  age_groups = age_groups,
  age_weights = w_age,
  cam_pop = cam_pop,  
  crs = 32633,
  pathogen_name = "ONNV"
)



cameroon_regions <- ne_states(country = "Cameroon", returnclass = "sf")
# rename 'name' to region
regions_sf <- cameroon_regions %>%
  dplyr::select(region = name) %>%   # standardise column name
  st_make_valid()

# same aggregation as used for model fitting (10x10 km)
cam_pop_agg <- terra::aggregate(cam_pop, fact = 10, fun = sum, na.rm = TRUE)

foi_region <- aggregate_predictions_by_region(
  pred_sf   = foi_onnv$foi_sf,
  regions_sf = regions_sf,
  cam_pop   = cam_pop_agg,
  value_col = "foi",
  agg_type  = "weighted_mean"
)
print(foi_region)


prev_region <- aggregate_predictions_by_region(
  pred_sf    = sero_onnv$prev_sf,
  regions_sf = regions_sf,
  cam_pop    = cam_pop_agg,
  value_col  = "prev",
  agg_type   = "weighted_mean"
)
print(prev_region)

infection_region <- aggregate_predictions_by_region(
  pred_sf    = infections_onnv$infections_sf,
  regions_sf = regions_sf,
  value_col  = "infections",
  agg_type   = "sum"
)
print(infection_region)

# save region level predictions
region_level_predictions <- list(
  foi = foi_region,
  prev = prev_region,
  infections = infection_region
)
saveRDS(region_level_predictions, here('Results/region_level_predictions.rds'))


# Discussion - prob of disease,  acute cases, arthralgic cases and  deaths per year occur in Cameroon each year
average_annual_infections <- metrics_mean$infections
ciL_annual_infections <- metrics_ciL$infections
ciU_annual_infections <- metrics_ciU$infections

prob_disease <- 0.5 # assuming same as CHIKV
prob_mild <- 0.88 # Given disease, probability of it being mild (from Gabrial paper)
prob_severe <- 0.12  # Given disease, probability of it being severe (from Gabrial paper)
prob_medically_attended <- 0.0113
prob_chronic_given_severe <- 0.44 #  Kang et al 

# acute cases 
acute_cases <- average_annual_infections * prob_disease 
acute_cases_ciL <- ciL_annual_infections * prob_disease 
acute_cases_ciU <- ciU_annual_infections * prob_disease 


cat("Estimated number of acute cases per year:", acute_cases)

# chronic arthlgic cases
#arthralgic_cases  <- average_annual_infections * prob_disease * prob_severe * prob_chronic_given_severe
#arthralgic_cases_ciL  <- ciL_annual_infections * prob_disease * prob_severe * prob_chronic_given_severe
#arthralgic_cases_ciU  <- ciU_annual_infections * prob_disease * prob_severe * prob_chronic_given_severe

arthralgic_cases  <- average_annual_infections  * prob_medically_attended/2 

cat("Estimated number of arthralgic cases per year:", arthralgic_cases)

prob_ifr <- 4.2 * 10^-5 #in Brazil, 

death <- average_annual_infections  * prob_ifr
death_ciL <- ciL_annual_infections * prob_disease * prob_ifr
death_ciU <- ciU_annual_infections * prob_disease * prob_ifr

cat("Estimated number of deaths per year:", death)



# supplementary data - region level prediction 

# supplementary data 
region_level_predictions <- readRDS(here("Results/region_level_predictions.rds"))
region_level_predictions


# comparison of malaria risk with ONNV risk 
# Number of newly diagnosed Plasmodium falciparum cases per 1,000 population (using 2024)

pf <- terra::rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/clippedlayers-4/202508_Global_Pf_Incidence_Rate_CMR_2024.tiff')
pf_rate <- pf[[1]]  # band 1 = incidence rate

global(pf_rate, fun = c("min", "max"), na.rm = TRUE)

points_vect <- terra::vect(data_points)
meta_data_with_coords$pf_incidence <- terra::extract(
  pf_rate,
  points_vect
)[,2]

meta_data_with_labels$pf_incidence <- meta_data_with_coords$pf_incidence
summary(meta_data_with_labels$pf_incidence)


breaks <- seq(
  min(meta_data_with_labels$pf_incidence, na.rm = TRUE),
  max(meta_data_with_labels$pf_incidence, na.rm = TRUE),
  length.out = 8
)

onnv_pf_incidence <- calculate_prop_by_variable(
  data = meta_data_with_labels,
  var_col = "pf_incidence",
  positive_col = "ONNV_pos",
  breaks_max = breaks[-1],
  breaks_min = breaks[-length(breaks)]
)

summary(onnv_pf_incidence$log_model)

prop_pf_onnv <- make_plot_onnv(
  onnv_pf_incidence$obs,
  meta_data_with_labels$pf_incidence,
  "Pf Incidence",
  color ="#16622b", pos_col =  "ONNV_pos"
)
print(prop_pf_onnv)


ggsave("Results/supplementary_fig4.png", 
       plot = prop_pf_onnv,
       width = 6, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")
