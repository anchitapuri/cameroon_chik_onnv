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
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/Spatial Analysis/Functions.R'))

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
meta_data_with_coords <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/meta_data_with_coords.rds')
colnames(meta_data_with_coords)
nrow(meta_data_with_coords)

# ----- Read labels data
meta_data_with_labels <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/final_meta_data_with_labels.csv')
colnames(meta_data_with_labels)
nrow(meta_data_with_labels)

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

model_data <- meta_data_with_labels
colnames(model_data)
nrow(model_data)

sum(is.na(model_data$AgeInYears))
sum(model_data$AgeInYears == 0, na.rm = TRUE)
sum(is.na(model_data$Easting) & !is.na(model_data$Northing))
sum(is.na(model_data$ONNV_pos))

# --- Run INLA model for ONNV (with historic year of intro 1900)
# --- Use population raster for prediction grid
onnv_results_pop_grid <- run_inla(
  year_intro = 1900,
  data = model_data,
  cam_pop = cam_pop,
  positive_col = "ONNV_pos")

# --- Save best model results 
saveRDS(onnv_results_pop_grid, '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/ONNV_INLAResults.rds')


# --- ESTIMATED LOCATION CALCULATIONS BY REGION ----
index_est <- inla.stack.index(onnv_results_pop_grid$stk.full, tag = "est")$data
index_pred <- inla.stack.index(onnv_results_pop_grid$stk.full, tag = "pred")$data
length(index_est)

# Extract the intercept
eta_est <- onnv_results_pop_grid$output$summary.linear.predictor[index_est, "mean"]
eta_pred <- onnv_results_pop_grid$output$summary.linear.predictor[index_pred, "mean"]

# account for age structure of data
age_est <- onnv_results_pop_grid$data_filtered$years_of_exposure[index_est]
log_lambda_est <- eta_est - log(age_est)
lambda_est <- exp(log_lambda_est)

# dont need to take into account age  
lambda_pred <- exp(eta_pred)

range(lambda_est)
range(lambda_pred)


# Add lamba est to data (to then group by region)
est_data <- onnv_results_pop_grid$data_filtered

# add region
est_data$region <- sapply(strsplit(est_data$IdNumber, "-"), `[`, 2)
# assuming ZST is mis-spelled and is EST
est_data <- est_data %>%
  mutate(region = recode(region, "ZST" = "EST"))
unique(est_data$region)

# add lambda at estimated locations
est_data$lambda_est <- lambda_est


pop_region <- est_data %>%
  distinct(region, district_lower, Total_Population) %>%
  group_by(region) %>%
  summarise(
    total_population = sum(Total_Population),
    .groups = "drop"
  )

sum(pop_region$total_population)

lambda_region <- est_data %>%
  group_by(region) %>%
  summarise(
    lambda_weighted = weighted.mean(lambda_est, Total_Population),
    lambda_unweighted = mean(lambda_est),
    .groups = "drop"
  )

region_lambda <- left_join(pop_region, lambda_region, by = "region")
print(region_lambda)
sum(region_lambda$total_population)
max(region_lambda$lambda_weighted)
min(region_lambda$lambda_weighted)



# --- SPATIAL PREDICTIONS: # Overall cameroon estimates ---- 
# overall FOI
foi_summary <- onnv_results_pop_grid$output$summary.fixed
est_cameroonwide_foi <- list(
  mean = exp(foi_summary$mean),
  ciL  = exp(foi_summary$`0.025quant`),
  ciU  = exp(foi_summary$`0.975quant`)
)
est_cameroonwide_foi

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
  individuals get infected each year and that in 2020, %s (95%% CI: %s-%s) individuals had a history of ONNV infection, 
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


# --- Cameroon wide maps 
foi_onnv <- plot_predicted_foi(onnv_results_pop_grid, onnv_results_pop_grid$coop, pathogen_name = "ONNV")
range(foi_onnv$foi_df$foi)

# --- Prob of seropositive proportion 
sero_onnv <- plot_predicted_seroprevalence(
  foi_result = foi_onnv,
  model = onnv_results_pop_grid,
  age_groups = age_groups,
  age_weights = w_age,
  crs = 32633,
  pathogen_name = "ONNV"
)


# --- Annual Infections 
infections_onnv <- plot_predicted_annual_infections(
  foi_result = foi_onnv,
  model = onnv_results_pop_grid,
  age_groups = age_groups,
  cam_pop = cam_pop,  
  age_weights = w_age, 
  crs = 32633,
  pathogen_name = "ONNV"
)


infections_onnv$total_infections  # Total infections across all locations
infections_onnv$susceptible_people  # Breakdown by age group
infections_onnv$seropositive_people  # Min and max by location

# combined plots 
maps <- foi_onnv$plot + sero_onnv$plot + infections_onnv$plot 

# --- Save Figure 4a
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig4b.png", 
       plot = maps,
       width = 12, 
       height = 15, 
       units = "in", 
       dpi = 300,
       bg = "white")




# --- Cameroon Wide prediction, aggregated by region 

cameroon_regions <- ne_states(country = "Cameroon", returnclass = "sf")
# rename 'name' to region
regions_sf <- cameroon_regions %>%
  dplyr::select(region = name) %>%   # standardise column name
  st_make_valid()

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
