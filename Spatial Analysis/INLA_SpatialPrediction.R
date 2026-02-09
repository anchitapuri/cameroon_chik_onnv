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

# --- Source functions
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/Spatial Analysis/Functions.R'))

# Population and mosquito rasts
anopheles_funestus <- rast('2010_Anopheles_funestus_CMR.tiff')
anopheles_gambiae <- rast('2010_Anopheles_gambiae_ss_CMR.tiff')
cam_pop <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_ppp_2020_UNadj.tif")
sum(values(cam_pop), na.rm = TRUE)

# Cameroon population by age
cameroon_age_2025 <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/CameroonAge2025.csv')
cameroon_age_2025 <- cameroon_age_2025 %>%
  mutate(total = M + F)

# ----- Read preprocessed data with coords 
meta_data_with_coords <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/meta_data_with_coords.rds')
colnames(meta_data_with_coords)
nrow(meta_data_with_coords)

# ----- Read labels data
meta_data_with_labels <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/final_meta_data_with_labels.csv')
nrow(meta_data_with_labels)

# Convert coords to Easting and Northing
sp_vill <- SpatialPoints(cbind(meta_data_with_coords$Longitude, meta_data_with_coords$Latitude))
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

# --- Run INLA model for ONNV (with historic year of intro 1900)
# --- Use population raster for prediction grid
onnv_results_pop_grid <- run_inla(
  year_intro = 1900,
  data = model_data,
  cam_pop = cam_pop,
  positive_col = "ONNV_pos")

                         
# --- Save best model results 
saveRDS(onnv_results_pop_grid, '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/ONNV_INLAResults.rds')


nrow(onnv_results_pop_grid$cooe)


# --- Extract and plot FOI
foi_onnv <- extract_and_plot_foi(onnv_results_pop_grid, onnv_results_pop_grid$coop, pathogen_name = "ONNV")
# overall FOI
est_cameroonwide_foi <-exp(onnv_results_pop_grid$output$summary.fixed$mean)
print(est_cameroonwide_foi)
range(foi_onnv$foi_sf$foi, na.rm = TRUE)

# --- Prob of seropositive proportion 
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

# Plot
sero_onnv <- plot_predicted_seroprevalence(
  foi_result = foi_onnv,
  model = onnv_results_pop_grid,
  age_groups = age_groups,
  age_weights = w_age,
  crs = 32633,
  pathogen_name = "ONNV"
)

sero_onnv$prev_range

# --- Annual Infections 
infections_onnv <- plot_predicted_annual_infections(
  foi_result = foi_onnv,
  model = onnv_results_pop_grid,
  age_groups = age_groups,
  cam_pop = cam_pop,  # Your spatial population data
  age_weights = w_age,  # Age distribution weights
  crs = 32633,
  pathogen_name = "ONNV"
)

# Access results
infections_onnv$total_infections  # Total infections across all locations
infections_onnv$susceptible_people  # Breakdown by age group
infections_onnv$seropositive_people  # Min and max by location

# combined plots 
maps <- foi_onnv$plot + sero_onnv$plot + infections_onnv$plot 


# --- Save Figure 1a
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig4b.png", 
       plot = maps,
       width = 12, 
       height = 15, 
       units = "in", 
       dpi = 300,
       bg = "white")


# --- Validate model----
avg_susceptible_by_age <- numeric(nrow(age_groups))

for (j in 1:nrow(age_groups)) {
  a_lower <- age_groups$age_lower[j]
  a_upper <- age_groups$age_upper[j]
  age_width <- a_upper - a_lower
  
  if (est_cameroonwide_foi > 1e-10) {
    avg_susceptible_by_age[j] <- (1/(est_cameroonwide_foi * age_width)) * 
                                  (exp(-est_cameroonwide_foi * a_lower) - 
                                   exp(-est_cameroonwide_foi * a_upper))
  } else {
    avg_susceptible_by_age[j] <- 1
  }
}


# Weight by age distribution to get overall average susceptible
cameroon_avg_susceptible <- sum(avg_susceptible_by_age * w_age)
cameroon_avg_susceptible
est_cameroonwide_foi

cameroon_avg_seroprev = 1 - cameroon_avg_susceptible
cameroon_avg_seroprev

total_cameroon_pop <- sum(values(cam_pop), na.rm = TRUE)  # ~26 million
expected_infections <- total_cameroon_pop * est_cameroonwide_foi * cameroon_avg_susceptible

sum(cameroon_age_2025$total)


# --- pop at predicted locations
pred_coords_sf <- st_as_sf(
  data.frame(X = onnv_results$coop[, "X"] * 1000, 
              Y = onnv_results$coop[, "Y"] * 1000),
  coords = c("X", "Y"),
  crs = 32633
)
# transform to match the cam_pop raster
pred_coords_transformed <- sf::st_transform(pred_coords_sf, crs = terra::crs(cam_pop))

# population raster has 13,710 × 9,233 = 126+ million cells
# 15,035 prediction points
# sampling less than 0.01% of the cells
# aggregate population to match your grid resolution
cam_pop_agg <- terra::aggregate(cam_pop, fact = 60, fun = sum, na.rm = TRUE)

pop_at_locations <- terra::extract(cam_pop_agg, pred_coords_transformed, ID = FALSE)

# 1. Check population coverage
captured_pop <- sum(pop_at_locations, na.rm = TRUE)
total_pop <- sum(values(cam_pop), na.rm = TRUE)
cat("Population captured: ", round(captured_pop/total_pop * 100, 1), "%\n")

# 2. Check weighted average FOI
weighted_avg_foi <- sum(lambda_pred * pop_at_locations, na.rm = TRUE) / sum(pop_at_locations, na.rm = TRUE)
cat("Cameroon-wide FOI: ", est_cameroonwide_foi, "\n")
cat("Population-weighted FOI: ", weighted_avg_foi, "\n")




