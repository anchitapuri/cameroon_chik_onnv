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

# --- Source functions
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/Spatial Analysis/Functions.R'))

# Get Cameroon boundary
cameroon <- ne_countries(country = "Cameroon", returnclass = "sf")

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
# ----- Read labels data
meta_data_with_labels <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/final_meta_data_with_labels.csv')
nrow(meta_data_with_labels)

# Convert coords to Easting and Northing
sp_vill <- SpatialPoints(cbind(meta_data_with_coords$Longitude, meta_data_with_coords$Latitude))
points_to_extract <- terra::vect(sp_vill)

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



# --- Run INLA model for ONNV (with historic year of intro 1900)
onnv_results <- run_inla(
  year_intro = 1900,
  data = model_data,
  cameroon = cameroon,
  positive_col = "ONNV_pos")



# --- Index of prediction and estimation stacks 
index_pred_onnv <- inla.stack.index(onnv_results$stk.full, "pred")$data
length(index_pred_onnv)
index_est_onnv <- inla.stack.index(onnv_results$stk.full, "est")$data
length(index_est_onnv)
str(cam_pop)
head(cam_pop)



# --- Extract and plot FOI
foi_onnv <- extract_and_plot_foi(onnv_results, onnv_results$coop, pathogen_name = "ONNV")


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


sero_onnv <- plot_predicted_seroprevalence(
  foi_result = foi_onnv,
  model = onnv_results,
  age_groups = age_groups,
  age_weights = w_age,
  crs = 32633,
  pathogen_name = "ONNV"
)


# --- Annual Infections 
infections_onnv <- plot_predicted_annual_infections(
  foi_result = foi_onnv,
  model = onnv_results,
  age_groups = age_groups,
  cam_pop = cam_pop,  # Your spatial population data
  age_weights = w_age,  # Age distribution weights
  crs = 32633,
  pathogen_name = "ONNV"
)

# Access results
infections_onnv$total_infections  # Total infections across all locations
infections_onnv$infections_by_age  # Breakdown by age group
infections_onnv$infections_range  # Min and max by location

# Check spacing between prediction points (in km)
pred_coords_df <- data.frame(
  X = onnv_results$coop[, "X"],
  Y = onnv_results$coop[, "Y"]
)

# Calculate typical spacing
x_spacing <- median(diff(sort(unique(pred_coords_df$X))))
y_spacing <- median(diff(sort(unique(pred_coords_df$Y))))

print(paste("Grid spacing: X =", x_spacing, "km, Y =", y_spacing, "km"))

# Use half the grid spacing as buffer radius
buffer_radius_km <- min(x_spacing, y_spacing) / 2
buffer_radius_deg <- buffer_radius_km / 111  # Convert km to degrees (approx)

# Then extract with this buffer
pred_coords_transformed <- st_transform(pred_coords_sf, crs = crs(cam_pop))


pred_buffers <- st_buffer(pred_coords_transformed, dist = buffer_radius_deg)
pop_extraction <- terra::extract(cam_pop, pred_buffers, fun = sum, na.rm = TRUE)

print(head(pop_extraction))
print(dim(pop_extraction))
print(names(pop_extraction))

st_crs(pred_coords_sf)
crs(cam_pop)

pop_at_locations <- pop_extraction[, 2]
range(pop_at_locations, na.rm = TRUE)
sum(pop_at_locations, na.rm = TRUE)





# --- Model fits 
plot_age_seroprevalence_model_fits(onnv_results$year, onnv_results$output, model_data, "ONNV_pos")


                         
# --- Save best model results 
saveRDS(onnv_results, 'ONNV_INLAResults.rds', compress = "gzip")



# --- Mosquito and population proportion vs proportion positive 