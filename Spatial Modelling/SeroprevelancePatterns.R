# 1) Extract population weighted Latitude and Longitude (from district centroids)
# 2) Patterns of observed seroprevelance (aggregated by district, mosquito proportion and age) 

# --- Source functions
source(here('Functions.R'))


# Load data with pop weighted coords (post preprocessing) 
cameroon_data <- readRDS('sf_meta_data_with_coords_pw.rds')
nrow(meta_data_with_coords)
length(unique(cameroon_data$Sample))
cameroon_data <- cameroon_data[!duplicated(cameroon_data$Sample), ]

cameroon_data$year_of_survey <- as.numeric(substr(cameroon_data$Sample, 1, 4))
unique(cameroon_data$year_of_survey)


# Load population rasters
cam_pop <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_ppp_2020_UNadj.tif")
cam_pop_den <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_pd_2020_1km_UNadj.tif")

# ---  Add UTM coordinates
sp_vill <- SpatialPoints(cbind(cameroon_data$Longitude, cameroon_data$Latitude))
points_to_extract <- terra::vect(sp_vill)

data_points <- cameroon_data %>%
  st_drop_geometry() %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
data_utm <- st_transform(data_points, crs = 32633)
coords_utm <- st_coordinates(data_utm) / 1000  # Convert to km
colnames(coords_utm) <- c("Easting", "Northing")

# Add Easting and Northing to dataframe
cameroon_data$Easting <- coords_utm[, "Easting"]
cameroon_data$Northing <- coords_utm[, "Northing"]

# Extract population density values at district lat / long 
density_values <- terra::extract(cam_pop_den, points_to_extract)
cameroon_data$pop_density <- density_values[, 2]
cameroon_data$logpopden <- log(cameroon_data$pop_density, 10)


N <- nrow(cameroon_data)


# ---- RUN FUNCTIONS + PLOTS ----

# 1) --- Prevelance by district 
# CHIK: very few infections (n=28) 
# Check where CHIK infection fall within cameroon (coloured by Aedes Aegypti proportion) 
chik_pos <- cameroon_data |>
  dplyr::filter(CHIK_pos == 1)

chik_pos_sf <- st_as_sf(
  chik_pos,
  coords = c("Longitude", "Latitude"),
  crs = 4326   # WGS84
)
# Highlight two biggest cities in Cameroon 
cities <- data.frame(
  city = c("Yaoundé", "Douala"),
  Latitude  = c(3.8617, 4.0511),
  Longitude  = c(11.5202, 9.7679)
)

cities_sf <- st_as_sf(
  cities,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

ggplot() +
  geom_sf(
    data = cameroon_data,
    fill = "grey95",
    colour = "grey40",
    linewidth = 0.3
  ) +
  geom_sf(
    data = st_centroid(chik_pos_sf),
    aes(colour = aeg_pw_district),
    size = 3,
    alpha = 0.8
  ) +
  scale_colour_gradient(
    low = "yellow",
    high = "red",
    name = "Aegypti"
  ) +
  geom_sf(
    data = cities_sf,
    aes(fill = city),
    shape = 21,
    colour = "black",
    size = 4,
    stroke = 0.6
  ) +
  scale_fill_manual(
    values = c("Yaoundé" = "lightblue", "Douala" = "darkblue"),
    name = "Cities"
  ) +
  theme_minimal() +
  labs(
    title = "CHIK-positive samples in Cameroon",
    subtitle = "Point color shows Aedes aegypti levels"
  )


# ONNV: Prevelance by district 
onnv_by_district <- calculate_district_prevalence(cameroon_data, "ONNV_pos")

ggplot(onnv_by_district) +
  geom_sf(aes(fill = mean_positive), colour = "grey20", linewidth = 0.3) +
  scale_fill_viridis_c(
    option = "magma",
    name = "ONNV seroprevelance \n by district",
    limits = c(0, max(onnv_by_district$mean_positive, na.rm = TRUE))
  ) +
  theme_minimal()




# --- 2) ODDS RATIO 
# Calculate for CHIK
or_chik <- calculate_odds_ratio(cameroon_data, "CHIK_pos", n_iter = 100)

# Calculate for ONNV
or_onnv <- calculate_odds_ratio(cameroon_data, "ONNV_pos", n_iter = 100)


# Plot CHIK
ggplot(or_chik$plot_data, aes(x = distance, y = odds_ratio)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10() +
  labs(x = "Distance (km)", y = "Odds Ratio", title = "CHIKV") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )

# Plot ONNV
ggplot(or_onnv$plot_data, aes(x = distance, y = odds_ratio)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10() +
  labs(x = "Distance (km)", y = "Odds Ratio", title = "ONNV") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )




# --- 3) Proportion positive by mosquito distributions (only ONNV) 

# Anopheles bins (adjust if needed based on your data range)
anoph_max <- seq(0, 1, 0.1)
anoph_min <- anoph_max - 0.5
anoph_min[which(anoph_min < 0)] <- 0


# Funestus
df_fun <- calculate_prop_by_variable(
  cameroon_data, "fun_pw_district", "ONNV_pos", anoph_max, anoph_min
)
df_fun$species <- "Funestus"

# Gambiae
df_gam <- calculate_prop_by_variable(
  cameroon_data, "gam_pw_district", "ONNV_pos", anoph_max, anoph_min
)
df_gam$species <- "Gambiae"

df_anopheles <- rbind(df_fun, df_gam)


# Funestus plot
prop_fun_prev <- ggplot(df_fun, aes(x = x, y = y)) +
  geom_point(color = "purple", size = 2) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "purple") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(x = "Proportion Anopheles funestus", y = "Proportion ONNV positive") +
  theme_classic()

# Gambiae plot
prop_gam_prev <- ggplot(df_gam, aes(x = x, y = y)) +
  geom_point(color = "orange", size = 2) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "orange") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(x = "Proportion Anopheles gambiae", y = "Proportion ONNV positive") +
  theme_classic()

# Combined plot
(prop_fun_prev + prop_gam_prev)



#--- 4) Proportion positive by age / sex 
# distributions by age 
plot_age_seroprevalence_by_year_obs(model_data, 'ONNV_pos')
plot_age_seroprevalence_by_year_obs(model_data, 'CHIK_pos')

# distributions by age statified bt sex 
plot_age_seroprevalence_by_year_gender_obs(model_data, 'ONNV_pos')
plot_age_seroprevalence_by_year_gender_obs(model_data, 'CHIK_pos')


