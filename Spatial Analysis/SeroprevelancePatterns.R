
# Prevelance patterns by age 

# --- Source functions
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/Spatial Analysis/Functions.R'))
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/MultiSeroModel/MultiSeroFunctions.R'))

# Import multisero fits 
fit <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/redone_final_model_fits.rds')
# Import INLA model fits 
onnv_results <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/ONNV_INLAResults.rds')
# Extract chains 
chains <- fit$draws(format='df')
chains_df <- as.data.frame(chains)

# Read preprocessed data
preprocessed_data <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/preprocessed_data_alpha.rds')

# --- Plot age seroprevalence model fits
# prepare data for stan
preprocessed_data$data$infM



quartz(width = 14, height = 14)
plot_age_seroprevalence_model_fits(
  year_intro = onnv_results$year,
  result = onnv_results, 
  data = model_data,
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = "a"
)






# --- Plot CHIK infection locations 
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





# --- Proportion positive by mosquito distributions (only ONNV) 
# Anopheles bins (adjust if needed based on your data range
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

