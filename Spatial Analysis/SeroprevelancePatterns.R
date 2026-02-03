
library(ggspatial)
library(cowplot)

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
chik_pos <- onnv_results$data_filtered |>
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

cameroon_sf <- ne_countries(
  country = "Cameroon",
  scale = "medium",
  returnclass = "sf"
)

ggplot() +
  geom_sf(
  data = cameroon_sf,
  fill = NA,          # no fill
  colour = "black",
  linewidth = 0.6
) +
  geom_sf(
    data = st_centroid(chik_pos_sf),
    aes(colour = aeg_pw_district),
    size = 3,
    alpha = 0.8
  ) +
  scale_colour_gradient(
    low = "#0157b2",
    high = "#8b3351",
    name = "Aegypti",
    guide = guide_colorbar(
      barheight = unit(2, "cm"),
      barwidth = unit(0.4, "cm"),
      ticks = TRUE,
      ticks.length = unit(0.15, "cm")
    )
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
  annotation_scale(
    location = "bl",        # bottom-left
    width_unit = "km",
    bar_cols = c("black", "white"),  # alternating black/white like the reference
    height = unit(0.2, "cm"),
    text_family = "sans"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )
 




# --- Proportion positive by mosquito distributions (only ONNV) 
# Anopheles bins (adjust if needed based on your data range
anoph_max <- seq(0, 1, 0.1)
anoph_min <- anoph_max - 0.5
anoph_min[which(anoph_min < 0)] <- 0


# Funestus
df_fun <- calculate_prop_by_variable(
  data = onnv_results$data_filtered, 
  var_col = "fun_pw_district", 
  chains_df = chains_df, 
  infM = preprocessed_data$data$infM, 
  pathogen_col = "a", 
  breaks_max = anoph_max, 
  breaks_min = anoph_min
)
df_fun$species <- "Funestus"

# Gambiae
df_gam <- calculate_prop_by_variable(
  data = onnv_results$data_filtered, 
  var_col = "gam_pw_district", 
  chains_df = chains_df, 
  infM = preprocessed_data$data$infM, 
  pathogen_col = "a", 
  breaks_max = anoph_max, 
  breaks_min = anoph_min
)
df_gam$species <- "Gambiae"

df_anopheles <- rbind(df_fun, df_gam)


df_fun <- old_calculate_prop_by_variable(
  data = onnv_results$data_filtered, 
  var_col = "fun_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = anoph_max, 
  breaks_min = anoph_min)
df_fun$species <- "Funestus"

# Gambiae
df_gam <- old_calculate_prop_by_variable(
  data = onnv_results$data_filtered, 
  var_col = "gam_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = anoph_max, 
  breaks_min = anoph_min)
df_gam$species <- "Gambiae"

df_anopheles <- rbind(df_fun, df_gam)

# Funestus plot
prop_fun_prev <- ggplot(df_fun, aes(x = x, y = y)) +
  geom_point(color = "#42026a", size = 2) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#42026a") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(x = "Proportion Anopheles funestus", y = "Proportion ONNV positive") +
  theme_classic()

# Gambiae plot
prop_gam_prev <- ggplot(df_gam, aes(x = x, y = y)) +
  geom_point(color = "#00a2ff", size = 2) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#00a2ff") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(x = "Proportion Anopheles gambiae", y = "Proportion ONNV positive") +
  theme_classic()

# Combined plot
combined <- prop_fun_prev + prop_gam_prev


ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/Figure4.png", plot = combined, width = 10, height = 4)
