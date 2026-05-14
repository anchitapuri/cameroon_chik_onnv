
library(ggspatial)
library(cowplot)
library(purrr)
library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(here)
library(patchwork)
library(cowplot)
library(ggpubr)
library(ggh4x)

# --- Source functions
source(here('R/Functions.R'))
source(here('R/MultiSeroFunctions.R'))


# general plot format 
base_theme <- theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm")
  )


# Import data 
cameroon_districts <- ne_states(country = "Cameroon", returnclass = "sf")
cameroon <- ne_countries(country = "Cameroon", returnclass = "sf")


# Import INLA model fits 
onnv_results_pop_grid <- readRDS(here('Results/ONNV_INLAResults.rds'))

# Multiset model +  chains + model data 
fit_full_model <- readRDS(here('Results/full_model_fits.rds'))
chains_full <- fit_full_model$draws(format='df')
chains_df_full <- as.data.frame(chains_full)
preprocessed_data_full_model <- readRDS('Results/preprocessed_data_full_model.rds')

# meta data with labels and coordinates
meta_data_with_coords <- readRDS('Results/meta_data_clean_with_coords.rds')
meta_data_with_labels <- read.csv('Results/meta_data_with_labels.csv')


# --- Proportion positive by mosquito distributions (only ONNV) (using multisero probs)
# Anopheles bins (adjust if needed based on your data range
anoph_max <- seq(0, 1, 0.1)
anoph_min <- anoph_max - 0.5
anoph_min[which(anoph_min < 0)] <- 0

aegmax <- seq(0,1,0.1)
aegmin <- aegmax - 0.5
aegmin[which(aegmin<0)] <- 0


df_fun_binary <- calculate_prop_by_variable (
  data = meta_data_with_labels,
  var_col = "fun_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = anoph_max, 
  breaks_min = anoph_min)

df_gam_binary <- calculate_prop_by_variable (
  data = meta_data_with_labels,
  var_col = "gam_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = anoph_max, 
  breaks_min = anoph_min)


df_aegypti_binary <- calculate_prop_by_variable (
  data = meta_data_with_labels,
  var_col = "aeg_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = aegmax, 
  breaks_min = aegmin)

df_albopictus_binary <- calculate_prop_by_variable (
  data = meta_data_with_labels,
  var_col = "alb_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = aegmax, 
  breaks_min = aegmin)

# Plot ONNV with AN. Gambia 
prop_gam_prev <- make_plot_onnv(
  df_gam_binary$obs,
  meta_data_with_labels$gam_pw_district,
  "Suitability Anopheles gambiae",
  color ="#165262", pos_col =  "ONNV_pos"
)


ggsave("Results/fig4d.png", 
       plot = prop_gam_prev,
       width = 6, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")


# beta 
print(summary(df_fun_binary$log_model))
print(summary(df_gam_binary$log_model))
print(summary(df_aegypti_binary$log_model))
print(summary(df_albopictus_binary$log_model))

# ---ONNV with population density 

meta_data_with_labels <- meta_data_with_labels %>%
  mutate(
    pop_density = Total_Population / area_km2,
    log_pop_density = log(pop_density + 1)  # +1 to avoid log(0) issues
  )

popdenmax <-seq(-1, 4.2, 0.1)
popdenmin <-popdenmax-2

df_onnv_pop <- calculate_prop_by_variable(
  data = meta_data_with_labels, 
  var_col = "log_pop_density", 
  positive_col = "ONNV_pos",
  breaks_max = popdenmax, 
  breaks_min = popdenmin)

summary(df_onnv_pop$log_model)


# CHIK with vectors 
# ---  CHIK 
df_fun_chik <- calculate_prop_by_variable(
  data = meta_data_with_labels,
  var_col = "fun_pw_district", 
  positive_col = "CHIK_pos",
  breaks_max = anoph_max, 
  breaks_min = anoph_min)

# Gambiae
df_gam_chik <- calculate_prop_by_variable(
  data = meta_data_with_labels,
  var_col = "gam_pw_district", 
  positive_col = "CHIK_pos",
  breaks_max = anoph_max, 
  breaks_min = anoph_min)

# Aegypti
df_aegypti_chik <- calculate_prop_by_variable(
  data = meta_data_with_labels, 
  var_col = "aeg_pw_district", 
  positive_col = "CHIK_pos",
  breaks_max = aegmax, 
  breaks_min = aegmin)

  # Albopictus
df_albopictus_chik <- calculate_prop_by_variable(
  data = meta_data_with_labels, 
  var_col = "alb_pw_district", 
  positive_col = "CHIK_pos",
  breaks_max = aegmax, 
  breaks_min = aegmin)



# beta 
print(summary(df_fun_chik$log_model))
print(summary(df_gam_chik$log_model))

print(summary(df_aegypti_chik$log_model))
print(summary(df_albopictus_chik$log_model))


# Plot where model inferred CHIKV pos cases fall 
chik_pos <- onnv_results_pop_grid$data_filtered |>
  dplyr::filter(CHIK_pos == 1)

chik_pos <- chik_pos %>%
  count(district_lower) %>%
  left_join(
    chik_pos %>% distinct(district_lower, Longitude, Latitude),
    by = "district_lower"
  ) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)



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


cities_buffer <- cities_sf %>%
  st_transform(crs = 32632) %>%       # project to metres (UTM zone 32N for Cameroon)
  st_buffer(dist = 100000) %>%         # 50km radius — adjust this
  st_transform(crs = 4326)            # back to WGS84


# 2. Find which chik_pos points fall within either buffer
chik_in_cities <- chik_pos %>%
  st_join(cities_buffer, join = st_within) %>%
  filter(!is.na(city))                # keep only those within a buffer

# 3. Calculate the proportion
total_cases    <- sum(chik_pos$n)
cases_in_cities <- sum(chik_in_cities$n)
prop <- cases_in_cities / total_cases * 100
print(prop)


cameroon_plot <- ggplot() +
  # Base Cameroon map
  geom_sf(
    data = cameroon,
    fill = "white",
    colour = "#252525",
    linewidth = 0.3
  ) +
  # All chik_pos districts as points
  geom_sf(
    data = chik_pos,
    colour = "grey60",
    size = 2,
    alpha = 0.7
  ) +
  # City label points
  geom_sf(
    data = cities_sf,
    shape = 18,   # diamond
    size = 4,
    colour = "black"
  ) +
  geom_sf(
    data = cities_buffer,
    fill = NA,
    colour = "black",
    linewidth = 0.5,
    linetype = "dashed"
  )


print(cameroon_plot)


