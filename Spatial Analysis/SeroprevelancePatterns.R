
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

age_prev_model_fits <- plot_age_seroprevalence_model_fits(
  year_intro = onnv_results$year,
  result = onnv_results, 
  data = model_data,
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = "a" # a == ONNV 
)
print(age_prev_model_fits)
class(age_prev_model_fits)
print(age_prev_model_fits[[1]])

ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig3.png", 
       plot = age_prev_model_fits[[1]],    # swap for your actual plot object name
       width = 18, 
       height = 12, 
       units = "in", 
       dpi = 300,
       bg = "white")


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
  
chik_infections <-ggplot() +
  geom_sf(
  data = cameroon_sf,
  fill = NA,          # no fill
  colour = "black",
  linewidth = 0.6
) +
  geom_sf(
    data = st_centroid(chik_pos_sf),
    aes(colour = aeg_pw_district),
    size = 8,
    alpha = 0.85
  ) +
  scale_colour_gradient(
   low = "#0298fc",
    high = "#d70048",
    name = "Aegypti",
    guide = guide_colorbar(
      direction = "horizontal",  
      barheight = unit(0.2, "cm"),      # thin (horizontal)
      barwidth = unit(6.5, "cm"),       # wide (horizontal)
      ticks = TRUE,
      title.position = "bottom",
      label.position = "bottom",
      ticks.length = unit(0.1, "cm"),
      title.vjust = 0.5)
    ) +
  geom_sf(
    data = cities_sf,
    aes(fill = city),
    shape = 24,
    colour = "#f1f3f4",
    size = 8,
    stroke = 0.6
  ) +
  scale_fill_manual(
    values = c("Yaoundé" = "#013018", "Douala" = "#013018"),
    name = ""
  ) +
  annotation_scale(
    bar_cols = c("black", "white"),  # alternating black/white like the reference
    height = unit(0.2, "cm"),
    text_family = "sans", 
    text_cex = 1.5
  ) +
  theme(
    legend.position = c(0.3, 0.6), 
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 20),        
    legend.title = element_text(size = 20),   
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
  )

print(chik_infections)

ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig4a.png", 
       plot = chik_infections,    # swap for your actual plot object name
       width = 18, 
       height = 12, 
       units = "in", 
       dpi = 300,
       bg = "white")



# --- Proportion positive by mosquito distributions (only ONNV) 
# Anopheles bins (adjust if needed based on your data range
anoph_max <- seq(0, 1, 0.1)
anoph_min <- anoph_max - 0.5
anoph_min[which(anoph_min < 0)] <- 0

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

df_fun_filtered <- df_fun %>% filter(n >= 10)  # Adjust threshold as needed
df_gam_filtered <- df_gam %>% filter(n >= 10)


df_anopheles <- rbind(df_fun_filtered, df_gam_filtered)


# Funestus plot
prop_fun_prev <- ggplot(df_fun_filtered, aes(x = x, y = y)) +
  geom_point(color = "#c1518b", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#c1518b") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
  labs(x = "Proportion Anopheles funestus", y = "Proportion ONNV positive") +
  theme_classic() + 
  theme(panel.grid = element_blank(),
    aspect.ratio = 0.75,
    axis.line = element_line(color = "black", linewidth = 0.7),  # Add x and y axis lines
    axis.title.x = element_text(size = 24),                             # X-axis label
    axis.title.y = element_text(size = 24),                             # Y-axis label
    axis.text.x = element_text(size = 20),                              # X-axis tick labels
    axis.text.y = element_text(size = 20),                              # Y-axis tick labels
    legend.title = element_text(size = 20),                             # Legend title
    legend.text = element_text(size = 20),                               # Legend text
    axis.ticks.x = element_line(color = "black", size = 0.5),  # X-axis ticks only
    axis.ticks.y = element_line(color = "black", size = 0.5),  # Y-axis ticks only
    axis.ticks.length = unit(0.2, "cm")

  )

# Gambiae plot
prop_gam_prev <- ggplot(df_gam_filtered, aes(x = x, y = y)) +
  geom_point(color = "#165262", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#165262") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05))+
  labs(x = "Proportion Anopheles gambiae", y = "Proportion ONNV positive") +
  theme_classic() + 
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 0.75,
    axis.line = element_line(color = "black", linewidth = 0.7),  # Add x and y axis lines
    axis.title.x = element_text(size = 24),                             # X-axis label
    axis.title.y = element_text(size = 24),                             # Y-axis label
    axis.text.x = element_text(size = 20),                              # X-axis tick labels
    axis.text.y = element_text(size = 20),                              # Y-axis tick labels
    legend.title = element_text(size = 20),                             # Legend title
    legend.text = element_text(size = 20),                               # Legend text
    axis.ticks.x = element_line(color = "black", size = 0.5),  # X-axis ticks only
    axis.ticks.y = element_line(color = "black", size = 0.5),  # Y-axis ticks only
    axis.ticks.length = unit(0.2, "cm")
  )

# Combined plot
combined_mosquito_plots <- prop_fun_prev + prop_gam_prev
print(combined_mosquito_plots)


# --- Save Figure 1a
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig4c.png", 
       plot = combined_mosquito_plots,
       width = 13, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")
