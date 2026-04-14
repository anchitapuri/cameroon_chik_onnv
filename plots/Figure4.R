
# --- Source functions
source(here('R/Functions.R'))

# --- Load data
cameroon <- ne_countries(country = "Cameroon", returnclass = "sf")


onnv_results_pop_grid <- readRDS(here('Results/ONNV_INLAResults.rds'))

# --- Cameroon wide maps 
foi_onnv <- predicted_foi(onnv_results_pop_grid, onnv_results_pop_grid$coop, pathogen_name = "ONNV")
range(foi_onnv$foi_df$foi)

# --- Prob of seropositive proportion 
sero_onnv <- predicted_seroprevalence(
  foi_result = foi_onnv,
  model = onnv_results_pop_grid,
  age_groups = age_groups,
  age_weights = w_age,
  crs = 32633,
  pathogen_name = "ONNV"
)


# --- Annual Infections 
infections_onnv <- predicted_annual_infections(
  foi_result = foi_onnv,
  model = onnv_results_pop_grid,
  age_groups = age_groups,
  cam_pop = cam_pop,  
  age_weights = w_age, 
  crs = 32633,
  pathogen_name = "ONNV"
)


# combined plots 
maps <- foi_onnv$plot + sero_onnv$plot + infections_onnv$plot  +  plot_layout(ncol = 3)

print(maps)
ggsave(here('Results/Fig4a.png'), 
       plot = maps,
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300,
       bg = "white")



# CHIK cases 
chik_pos <- onnv_results_pop_grid$data_filtered |>
  dplyr::filter(CHIK_pos == 1)

chik_pos <- chik_pos %>%
  count(district_lower) %>%
  left_join(
    chik_pos %>% distinct(district_lower, Longitude, Latitude),
    by = "district_lower"
  ) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)



set.seed(42)  # for reproducibility
chik_centroids_jittered <- chik_pos |>
  dplyr::mutate(geometry = geometry + sf::st_sfc(lapply(seq_len(dplyr::n()), function(i) {
    sf::st_point(c(runif(1, -0.05, 0.05), runif(1, -0.05, 0.05)))
  })))
sf::st_crs(chik_centroids_jittered) <- sf::st_crs(chik_pos)

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

chik_infections <-ggplot()  +
  geom_sf(
    data = cameroon,
    fill = "white",
    colour = "#252525",
    linewidth = 0.3
  ) +
  geom_sf(
    data = chik_centroids_jittered,
    colour = '#004E66', 
    size = 8,
    alpha = 0.85)  + 
  scale_size(
  range = c(4, 11),
  breaks = c(1, 2, 3, 4),
  name = "CHIK+ve samples"
) +
  geom_sf(
    data = cities_sf,
    aes(fill = city),
    shape = 23,
    colour = "#00060f",
    size = 15,
    stroke = 0.6
  ) +
  scale_fill_manual(
    values = c("Yaoundé" = "#F0D3F7", "Douala" = "#B6C8A9"),
    name = ""
  ) +
  annotation_scale(
    bar_cols = c("black", "white"),  # alternating black/white like the reference
    height = unit(0.2, "cm"),
    text_family = "sans", 
    text_cex = 1.5
  ) +
  theme(
    legend.position = c(0.32, 0.66), 
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 20),        
    legend.title = element_text(size = 20),   
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
  )
print(chik_infections)

ggsave(here('Results/fig4f.png'), 
       plot = chik_infections,    # swap for your actual plot object name
       width = 10, 
       height = 10, 
       units = "in", 
       dpi = 300,
       bg = "white")



# Proportion of Anopheles positive vs proportion of ONNV positive
anoph_max <- seq(0, 1, 0.1)
anoph_min <- anoph_max - 0.5
anoph_min[which(anoph_min < 0)] <- 0


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

make_plot <- function(df_obs, raw_data, xlab, color, pos_col = "ONNV_pos") {

  ylab <- paste0("Proportion ", gsub("_pos", "", pos_col), "positive")

  obs_clean <- df_obs[!is.nan(df_obs$x), ]

  # truncate CIs to [0, 0.5]
  obs_clean$ymin <- pmax(obs_clean$ymin, 0)
  obs_clean$ymax <- pmin(obs_clean$ymax, 0.5)

  hist_df <- data.frame(x = as.numeric(raw_data))
  hist_df <- hist_df[!is.na(hist_df$x), , drop = FALSE]

    x_scale <- scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    expand = c(0, 0)
  )

  plot_hist <- ggplot(hist_df, aes(x = x)) +
    geom_histogram(fill = color, alpha = 0.5, bins = 30, color = NA) +
    x_scale +
    labs(x = NULL, y = "Count") +
    base_theme +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(t = 6, r = 14, b = 10, l = 14)
    )

  plot_scatter <- ggplot(obs_clean, aes(x = x, y = y)) +
    geom_point(color = color, size = 4, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0, color = color, alpha = 0.6, linewidth = 0.6
    ) +
    x_scale +
    coord_cartesian(ylim = c(0, 0.4), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, 0.4, 0.1)) +
    labs(x = xlab, y = ylab) +
    base_theme +
    theme(
      plot.margin = margin(t = 10, r = 14, b = 12, l = 14)
    )

  plot_hist / plot_scatter +
    patchwork::plot_layout(heights = c(2, 4))
  
}


df_gam_binary <- calculate_prop_by_variable (
  data = meta_data_with_labels,
  var_col = "gam_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = anoph_max, 
  breaks_min = anoph_min)


prop_gam_prev <- make_plot(
  df_gam_binary$obs,
  meta_data_with_labels$gam_pw_district,
  "Proportion Anopheles gambiae",
  color ="#165262", pos_col =  "ONNV_pos"
)

# --- Save Figures
ggsave(here('Results/fig4d.png'), 
       plot = prop_gam_prev,
       width = 8, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")
