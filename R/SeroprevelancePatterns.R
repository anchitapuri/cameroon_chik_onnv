
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
onnv_results_pop_grid <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/ONNV_INLAResults.rds')

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
# First create the density column in your data
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










# vectors + CHIK and ONNV all in one plot 
virus_colors <- c(
  CHIK = "#2e86ab",
  ONNV = "#b31459"
)

# ── Binning for all 4 vectors ─────────────────────────────────────────────────
aegmax <- seq(0, 1, 0.1); aegmin <- pmax(aegmax - 0.5, 0)
anoph_max <- seq(0, 1, 0.1); anoph_min <- pmax(anoph_max - 0.5, 0)

# Aegypti
df_aeg_chik_multisero <- calculate_prop_by_variable(meta_data_with_labels, "aeg_pw_district", "CHIK_pos", aegmax, aegmin)
df_aeg_onnv_multisero <- calculate_prop_by_variable(meta_data_with_labels, "aeg_pw_district", "ONNV_pos", aegmax, aegmin)

# Albopictus
df_alb_chik_multisero <- calculate_prop_by_variable(meta_data_with_labels, "alb_pw_district", "CHIK_pos", aegmax, aegmin)
df_alb_onnv_multisero <- calculate_prop_by_variable(meta_data_with_labels, "alb_pw_district", "ONNV_pos", aegmax, aegmin)

# Funestus
df_fun_chik_multisero <- calculate_prop_by_variable(meta_data_with_labels, "fun_pw_district", "CHIK_pos", anoph_max, anoph_min)
df_fun_onnv_multisero <- calculate_prop_by_variable(meta_data_with_labels, "fun_pw_district", "ONNV_pos", anoph_max, anoph_min)

# Gambiae
df_gam_chik_multisero <- calculate_prop_by_variable(meta_data_with_labels, "gam_pw_district", "CHIK_pos", anoph_max, anoph_min)
df_gam_onnv_multisero <- calculate_prop_by_variable(meta_data_with_labels, "gam_pw_district", "ONNV_pos", anoph_max, anoph_min)


make_plot_multi_virus <- function(obs_chik, obs_onnv,
                                  raw_data,
                                  xlab,
                                  ylim_upper = 0.4,
                                  colors = virus_colors) {

  combined <- dplyr::bind_rows(
    obs_chik %>%
      dplyr::filter(!is.nan(x)) %>%
      dplyr::mutate(
        virus = "CHIK",
        ymin  = pmax(ymin, 0),
        ymax  = pmin(ymax, ylim_upper)
      ),
    obs_onnv %>%
      dplyr::filter(!is.nan(x)) %>%
      dplyr::mutate(
        virus = "ONNV",
        ymin  = pmax(ymin, 0),
        ymax  = pmin(ymax, ylim_upper)
      )
  )
  combined$virus <- factor(combined$virus, levels = names(colors))

  hist_df <- data.frame(x = as.numeric(raw_data))
  hist_df <- hist_df[!is.na(hist_df$x), , drop = FALSE]

  x_scale <- scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    expand = c(0, 0)
  )

  # Histogram (top panel, shared across both virus rows)
  plot_hist <- ggplot(hist_df, aes(x = x)) +
    geom_histogram(fill = "grey60", alpha = 0.6, bins = 30, color = NA) +
    x_scale +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y      = element_text(size = 14),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      plot.margin  = margin(t = 0, r = 0, b = 0, l = 0)
    )

  # Scatter — faceted by virus, each with free y scale
  plot_scatter <- ggplot(combined, aes(x = x, y = y, color = virus)) +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0, alpha = 0.6, linewidth = 0.6
    ) +
    scale_color_manual(values = colors, name = NULL) +
    x_scale +
    # Free y scale so CHIK gets its own range rather than being squashed
    facet_wrap(
      ~virus,
      ncol   = 1,
      scales = "free_y") +
     ggh4x::facetted_pos_scales(
        y = list(
        virus == "CHIK" ~ scale_y_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)),
        virus == "ONNV" ~ scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1))
      )
     ) +
    labs(x = xlab, y = NULL) +
    theme(
      legend.position = "none",
      plot.margin      = margin(t = 0, r = 4, b = 8, l = 4),
      strip.background = element_blank(),
      strip.text       = element_blank(),
      strip.text.y     = element_blank(),
      axis.text.x      = element_text(size = 16),
      axis.text.y      = element_text(size = 14),
      axis.title.x     = element_text(size = 16),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.spacing    = unit(4, "pt")    # tighter gap between CHIK / ONNV panels
    )

  plot_hist / plot_scatter +
    patchwork::plot_layout(heights = c(1, 4))   # histogram smaller relative to 2 scatter rows
}


# ── Build 4 plots (unchanged calls) ──────────────────────────────────────────
prop_aeg_prev <- make_plot_multi_virus(
  df_aeg_chik_multisero$obs, df_aeg_onnv_multisero$obs,
  meta_data_with_labels$aeg_pw_district,
  xlab = "Proportion\nAedes aegypti"
)

prop_alb_prev <- make_plot_multi_virus(
  df_alb_chik_multisero$obs, df_alb_onnv_multisero$obs,
  meta_data_with_labels$alb_pw_district,
  xlab = "Proportion\nAedes albopictus"
)

prop_fun_prev <- make_plot_multi_virus(
  df_fun_chik_multisero$obs, df_fun_onnv_multisero$obs,
  meta_data_with_labels$fun_pw_district,
  xlab = "Proportion\nAnopheles funestus"
)

prop_gam_prev <- make_plot_multi_virus(
  df_gam_chik_multisero$obs, df_gam_onnv_multisero$obs,
  meta_data_with_labels$gam_pw_district,
  xlab = "Proportion\nAnopheles gambiae"
)


# ── Arrange 2×2 with shared y-axis label ─────────────────────────────────────
multisero_mosquito_pos_plots <- patchwork::wrap_plots(
  prop_aeg_prev, prop_alb_prev,
  prop_fun_prev, prop_gam_prev,
  ncol = 4
) 

print(multisero_mosquito_pos_plots)



# save plots 
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/xStarPres/multisero_mosquito_pos_plots.png", 
       plot = multisero_mosquito_pos_plots,
       width = 10.5, 
       height = 6.5, 
       units = "in", 
       dpi = 300,
       bg = "white")


