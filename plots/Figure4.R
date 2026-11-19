
# --- Source functions
source(here('R/Functions.R'))

# --- Load data
cameroon <- ne_countries(country = "Cameroon", returnclass = "sf")


onnv_results_pop_grid <- readRDS(here('Results/ONNV_INLAResults.rds'))
meta_data_with_labels <- read.csv(here('Results/meta_data_with_labels.csv'))

# --- Cameroon wide maps 
foi_onnv <- predicted_foi_diff_col(onnv_results_pop_grid, onnv_results_pop_grid$coop, pathogen_name = "ONNV")
range(foi_onnv$foi_df$foi)

print(foi_onnv$plot)

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



# ONNV +ve iwht all vectors 
# Sliding windows: width = `frac` of the observed range, placed at `n` positions
# across that range. Windows overlap, so an individual falls in several.
make_windows <- function(x, frac = 0.3, n = 11) {
  rng   <- range(x, na.rm = TRUE)
  width <- frac * diff(rng)
  upper <- seq(rng[1] + width, rng[2], length.out = n)
  data.frame(min = upper - width, max = upper)
}

w_fun <- make_windows(model_data$fun_pw_district)
w_gam <- make_windows(model_data$gam_pw_district)
w_aeg <- make_windows(model_data$aeg_pw_district)
w_alb <- make_windows(model_data$alb_pw_district)


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



df_fun_binary <- calculate_prop_by_variable (
  data = meta_data_with_labels,
  var_col = "fun_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = w_fun$max, 
  breaks_min = w_fun$min)

df_gam_binary <- calculate_prop_by_variable (
  data = meta_data_with_labels,
  var_col = "gam_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = w_gam$max, 
  breaks_min = w_gam$min)


df_aegypti_binary <- calculate_prop_by_variable (
  data = meta_data_with_labels,
  var_col = "aeg_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = w_aeg$max, 
  breaks_min = w_aeg$min)

df_albopictus_binary <- calculate_prop_by_variable (
  data = meta_data_with_labels,
  var_col = "alb_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = w_alb$max, 
  breaks_min = w_alb$min)


# --- Plots
prop_fun_prev <- make_plot_onnv(
  df_fun_binary$obs,
  meta_data_with_labels$fun_pw_district,
   "Suitability<br>*Anopheles funestus*",
  color ="#023e8a", pos_col =  "ONNV_pos"
)

prop_gam_prev <- make_plot_onnv(
  df_gam_binary$obs,
  meta_data_with_labels$gam_pw_district,
  "Suitability<br>*Anopheles gambiae*",
  color ="#165262", pos_col =  "ONNV_pos"
)

prop_aeg_prev <- make_plot_onnv(
  df_aegypti_binary$obs,
  meta_data_with_labels$aeg_pw_district,
  "Suitability<br>*Aedes aegypti*",
  color ="#c1518b", pos_col =  "ONNV_pos"
)

prop_albo_prev <- make_plot_onnv(
  df_albopictus_binary$obs,
  meta_data_with_labels$alb_pw_district,
  "Suitability<br>*Aedes albopictus*",
  color = "#430726", pos_col =  "ONNV_pos"
)


onnv_with_vectors <- (prop_fun_prev | prop_gam_prev | prop_aeg_prev | prop_albo_prev) 


# --- Save Figures
ggsave(here('Results/fig4d_e.png'), 
       plot = onnv_with_vectors,
       width = 20, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")




# Proportion of Anopheles positive vs proportion of ONNV positive
#anoph_max <- seq(0, 1, 0.1)
#anoph_min <- anoph_max - 0.5
#anoph_min[which(anoph_min < 0)] <- 0

# Proportion of Anopheles positive vs proportion of ONNV positive
anoph_max <- seq(0, 1, 0.1)
anoph_min <- anoph_max - 0.3
anoph_min[which(anoph_min < 0)] <- 0










# --- old code -----
df_gam_binary <- calculate_prop_by_variable (
  data = meta_data_with_labels,
  var_col = "gam_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = anoph_max, 
  breaks_min = anoph_min)

summary(df_gam_binary$log_model)
prop_gam_prev <- make_plot(
  df_gam_binary$obs,
  meta_data_with_labels$gam_pw_district,
  "Proportion Anopheles gambiae",
  color ="#165262", pos_col =  "ONNV_pos"
)
print(prop_gam_prev)