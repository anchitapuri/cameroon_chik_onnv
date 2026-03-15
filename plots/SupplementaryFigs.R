# Libraries
library(forcats)
library(ggtext)
# --- Source functions
source(here('R/Functions.R'))
source(here('R/MultiSeroFunctions.R'))


# --- Plot functions 
make_plot_onnv <- function(df_obs, raw_data, xlab, color, pos_col = "ONNV_pos") {

  ylab <- NULL

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
    labs(x = NULL, y = NULL) +
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


make_plot_chik <- function(df_obs, raw_data, xlab, color, pos_col = "ONNV_pos") {

  ylab <- NULL

  obs_clean <- df_obs[!is.nan(df_obs$x), ]

  # truncate CIs to [0, 0.5]
  obs_clean$ymin <- pmax(obs_clean$ymin, 0)
  obs_clean$ymax <- pmin(obs_clean$ymax, 0.1)

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
    labs(x = NULL, y = NULL) +
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
    coord_cartesian(ylim = c(0, 0.05), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, 0.05, 0.01)) +
    labs(x = xlab, y = ylab) +
    base_theme +
    theme(
      plot.margin = margin(t = 10, r = 14, b = 12, l = 14)
    )

  plot_hist / plot_scatter +
    patchwork::plot_layout(heights = c(2, 4))
  
}

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


# --- Supplementary figure 1: log-likelihood comparison of multisero models
# Read Multisero comparison log likelihoods
lli <- read.csv(here("Results/loglik_model_comparison.csv"))

lli <- lli |>
  dplyr::mutate(model = fct_reorder(model, med))  # order by descending log-likelihood

multisero_loglik_plot <- ggplot(lli, aes(x = model, y = med, group = 1)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_discrete(labels = c(
  "CHIK_only_model" = "CHIK only model",
  "ONNV_only_model" = "ONNV only model",
  "Full_Model"      = "ONNV+CHIK model"
  )) +
  labs(
    x = "Model",
    y = "Log Likelihood"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 0.75,
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 18,hjust = 0.5),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(t = 10, r = 40, b = 10, l = 10, unit = "pt")
  )
print(multisero_loglik_plot)

ggsave(here("Results/supplementary_fig1.png"), 
       plot = multisero_loglik_plot,
       width = 10, 
       height = 7,
       units = "in", 
       dpi = 300,
       bg = "white")



# --- Supplementary Figure 2: ONNV seroprevalence with all vecotrs 
# Anopheles and Aedes bins (adjust if needed based on your data range
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


# --- Plots
prop_fun_prev <- make_plot_onnv(
  df_fun_binary$obs,
  meta_data_with_labels$fun_pw_district,
  "Proportion \nAnopheles funestus",
  color ="#023e8a", pos_col =  "ONNV_pos"
)

prop_gam_prev <- make_plot_onnv(
  df_gam_binary$obs,
  meta_data_with_labels$gam_pw_district,
  "Proportion \nAnopheles gambiae",
  color ="#165262", pos_col =  "ONNV_pos"
)

prop_aeg_prev <- make_plot_onnv(
  df_aegypti_binary$obs,
  meta_data_with_labels$aeg_pw_district,
  "Proportion \nAedes aegypti",
  color ="#c1518b", pos_col =  "ONNV_pos"
)

prop_albo_prev <- make_plot_onnv(
  df_albopictus_binary$obs,
  meta_data_with_labels$alb_pw_district,
  "Proportion \nAedes albopictus",
  color = "#430726", pos_col =  "ONNV_pos"
)

# Model Comparison Plot using AIC and log likelihood
metrics_df <- data.frame(
  species = c(
    "Anopheles \n funestus",
    "Anopheles \n gambiae",
    "Aedes \n aegypti",
    "Aedes \n albopictus"
  ),
  AIC = c(
    AIC(df_fun_binary$log_model),
    AIC(df_gam_binary$log_model),
    AIC(df_aegypti_binary$log_model),
    AIC(df_albopictus_binary$log_model)),
  
    logLik = c(
    as.numeric(logLik(df_fun_binary$log_model)),
    as.numeric(logLik(df_gam_binary$log_model)),
    as.numeric(logLik(df_aegypti_binary$log_model)),
    as.numeric(logLik(df_albopictus_binary$log_model))
  )
)

#metrics_df <- metrics_df[order(metrics_df$logLik), ]
metrics_df$k <- 1:nrow(metrics_df)  

loglik_onnv_mosquito_vectors <- ggplot(metrics_df,
                      aes(x = k, y = logLik)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_continuous(
  breaks = metrics_df$k,
  labels = paste0("*", metrics_df$species, "*")
  )+
  labs(
    x = "Vector Species",
    y = "Log \nLikelihood"
  ) +
  theme_minimal()  + 
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_markdown(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(t = 10, r = 40, b = 10, l = 10, unit = "pt")
  )


anopheles_and_aedes_onnv <- 
  (prop_fun_prev | prop_gam_prev | prop_aeg_prev | prop_albo_prev)  / loglik_onnv_mosquito_vectors


ggsave("Results/supplementary_fig2.png", 
       plot = anopheles_and_aedes_onnv,
       width = 18, 
       height = 10, 
       units = "in", 
       dpi = 300,
       bg = "white")





# --- Supplementary figure 3: CHIK seroprevalence with all vecotrs 
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


# --- Plots
prop_fun_prev_chik <- make_plot_chik(
  df_fun_chik$obs,
  meta_data_with_labels$fun_pw_district,
  "Proportion \nAnopheles funestus",
  color ="#023e8a", pos_col =  "CHIK_pos"
)

prop_gam_prev_chik <- make_plot_chik(
  df_gam_chik$obs,
  meta_data_with_labels$gam_pw_district,
  "Proportion \nAnopheles gambiae",
  color ="#165262", pos_col =  "CHIK_pos"
)

prop_aeg_prev_chik <- make_plot_chik(
  df_aegypti_chik$obs,
  meta_data_with_labels$aeg_pw_district,
  "Proportion \nAedes aegypti",
  color ="#c1518b", pos_col =  "CHIK_pos"
)

prop_albo_prev_chik <- make_plot_chik(
  df_albopictus_chik$obs,
  meta_data_with_labels$alb_pw_district,
  "Proportion \nAedes albopictus",
  color = "#430726", pos_col =  "CHIK_pos"
)

anopheles_and_aedes_chik <- prop_fun_prev_chik | prop_gam_prev_chik | prop_aeg_prev_chik | prop_albo_prev_chik

ggsave("Results/supplementary_fig3.png", 
       plot = anopheles_and_aedes_chik,
       width = 20, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")






