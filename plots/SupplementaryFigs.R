# Libraries
library(forcats)
library(ggtext)

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
#anoph_max <- seq(0, 1, 0.1)
#anoph_min <- anoph_max - 0.3
anoph_min[which(anoph_min < 0)] <- 0

#aegmax <- seq(0,1,0.1)
#aegmin <- aegmax - 0.3
#aegmin[which(aegmin<0)] <- 0


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

summary(df_gam_binary$log_model) 

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


metrics_df <- data.frame(
  species = c(
    "Anopheles\nfunestus",
    "Anopheles\ngambiae",
    "Aedes\naegypti",
    "Aedes\nalbopictus"
  ),
  AIC = c(
    AIC(df_fun_binary$log_model),
    AIC(df_gam_binary$log_model),
    AIC(df_aegypti_binary$log_model),
    AIC(df_albopictus_binary$log_model)
  ),
  logLik = c(
    as.numeric(logLik(df_fun_binary$log_model)),
    as.numeric(logLik(df_gam_binary$log_model)),
    as.numeric(logLik(df_aegypti_binary$log_model)),
    as.numeric(logLik(df_albopictus_binary$log_model))
  )
)

metrics_df$k <- 1:nrow(metrics_df)

loglik_onnv_mosquito_vectors <- ggplot(metrics_df, aes(x = k, y = logLik)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = metrics_df$k,
    labels = metrics_df$species
  ) +
  labs(
    x = "Vector Species",
    y = "Log \nLikelihood"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, face = "italic"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(t = 10, r = 40, b = 10, l = 10, unit = "pt")
  )


anopheles_and_aedes_onnv <- 
  (prop_fun_prev | prop_gam_prev | prop_aeg_prev | prop_albo_prev)  / loglik_onnv_mosquito_vectors
print(anopheles_and_aedes_onnv)

ggsave("Results/NEW_supplementary_fig2.png", 
       plot = anopheles_and_aedes_onnv,
       width = 16, 
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
  breaks_max = w_fun$max, 
  breaks_min = w_fun$min)

# Gambiae
df_gam_chik <- calculate_prop_by_variable(
  data = meta_data_with_labels,
  var_col = "gam_pw_district", 
  positive_col = "CHIK_pos",
  breaks_max = w_gam$max, 
  breaks_min = w_gam$min)

# Aegypti
df_aegypti_chik <- calculate_prop_by_variable(
  data = meta_data_with_labels, 
  var_col = "aeg_pw_district", 
  positive_col = "CHIK_pos",
  breaks_max = w_aeg$max, 
  breaks_min = w_aeg$min)

  # Albopictus
df_albopictus_chik <- calculate_prop_by_variable(
  data = meta_data_with_labels, 
  var_col = "alb_pw_district", 
  positive_col = "CHIK_pos",
  breaks_max = w_alb$max, 
  breaks_min = w_alb$min)


# --- Plots
prop_fun_prev_chik <- make_plot_chik(
  df_fun_chik$obs,
  meta_data_with_labels$fun_pw_district,
  "Suitability<br>*Anopheles funestus*",
  color ="#023e8a", pos_col =  "CHIK_pos"
)

prop_gam_prev_chik <- make_plot_chik(
  df_gam_chik$obs,
  meta_data_with_labels$gam_pw_district,
  "Suitability<br>*Anopheles gambiae*",
  color ="#165262", pos_col =  "CHIK_pos"
)

prop_aeg_prev_chik <- make_plot_chik(
  df_aegypti_chik$obs,
  meta_data_with_labels$aeg_pw_district,
  "Suitability<br>*Aedes aegypti*",
  color ="#c1518b", pos_col =  "CHIK_pos"
)

prop_albo_prev_chik <- make_plot_chik(
  df_albopictus_chik$obs,
  meta_data_with_labels$alb_pw_district,
  "Suitability<br>*Aedes albopictus*",
  color = "#430726", pos_col =  "CHIK_pos"
)

anopheles_and_aedes_chik <- prop_fun_prev_chik | prop_gam_prev_chik | prop_aeg_prev_chik | prop_albo_prev_chik

ggsave("Results/NEW_supplementary_fig3.png", 
       plot = anopheles_and_aedes_chik,
       width = 20, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")




# Supplementary figure 4: An. Gambiae suitability and malaria 
# 2) Pf_Incidence_Rate vs ONNV seroprevalence 
pf <- terra::rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/clippedlayers-4/202508_Global_Pf_Incidence_Rate_CMR_2024.tiff')
pf_rate <- pf[[1]]  # band 1 = incidence rate
  
global(pf_rate, fun = c("min", "max"), na.rm = TRUE)

points_vect <- terra::vect(data_points)
meta_data_with_coords$pf_incidence <- terra::extract(
  pf_rate,
  points_vect
)[,2]

meta_data_with_labels$pf_incidence <- meta_data_with_coords$pf_incidence
summary(meta_data_with_labels$pf_incidence)

w_pf <- make_windows(meta_data_with_labels$pf_incidence)


onnv_pf_incidence <- calculate_prop_by_variable(
  data = meta_data_with_labels,
  var_col = "pf_incidence",
  positive_col = "ONNV_pos",
  breaks_max = w_pf$max,
  breaks_min = w_pf$min
)

summary(onnv_pf_incidence$log_model)

global(pf, fun = c("min", "max"), na.rm = TRUE)     # all 4 bands, so we can see if 1 vs 2-4 differ in scale
summary(meta_data_with_labels$pf_incidence)          # range of the actual extracted covariate
summary(onnv_pf_incidence$log_model)                 # or at least exp(coef(onnv_pf_incidence$log_model))

prop_pf_onnv <- make_plot_onnv(
  onnv_pf_incidence$obs,
  meta_data_with_labels$pf_incidence,
  "<br>*P. falciparum incidence* rate",  
  color ="#16622b", pos_col =  "ONNV_pos"
)
print(prop_pf_onnv)


ggsave("Results/NEW_supplementary_fig5.png", 
       plot = prop_pf_onnv,
       width = 6, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")



