
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



# --- Source functions
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/Spatial Analysis/Functions.R'))
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/MultiSeroModel/MultiSeroFunctions.R'))

cameroon_districts <- ne_states(country = "Cameroon", returnclass = "sf")
cameroon <- ne_countries(country = "Cameroon", returnclass = "sf")
# Import multisero fits 
fit <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/adapted_full_model_fits.rds')
# Import INLA model fits 
onnv_results_pop_grid <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/ONNV_INLAResults.rds')
# Extract chains 
chains <- fit$draws(format='df')
chains_df <- as.data.frame(chains)

meta_data_with_coords <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/meta_data_with_coords.rds')

# Read preprocessed data
preprocessed_data <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/preprocessed_data_full_model.rds')

# shapefile with district geometries 
cam_shapefile_districts_merged <- readRDS("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cam_shapefile_districts_merged.rds")


meta_data_with_labels <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/final_meta_data_with_labels.csv')
nrow(meta_data_with_labels)


nrow(meta_data_with_labels)
nrow(onnv_results_pop_grid$data_filtered)

# --- Plot age seroprevalence model fits
# prepare data for stan
age_prev_model_fits <- plot_age_seroprevalence_model_fits(
  year_intro = onnv_results_pop_grid$year,
  result = onnv_results_pop_grid, 
  data = onnv_results_pop_grid$data_filtered,
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = "a" # a == ONNV 
)
print(age_prev_model_fits)

stack_data <- inla.stack.data(onnv_results_pop_grid$stk.full)
idx_est <- inla.stack.index(onnv_results_pop_grid$stk.full, tag = "est")$data

head(stack_data$id[idx_est])
head(data_plot$id)


ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig3.png", 
       plot = age_prev_model_fits[[1]],    # swap for your actual plot object name
       width = 18, 
       height = 12, 
       units = "in", 
       dpi = 300,
       bg = "white")




chik_pos <- onnv_results_pop_grid$data_filtered |>
  dplyr::filter(CHIK_pos == 1)
chik_counts <- chik_pos |>
  count(district_lower, name = "n")
sum(chik_counts$n)


chik_districts_sf <- cam_shapefile_districts_merged |>
  dplyr::left_join(chik_counts,
                   by = c("shapefile_district_lower" = "district_lower"))
chik_districts_sf$n[(chik_districts_sf$n)] <- 0
nrow(chik_districts_sf)

chik_centroids <- chik_districts_sf |>
  st_centroid()

chik_centroids_jittered <- chik_centroids |>
  dplyr::filter(n > 0) |>
  dplyr::mutate(geometry = geometry + sf::st_sfc(lapply(seq_len(dplyr::n()), function(i) {
    sf::st_point(c(runif(1, -0.05, 0.05), runif(1, -0.05, 0.05)))
  })))
sf::st_crs(chik_centroids_jittered) <- sf::st_crs(chik_centroids)

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
    colour = '#065a82', 
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
    shape = 24,
    colour = "#1b3b6f",
    size = 10,
    stroke = 0.6
  ) +
  scale_fill_manual(
    values = c("Yaoundé" = "#db2e6e", "Douala" = "#6a4c93"),
    name = ""
  ) +
  annotation_scale(
    bar_cols = c("black", "white"),  # alternating black/white like the reference
    height = unit(0.2, "cm"),
    text_family = "sans", 
    text_cex = 1.5
  ) +
  theme(
    legend.position = c(0.37, 0.66), 
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 20),        
    legend.title = element_text(size = 20),   
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
  )
quartz()
print(chik_infections)

# proportion of chik by region 
chik_pos_by_district <- onnv_results_pop_grid$data_filtered |>
  dplyr::group_by(district_lower) |>  # replace 'district' with your actual district column name
  dplyr::summarise(
    chik_pos = sum(CHIK_pos == 1, na.rm = TRUE),
    total = n(),
    proportion = chik_pos / total
  ) |>
  dplyr::arrange(dplyr::desc(proportion))


ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig4a.png", 
       plot = chik_infections,    # swap for your actual plot object name
       width = 10, 
       height = 10, 
       units = "in", 
       dpi = 300,
       bg = "white")



# --- Proportion positive by mosquito distributions (only ONNV) (using multisero probs)
# Anopheles bins (adjust if needed based on your data range
anoph_max <- seq(0, 1, 0.1)
anoph_min <- anoph_max - 0.5
anoph_min[which(anoph_min < 0)] <- 0



# general plot format 
base_theme <- theme_classic() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 0.75,
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
  
  ylab <- paste0("Proportion ", gsub("_pos", "", pos_col), " positive")
  
  obs_clean <- df_obs[!is.nan(df_obs$x), ]
  hist_df <- data.frame(x = as.numeric(raw_data))
  hist_df <- hist_df[!is.na(hist_df$x), , drop = FALSE]
  hist_breaks <- seq(0, 1, length.out = 31)
  hist_counts <- hist(hist_df$x, breaks = hist_breaks, plot = FALSE)$counts
  max_count <- if (length(hist_counts) > 0) max(hist_counts, na.rm = TRUE) else 0
  
  plot1 <- ggplot(obs_clean, aes(x = x, y = y)) +
    geom_point(color = color, size = 4, alpha = 0.9) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0, color = color, alpha = 0.6, linewidth = 0.6) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous() +
    labs(x = xlab, y = ylab) +
    theme_classic() + base_theme
  
  plot2 <- ggplot(hist_df, aes(x = x)) +
    geom_histogram(fill = color, alpha = 0.5, bins = 30, color = NA) +
    labs(x = xlab, y = "Count") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, max_count)) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      aspect.ratio = 0.75,
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      axis.ticks.x = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      axis.ticks.length = unit(0.2, "cm")
    )
  
  inset_pos <- if (grepl("CHIK", pos_col)) {
    inset_element(plot2, 0.6, 0.6, 1, 1)   # top right
  } else {
    inset_element(plot2, 0.6, 0, 1, 0.4)   # bottom right
  }
  
  plot1 + inset_pos
}



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


# ----onnv correlation with aedes 
aegmax <- seq(0,1,0.1)
aegmin <- aegmax - 0.5
aegmin[which(aegmin<0)] <- 0


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
prop_fun_prev <- make_plot(
  df_fun_binary$obs,
  meta_data_with_labels$fun_pw_district,
  "Proportion Anopheles funestus",
  color ="#023e8a", pos_col =  "ONNV_pos"
)


prop_gam_prev <- make_plot(
  df_gam_binary$obs,
  meta_data_with_labels$gam_pw_district,
  "Proportion Anopheles gambiae",
  color ="#165262", pos_col =  "ONNV_pos"
)
quartz()
print(prop_gam_prev)

prop_aeg_prev <- make_plot(
  df_aegypti_binary$obs,
  meta_data_with_labels$aeg_pw_district,
  "Proportion Aedes aegypti",
  color ="#c1518b", pos_col =  "ONNV_pos"
)

prop_albo_prev <- make_plot(
  df_albopictus_binary$obs,
  meta_data_with_labels$alb_pw_district,
  "Proportion Aedes albopictus",
  color = "#430726", pos_col =  "ONNV_pos"
)

anopheles_and_aedes_onnv <- 
  (prop_fun_prev + prop_gam_prev) /
  (prop_aeg_prev + prop_albo_prev)


dfun <- function(model) {
  summary(model)$dispersion
}

# compare models using AIC and log likelihood
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

metrics_df$delta_AIC <- metrics_df$AIC - min(metrics_df$AIC)
metrics_df <- metrics_df[order(metrics_df$logLik), ]
metrics_df$k <- 1:nrow(metrics_df)  

loglik_plot <- ggplot(metrics_df,
                      aes(x = k, y = logLik)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = metrics_df$k, labels = metrics_df$species) +
  labs(
    x = "Species",
    y = "Log Likelihood"
  ) +
  theme_minimal()  +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 0.75,
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(t = 10, r = 40, b = 10, l = 10, unit = "pt")
  )

  
quartz()
print(loglik_plot)


# beta 
print(summary(df_fun_binary$log_model))
print(summary(df_gam_binary$log_model))
print(summary(df_aegypti_binary$log_model))
print(summary(df_albopictus_binary$log_model))

# --- Save Figures
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig4c_new.png", 
       plot = prop_gam_prev,
       width = 8, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")

ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/sup_Fig1.png", 
       plot = anopheles_and_aedes_onnv,
       width = 16, 
       height = 11, 
       units = "in", 
       dpi = 300,
       bg = "white")


ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/sup_Fig2.png", 
       plot = loglik_plot,
       width = 8, 
       height = 6,
       units = "in", 
       dpi = 300,
       bg = "white")




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
prop_fun_prev_chik <- make_plot(
  df_fun_chik$obs,
  meta_data_with_labels$fun_pw_district,
  "Proportion Anopheles funestus",
  color ="#023e8a", pos_col =  "CHIK_pos"
)


prop_gam_prev_chik <- make_plot(
  df_gam_chik$obs,
  meta_data_with_labels$gam_pw_district,
  "Proportion Anopheles gambiae",
  color ="#165262", pos_col =  "CHIK_pos"
)


prop_aeg_prev_chik <- make_plot(
  df_aegypti_chik$obs,
  meta_data_with_labels$aeg_pw_district,
  "Proportion Aedes aegypti",
  color ="#c1518b", pos_col =  "CHIK_pos"
)

prop_albo_prev_chik <- make_plot(
  df_albopictus_chik$obs,
  meta_data_with_labels$alb_pw_district,
  "Proportion Aedes albopictus",
  color = "#430726", pos_col =  "CHIK_pos"
)

anopheles_and_aedes_chik <- 
  (prop_fun_prev_chik + prop_gam_prev_chik) /
  (prop_aeg_prev_chik + prop_albo_prev_chik)


ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/sup_Fig3.png", 
       plot = anopheles_and_aedes_chik,
       width = 16, 
       height = 11, 
       units = "in", 
       dpi = 300,
       bg = "white")


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





# --- Multisero probs with mosquito distributions 

# Funestus
df_fun <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels,
  var_col = "fun_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a', # a = ONNV
  breaks_max = anoph_max, 
  breaks_min = anoph_min)

# Gambiae
df_gam <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels, 
  var_col = "gam_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a',  # a = ONNV
  breaks_max = anoph_max, 
  breaks_min = anoph_min)

# Aegypti
df_aegypti <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels,
  var_col = "aeg_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a',  # a = ONNV
  breaks_max = aegmax, 
  breaks_min = aegmin)
# Albopictus
df_albopictus <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels,
  var_col = "alb_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a',  # a = ONNV
  breaks_max = aegmax, 
  breaks_min = aegmin)



