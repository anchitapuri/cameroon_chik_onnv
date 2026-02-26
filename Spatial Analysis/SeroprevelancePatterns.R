
library(ggspatial)
library(cowplot)
library(purrr)

# --- Source functions
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/Spatial Analysis/Functions.R'))
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/MultiSeroModel/MultiSeroFunctions.R'))

cameroon_districts <- ne_states(country = "Cameroon", returnclass = "sf")

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


chik_infections <-ggplot() +
  geom_sf(
    data = cam_shapefile_districts_merged,
    fill = "white",
    colour = "black",
    linewidth = 0.3
) +
  geom_sf(
    data = chik_centroids |> dplyr::filter(n > 0),
    colour = '#065a82', 
    aes(size = n),
    alpha = 0.85
  )  + 
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

# How many individuals per bin?
table(cut(meta_data_with_labels$fun_pw_district, 
          breaks = unique(c(anoph_min, anoph_max)), 
          include.lowest = TRUE))

df_fun <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels,
  var_col = "fun_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a', # a = ONNV
  breaks_max = anoph_max, 
  breaks_min = anoph_min)
df_fun$log_model

cor(df_fun$x, df_fun$y, use = "complete.obs", method = "spearman")


# Gambiae
df_gam <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels, 
  var_col = "gam_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a',  # a = ONNV
  breaks_max = anoph_max, 
  breaks_min = anoph_min)

df_obs_gam <- df_gam$obs



# supplementary figure 
# plot label 
make_label <- function(model, pathogen, vector) {
  beta <- round(coef(model)[2], 3)
  ci   <- round(confint(model)[2, ], 3)
  df$strip_label <- paste0(β = ", beta",
                           " (", ci[1], ", ", ci[2], ")")
}


label_fun  <- make_label(df_fun$log_model,   "ONNV", "An. funestus")
label_gam  <- make_label(df_gam$log_model,    "ONNV", "An. gambiae")


# --- Plots 
# Funestus plot
prop_fun_prev <- ggplot(df_fun$obs, aes(x = x, y = y)) +
  geom_point(color = "#023e8a", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#c1518b") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0.10, 0.30), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05)) +
  labs(x = "Proportion Anopheles funestus", y = "Proportion ONNV positive",
  title = label_fun) +
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
print(prop_fun_prev)


# Gambiae plot
prop_gam_prev <- ggplot(df_obs_gam, aes(x = x, y = y)) +
  geom_point(color = "#165262", size = 5)  +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#165262") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0.10, 0.25), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05))+
  labs(x = "Proportion Anopheles gambiae", y = "Proportion ONNV positive",
  title = label_gam) +
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
print(prop_gam_prev)


(prop_fun_prev + prop_gam_prev)

# --- Save Figure 1a
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig4c_new.png", 
       plot = prop_gam_prev,
       width = 10, 
       height = 10, 
       units = "in", 
       dpi = 300,
       bg = "white")


# Aedes with ONNV 

# ----onnv correlation with aedes 
aegmax <- seq(0,1,0.1)
aegmin <- aegmax - 0.5
aegmin[which(aegmin<0)] <- 0

# Aegypti
df_aegypti <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels,
  var_col = "aeg_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a',  # a = ONNV
  breaks_max = aegmax, 
  breaks_min = aegmin)

df_obs_aeg <- df_aegypti$obs

# Albopictus
df_albopictus <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels,
  var_col = "alb_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a',  # a = ONNV
  breaks_max = aegmax, 
  breaks_min = aegmin)

df_obs_albo <- df_albopictus$obs


label_aeg  <- make_label(df_aegypti$log_model,    "ONNV", "Ae. aegypti")
label_albo <- make_label(df_albopictus$log_model, "ONNV", "Ae. albopictus")

prop_aeg_prev <- ggplot(df_obs_aeg, aes(x = x, y = y)) +
  geom_point(color = "#c1518b", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#c1518b") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05))+
  labs(x = "Proportion Aedes Aegypti", y = "Proportion ONNV positive", 
      title = label_aeg) +
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


prop_albo_prev <- ggplot(df_obs_albo, aes(x = x, y = y)) +
  geom_point(color = "#c1518b", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#c1518b") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05))+
  labs(x = "Proportion Aedes Albopictus", y = "Proportion ONNV positive",
      title = label_albo) +
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

anopheles_and_aedes_onnv <- (prop_fun_prev + prop_gam_prev) / (prop_aeg_prev + prop_albo_prev)
print(anopheles_and_aedes_onnv)

#
#


# ---  CHIK 
df_fun_chik <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels,
  var_col = "fun_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'b', # b = CHIK
  breaks_max = anoph_max, 
  breaks_min = anoph_min)
df_fun$log_model



# Gambiae
df_gam_chik <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels, 
  var_col = "gam_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'b',  # b = CHIK
  breaks_max = anoph_max, 
  breaks_min = anoph_min)


# --- Plots 
# Funestus plot
chik_prop_fun_prev <- ggplot(df_fun_chik$obs, aes(x = x, y = y)) +
  geom_point(color = "#023e8a", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#023e8a") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0.10, 0.30), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05)) +
  labs(x = "Proportion Anopheles funestus", y = "Proportion CHIK positive") +
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
print(chik_prop_fun_prev)


# Gambiae plot
chik_prop_gam_prev <- ggplot(df_gam_chik$obs, aes(x = x, y = y)) +
  geom_point(color = "#165262", size = 5)  +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#165262") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0.10, 0.25), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05))+
  labs(x = "Proportion Anopheles gambiae", y = "Proportion CHIK positive") +
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
print(chik_prop_gam_prev)
















# Aegypti
df_aegypti_chik <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels, 
  var_col = "aeg_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'b',  # b = CHIK 
  breaks_max = aegmax, 
  breaks_min = aegmin)

chik_prop_aeg_prev <- ggplot(df_aegypti_chik$obs, aes(x = x, y = y)) +
  geom_point(color = "#c1518b", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#c1518b") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05))+
  labs(x = "Proportion Aedes Aegypti", y = "Proportion CHIK positive") +
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


# Albopictus
df_albopictus_chik <- calculate_prop_by_variable_multisero_probs(
  data = meta_data_with_labels, 
  var_col = "alb_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'b',  # b = CHIK 
  breaks_max = aegmax, 
  breaks_min = aegmin)

chik_prop_albo_prev <- ggplot(df_albopictus_chik$obs, aes(x = x, y = y)) +
  geom_point(color = "#c1518b", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#c1518b") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05))+
  labs(x = "Proportion Aedes Albopictus", y = "Proportion CHIK positive") +
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





# ---ONNV with population density 
# First create the density column in your data
onnv_results_pop_grid$data_filtered <- onnv_results_pop_grid$data_filtered %>%
  mutate(
    pop_density = Total_Population / area_km2,
    log_pop_density = log(pop_density + 1)  # +1 to avoid log(0) issues
  )

popdenmax <-seq(-1, 4.2, 0.1)
popdenmin <-popdenmax-2

df_pop <- calculate_prop_by_variable_multisero_probs(
  data = onnv_results_pop_grid$data_filtered,
  var_col = "log_pop_density",
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a',  # a = ONNV
  breaks_max = popdenmax,
  breaks_min = popdenmin
)

cor(df_pop$x, df_pop$y, use = "complete.obs", method = "pearson")
# Plot
prop_pop_prev <- ggplot(df_pop, aes(x = x, y = y)) +
  geom_point(color = "#8b145b", size = 4) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#8B6914") +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05)) +
  labs(x = "Log population density (per km²)", y = "Proportion ONNV positive") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 0.75,
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm")
  )
print(prop_pop_prev)
