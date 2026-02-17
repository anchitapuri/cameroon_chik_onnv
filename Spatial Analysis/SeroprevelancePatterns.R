
library(ggspatial)
library(cowplot)
library(purrr)

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
chik_pos <- onnv_results_pop_grid$data_filtered |>
  dplyr::filter(CHIK_pos == 1)
nrow(chik_pos)

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
    colour = '#00153e', 
    size = 8,
    alpha = 0.85
  )  +
  geom_sf(
    data = cities_sf,
    aes(fill = city),
    shape = 24,
    colour = "#f1f3f4",
    size = 8,
    stroke = 0.6
  ) +
  scale_fill_manual(
    values = c("Yaoundé" = "#db2e6e", "Douala" = "#db2e6e"),
    name = ""
  ) +
  annotation_scale(
    bar_cols = c("black", "white"),  # alternating black/white like the reference
    height = unit(0.2, "cm"),
    text_family = "sans", 
    text_cex = 1.5
  ) +
  theme(
    legend.position = c(0.1, 0.6), 
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


df_fun_old <- calculate_prop_by_variable(
  data = onnv_results_pop_grid$data_filtered, 
  var_col = "fun_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = anoph_max, 
  breaks_min = anoph_min)
df_fun_old$species <- "Funestus"
cor(df_fun_old$x, df_fun_old$y, use = "complete.obs", method = "pearson")

# Gambiae
df_gam_old <- calculate_prop_by_variable(
  data = onnv_results_pop_grid$data_filtered, 
  var_col = "gam_pw_district", 
  positive_col = "ONNV_pos",
  breaks_max = anoph_max, 
  breaks_min = anoph_min)
df_gam_old$species <- "Gambiae"
cor(df_gam_old$x, df_gam_old$y, use = "complete.obs", method = "pearson")



# Funestus plot
prop_fun_prev_old <- ggplot(df_fun_old, aes(x = x, y = y)) +
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
prop_gam_prev_old <- ggplot(df_gam_old, aes(x = x, y = y)) +
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
combined_mosquito_plots_old <- prop_fun_prev_old + prop_gam_prev_old
print(combined_mosquito_plots_old)

# --- Save Figure 1a
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig4c.png", 
       plot = combined_mosquito_plots,
       width = 13, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")





# --- USING MULTISERO MODEL PROBABILITIES
# Anopheles bins (adjust if needed based on your data range
preprocessed_data$data$y$model_row_id <- 1:nrow(preprocessed_data$data$y)
data_filtered <- preprocessed_data$data$y %>%
  filter(!is.na(AgeIn) & age != 0)
head(preprocessed_data$data$y)

range(onnv_results_pop_grid$data_filtered$model_row_id)
max(onnv_results_pop_grid$data_filtered$model_row_id)
nrow(preprocessed_data$data$y)

head(onnv_results_pop_grid$data_filtered$model_row_id)
head(colnames(chains_df)[grep("post_prob", colnames(chains_df))][1:10])

all.equal(
  preprocessed_data$data$y$ONNV_VLP_log,
  onnv_results_pop_grid$data_filtered$ONNV_VLP_log
)


anoph_max <- seq(0, 1, 0.1)
anoph_min <- anoph_max - 0.5
anoph_min[which(anoph_min < 0)] <- 0


df_fun <- calculate_prop_by_variable_multisero_probs(
  data = onnv_results_pop_grid$data_filtered, 
  var_col = "fun_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a', # a = ONNV
  breaks_max = anoph_max, 
  breaks_min = anoph_min)
df_fun$species <- "Funestus"



# Gambiae
df_gam <- calculate_prop_by_variable_multisero_probs(
  data = onnv_results_pop_grid$data_filtered, 
  var_col = "gam_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a',  # a = ONNV
  breaks_max = anoph_max, 
  breaks_min = anoph_min)
df_gam$species <- "Gambiae"
cor(df_gam$x, df_gam$y, use = "complete.obs", method = "pearson")
df_gam

dim(chains_df)
head(colnames(chains_df)[1:20])

nrow(preprocessed_data$data$y)
nrow(onnv_results_pop_grid$data_filtered)
# --- Plots 
# Funestus plot
prop_fun_prev <- ggplot(df_fun, aes(x = x, y = y)) +
  geom_point(color = "#c1518b", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#c1518b") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0.10, 0.30), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05)) +
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
print(prop_fun_prev)
# Gambiae plot
prop_gam_prev <- ggplot(df_gam, aes(x = x, y = y)) +
  geom_point(color = "#165262", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#165262") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0.10, 0.25), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05))+
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
print(prop_gam_prev)
(prop_fun_prev + prop_gam_prev)
combined_mosquito_plots_new <- (prop_fun_prev + prop_gam_prev)
# --- Save Figure 1a
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig4c_new.png", 
       plot = combined_mosquito_plots_new,
       width = 13, 
       height = 8, 
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
  data = onnv_results_pop_grid$data_filtered, 
  var_col = "aeg_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a',  # a = ONNV
  breaks_max = aegmax, 
  breaks_min = aegmin)
df_aegypti$species <- "Aegypti"
cor(df_aegypti$x, df_aegypti$y, use = "complete.obs", method = "pearson")

ggplot(df_aegypti, aes(x = x, y = y)) +
  geom_point(color = "#c1518b", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#c1518b") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05))+
  labs(x = "Proportion Aedes Aegypti", y = "Proportion ONNV positive") +
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
df_albopictus <- calculate_prop_by_variable_multisero_probs(
  data = onnv_results_pop_grid$data_filtered, 
  var_col = "alb_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'a',  # a = ONNV
  breaks_max = aegmax, 
  breaks_min = aegmin)
df_albopictus$species <- "Albopictus"
cor(df_albopictus$x, df_albopictus$y, use = "complete.obs", method = "pearson")


ggplot(df_albopictus, aes(x = x, y = y)) +
  geom_point(color = "#c1518b", size = 5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "#c1518b") +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
  scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 0.35, 0.05))+
  labs(x = "Proportion Aedes Albopictus", y = "Proportion ONNV positive") +
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




# --- Aedes with CHIK 
# Aegypti
df_aegypti_chik <- calculate_prop_by_variable_multisero_probs(
  data = onnv_results_pop_grid$data_filtered, 
  var_col = "aeg_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'b',  # b = CHIK 
  breaks_max = aegmax, 
  breaks_min = aegmin)
df_aegypti_chik$species <- "Aegypti"
cor(df_aegypti_chik$x, df_aegypti_chik$y, use = "complete.obs", method = "pearson")


# Albopictus
df_albopictus_chik <- calculate_prop_by_variable_multisero_probs(
  data = onnv_results_pop_grid$data_filtered, 
  var_col = "alb_pw_district", 
  chains_df = chains_df,
  infM = preprocessed_data$data$infM,
  pathogen_col = 'b',  # b = CHIK 
  breaks_max = aegmax, 
  breaks_min = aegmin)
df_albopictus_chik$species <- "Albopictus"
cor(df_albopictus_chik$x, df_albopictus_chik$y, use = "complete.obs", method = "pearson")






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
