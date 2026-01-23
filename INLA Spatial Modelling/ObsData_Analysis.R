# 1) Extract population weighted Latitude and Longitude (from district centroids)
# 2) Patterns of observed seroprevelance (aggregated by district, mosquito proportion and age) 

# --- Source functions
source(here('SpatialModellingFunctions.R'))


# Load data with coords (post preprocessing) 
meta_data_with_coords <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/meta_data_with_coords.rds')
nrow(meta_data_with_coords)

# Load population rasters
cam_pop <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_ppp_2020_UNadj.tif")
cam_pop_den <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_pd_2020_1km_UNadj.tif")
plot(cam_pop_den)
plot(cam_pop)

# Load mosquito proportions in Cameroon
aegypti <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/Aedes_maps_public/aegypti.tif')
albopictus <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/Aedes_maps_public/albopictus.tif')
anopheles_funestus <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/2010_Anopheles_funestus_CMR.tiff')
anopheles_gambiae <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/2010_Anopheles_gambiae_ss_CMR.tiff')

# Load population statified by age group data 
cameroon_age_2025 <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/CameroonAge2025.csv')
cameroon_age_2025 <- cameroon_age_2025 %>%
  mutate(total = M + F)

# --- Population weighted centroids ----
# Prepare spatial data
sf_meta_data_with_coords <- st_as_sf(meta_data_with_coords) %>%
  filter(!st_is_empty(geometry)) %>%
  filter(!is.na(st_dimension(geometry)))
nrow(sf_meta_data_with_coords)

# Create district polygons with area
districts <- sf_meta_data_with_coords %>%
  st_transform(32633) %>%  # UTM Zone 33N
  group_by(district_lower) %>%
  summarise(
    geometry = st_union(geometry),
    .groups = "drop"
  ) %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_transform(crs = crs(cam_pop)) %>%
  mutate(district_id = dplyr::row_number())


# Create population stack
area_rast_pop <- terra::cellSize(cam_pop, unit = "km")
cam_stack_pop <- c(cam_pop, area_rast_pop)
names(cam_stack_pop) <- c("population_per_pixel", "cell_area")

# Function to  Calculate population-weighted centroid
pop_weighted_centroid <- function(df, ...) {
  pop_per_cell <- df$population_per_pixel * df$coverage_fraction
  total_pop <- sum(pop_per_cell, na.rm = TRUE)
  
  if (is.na(total_pop) || total_pop == 0) {
    return(data.frame(
      Longitude = NA_real_,
      Latitude = NA_real_,
      Total_Population = 0
    ))
  }
  
  cx <- sum(df$x * pop_per_cell, na.rm = TRUE) / total_pop
  cy <- sum(df$y * pop_per_cell, na.rm = TRUE) / total_pop
  
  data.frame(
    Longitude = cx,
    Latitude = cy,
    Total_Population = total_pop
  )
}

# Calculate centroids
centroids_df <- exact_extract(
  cam_stack_pop,
  districts,
  fun = pop_weighted_centroid,
  include_xy = TRUE,
  summarize_df = TRUE,
  progress = TRUE
) %>%
  dplyr::bind_cols(districts %>% st_drop_geometry() %>% dplyr::select(district_lower, district_id)) %>%
  left_join(
    districts %>% st_drop_geometry() %>% select(district_lower, area_km2),  # Drop geometry here
    by = "district_lower"
  )

# valide: sum(total population) == approx population of Cameroon
sum(centroids_df$Total_Population)


# Resample mosquito maps to population grid
# -- Aedes aegypti 
aegypti_cam_pop <- resample(aegypti, cam_pop, method = "bilinear")
stack_pop_aeg <- c(cam_pop, area_rast_pop, aegypti_cam_pop)
names(stack_pop_aeg) <- c("population_per_pixel", "cell_area", "aeg")

# -- Aedes albopictus
albopictus_cam_pop <- resample(albopictus, cam_pop, method = "bilinear")
stack_pop_alb <- c(cam_pop, area_rast_pop, albopictus_cam_pop)
names(stack_pop_alb) <- c("population_per_pixel", "cell_area", "alb")

# --- Anopheles funestus 
funestus_cam_pop <- resample(anopheles_funestus, cam_pop, method = "bilinear")
stack_pop_fun <- c(cam_pop, area_rast_pop, funestus_cam_pop)
names(stack_pop_fun) <- c("population_per_pixel", "cell_area", "fun")

# --- Anopheles funestus 
gambiae_cam_pop <- resample(anopheles_gambiae, cam_pop, method = "bilinear")
stack_pop_gam <- c(cam_pop, area_rast_pop, gambiae_cam_pop)
names(stack_pop_gam) <- c("population_per_pixel", "cell_area", "gam")



# Function for population-weighted mean of mosquito distribution per district
pop_weighted_mean_fun <- function(df, var_name, ...) {
  pop_count <- df$population_per_pixel * df$coverage_fraction
  total_pop <- sum(pop_count, na.rm = TRUE)
  
  if (is.na(total_pop) || total_pop == 0) {
    return(setNames(data.frame(NA_real_), paste0(var_name, "_pw_district")))
  }
  
  weighted_mean <- sum(df[[var_name]] * pop_count, na.rm = TRUE) / total_pop
  setNames(data.frame(weighted_mean), paste0(var_name, "_pw_district"))
}

# Calculate mosquito distributions
aeg_pw_df <- exact_extract(
  stack_pop_aeg,
  districts,
  fun = function(df, ...) pop_weighted_mean_fun(df, "aeg", ...),
  summarize_df = TRUE,
  include_xy = TRUE,
  progress = TRUE
) %>%
  mutate(district_id = districts$district_id, district_lower = districts$district_lower)

alb_pw_df <- exact_extract(
  stack_pop_alb,
  districts,
  fun = function(df, ...) pop_weighted_mean_fun(df, "alb", ...),
  summarize_df = TRUE,
  include_xy = TRUE,
  progress = TRUE
) %>%
  mutate(district_id = districts$district_id, district_lower = districts$district_lower)

# Calculate Anopheles distributions (population-weighted)
fun_pw_df <- exact_extract(
  stack_pop_fun,
  districts,
  fun = function(df, ...) pop_weighted_mean_fun(df, "fun", ...),
  summarize_df = TRUE,
  include_xy = TRUE,
  progress = TRUE
) %>%
  mutate(district_id = districts$district_id, district_lower = districts$district_lower)

gam_pw_df <- exact_extract(
  stack_pop_gam,
  districts,
  fun = function(df, ...) pop_weighted_mean_fun(df, "gam", ...),
  summarize_df = TRUE,
  include_xy = TRUE,
  progress = TRUE
) %>%
  mutate(district_id = districts$district_id, district_lower = districts$district_lower)


# --- Combine all into one dataframe 
sf_meta_data_with_coords_pw <- sf_meta_data_with_coords %>%
  left_join(
    centroids_df %>% 
      st_drop_geometry() %>%  # Drop geometry to avoid duplication
      select(-district_id), 
    by = "district_lower"
  ) %>%
  left_join(aeg_pw_df %>% select(district_lower, aeg_pw_district), by = "district_lower") %>%
  left_join(alb_pw_df %>% select(district_lower, alb_pw_district), by = "district_lower") %>%
  left_join(fun_pw_df %>% select(district_lower, fun_pw_district), by = "district_lower") %>%
  left_join(gam_pw_df %>% select(district_lower, gam_pw_district), by = "district_lower")


# --- Validate: Plot Population weighted vs unweighted coords 
sf_use_s2(FALSE)
unweighted_centroids <- st_centroid(sf_meta_data_with_coords)
unweighted_centroids <- st_coordinates(unweighted_centroids)

ggplot() +
  # Plot the polygons
  geom_sf(data = sf_meta_data_with_coords, fill = "white", color = "black", linewidth = 0.3) +
  # Add the weighted centroids as points using lon/lat coordinates
  geom_point(data = centroids_df, aes(x = Longitude, y = Latitude), 
             color = "red", size = 2, shape = 16) +
  geom_point(data = unweighted_centroids, aes(x = X, y = Y), 
             color = "darkblue", size = 2, shape = 16) +
  
  theme_minimal() +
  labs(title = "Cameroon: Population-Weighted Centroids",
       x = "Longitude", y = "Latitude")




# --- Seroprevalence Patterns 
cameroon_data <- sf_meta_data_with_coords_pw 


# ---  Add UTM coordinates
sp_vill <- SpatialPoints(cbind(cameroon_data$Longitude, cameroon_data$Latitude))
points_to_extract <- terra::vect(sp_vill)

data_points <- cameroon_data %>%
  st_drop_geometry() %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
data_utm <- st_transform(data_points, crs = 32633)
coords_utm <- st_coordinates(data_utm) / 1000  # Convert to km
colnames(coords_utm) <- c("Easting", "Northing")

# Add Easting and Northing to dataframe
cameroon_data$Easting <- coords_utm[, "Easting"]
cameroon_data$Northing <- coords_utm[, "Northing"]

# Extract population density values at district lat / long 
density_values <- terra::extract(cam_pop_den, points_to_extract)
cameroon_data$pop_density <- density_values[, 2]
cameroon_data$logpopden <- log(cameroon_data$pop_density, 10)


N <- nrow(cameroon_data)


# --- FUNCTIONS ----

# --- SEROPOSITIVITY BY DISTRICT
calculate_district_prevalence <- function(data, positive_col) {
  data %>%
    group_by(district_lower) %>%
    summarise(
      mean_positive = mean(.data[[positive_col]]),
      n_samples = dplyr::n(),
      n_positive = sum(.data[[positive_col]] == 1, na.rm = TRUE),
      .groups = "drop"
    )
}

# ---ODDS RATIO ANALYSIS (FUNCTION)
calculate_odds_ratio <- function(data, positive_col, n_iter = 100) {
  
  N <- nrow(data)
  
  # Distance matrix
  dmat <- as.matrix(dist(cbind(data$Easting, data$Northing)))
  diag(dmat) <- NA
  
  # Get positive indices
  ind <- which(data[[positive_col]] == 1)
  dmat2_neg <- dmat2_pos <- dmat2 <- dmat[ind, ]
  dmat2_pos[, -ind] <- NA
  dmat2_neg[, ind] <- NA
  
  # Distance bins
  distmax <- seq(0, 500, 20)
  distmin <- distmax - 50
  distmin[which(distmin < 0)] <- 10
  distmin[1] <- 0
  distmid <- (distmax + distmin) / 2
  
  # Calculate OR
  counts_pos_1 <- cumsum(hist(dmat2_pos, breaks = c(distmax, 1e10), plot = FALSE)$counts)
  counts_neg_1 <- cumsum(hist(dmat2_neg, breaks = c(distmax, 1e10), plot = FALSE)$counts)
  counts_pos_2 <- cumsum(hist(dmat2_pos, breaks = c(distmin, 1e10), plot = FALSE)$counts)
  counts_neg_2 <- cumsum(hist(dmat2_neg, breaks = c(distmin, 1e10), plot = FALSE)$counts)
  
  counts_pos_win <- counts_pos_1 - counts_pos_2
  counts_neg_win <- counts_neg_1 - counts_neg_2
  
  allPos <- sum(dmat2_pos >= 0, na.rm = TRUE)
  allNeg <- sum(dmat2_neg >= 0, na.rm = TRUE)
  overallPropPos <- allPos / allNeg
  
  OR <- (counts_pos_win / counts_neg_win) / overallPropPos
  
  # Bootstrap
  bs_out <- matrix(NaN, length(distmax), n_iter)
  for (i in 1:n_iter) {
    a <- sample(N, replace = TRUE)
    ind_bs <- which(data[[positive_col]][a] == 1)
    dmat_bs <- dmat[a, a]
    
    dmat2_neg_bs <- dmat2_pos_bs <- dmat_bs[ind_bs, ]
    dmat2_pos_bs[, -ind_bs] <- NA
    dmat2_neg_bs[, ind_bs] <- NA
    
    counts_pos_1 <- cumsum(hist(dmat2_pos_bs, breaks = c(distmax, 1e10), plot = FALSE)$counts)
    counts_neg_1 <- cumsum(hist(dmat2_neg_bs, breaks = c(distmax, 1e10), plot = FALSE)$counts)
    counts_pos_2 <- cumsum(hist(dmat2_pos_bs, breaks = c(distmin, 1e10), plot = FALSE)$counts)
    counts_neg_2 <- cumsum(hist(dmat2_neg_bs, breaks = c(distmin, 1e10), plot = FALSE)$counts)
    
    counts_pos_win <- counts_pos_1 - counts_pos_2
    counts_neg_win <- counts_neg_1 - counts_neg_2
    
    allPos <- sum(dmat2_pos_bs >= 0, na.rm = TRUE)
    allNeg <- sum(dmat2_neg_bs >= 0, na.rm = TRUE)
    overallPropPos <- allPos / allNeg
    
    bs_out[, i] <- (counts_pos_win / counts_neg_win) / overallPropPos
  }
  
  CI <- apply(bs_out, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Return results
  list(
    distance = distmid,
    OR = OR,
    CI_lower = CI[1, ],
    CI_upper = CI[2, ],
    plot_data = data.frame(
      distance = distmid,
      odds_ratio = OR,
      ci_lower = CI[1, ],
      ci_upper = CI[2, ]
    )
  )
}

# --- PROP POS BY MOSQUITO + LOG POP DENSITY 
calculate_prop_by_variable <- function(data, var_col, positive_col, breaks_max, breaks_min) {
  var_mid <- rep(NaN, length(breaks_max))
  prop_pos <- matrix(NaN, length(breaks_max), 3)
  
  for (i in 1:length(breaks_max)) {
    tmp <- which(data[[var_col]] < breaks_max[i] & data[[var_col]] >= breaks_min[i])
    if (length(tmp) > 5) {
      prop_pos[i, 1] <- mean(data[[positive_col]][tmp], na.rm = TRUE)
      a <- prop.test(sum(data[[positive_col]][tmp]), length(tmp))
      prop_pos[i, 2:3] <- a$conf.int
      var_mid[i] <- mean(data[[var_col]][tmp], na.rm = TRUE)
    }
  }
  
  data.frame(
    x = var_mid,
    y = prop_pos[, 1],
    ymin = prop_pos[, 2],
    ymax = prop_pos[, 3]
  )
}

table(cameroon_data$CHIK_pos)

# ---- RUN FUNCTIONS + PLOTS ----

# 1) --- Prevelance by district 
# CHIK: very few infections (n=28) 
# Check where CHIK infection fall within cameroon (coloured by Aedes Aegypti proportion) 
chik_pos <- cameroon_data |>
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

ggplot() +
  geom_sf(
    data = cameroon_data,
    fill = "grey95",
    colour = "grey40",
    linewidth = 0.3
  ) +
  geom_sf(
    data = st_centroid(chik_pos_sf),
    aes(colour = aeg_pw_district),
    size = 3,
    alpha = 0.8
  ) +
  scale_colour_gradient(
    low = "yellow",
    high = "red",
    name = "Aegypti"
  ) +
  geom_sf(
    data = cities_sf,
    aes(fill = city),
    shape = 21,
    colour = "black",
    size = 4,
    stroke = 0.6
  ) +
  scale_fill_manual(
    values = c("Yaoundé" = "lightblue", "Douala" = "darkblue"),
    name = "Cities"
  ) +
  theme_minimal() +
  labs(
    title = "CHIK-positive samples in Cameroon",
    subtitle = "Point color shows Aedes aegypti levels"
  )


# ONNV: Prevelance by district 
onnv_by_district <- calculate_district_prevalence(cameroon_data, "ONNV_pos")

ggplot(onnv_by_district) +
  geom_sf(aes(fill = mean_positive), colour = "grey20", linewidth = 0.3) +
  scale_fill_viridis_c(
    option = "magma",
    name = "ONNV seroprevelance \n by district",
    limits = c(0, max(onnv_by_district$mean_positive, na.rm = TRUE))
  ) +
  theme_minimal()




# --- 2) ODDS RATIO 
# Calculate for CHIK
or_chik <- calculate_odds_ratio(cameroon_data, "CHIK_pos", n_iter = 100)

# Calculate for ONNV
or_onnv <- calculate_odds_ratio(cameroon_data, "ONNV_pos", n_iter = 100)


# Plot CHIK
ggplot(or_chik$plot_data, aes(x = distance, y = odds_ratio)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10() +
  labs(x = "Distance (km)", y = "Odds Ratio", title = "CHIKV") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )

# Plot ONNV
ggplot(or_onnv$plot_data, aes(x = distance, y = odds_ratio)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10() +
  labs(x = "Distance (km)", y = "Odds Ratio", title = "ONNV") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15)
  )




# --- 3) Proportion positive by mosquito distributions (only ONNV) 

# Anopheles bins (adjust if needed based on your data range)
anoph_max <- seq(0, 1, 0.1)
anoph_min <- anoph_max - 0.5
anoph_min[which(anoph_min < 0)] <- 0


# Funestus
df_fun <- calculate_prop_by_variable(
  cameroon_data, "fun_pw_district", "ONNV_pos", anoph_max, anoph_min
)
df_fun$species <- "Funestus"

# Gambiae
df_gam <- calculate_prop_by_variable(
  cameroon_data, "gam_pw_district", "ONNV_pos", anoph_max, anoph_min
)
df_gam$species <- "Gambiae"

df_anopheles <- rbind(df_fun, df_gam)


# Funestus plot
prop_fun_prev <- ggplot(df_fun, aes(x = x, y = y)) +
  geom_point(color = "purple", size = 2) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "purple") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(x = "Proportion Anopheles funestus", y = "Proportion ONNV positive") +
  theme_classic()

# Gambiae plot
prop_gam_prev <- ggplot(df_gam, aes(x = x, y = y)) +
  geom_point(color = "orange", size = 2) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, color = "orange") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(x = "Proportion Anopheles gambiae", y = "Proportion ONNV positive") +
  theme_classic()

# Combined plot
(prop_fun_prev + prop_gam_prev)



# save RDS 
saveRDS(cameroon_data, '/Users/ap2488/Desktop/Cameroon_Analysis_2025/16thJan2026_cameroon_data.rds')
