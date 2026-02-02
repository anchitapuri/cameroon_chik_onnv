
# Import libraries
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(matrixStats)
library(stringr)
library(data.table)
library(dplyr)
library(scales)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(here)
library(terra)
library(exactextractr)
library(raster)

# read data files with labels 
meta_data <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/base_complete_MFI_meta.csv')
nrow(meta_data)


# Load shapefile
cam_shapefile_districts <- read_sf('/Users/ap2488/Desktop/Cameroon_Analysis_2025/S4_Cameroon_health_districts_files/Caedistricts179_region.shp')
# Second shapefile used (to find remaining mismatched districts)
cam_shapefile_districts2 <- read_sf('/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_admin_boundaries/cmr_admin3.shp')


# Load population rasters
cam_pop <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_ppp_2020_UNadj.tif")
cam_pop_den <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_pd_2020_1km_UNadj.tif")

# Load mosquito maps
aegypti <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/Aedes_maps_public/aegypti.tif')
albopictus <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/Aedes_maps_public/albopictus.tif')
anopheles_funestus <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/2010_Anopheles_funestus_CMR.tiff')
anopheles_gambiae <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/2010_Anopheles_gambiae_ss_CMR.tiff')

# Load population by gender / age data 
cameroon_age_2025 <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/CameroonAge2025.csv')
cameroon_age_2025 <- cameroon_age_2025 %>%
  mutate(total = M + F)

# ---1) Match district names in data with shapefiles to extract geometry for each district 
# Clean names + create lower case district column
cam_shapefile_districts$NAME2 <- gsub("DS_", "", cam_shapefile_districts$NAME2)
cam_shapefile_districts <- cam_shapefile_districts %>%
  mutate(shapefile_district_lower = tolower(NAME2))

# Check how many names and districts 
length(unique(cam_shapefile_districts$NAME2))
length(unique(cam_shapefile_districts$geometry))

# Find districts with multiple geometries
districts_with_multiple_geoms <- cam_shapefile_districts %>%
  group_by(NAME2) %>%
  mutate(n_geoms = n()) %>%
  filter(n_geoms > 1) %>%
  arrange(NAME2)

# manoka == 5 geometires - merge these 
manoka_merged <- cam_shapefile_districts %>%
  filter(shapefile_district_lower == "manoka") %>%
  st_make_valid() %>%
  summarise(
    NAME2 = first(NAME2),
    NAME1 = first(NAME1),
    COUNTRY3 = first(COUNTRY3),
    shapefile_district_lower = first(shapefile_district_lower),
    geometry = st_union(geometry)
  )
# Keep all other districts 
other_districts <- cam_shapefile_districts %>%
  filter(shapefile_district_lower != "manoka")
# Combine
cam_shapefile_districts_merged <- bind_rows(other_districts, manoka_merged)



# --- Meta Data ---
# Create lowercase district column for meta_data
meta_data <- meta_data %>%
  mutate(district_lower = tolower(DistrictOfresidence))
length(unique(meta_data$district_lower))


# --- Non-matching districts between shapefile 1 and data 
non_matching_meta <- meta_data %>%
  filter(!district_lower %in% cam_shapefile_districts_merged$shapefile_district_lower)

# Print unique non-matching district names
cat("\nNon-matching districts:\n")
print(unique(non_matching_meta$district_lower))
length(unique(non_matching_meta$district_lower))

# Count rows for each non-matching district
non_matching_counts <- non_matching_meta %>%
  count(district_lower) %>%
  arrange(desc(n))
cat("\nCount of non-matching districts:\n")
print(non_matching_counts)
cat("\nTotal non-matching rows:", nrow(non_matching_meta), "\n")


# --- Manual mapping ---- 
meta_data_cleaned <- meta_data %>%
  mutate(district_lower = case_when(
    district_lower == "njombe penja" ~ "njombe-penja",
    district_lower == "tchollire" ~"tcholire",
    district_lower == "cite verte" ~ "cite vert",
    district_lower == "malentouen" ~"malantouen",
    district_lower == "nkongsamba" ~ "nkonsamba",
    district_lower == "guidiguis" ~ "guidiguise",
    district_lower == "maroua 3" ~ "maroua rural",
    district_lower == "maroua 1" ~ "maroua urbain",
    district_lower == "kumba-north" ~ "kumba",
    district_lower == "bamenda 3" ~ "bamenda",
    district_lower == "nkongsamba" ~ "nkonsamba",
    district_lower == "garoua 1" ~ "garoua i",
    district_lower == "ngaoundal" ~ "ngaoundere rural",
    district_lower == "garoua urbain" ~ "garoua boulai",
    district_lower == "eyumodjock" ~ "eyumojock",
    district_lower == "garoua 2" ~ "garoua ii",
    district_lower == "ndikinimeki" ~ "ndikinimiki",
    district_lower == "mbandjock" ~ "mbanjock",
    district_lower == 'bandjoun' ~ "banjoun", 
    district_lower == 'bangangte' ~ "bangante",
    district_lower == 'bangourain'~ "bangorain",
    
    TRUE ~ district_lower
  ))


# Use hapefile 2 since there are still districts in data missing from shapefile 1
# Merge geometries in the second shapefile (in case it has duplicates too)
cam_shapefile_districts2 <- cam_shapefile_districts2 %>%
  mutate(shapefile_district_lower2 = tolower(adm3_name1))

cam_shapefile_districts2_merged <- cam_shapefile_districts2 %>%
  group_by(shapefile_district_lower2) %>%
  summarise(
    adm3_name = first(adm3_name),
    geometry = st_union(geometry),
    .groups = "drop"
  )
length(unique(cam_shapefile_districts2_merged$shapefile_district_lower2))
length(unique(cam_shapefile_districts2_merged$geometry))

# --- manually selected these districts that are in shapefile 2 and data
districts_new_shapefile <- c('belabo', 'belel', 'evodoula', 'fotokol',
                             'gazawa', 'goulfey', 'nguelemendouka', 'njombe-penja',
                             'oku')
rows_from_shapefile2 <- cam_shapefile_districts2 %>%
  filter(shapefile_district_lower2 %in% districts_new_shapefile)
rows_to_add <- rows_from_shapefile2 %>%
  st_transform(st_crs(cam_shapefile_districts)) %>%
  transmute(
    shapefile_district_lower = shapefile_district_lower2,
    geometry = geometry
    # Map other columns as needed
  )
cam_shapefile_districts_merged <- cam_shapefile_districts_merged %>%
  bind_rows(rows_to_add)
cam_shapefile_districts_merged <- cam_shapefile_districts_merged %>%
  st_make_valid()


# One geometry per district 
cam_shapefile_districts_unique <- cam_shapefile_districts_merged %>%
  group_by(shapefile_district_lower) %>%
  slice(1) %>%  # Just take the first geometry for each district
  ungroup()


# Join 
meta_data_with_coords <- meta_data_cleaned %>%
  left_join(cam_shapefile_districts_unique, 
            by = c("district_lower" = "shapefile_district_lower"))

# Check - should be 6336 rows
nrow(meta_data_with_coords)


# --- Check the remianing istricts in meta_data but NOT in shapefile
# Rows lost == 218
unmatched_districts <- meta_data_with_coords %>%
  filter(!district_lower %in% cam_shapefile_districts_unique$shapefile_district_lower) %>%
  count(district_lower, sort = TRUE)
cat("\nMerge summary:\n")
cat("Total rows in meta_data:", nrow(meta_data), "\n")
cat("Total rows after merge:", nrow(meta_data_with_coords), "\n")
cat("Rows with geometry:", sum(!is.na(st_dimension(meta_data_with_coords$geometry))), "\n")

print(unmatched_districts)
sum(unmatched_districts$n)

# Plot to validate districts 
sf_meta_data_with_coords <- st_as_sf(meta_data_with_coords)
ggplot(sf_meta_data_with_coords) +
  geom_sf() +
  geom_sf_text(aes(label = district_lower), size = 2, check_overlap = TRUE) +
  theme_minimal()




# --- 2) Population-Weighted Centroids 
# Calculates the geographic center of each district weighted by where people actually live,
# rather than the simple geometric cente
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

# Function to calculate population-weighted centroid
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




# --- 3) Population-Weighted Mosquito Value
# Calculates the average mosquito density (Aedes aegypti, Aedes albopictus, Anopheles funestus, Anopheles gambiae) experienced by the population in each district. 
# Each raster cell's mosquito value is weighted by the number of people living in that cell, then averaged across the district. 
# This represents the district-wide mosquito exposure of the population, accounting for both spatial variation in mosquito density and population distribution.
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


# Drop Nas + remove duplicates
sf_meta_data_with_coords_pw <- sf_meta_data_with_coords_pw %>%
  drop_na(CHIKV_sE2, MAYV_E2, ONNV_VLP)

sf_meta_data_with_coords_pw <- sf_meta_data_with_coords_pw[!duplicated(sf_meta_data_with_coords_pw$Sample), ]


# Add year of survey column
sf_meta_data_with_coords_pw$year_of_survey <- as.numeric(substr(sf_meta_data_with_coords_pw$Sample, 1, 4))

nrow(sf_meta_data_with_coords_pw)

# --- Save file with pop weighted coords and mosquito proportions for spatial analysis 
saveRDS(sf_meta_data_with_coords_pw, '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/meta_data_with_coords.rds')
# Also save dataframe without geometry for Stan Multisero model
preprocessed_meta_data_without_coords <- sf_meta_data_with_coords_pw %>%
  sf::st_drop_geometry()
write.csv(preprocessed_meta_data_without_coords, 
          '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/meta_data_without_coords.csv', 
          row.names = FALSE)


sf_meta_data_with_coords_pw <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/meta_data_with_coords.rds')


# --- Figure 1 ----
location_counts <- sf_meta_data_with_coords_pw %>%
  st_drop_geometry() %>%  # Remove spatial features
  group_by(district_lower, Longitude, Latitude) %>%
  summarise(n_samples = n(), .groups = 'drop')

# Convert raster to data frame for ggplot
cam_pop_df <- as.data.frame(cam_pop_den, xy = TRUE)
colnames(cam_pop_df) <- c("x", "y", "pop_density")

# Create the population density inset map
inset_map <- ggplot(cam_pop_df, aes(x = x, y = y, fill = pop_density)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Population\nDensity", 
                       option = "plasma",
                       trans = "log10",
                       na.value = "transparent") +
  coord_equal() +
  labs(title = "Log Population Density") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    plot.margin = margin(3, 3, 3, 3)
  )

# Figure 1a: Map of Cameroon with sample collection locations
fig1a <- ggplot() +
  geom_sf(data = sf_meta_data_with_coords_pw, fill = "#ffffff", color = "#6d7275") +
  geom_point(data = location_counts, 
             aes(x = Longitude, y = Latitude, size = n_samples),
             color = "#00798c", alpha = 0.9) +
  scale_size_continuous(name = "Number of Samples", range = c(2, 10),
                        breaks = seq(0, max(location_counts$n_samples), by = 30), limits = c(0, max(location_counts$n_samples)))  +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 24),  
    axis.title.x = element_text(size = 24),                             # X-axis label
    axis.title.y = element_text(size = 24),                             # Y-axis label
    axis.text.x = element_text(size = 20),                              # X-axis tick labels
    axis.text.y = element_text(size = 20),                              # Y-axis tick labels
    legend.title = element_text(size = 20),                             # Legend title
    legend.text = element_text(size = 20)                               # Legend text
  )

# Add the inset map using patchwork
fig1a_with_inset <- fig1a +
  inset_element(inset_map, 
                left = 0.05, bottom = 0.60, 
                right = 0.45, top = 1,
                align_to = "panel")

print(fig1a_with_inset)


# Figure 1b: Number of samples by year of survey
fig1b <- sf_meta_data_with_coords_pw %>%
  st_drop_geometry() %>%  # Remove geometry for plotting
  group_by(year_of_survey) %>%
  summarise(n_samples = n()) %>%
  ggplot(aes(x = factor(year_of_survey), y = n_samples)) +
  geom_bar(stat = "identity", fill = "#00798c") +
  geom_text(size = 8, aes(label = n_samples), vjust = -0.5) +
  theme_minimal() +
  labs(x = "Year of Survey",
       y = "Number of Samples") +
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

print(fig1b)

# Figure 1c: Male vs Female by Age
# Recode Sex variable (1 = Male, 2 = Female)
census_totals <- cameroon_age_2025 %>%
summarise(
  total_M = sum(M),
  total_F = sum(F)
)

pyramid_data <- sf_meta_data_with_coords_pw %>%
  st_drop_geometry() %>%
  filter(!is.na(Sex), !is.na(AgeInYears)) %>%
  mutate(
    Sex_label = case_when(
      Sex == 1 ~ "Male",
      Sex == 2 ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Sex_label)) %>%
  mutate(
    Sex_label = factor(Sex_label, levels = c("Male", "Female")),
    age_group = cut(
      AgeInYears,
      breaks = seq(0, 110, by = 10),
      include.lowest = TRUE,
      right = FALSE,
      labels = c(
        "0-9", "10-19", "20-29", "30-39", "40-49",
        "50-59", "60-69", "70-79", "80-89", "90-99", "100+"
      )
    )
  ) %>%
  group_by(age_group, Sex_label) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    Sex_label = droplevels(Sex_label),
    count = ifelse(Sex_label == "Female", -count, count)
  )

sample_totals <- pyramid_data %>%
  group_by(Sex_label) %>%
  summarise(total_samples = sum(abs(count)), .groups = 'drop')

# Calculate expected samples based on census proportions
# You'll need to combine your census data into 10-year age groups too
expected_data <- cameroon_age_2025 %>%
  mutate(
    age_group_10 = case_when(
      Age %in% c("0-4", "5-9") ~ "0-9",
      Age %in% c("10-14", "15-19") ~ "10-19",
      Age %in% c("20-24", "25-29") ~ "20-29",
      Age %in% c("30-34", "35-39") ~ "30-39",
      Age %in% c("40-44", "45-49") ~ "40-49",
      Age %in% c("50-54", "55-59") ~ "50-59",
      Age %in% c("60-64", "65-69") ~ "60-69",
      Age %in% c("70-74", "75-79") ~ "70-79",
      Age %in% c("80-84", "85-89") ~ "80-89",
      Age %in% c("90-94", "95-99") ~ "90-99",
      Age == "100+" ~ "100+"
    )
  ) %>%
  group_by(age_group_10) %>%
  summarise(M = sum(M), F = sum(F), .groups = "drop") %>%
  pivot_longer(cols = c(M, F), names_to = "Sex_label", values_to = "census_count") %>%
  mutate(
    Sex_label = case_when(
      Sex_label == "M" ~ "Male",
      Sex_label == "F" ~ "Female"
    ),
    Sex_label = factor(Sex_label, levels = c("Male", "Female"))
  ) %>%
  left_join(sample_totals, by = "Sex_label") %>%
  mutate(
    # Calculate total census population by sex
    total_census = ifelse(Sex_label == "Male", census_totals$total_M, census_totals$total_F),
    # Calculate proportion of total population in this age group
    proportion = census_count / total_census,
    # Expected samples = total samples for this sex * proportion in this age group
    expected_count = total_samples * proportion,
    # Make female counts negative for pyramid
    expected_count = ifelse(Sex_label == "Female", -expected_count, expected_count),
    age_group = factor(age_group_10, levels = c(
      "0-9", "10-19", "20-29", "30-39", "40-49",
      "50-59", "60-69", "70-79", "80-89", "90-99", "100+"
    ))
  )

fig1c <- ggplot(pyramid_data, aes(x = age_group, y = count, fill = Sex_label)) +
  geom_bar(stat = "identity", width = 0.9) +
  # Add expected distribution as lines
  geom_line(data = expected_data, 
            aes(x = age_group, y = expected_count, color = Sex_label, group = Sex_label),
            linewidth = 1.2, linetype = "solid") +
  geom_point(data = expected_data,
             aes(x = age_group, y = expected_count, color = Sex_label),
             size = 3) +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("Male" = "#b66577", "Female" = "#379392"),
                    name = "Observed") +
  scale_color_manual(values = c("Male" = "#8b4456", "Female" = "#2a6e6d"),
                     name = "Expected (Census)") +
  guides(
    fill = guide_legend(override.aes = list(shape = NA)),
    color = guide_legend(override.aes = list(linetype = 1, shape = 16))
  ) +
  theme_minimal() +
  labs(x = "Age Group",
       y = "Number of Samples") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.line = element_line(color = "black", linewidth = 0.7),
        axis.ticks.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        legend.position.inside = c(0.95, 0.5),
        panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.title = element_text(size = 24),
        aspect.ratio = 0.75,
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 20))

print(fig1c)

# use patchwork to combine figures
(fig1a_with_inset | (fig1b / fig1c)) +
  plot_layout(widths = c(3, 2), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    plot.tag = element_text(size = 20, face = "bold", hjust = 0, vjust = 1),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/Figure1.png", 
       width = 16, height = 8, units = "in", dpi = 300)
