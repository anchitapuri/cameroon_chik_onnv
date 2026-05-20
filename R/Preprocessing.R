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
library(readxl)
library(geodata)
library(ggspatial)
library(patchwork)

# --- Read data files 
# original data
meta_data <- read.csv(here(data, 'MFI_meta.csv'))
                   

# shapefile #1
cam_shapefile_districts <- read_sf(here(data, 'Caedistricts179_region.shp'))
# Second shapefile used (to find remaining mismatched districts)
cam_shapefile_districts2 <- read_sf(here(data, 'cmr_admin3.shp'))



# Load population rasters
cam_pop <- rast(here(data, 'cmr_ppp_2020_UNadj.tif'))
cam_pop_den <- rast(here(data, 'cmr_pd_2020_1km_UNadj.tif'))

# Load mosquito maps
aegypti <- rast(here(data, 'Aedes_maps_public', 'aegypti.tif'))
albopictus <- rast(here(data, 'Aedes_maps_public', 'albopictus.tif'))
anopheles_funestus <- rast(here(data, 'Anopheles_maps', '2010_Anopheles_funestus_CMR.tiff'))
anopheles_gambiae <- rast(here(data, 'Anopheles_maps', '2010_Anopheles_gambiae_ss_CMR.tiff'))



# Load population by gender / age data 
cameroon_age_2025 <- read.csv(here(data, 'CameroonAge2025.csv'))
cameroon_age_2025 <- cameroon_age_2025 %>%
  mutate(total = M + F)


# excel sheet with additional district geometeries - these were missing from both the shapefiles
missing_districts_geometeries <- read_excel(here(data, 'Districts_sante_2021.xls'), sheet = "Sheet2")


# drop NAs 
nrow(meta_data) #6336
length(unique(tolower(meta_data$DistrictOfresidence))) #208


sum(is.na(meta_data$CHIKV_sE2)) #920
sum(is.na(meta_data$ONNV_VLP)) #11
sum(is.na(meta_data$MAYV_E2)) #15
sum(is.na(meta_data$AgeInYears)) #7
sum(meta_data$AgeInYears == 0, na.rm = TRUE) #117 
sum(is.na(meta_data$Sex))  #22
unique(meta_data$Sex)
sum(meta_data$Sex == 9, na.rm = TRUE) #9
sum(is.na(meta_data$DistrictOfresidence)) #3


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
length(unique(cam_shapefile_districts_merged$geometry)) #179 unique geometeries and districts 
View(cam_shapefile_districts_merged)


# --- Meta Data ---
# Create lowercase district column for meta_data
meta_data <- meta_data %>%
  mutate(district_lower = tolower(DistrictOfresidence))
length(unique(meta_data$district_lower)) #208 unique districts 




# --- Non-matching districts between shapefile 1 and data 
non_matching_meta <- meta_data %>%
  filter(!district_lower %in% cam_shapefile_districts_merged$shapefile_district_lower)

# Print unique non-matching district names
cat("\nNon-matching districts:\n")
print(unique(non_matching_meta$district_lower))
length(unique(non_matching_meta$district_lower)) #42 non matching districts 

# Count rows for each non-matching district
non_matching_counts <- non_matching_meta %>%
  count(district_lower) %>%
  arrange(desc(n))
cat("\nCount of non-matching districts:\n")
print(non_matching_counts)
cat("\nTotal non-matching rows:", nrow(non_matching_meta), "\n") #this is 976 rows


# --- Manual mapping ---- 
meta_data_districts_added <- meta_data %>%
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
    district_lower == "garoua urbain" ~ "garoua i", # Garoua urban, the shape file should be Garoua I or 1
    district_lower == "eyumodjock" ~ "eyumojock",
    district_lower == "garoua 2" ~ "garoua ii",
    district_lower == "ndikinimeki" ~ "ndikinimiki",
    district_lower == "mbandjock" ~ "mbanjock",
    district_lower == 'bandjoun' ~ "banjoun", 
    district_lower == 'bangangte' ~ "bangante",
    district_lower == 'bangourain'~ "bangorain",
    district_lower == 'garoua rural' ~ 'garoua ii' , # garoua rural == garoua ii
    district_lower == 'maroua 2' ~ 'maroua rural', # same logic as before 2 == rural
    TRUE ~ district_lower
  ))
length(unique(meta_data_districts_added$district_lower)) #198 districts - some names merged 

# Use shapefile 2 since there are still districts in data missing from shapefile 1
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
                             'oku', 'ngaoundal')
rows_from_shapefile2 <- cam_shapefile_districts2 %>%
  filter(shapefile_district_lower2 %in% districts_new_shapefile)
rows_to_add <- rows_from_shapefile2 %>%
  st_transform(st_crs(cam_shapefile_districts)) %>%
  transmute(
    shapefile_district_lower = shapefile_district_lower2,
    geometry = geometry
    # Map other columns as needed
  )
# add rows from shapefile #2 to original file
cam_shapefile_districts_merged <- cam_shapefile_districts_merged %>%
  bind_rows(rows_to_add)
cam_shapefile_districts_merged <- cam_shapefile_districts_merged %>%
  st_make_valid()
# One geometry per district 
cam_shapefile_districts_unique <- cam_shapefile_districts_merged %>%
  group_by(shapefile_district_lower) %>%
  slice(1) %>%  # Just take the first geometry for each district
  ungroup()
length(unique(meta_data_districts_added$district_lower))

# --- Check the remianing istricts in meta_data but NOT in shapefile
# Rows lost currently == 192
unmatched_districts <- meta_data_districts_added %>%
  filter(!district_lower %in% cam_shapefile_districts_unique$shapefile_district_lower) %>%
  count(district_lower, sort = TRUE)
cat("\nMerge summary:\n")
cat("Total rows in meta_data:", nrow(meta_data), "\n")
cat("Total rows after merge:", nrow(meta_data_districts_added), "\n")
 

print(unmatched_districts)
sum(unmatched_districts$n)


# Drop NAs and Chad rows
meta_data_districts_added <- meta_data_districts_added |>
  filter(!is.na(district_lower),
         !district_lower %in% c("abeche", "biltine"))
nrow(meta_data_districts_added) #6331

length(unique(meta_data_districts_added$district_lower)) # 195

# districts still unmatched 
remaining_unmatched_districts <- c(
  "boko",
  "dang",
  "mozogo",
  "bangue",
  "japoma",
  "odza",
  "abo",
  "nkolbisson",
  "mvog-ada" 
)
subset_missing_districts_geometeries <- missing_districts_geometeries %>%
  filter(tolower(District) %in% remaining_unmatched_districts)
nrow(subset_missing_districts_geometeries)

# visualise where the missing geometeries fall 
missing_districts_geometeries_sf <- subset_missing_districts_geometeries %>%
  mutate(
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  ) %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Make sure shapefile CRS matches (reproject if needed)
cam_shapefile_districts_merged <- st_transform(cam_shapefile_districts_merged, crs = 4326)

# Plot
ggplot() +
  geom_sf(data = cam_shapefile_districts_merged, fill = "lightgrey", color = "white", linewidth = 0.3) +
  geom_sf(data = missing_districts_geometeries_sf, aes(color = Region), size = 3, alpha = 0.8) +
  labs(
    title = "Cameroon Districts with Survey Points",
    color = "Region"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )


sf_use_s2(FALSE)
# Spatial join: finds which polygon each point falls within
mapped_districts <- st_join(
  missing_districts_geometeries_sf,
  cam_shapefile_districts_merged[, c("shapefile_district_lower")],
  join = st_intersects
)

district_lookup <- mapped_districts %>%
  st_drop_geometry() %>%
  dplyr::select(District, shapefile_district_lower)
district_lookup
# replace meta data missing districs with mapped districts 
meta_data_districts_added <- meta_data_districts_added %>%
  left_join(district_lookup,
            by = c("district_lower" = "District")) %>%
  mutate(district_lower = coalesce(shapefile_district_lower, district_lower)) %>%
  dplyr::select(-shapefile_district_lower)


# Join 
meta_data_with_coords <- meta_data_districts_added %>%
  left_join(cam_shapefile_districts_unique, 
            by = c("district_lower" = "shapefile_district_lower"))


cat("Rows with geometry:", sum(!is.na(st_is_empty(meta_data_with_coords$geometry))), "\n")


# Plot to validate districts 
sf_meta_data_with_coords <- st_as_sf(meta_data_with_coords)
length(unique(sf_meta_data_with_coords$geometry))
length(unique(sf_meta_data_with_coords$district_lower))
quartz() 
ggplot(sf_meta_data_with_coords) +
  geom_sf() +
  geom_sf_text(aes(label = district_lower), size = 2, check_overlap = TRUE) +
  theme_minimal()


# total 6331 rows with valid geometry!


# --- 2) Population-Weighted Centroids 
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
    districts %>% st_drop_geometry() %>% dplyr::select(district_lower, area_km2),  # Drop geometry here
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
      dplyr::select(-district_id), 
    by = "district_lower"
  ) %>%
  left_join(aeg_pw_df %>% dplyr::select(district_lower, aeg_pw_district), by = "district_lower") %>%
  left_join(alb_pw_df %>% dplyr::select(district_lower, alb_pw_district), by = "district_lower") %>%
  left_join(fun_pw_df %>% dplyr::select(district_lower, fun_pw_district), by = "district_lower") %>%
  left_join(gam_pw_df %>% dplyr::select(district_lower, gam_pw_district), by = "district_lower")


# --- Validate: Plot Population weighted vs unweighted coords 
sf_use_s2(FALSE)
unweighted_centroids <- st_centroid(sf_meta_data_with_coords)
unweighted_centroids <- st_coordinates(unweighted_centroids)
quartz()
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

# --- Save file with pop weighted coords and mosquito proportions for spatial analysis 

# -- drop duplicated ids (use ids as identification of each sample in later analysis)
length(unique(sf_meta_data_with_coords_pw$id[duplicated(sf_meta_data_with_coords_pw$id)]))
sf_meta_data_with_coords_pw[duplicated(sf_meta_data_with_coords_pw$id) | duplicated(sf_meta_data_with_coords_pw$id, fromLast = TRUE), ]

sf_meta_data_with_coords_pw_filtered <- sf_meta_data_with_coords_pw |>
  distinct(id, .keep_all = TRUE)
nrow(sf_meta_data_with_coords_pw_filtered) #6324

# Verify
nrow(sf_meta_data_with_coords_pw_filtered)  # should be 6324
length(unique(sf_meta_data_with_coords_pw_filtered$geometry))
length(unique(sf_meta_data_with_coords_pw_filtered$district_lower))

# Add year of survey column
sf_meta_data_with_coords_pw_filtered$year_of_survey <- as.numeric(substr(sf_meta_data_with_coords_pw_filtered$Sample, 1, 4))
unique(sf_meta_data_with_coords_pw_filtered$year_of_survey)
nrow(sf_meta_data_with_coords_pw_filtered)

#add population density to dataframe
sf_meta_data_with_coords_pw_filtered <- sf_meta_data_with_coords_pw_filtered %>%
  mutate(
    pop_density = Total_Population / area_km2,
    log_pop_density = log(pop_density + 1)  # +1 to avoid log(0) issues
  )
  nrow(sf_meta_data_with_coords_pw_filtered)


saveRDS(sf_meta_data_with_coords_pw_filtered, here('Results/meta_data_with_coords.rds'))

# Also save dataframe without geometry for Stan Multisero model
preprocessed_meta_data_without_coords <- sf_meta_data_with_coords_pw_filtered %>%
  sf::st_drop_geometry()
write.csv(preprocessed_meta_data_without_coords, 
          here('Results/meta_data_without_coords.csv'), 
          row.names = FALSE)



sf_meta_data_with_coords_pw_filtered <- readRDS('Results/meta_data_with_coords.rds')
preprocessed_meta_data_without_coords <- read.csv('Results/meta_data_without_coords.csv')
nrow(preprocessed_meta_data_without_coords)

# # Drop NAs and age = 0 and sex = 9 
meta_data_clean_with_coords <- subset(
  sf_meta_data_with_coords_pw_filtered,
  !is.na(CHIKV_sE2) &
  !is.na(ONNV_VLP) &
  !is.na(MAYV_E2) &
  !is.na(AgeInYears) &
  AgeInYears != 0 &
  !is.na(Sex) &
  Sex != 9
)

# # Drop NAs and age = 0 and sex = 9 
meta_data_clean_without_coords <- subset(
  preprocessed_meta_data_without_coords,
  !is.na(CHIKV_sE2) &
  !is.na(ONNV_VLP) &
  !is.na(MAYV_E2) &
  !is.na(AgeInYears) &
  AgeInYears != 0 &
  !is.na(Sex) &
  Sex != 9
)


nrow(meta_data_clean_with_coords) #5272
nrow(meta_data_clean_without_coords) #5272

# save RDS  - these are used for ALL downstream analysis 
saveRDS(meta_data_clean_with_coords, here('Results/meta_data_clean_with_coords.rds'))
saveRDS(meta_data_clean_without_coords, here('Results/meta_data_clean_without_coords.rds'))

 # save another version with CHIK,ONNV and MAYV NAs 
 # this is for supplementary materials - to show that results are similar when including all samples vs only those with complete data
meta_data_without_coords_supp_materials <- subset(
  preprocessed_meta_data_without_coords,
  !is.na(AgeInYears) &
  AgeInYears != 0 &
  !is.na(Sex) &
  Sex != 9
)


saveRDS(meta_data_without_coords_supp_materials, here('Results/meta_data_without_coords_supp_materials.rds'))
