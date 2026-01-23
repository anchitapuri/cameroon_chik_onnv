source(here('~/Desktop/Cameroon_Analysis_2025/Functions.R'))

# read data files with labels 
meta_data <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/final_meta_data_with_labels.csv')
nrow(meta_data)

# Load shapefile
cam_shapefile_districts <- read_sf('Caedistricts179_region.shp')
# clean names + create lower case district column
cam_shapefile_districts$NAME2 <- gsub("DS_", "", cam_shapefile_districts$NAME2)
cam_shapefile_districts <- cam_shapefile_districts %>%
  mutate(shapefile_district_lower = tolower(NAME2))

# Second shapefile used (to find remaining mismatched districts)
cam_shapefile_districts2 <- read_sf('cmr_admin3.shp')
cam_shapefile_districts2 <- cam_shapefile_districts2 %>%
  mutate(shapefile_district_lower2 = tolower(adm3_name1))

# check how many names and districts 
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

# Keep all other districts as-is
other_districts <- cam_shapefile_districts %>%
  filter(shapefile_district_lower != "manoka")

# Combine them
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


# Merge geometries in the second shapefile (in case it has duplicates too)
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

# Check - should be 5407 rows
nrow(meta_data_with_coords)


# --- Check the remianing istricts in meta_data but NOT in shapefile
# Rows lost == 213
unmatched_districts <- meta_data_with_coords %>%
  filter(!district_lower %in% cam_shapefile_districts$shapefile_district_lower) %>%
  count(district_lower, sort = TRUE)
cat("\nMerge summary:\n")
cat("Total rows in meta_data:", nrow(meta_data), "\n")
cat("Total rows after merge:", nrow(meta_data_with_coords), "\n")
cat("Rows with geometry:", sum(!is.na(st_dimension(meta_data_with_coords$geometry))), "\n")


# Plot to validate districts 
sf_meta_data_with_coords <- st_as_sf(meta_data_with_coords)
ggplot(sf_meta_data_with_coords) +
  geom_sf() +
  geom_sf_text(aes(label = district_lower), size = 2, check_overlap = TRUE) +
  theme_minimal()

# Save files with coords for downstream analysis
saveRDS(meta_data_with_coords, 'meta_data_with_coords.rds')
