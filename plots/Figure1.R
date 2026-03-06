

# read preprocessed data 
sf_meta_data_with_coords_pw_filtered <- read.RDS(here('Results/sf_meta_data_with_coords_pw_filtered.rds'))


#demongraphy data 
cam_pop <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_ppp_2020_UNadj.tif")
cam_pop_den <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_pd_2020_1km_UNadj.tif")

# --- Figure 1 ----
location_counts <- sf_meta_data_with_coords_pw_filtered %>%
  st_drop_geometry() %>%  # Remove spatial features
  group_by(district_lower, Longitude, Latitude) %>%
  summarise(n_samples = n(), .groups = 'drop')


# Convert raster to data frame for ggplot
cam_pop_df <- as.data.frame(cam_pop_den, xy = TRUE)
colnames(cam_pop_df) <- c("x", "y", "pop_density")

# Create the population density inset map
inset_map <- ggplot(cam_pop_df, aes(x = x, y = y, fill = pop_density)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Log Population Density (per km²)", 
                       option = "plasma",
                       trans = "log10",
                       na.value = "transparent",
                       breaks = c(1, 10, 100, 1000, 10000),          # set the tick positions
                       labels = c("1", "10", "100", "1,000", "10,000"),
                       guide = guide_colorbar(
                       barheight = unit(0.2, "cm"),      # thin (horizontal)
                       barwidth = unit(6.5, "cm"),       # wide (horizontal)
                       ticks = TRUE,
                       title.position = "top",
                       label.position = "bottom",
                       ticks.length = unit(0.1, "cm"),
                       title.vjust = 0.5)) +
  coord_equal() +
  theme_void() +
  theme(
     legend.position  = c(0.5, 0.05),   
     legend.text = element_text(size = 11),
     legend.title = element_text(size = 11),
     legend.direction = "horizontal",    
     legend.margin = margin(20, 0, 0, 0),  # adds space above the legend
     #plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
     plot.margin = margin(10, 5, 15, 5)
  )

# Figure 1a: Map of Cameroon with sample collection locations
fig1a <- ggplot() +
  geom_sf(data = sf_meta_data_with_coords_pw_filtered, fill = "#ffffff", color = "#6d7275") +
  geom_point(data = location_counts, 
             aes(x = Longitude, y = Latitude, size = n_samples),
             shape = 21, fill = "#015b69", colour = "white", alpha = 0.85) +
  scale_size_continuous(name = "Number of \nSamples", range = c(2, 10),
                        breaks = seq(0, max(location_counts$n_samples), by = 30)) +
  annotation_scale(
    plot_unit = "km",
    bar_cols = c("black", "white"),  # alternating black/white like the reference
    height = unit(0.2, "cm"),
    text_family = "sans",
    pad_y = unit(0.8, "cm"),
    text_cex = 1.5      
  ) +
  theme_minimal(base_size = 11)  +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(size = 20),                             # Legend title
    legend.text = element_text(size = 20) ,                              # Legend text
    legend.position = c(1.05, 0.4),
    legend.key.height = unit(0.4, "cm"),
    legend.spacing.y  = unit(0.2, "cm")
  )
quartz()
fig1a_with_inset <- fig1a +
  inset_element(
    inset_map, left = -1, bottom = 0.5, right = 1, top = 1, align_to = 'plot')

print(fig1a_with_inset)



# --- Save Figure 1a
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig1a.png", 
       plot = fig1a_with_inset,  
       width = 10, 
       height = 10, 
       units = "in", 
       dpi = 300,
       bg = "white")

# Figure 1b: Number of samples by year of survey
fig1b <- sf_meta_data_with_coords_pw_filtered %>%
  st_drop_geometry() %>%  # Remove geometry for plotting
  group_by(year_of_survey) %>%
  summarise(n_samples = n()) %>%
  ggplot(aes(x = factor(year_of_survey), y = n_samples)) +
    scale_y_continuous(limits = c(0, 1650)) +   
  geom_bar(stat = "identity", fill = "#015b69") +
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
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig1b.png", 
       plot = fig1b,    # swap this for your actual plot object name
       width = 7, 
       height = 10, 
       units = "in", 
       dpi = 300,
       bg = "white")

# Figure 1c: Male vs Female by Age
# Recode Sex variable (1 = Male, 2 = Female)
census_totals <- cameroon_age_2025 %>%
summarise(
  total_M = sum(M),
  total_F = sum(F)
)
nrow(sf_meta_data_with_coords_pw_filtered) #6324
sum(is.na(sf_meta_data_with_coords_pw_filtered$Sex)) #21
sum(sf_meta_data_with_coords_pw_filtered$Sex == 9, na.rm = TRUE) #9 
sum(is.na(sf_meta_data_with_coords_pw_filtered$AgeInYears)) # 6 
table(sf_meta_data_with_coords_pw_filtered$Sex)

# Mean age = 18
mean(sf_meta_data_with_coords_pw_filtered$AgeInYears, na.rm = TRUE)


pyramid_data <- sf_meta_data_with_coords_pw_filtered %>%
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
sample_totals
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
  scale_fill_manual(values = c("Male" = "#b84f74", "Female" = "#00798c"),
                    name = "Observed") +
  scale_color_manual(values = c("Male" = "#7c334d", "Female" = "#014751"),
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

fig1 <- (fig1a_with_inset | (fig1b / fig1c)) + plot_layout(widths = c(2, 1))



ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/fig1.png", 
       plot = fig1,    # swap for your actual plot object name
       width = 19.5, 
       height = 12, 
       units = "in", 
       dpi = 300,
       bg = "white")

