# Fit geostatistical models to predict ONNV FOI and prevelance across Cameroon
# using the stochastic partial differential equation (SPDE) approach and the R-INLA package.
library(ggraph)
library(igraph)
library(lhs)
library(matrixStats)
library(mvtnorm)
library(matrixcalc)
library(here)
library(mixR)
library(raster)
library(INLA)
library(PBSmapping)  
library(sp)
library(gstat)
library(sf)
library(rnaturalearth)
library(viridis)
library(patchwork) 
library(spatstat)
library(dplyr)
library(tidyr)
library(rworldmap)
library(exactextractr)
library(rworldxtra)
library(centr)
library(terra)
library(colorspace)
library(ggspatial)


# --- Source functions
source(here('R/Functions.R'))

# Population and mosquito rasts
# Load mosquito maps
data_dir <- 'Data'
aegypti <- rast(here(data_dir, 'Aedes_maps_public', 'aegypti.tif'))
albopictus <- rast(here(data_dir, 'Aedes_maps_public', 'albopictus.tif'))
anopheles_funestus <- rast(here(data_dir, '2010_Anopheles_funestus_CMR.tiff'))
anopheles_gambiae <- rast(here(data_dir, '2010_Anopheles_gambiae_ss_CMR.tiff'))


cam_pop <- rast(here(data_dir, 'cmr_ppp_2020_UNadj.tif'))
sum(values(cam_pop), na.rm = TRUE)
cam_pop_agg <- terra::aggregate(cam_pop, fact = 10, fun = sum, na.rm = TRUE) # (~1 km, from aggregating the ~100 m WorldPop grid by a factor of 10)".

# Cameroon population by age
cameroon_age_2025 <- read.csv(here(data_dir, 'CameroonAge2025.csv'))
cameroon_age_2025 <- cameroon_age_2025 %>%
  mutate(total = M + F)

cameroon_age_2025$total
w_age <- cameroon_age_2025$total / sum(cameroon_age_2025$total)
sum(cameroon_age_2025$total)

age_groups <- data.frame(
  age_string = cameroon_age_2025$Age,
  stringsAsFactors = FALSE
)
# Extract lower and upper bounds from strings like "0-4", "5-9", etc.
age_groups$age_lower <- as.numeric(sub("-.*", "", age_groups$age_string))
age_groups$age_upper <- age_groups$age_lower + 5                              # contiguous 5-yr bands


# Handle the last age group if it's something like "80+" 
if (grepl("\\+", age_groups$age_string[nrow(age_groups)])) {
  age_groups$age_lower[nrow(age_groups)] <- as.numeric(sub("\\+", "", age_groups$age_string[nrow(age_groups)]))
  age_groups$age_upper[nrow(age_groups)] <- 120  # or whatever max age you want
}


# regions for foi, prev and infections stratified by region
cameroon_regions <- ne_states(country = "Cameroon", returnclass = "sf")
# rename 'name' to region
regions_sf <- cameroon_regions %>%
  dplyr::select(region = name) %>%   # standardise column name
  st_make_valid()



# ----- Read preprocessed data with coords 
meta_data_with_coords <- readRDS(here('Results/meta_data_clean_with_coords.rds'))
colnames(meta_data_with_coords)
nrow(meta_data_with_coords)

meta_data_with_coords_supp_materials <- readRDS(here('Results/meta_data_with_coords_supp_materials.rds'))
nrow(meta_data_with_coords_supp_materials)

meta_data_without_coords_supp_materials <- subset(
  meta_data_with_coords_supp_materials,
  !is.na(AgeInYears) &
  AgeInYears != 0 &
  !is.na(Sex) &
  Sex != 9
)
meta_data_with_coords_supp_materials <- subset(
  meta_data_with_coords_supp_materials,
  !is.na(AgeInYears) &
  AgeInYears != 0 &
  !is.na(Sex) &
  Sex != 9
)
nrow(meta_data_with_coords_supp_materials)



# ----- Read labelled data
meta_data_with_labels <- read.csv(here('Results/meta_data_with_labels.csv'))
colnames(meta_data_with_labels)
nrow(meta_data_with_labels)

unique(meta_data_with_labels$year_of_survey)


meta_data_onnv_samples <- read.csv(here('Results/meta_data_onnv_samples_with_labels.csv'))
nrow(meta_data_onnv_samples)

# Convert coords to Easting and Northing
sp_vill <- SpatialPoints(cbind(meta_data_with_labels$Longitude, meta_data_with_labels$Latitude))
points_to_extract <- terra::vect(sp_vill)
length(points_to_extract)

data_points <- meta_data_with_coords %>%
  st_drop_geometry() %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

data_points_supp_materials <- meta_data_with_coords_supp_materials %>%
  st_drop_geometry() %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

data_utm <- st_transform(data_points, crs = 32633)
data_utm_supp_materials <- st_transform(data_points_supp_materials, crs = 32633)
coords_utm <- st_coordinates(data_utm) / 1000  # Convert to km
coords_utm_supp_materials <- st_coordinates(data_utm_supp_materials) / 1000  # Convert to km
colnames(coords_utm) <- c("Easting", "Northing")
colnames(coords_utm_supp_materials) <- c("Easting", "Northing")

# Add Easting and Northing to dataframe
meta_data_with_coords$Easting <- coords_utm[, "Easting"]
meta_data_with_coords$Northing <- coords_utm[, "Northing"]

meta_data_with_coords_supp_materials$Easting <- coords_utm_supp_materials[, "Easting"]
meta_data_with_coords_supp_materials$Northing <- coords_utm_supp_materials[, "Northing"]

meta_data_with_labels$Easting <- meta_data_with_coords$Easting
meta_data_with_labels$Northing <- meta_data_with_coords$Northing

meta_data_onnv_samples$Easting <- meta_data_with_coords_supp_materials$Easting
meta_data_onnv_samples$Northing <- meta_data_with_coords_supp_materials$Northing

colnames(meta_data_with_labels)

# check 
sum(is.na(meta_data_with_coords$Longitude))                      
nrow(meta_data_with_coords) == nrow(coords_utm)                
all(meta_data_with_labels$Sample == meta_data_with_coords$Sample) # same rows, same order


# Rename
model_data <- meta_data_with_labels
model_data_onnv_samples <- meta_data_onnv_samples
write.csv(model_data, here('Results/model_data.csv'), row.names = FALSE)

model_data <- read.csv(here('Results/model_data.csv'))
mean(model_data$AgeInYears)
table(model_data$Sex)


# Check for missing data / NA
sum(is.na(model_data$AgeInYears))
sum(model_data$AgeInYears == 0, na.rm = TRUE)
sum(is.na(model_data$Easting) & !is.na(model_data$Northing))
sum(is.na(model_data$ONNV_pos))



# --- Run INLA model for ONNV (with historic year of intro 1900)
onnv_results_pop_grid <- run_inla(
  year_intro = 1900,
  data = model_data,
  cam_pop = cam_pop,
  positive_col = "ONNV_pos")

# --- Run INLA model for ONNV (with historic year of intro 1900)
onnv_results_pop_grid_onnv_samples <- run_inla(
  year_intro = 1900,
  data = model_data_onnv_samples,
  cam_pop = cam_pop,
  positive_col = "ONNV_pos")


# --- Run INLA model with covariates - An. Gam and log pop density
# extract covariates for each grid cell
cam_pop_agg <- terra::aggregate(cam_pop, fact = 10, fun = sum, na.rm = TRUE)   # ~1 km cells
dens_agg    <- cam_pop_agg / terra::cellSize(cam_pop_agg, unit = "km")         # people / km²

# Covariates at observation points for fitting
obs_pts <- terra::vect(model_data, geom = c("Longitude","Latitude"), crs = "EPSG:4326")

model_data$gam_ras      <- terra::extract(anopheles_gambiae, obs_pts)[, 2] # [,2] = values, [,1] = ID
model_data$log_dens_ras <- log(terra::extract(dens_agg,      obs_pts)[, 2] + 1)

# Covariates at prediction grid points -> for the MAP ---------------
# align with coop (prediction grid)
grid_pts <- terra::as.points(cam_pop_agg, values = TRUE, na.rm = TRUE)

covar_grid <- data.frame(
  gam_ras      = terra::extract(anopheles_gambiae, grid_pts)[, 2],
  log_dens_ras = log(terra::extract(dens_agg,      grid_pts)[, 2] + 1)
)

# fill any grid NAs (cells where a raster has no data) so INLA give error
covar_grid$gam_ras[is.na(covar_grid$gam_ras)]  <- mean(covar_grid$gam_ras, na.rm = TRUE)
covar_grid$log_dens_ras[is.na(covar_grid$log_dens_ras)] <- mean(covar_grid$log_dens_ras, na.rm = TRUE)

sum(is.na(covar_grid$gam_ras))         # should be 0
sum(is.na(covar_grid$log_dens_ras))    # should be 0


onnv_results_pop_grid_multivariable <- run_inla_multivariable(
  year_intro   = 1900,
  data         = model_data,
  cam_pop      = cam_pop,
  positive_col = "ONNV_pos",
  covars       = c("gam_ras", "log_dens_ras"),
  covar_grid   = covar_grid
) 


# save results 
# --- Save prediction results 
saveRDS(onnv_results_pop_grid, here('Results/ONNV_INLAResults.rds'))
saveRDS(onnv_results_pop_grid_onnv_samples, here('Results/ONNV_INLAResults_ONNV_samples.rds'))
saveRDS(onnv_results_pop_grid_multivariable, here('Results/ONNV_INLAResults_multivariable.rds'))

# Read saved results
onnv_results_pop_grid <- readRDS(here('Results/ONNV_INLAResults.rds'))
onnv_results_pop_grid_onnv_samples <- readRDS(here('Results/ONNV_INLAResults_ONNV_samples.rds'))
onnv_results_pop_grid_multivariable <- readRDS(here('Results/ONNV_INLAResults_multivariable.rds'))


# Cameroon wide estimates (FOI, seroprevalence, infections) using the posterior samples from the INLA model + pop weighted 
national_foi_spatial <- national_foi(onnv_results_pop_grid, cam_pop)
national_foi_onnv_samples   <- national_foi(onnv_results_pop_grid_onnv_samples, cam_pop)
national_foi_mv <- national_foi_covars(onnv_results_pop_grid_multivariable, cam_pop,
                                       covar_grid = covar_grid)

                                 

# save 
saveRDS(national_foi_spatial, here('Results/national_foi_spatial.rds'))
saveRDS(national_foi_mv, here('Results/national_foi_mv.rds'))
saveRDS(national_foi_onnv_samples, here('Results/national_foi_onnv_samples.rds'))



national_foi_spatial <- readRDS(here('Results/national_foi_spatial.rds'))
national_foi_mv <- readRDS(here('Results/national_foi_mv.rds'))
national_foi_onnv_samples <- readRDS(here('Results/national_foi_onnv_samples.rds'))

# print values 
national_foi_spatial[c("mean","ciL","ciU")] 
national_foi_mv[c("mean","ciL","ciU")] 
national_foi_onnv_samples[c("mean","ciL","ciU")] 


assign_regions <- function(coop, regions_sf, region_col = "region", crs = 4326) {
  pts <- sf::st_as_sf(as.data.frame(coop), coords = 1:2, crs = crs)
  regions_sf <- sf::st_transform(regions_sf, crs)
  hit   <- sf::st_intersects(pts, regions_sf)
  first <- vapply(hit, function(i) if (length(i)) i[1] else NA_integer_, integer(1))
  miss  <- is.na(first)
  if (any(miss)) {
    message(sum(miss), " cells outside all polygons; assigned to nearest region.")
    first[miss] <- sf::st_nearest_feature(pts[miss, ], regions_sf)
  }
  regions_sf[[region_col]][first]
}

cell_region <- assign_regions(onnv_results_pop_grid$coop, regions_sf, "region")

national_metrics_summary_new <- function(res, cam_pop, age_groups, w_age,
                                     cell_region = NULL,
                                     agg_factor = 10, n = 1000, seed = 1) {
  set.seed(seed)
  stopifnot(length(w_age) == nrow(age_groups))

  cam_pop_agg <- terra::aggregate(cam_pop, fact = agg_factor, fun = sum, na.rm = TRUE)
  pop <- as.data.frame(terra::as.points(cam_pop_agg, values = TRUE, na.rm = TRUE))[, 1]
  stopifnot(length(pop) == nrow(res$coop))
  ok <- !is.na(pop); pop <- pop[ok]; w <- pop / sum(pop)

  sel   <- list(Intercept = 1, spatial.field = seq_len(res$mesh$n))
  samp  <- inla.posterior.sample(n, res$output, selection = sel)
  Ap_ok <- inla.spde.make.A(mesh = res$mesh, loc = as.matrix(res$coop))[ok, , drop = FALSE]

  L <- age_groups$age_lower; U <- age_groups$age_upper; dW <- U - L
  suscept <- function(lam) {
    s <- sapply(seq_along(L), function(j)
      ifelse(lam > 1e-10, (1/(lam*dW[j])) * (exp(-lam*L[j]) - exp(-lam*U[j])), 1))
    as.numeric(s %*% w_age)
  }

  # --- regional setup ---
  do_reg <- !is.null(cell_region)
  if (do_reg) {
    stopifnot(length(cell_region) == nrow(res$coop))
    grp     <- factor(cell_region[ok])
    gi      <- as.integer(grp)
    reg_pop <- as.numeric(rowsum(pop, gi))
    reg_foi <- reg_sero <- reg_inf <-
      matrix(NA_real_, nlevels(grp), n, dimnames = list(levels(grp), NULL))
  }

  # running sums for per-cell posterior means
  cell_foi <- cell_sero <- cell_inf <- numeric(length(pop))

  foi_n <- seroprev_n <- infect_n <- numeric(n)
  for (d in seq_len(n)) {
    lat <- samp[[d]]$latent
    b0  <- lat[grep("^Intercept",     rownames(lat)), 1]
    fld <- lat[grep("^spatial.field", rownames(lat)), 1]
    lam <- exp(b0 + as.numeric(Ap_ok %*% fld))
    S   <- suscept(lam)
    inf <- pop * lam * S

    foi_n[d]      <- sum(w * lam)
    seroprev_n[d] <- sum(w * (1 - S))
    infect_n[d]   <- sum(inf)

    cell_foi  <- cell_foi  + lam
    cell_sero <- cell_sero + (1 - S)
    cell_inf  <- cell_inf  + inf

    if (do_reg) {
      reg_foi[, d]  <- rowsum(pop * lam,     gi)[, 1] / reg_pop
      reg_sero[, d] <- rowsum(pop * (1 - S), gi)[, 1] / reg_pop
      reg_inf[, d]  <- rowsum(inf,           gi)[, 1]
    }
  }

  summ <- function(x) c(mean = mean(x),
                        ciL  = unname(quantile(x, 0.025)),
                        ciU  = unname(quantile(x, 0.975)))
  summ_rows <- function(m) {
    out <- t(apply(m, 1, summ))
    data.frame(region = rownames(m), out, row.names = NULL)
  }

  out <- list(
    foi = summ(foi_n), seroprev = summ(seroprev_n),
    infections = summ(infect_n), pop_covered = sum(pop),
    cells = data.frame(res$coop[ok, , drop = FALSE], population = pop,
                       foi = cell_foi / n, seroprev = cell_sero / n,
                       infections = cell_inf / n)
  )

  if (do_reg) {
    stopifnot(isTRUE(all.equal(colSums(reg_inf), infect_n)))   # regions add to nation, per draw
    out$regional <- list(
      foi        = summ_rows(reg_foi),
      seroprev   = summ_rows(reg_sero),
      infections = summ_rows(reg_inf),
      population = data.frame(region = levels(grp), population = reg_pop)
    )
  }
  out
} 

ns_spatial <- national_metrics_summary_new(onnv_results_pop_grid, cam_pop, age_groups, w_age, n = 1000, cell_region = cell_region)
# save national summary
saveRDS(ns_spatial, here('Results/national_summary_spatial.rds'))

ns_spatial <- readRDS(here('Results/national_summary_spatial.rds'))
# print
ns_spatial$foi
ns_spatial$seroprev
ns_spatial$infections

seropos_mean <- ns_spatial$seroprev[['mean']] * ns_spatial$pop_covered
seropos_ciL  <- ns_spatial$seroprev[['ciL']]  * ns_spatial$pop_covered
seropos_ciU  <- ns_spatial$seroprev[['ciU']]  * ns_spatial$pop_covered

cat(sprintf(
  "We estimated that an average of %s (95%% CI: %s-%s)
  individuals get infected each year and that in 2020,
  %s (95%% CI: %s-%s) individuals had a history of ONNV infection,
  representing %.1f%% (95%% CI: %.1f%%-%.1f%%) of the population.",
  formatC(ns_spatial$infections[['mean']], format = "f", digits = 0),
  formatC(ns_spatial$infections[['ciL']],  format = "f", digits = 0),
  formatC(ns_spatial$infections[['ciU']],  format = "f", digits = 0),
  formatC(seropos_mean, format = "f", digits = 0),
  formatC(seropos_ciL,  format = "f", digits = 0),
  formatC(seropos_ciU,  format = "f", digits = 0),
  ns_spatial$seroprev[['mean']] * 100,
  ns_spatial$seroprev[['ciL']]  * 100,
  ns_spatial$seroprev[['ciU']]  * 100
))

# By region
reg <- ns_spatial$regional
region_table <- reg$population |>
  left_join(rename_with(reg$foi,        ~ paste0("foi_", .x),  -region), by = "region") |>
  left_join(rename_with(reg$seroprev,   ~ paste0("sero_", .x), -region), by = "region") |>
  left_join(rename_with(reg$infections, ~ paste0("inf_", .x),  -region), by = "region")


means_table <- region_table |>
  dplyr::select(region, population, foi_mean, sero_mean, inf_mean)

write.csv(means_table, here('Results/inla_pred_by_region.csv'), row.names = FALSE)

# Cameroon wide prediction (maps + regional foi estimates)
foi_onnv <- predicted_foi(onnv_results_pop_grid, onnv_results_pop_grid$coop, pathogen_name = "ONNV")
sero_onnv <- predicted_seroprevalence( foi_result = foi_onnv, model = onnv_results_pop_grid,
  age_groups = age_groups,
  age_weights = w_age,
  crs = 32633,
  pathogen_name = "ONNV"
)
infections_onnv <- predicted_annual_infections(
  foi_result = foi_onnv,
  model = onnv_results_pop_grid,
  age_groups = age_groups,
  age_weights = w_age,
  cam_pop = cam_pop,  
  crs = 32633,
  pathogen_name = "ONNV"
)






# same aggregation as used for model fitting (10x10 km)
cam_pop_agg <- terra::aggregate(cam_pop, fact = 10, fun = sum, na.rm = TRUE)

foi_region <- aggregate_predictions_by_region(
  pred_sf   = foi_onnv$foi_sf,
  regions_sf = regions_sf,
  value_col = "foi",
  agg_type  = "weighted_mean"
)
print(foi_region)

prev_region <- aggregate_predictions_by_region(
  pred_sf    = sero_onnv$prev_sf,
  regions_sf = regions_sf,
  cam_pop    = cam_pop_agg,
  value_col  = "prev",
  agg_type   = "weighted_mean"
)
print(prev_region)
 
infection_region <- aggregate_predictions_by_region(
  pred_sf    = infections_onnv$infections_sf,
  regions_sf = regions_sf,
  value_col  = "infections",
  agg_type   = "sum"
)

print(infection_region)
print(infection_region$region)
sum(infection_region$value_weighted)

# save region level predictions
region_level_predictions <- list(
  foi = foi_region,
  prev = prev_region,
  infections = infection_region
)
saveRDS(region_level_predictions, here('Results/region_level_predictions.rds'))
region_level_predictions <- readRDS('Results/region_level_predictions.rds')


sum(infection_region$total_value)              # ~299,556
sum(infection_region$n_points)                 # 542,883
nrow(infections_onnv$infections_sf)                                  # if larger -> points dropped by the join
sum(infections_onnv$infections_sf$infections, na.rm = TRUE)        # should be ~307,924
sum(is.na(infections_onnv$infections_sf$infections))               # NAs being dropped?



# ---  prob of disease,  acute cases, arthralgic cases and  deaths per year occur in Cameroon each year
# --- Using CHIKV estimates 
average_annual_infections <- ns_spatial$infections[['mean']]
ciL_annual_infections <- ns_spatial$infections[['ciL']]
ciU_annual_infections <- ns_spatial$infections[['ciU']]

severe_acute_cases <- 1.13 / 100 #severe_acute_cases = cases that were detected by surveillance systems (Oscar CHIK paper) - 1.13%
prob_disease <- 0.5 # assuming same as CHIKV
prob_mild <- 0.88 
prob_severe <- 0.12 
prob_medically_attended <- 0.0113
prob_chronic_given_severe <- 0.44 

# 1) acute cases 
acute_cases <- average_annual_infections * prob_disease 
acute_cases_ciL <- ciL_annual_infections * prob_disease 
acute_cases_ciU <- ciU_annual_infections * prob_disease 

cat("Estimated number of acute cases per year:", acute_cases)


# 2) severe acute cases
severe_acute_cases <- average_annual_infections * severe_acute_cases 
severe_acute_cases_ciL <- ciL_annual_infections * severe_acute_cases 
severe_acute_cases_ciU <- ciU_annual_infections * severe_acute_cases 

cat("Estimated number of severe acute cases per year:", severe_acute_cases)


# 3) arthralgic cases
arthralgic_cases  <- average_annual_infections  * prob_medically_attended/2 

cat("Estimated number of arthralgic cases per year:", arthralgic_cases)

# 4) deathts
prob_ifr <- 4.2 * 10^-5 #in Brazil, 

death <- average_annual_infections  * prob_ifr
death_ciL <- ciL_annual_infections * prob_disease * prob_ifr
death_ciU <- ciU_annual_infections * prob_disease * prob_ifr

cat("Estimated number of deaths per year:", death)

# chronic arthlgic cases
#arthralgic_cases  <- average_annual_infections * prob_disease * prob_severe * prob_chronic_given_severe
#arthralgic_cases_ciL  <- ciL_annual_infections * prob_disease * prob_severe * prob_chronic_given_severe
#arthralgic_cases_ciU  <- ciU_annual_infections * prob_disease * prob_severe * prob_chronic_given_severe


# INLA related Supplementary material
# 1) Region level prediction 
region_level_predictions <- readRDS(here("Results/region_level_predictions.rds"))


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


breaks <- seq(
  min(meta_data_with_labels$pf_incidence, na.rm = TRUE),
  max(meta_data_with_labels$pf_incidence, na.rm = TRUE),
  length.out = 8
)

onnv_pf_incidence <- calculate_prop_by_variable(
  data = meta_data_with_labels,
  var_col = "pf_incidence",
  positive_col = "ONNV_pos",
  breaks_max = breaks[-1],
  breaks_min = breaks[-length(breaks)]
)

summary(onnv_pf_incidence$log_model)

global(pf, fun = c("min", "max"), na.rm = TRUE)     # all 4 bands, so we can see if 1 vs 2-4 differ in scale
summary(meta_data_with_labels$pf_incidence)          # range of the actual extracted covariate
summary(onnv_pf_incidence$log_model)                 # or at least exp(coef(onnv_pf_incidence$log_model))

prop_pf_onnv <- make_plot_onnv(
  onnv_pf_incidence$obs,
  meta_data_with_labels$pf_incidence,
  "P. falciparum incidence rate",
  color ="#16622b", pos_col =  "ONNV_pos"
)
print(prop_pf_onnv)


ggsave("Results/supplementary_fig5.png", 
       plot = prop_pf_onnv,
       width = 6, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")




# Covariate INLA model - FOI map 
foi_onnv_covariate_model <- predicted_foi(onnv_results_pop_grid_multivariable, onnv_results_pop_grid_multivariable$coop, pathogen_name = "ONNV")

foi_onnv_covariate_model$plot

# save plot 
ggsave(here('Results/Supplementary_MultivariableSpatialAnalysis_FOImap.png'), 
       plot = foi_onnv_covariate_model$plot,
       width = 10, 
       height = 7, 
       units = "in", 
       dpi = 300,
       bg = "white")


# interpret covariate model results
interpret_covars <- function(res, digits = 3) {
  sf <- res$output$summary.fixed
  sf <- sf[rownames(sf) != "Intercept", , drop = FALSE]      # drop intercept

  covars <- res$covars
  sds    <- res$covar_sds[covars]                            # SDs used to standardise

  out <- data.frame(
    covariate      = covars,
    # per 1 SD (the standardised coefficient) — FOI rate ratio
    RR_perSD       = round(exp(sf$mean), digits),
    RR_perSD_ciL   = round(exp(sf$`0.025quant`), digits),
    RR_perSD_ciU   = round(exp(sf$`0.975quant`), digits),
    # per 1 natural unit of the covariate (divide coef by its SD, then exp)
    RR_perUnit     = round(exp(sf$mean / sds), digits),
    RR_perUnit_ciL = round(exp(sf$`0.025quant` / sds), digits),
    RR_perUnit_ciU = round(exp(sf$`0.975quant` / sds), digits),
    row.names = NULL
  )

  # readable sentences
  for (i in seq_len(nrow(out))) {
    incon <- out$RR_perSD_ciL[i] < 1 & out$RR_perSD_ciU[i] > 1
    cat(sprintf(
      "%s: a 1-SD increase is associated with a %.0f%% %s in FOI (rate ratio %.2f, 95%% CrI %.2f-%.2f)%s\n",
      out$covariate[i],
      abs(out$RR_perSD[i] - 1) * 100,
      ifelse(out$RR_perSD[i] >= 1, "increase", "decrease"),
      out$RR_perSD[i], out$RR_perSD_ciL[i], out$RR_perSD_ciU[i],
      ifelse(incon, " — CrI crosses 1, not clearly distinguishable from no effect", "")
    ))
  }
  invisible(out)
}

covar_effects <- interpret_covars(onnv_results_pop_grid_multivariable)
covar_effects
