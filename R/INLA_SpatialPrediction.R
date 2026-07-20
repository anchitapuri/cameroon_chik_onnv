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
anopheles_funestus <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/2010_Anopheles_funestus_CMR.tiff')
anopheles_gambiae <- rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/2010_Anopheles_gambiae_ss_CMR.tiff')
cam_pop <- rast("/Users/ap2488/Desktop/Cameroon_Analysis_2025/cmr_ppp_2020_UNadj.tif")
sum(values(cam_pop), na.rm = TRUE)
cam_pop_agg <- terra::aggregate(cam_pop, fact = 10, fun = sum, na.rm = TRUE) # (~1 km, from aggregating the ~100 m WorldPop grid by a factor of 10)".

# Cameroon population by age
cameroon_age_2025 <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/CameroonAge2025.csv')
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

# Rename
model_data <- meta_data_with_labels
model_data_onnv_samples <- meta_data_onnv_samples


colnames(model_data)

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

## 2. Covariates at OBSERVATION points  -> for FITTING ------------------
##    (one value per point = your district locations; overwrites the cols you already have)
obs_pts <- terra::vect(model_data, geom = c("Longitude","Latitude"), crs = "EPSG:4326")

model_data$gam_ras      <- terra::extract(anopheles_gambiae, obs_pts)[, 2]
model_data$log_dens_ras <- log(terra::extract(dens_agg,      obs_pts)[, 2] + 1)

## 3. Covariates at PREDICTION grid points -> for the MAP ---------------
##    built from the SAME aggregate + as.points(na.rm=TRUE) the function uses,
##    so rows align 1:1 with coop
grid_pts <- terra::as.points(cam_pop_agg, values = TRUE, na.rm = TRUE)

covar_grid <- data.frame(
  gam_ras      = terra::extract(anopheles_gambiae, grid_pts)[, 2],
  log_dens_ras = log(terra::extract(dens_agg,      grid_pts)[, 2] + 1)
)

# fill any grid NAs (cells where a raster has no data) so INLA doesn't choke
covar_grid$gam_ras[is.na(covar_grid$gam_ras)]           <- mean(covar_grid$gam_ras, na.rm = TRUE)
covar_grid$log_dens_ras[is.na(covar_grid$log_dens_ras)] <- mean(covar_grid$log_dens_ras, na.rm = TRUE)


onnv_results_pop_grid_multivariable <- run_inla_multivariable(
  year_intro   = 1900,
  data         = model_data,
  cam_pop      = cam_pop,
  positive_col = "ONNV_pos",
  covars       = c("gam_ras", "log_dens_ras"),
  covar_grid   = covar_grid
) 
foi_mv <- predicted_foi(onnv_results_pop_grid_multivariable,
                        coop = onnv_results_pop_grid_multivariable$coop)
foi_mv$plot

covar_table <- function(res, digits = 2) {
  sf <- res$output$summary.fixed
  sf <- sf[rownames(sf) != "Intercept", , drop = FALSE]   # drop intercept

  data.frame(
    Covariate = sub("_z$", " (per 1 SD)", rownames(sf)),
    RateRatio = round(exp(sf$mean), digits),
    ciL       = round(exp(sf$`0.025quant`), digits),
    ciU       = round(exp(sf$`0.975quant`), digits),
    row.names = NULL
  )
}
covar_table(onnv_results_pop_grid_multivariable)

# --- Save prediction results 
saveRDS(onnv_results_pop_grid, here('Results/ONNV_INLAResults.rds'))
saveRDS(onnv_results_pop_grid_onnv_samples, here('Results/ONNV_INLAResults_ONNV_samples.rds'))
saveRDS(onnv_results_pop_grid_multivariable, here('Results/ONNV_INLAResults_multivariable.rds'))

# Read saved results
onnv_results_pop_grid <- readRDS(here('Results/ONNV_INLAResults.rds'))
onnv_results_pop_grid_multivariable <- readRDS(here('Results/ONNV_INLAResults_multivariable.rds'))

# --- SPATIAL PREDICTIONS: # Overall cameroon estimates ---- 
# population weighted national FOI estimate, using the posterior samples from the INLA model
national_foi <- function(res, cam_pop, agg_factor = 10, n = 1000) {

  # population per grid cell, in the SAME order as res$coop
  cam_pop_agg    <- terra::aggregate(cam_pop, fact = agg_factor, fun = sum, na.rm = TRUE)
  cam_pop_points <- terra::as.points(cam_pop_agg, values = TRUE, na.rm = TRUE)
  pop_grid       <- as.data.frame(cam_pop_points)[, 1]

  stopifnot(length(pop_grid) == nrow(res$coop))   # alignment check

  # posterior FOI surface: exp(Intercept + spatial.field), covariates = 0 at grid
  samp   <- inla.posterior.sample(n, res$output)
  b0_idx <- grep("^Intercept",     rownames(samp[[1]]$latent))
  sf_idx <- grep("^spatial.field", rownames(samp[[1]]$latent))
  Ap     <- inla.spde.make.A(mesh = res$mesh, loc = as.matrix(res$coop))

  foi_draws <- sapply(samp, function(s) {
    field <- as.numeric(Ap %*% s$latent[sf_idx, 1])
    exp(s$latent[b0_idx, 1] + field)
  })                                               # grid × n

  # population-weighted national FOI, one value per draw
  ok <- !is.na(pop_grid)
  w  <- pop_grid[ok] / sum(pop_grid[ok])
  foi_natl <- colSums(foi_draws[ok, ] * w)

  list(
    mean = mean(foi_natl),
    ciL  = unname(quantile(foi_natl, 0.025)),
    ciU  = unname(quantile(foi_natl, 0.975)),
    surface_mean = rowMeans(foi_draws)             # per-cell FOI map if you want it
  )
}

national_foi_covars <- function(res, cam_pop, covar_grid = NULL, agg_factor = 10, n = 1000) {

  # population per grid cell, aligned to res$coop
  cam_pop_agg    <- terra::aggregate(cam_pop, fact = agg_factor, fun = sum, na.rm = TRUE)
  pop_grid       <- as.data.frame(terra::as.points(cam_pop_agg, values = TRUE, na.rm = TRUE))[, 1]
  stopifnot(length(pop_grid) == nrow(res$coop))
  ok <- !is.na(pop_grid); w <- pop_grid[ok] / sum(pop_grid[ok])

  # standardise grid covariates with the SAME centre/scale used at fitting
  has_cov <- !is.null(covar_grid) && length(res$covars) > 0
  if (has_cov) {
    stopifnot(nrow(covar_grid) == nrow(res$coop))
    covars   <- res$covars
    covars_z <- paste0(covars, "_z")
    Xg <- sapply(seq_along(covars), function(i)
            (covar_grid[[covars[i]]] - res$covar_means[i]) / res$covar_sds[i])
    Xg <- Xg[ok, , drop = FALSE]
  }

  # sample only what we need: intercept, spatial field, and covariate coefficients
  sel <- list(Intercept = 1, spatial.field = seq_len(res$mesh$n))
  if (has_cov) for (cz in covars_z) sel[[cz]] <- 1
  samp  <- inla.posterior.sample(n, res$output, selection = sel)
  Ap_ok <- inla.spde.make.A(mesh = res$mesh, loc = as.matrix(res$coop))[ok, , drop = FALSE]

  foi_natl <- numeric(n); surf_sum <- numeric(sum(ok))
  for (d in seq_len(n)) {
    lat <- samp[[d]]$latent
    b0  <- lat[grep("^Intercept",     rownames(lat)), 1]
    fld <- lat[grep("^spatial.field", rownames(lat)), 1]
    eta <- b0 + as.numeric(Ap_ok %*% fld)
    if (has_cov) {
      beta <- sapply(covars_z, function(cz) lat[grep(paste0("^", cz), rownames(lat)), 1])
      eta  <- eta + as.numeric(Xg %*% beta)               # + covariate contribution
    }
    lam <- exp(eta)
    foi_natl[d] <- sum(w * lam)
    surf_sum    <- surf_sum + lam
  }

  list(mean = mean(foi_natl),
       ciL  = unname(quantile(foi_natl, 0.025)),
       ciU  = unname(quantile(foi_natl, 0.975)),
       surface_mean = surf_sum / n)                        # length = sum(ok)
}


national_foi_spatial <- national_foi(onnv_results_pop_grid, cam_pop)
national_foi_onnv_samples   <- national_foi(onnv_results_pop_grid_onnv_samples, cam_pop)
national_foi_mv   <- national_foi_covars(onnv_results_pop_grid_multivariable, cam_pop)

# save 
saveRDS(national_foi_spatial, here('Results/national_foi_spatial.rds'))
saveRDS(national_foi_mv, here('Results/national_foi_mv.rds'))
saveRDS(national_foi_onnv_samples, here('Results/national_foi_onnv_samples.rds'))


# print values 
national_foi_spatial[c("mean","ciL","ciU")] 
national_foi_mv[c("mean","ciL","ciU")] 
national_foi_onnv_samples[c("mean","ciL","ciU")] 


# -- additional national metrics
national_summary <- function(res, cam_pop, age_groups, w_age,
                             agg_factor = 10, n = 1000, seed = 1) {
  set.seed(seed)
  stopifnot(length(w_age) == nrow(age_groups))       # guard: weights must match age groups

  # per-cell population aligned to res$coop
  cam_pop_agg <- terra::aggregate(cam_pop, fact = agg_factor, fun = sum, na.rm = TRUE)
  pop <- as.data.frame(terra::as.points(cam_pop_agg, values = TRUE, na.rm = TRUE))[, 1]
  stopifnot(length(pop) == nrow(res$coop))
  ok <- !is.na(pop); pop <- pop[ok]; w <- pop / sum(pop)

  # sample ONLY Intercept + spatial.field (this is what stops the crash)
  sel   <- list(Intercept = 1, spatial.field = seq_len(res$mesh$n))
  samp  <- inla.posterior.sample(n, res$output, selection = sel)
  Ap_ok <- inla.spde.make.A(mesh = res$mesh, loc = as.matrix(res$coop))[ok, , drop = FALSE]

  # age-averaged susceptible proportion for a vector of FOIs
  L <- age_groups$age_lower; U <- age_groups$age_upper; dW <- U - L
  suscept <- function(lam) {
    s <- sapply(seq_along(L), function(j)
      ifelse(lam > 1e-10, (1/(lam*dW[j])) * (exp(-lam*L[j]) - exp(-lam*U[j])), 1))
    as.numeric(s %*% w_age)
  }

  # loop over draws — no giant matrix
  foi_n <- seroprev_n <- infect_n <- numeric(n)
  for (d in seq_len(n)) {
    lat <- samp[[d]]$latent
    b0  <- lat[grep("^Intercept",     rownames(lat)), 1]
    fld <- lat[grep("^spatial.field", rownames(lat)), 1]
    lam <- exp(b0 + as.numeric(Ap_ok %*% fld))
    S   <- suscept(lam)
    foi_n[d]      <- sum(w * lam)
    seroprev_n[d] <- sum(w * (1 - S))
    infect_n[d]   <- sum(pop * lam * S)
  }

  summ <- function(x) c(mean = mean(x),
                        ciL = unname(quantile(x, 0.025)),
                        ciU = unname(quantile(x, 0.975)))
  list(foi = summ(foi_n), seroprev = summ(seroprev_n),
       infections = summ(infect_n), pop_covered = sum(pop))
}


ns_spatial <- national_summary(onnv_results_pop_grid, cam_pop, age_groups, w_age, n = 1000)

# save national summary
saveRDS(ns_spatial, here('Results/national_summary_spatial.rds'))

# print
ns_spatial$foi
ns_spatial$seroprev
ns_spatial$infections


# --- Cameroon Wide prediction, aggregated by region 
foi_onnv <- predicted_foi(onnv_results_pop_grid, onnv_results_pop_grid$coop, pathogen_name = "ONNV")
foi_onnv$plot
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


cameroon_regions <- ne_states(country = "Cameroon", returnclass = "sf")
# rename 'name' to region
regions_sf <- cameroon_regions %>%
  dplyr::select(region = name) %>%   # standardise column name
  st_make_valid()

# same aggregation as used for model fitting (10x10 km)
cam_pop_agg <- terra::aggregate(cam_pop, fact = 10, fun = sum, na.rm = TRUE)

foi_region <- aggregate_predictions_by_region(
  pred_sf   = foi_onnv$foi_sf,
  regions_sf = regions_sf,
  cam_pop   = cam_pop_agg,
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

# save region level predictions
region_level_predictions <- list(
  foi = foi_region,
  prev = prev_region,
  infections = infection_region
)
saveRDS(region_level_predictions, here('Results/region_level_predictions.rds'))
print(region_level_predictions)


ns_spatial$infections[['ciL']]

# prob of disease,  acute cases, arthralgic cases and  deaths per year occur in Cameroon each year
average_annual_infections <- ns_spatial$infections[['mean']]
ciL_annual_infections <- ns_spatial$infections[['ciL']]
ciU_annual_infections <- ns_spatial$infections[['ciU']]

severe_acute_cases <- 1.13 / 100 #severe_acute_cases = cases that were detected by surveillance systems (Oscar CHIK paper) - 1.13%
prob_disease <- 0.5 # assuming same as CHIKV
prob_mild <- 0.88 
prob_severe <- 0.12 
prob_medically_attended <- 0.0113
prob_chronic_given_severe <- 0.44 

# acute cases 
acute_cases <- average_annual_infections * prob_disease 
acute_cases_ciL <- ciL_annual_infections * prob_disease 
acute_cases_ciU <- ciU_annual_infections * prob_disease 


cat("Estimated number of acute cases per year:", acute_cases)

# chronic arthlgic cases
#arthralgic_cases  <- average_annual_infections * prob_disease * prob_severe * prob_chronic_given_severe
#arthralgic_cases_ciL  <- ciL_annual_infections * prob_disease * prob_severe * prob_chronic_given_severe
#arthralgic_cases_ciU  <- ciU_annual_infections * prob_disease * prob_severe * prob_chronic_given_severe

arthralgic_cases  <- average_annual_infections  * prob_medically_attended/2 

cat("Estimated number of arthralgic cases per year:", arthralgic_cases)

prob_ifr <- 4.2 * 10^-5 #in Brazil, 

death <- average_annual_infections  * prob_ifr
death_ciL <- ciL_annual_infections * prob_disease * prob_ifr
death_ciU <- ciU_annual_infections * prob_disease * prob_ifr

cat("Estimated number of deaths per year:", death)



# supplementary data - region level prediction 
# supplementary data 
region_level_predictions <- readRDS(here("Results/region_level_predictions.rds"))
region_level_predictions

# comparison of malaria risk with ONNV risk 
# Number of newly diagnosed Plasmodium falciparum cases per 1,000 population (using 2024)

pf <- terra::rast('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/clippedlayers-4/202508_Global_Pf_Incidence_Rate_CMR_2024.tiff')
names(pf)
pf[[1]]
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





# old FOI estimates

# overall FOI
foi_summary <- onnv_results_pop_grid$output$summary.fixed
est_cameroonwide_foi <- list(
  mean = exp(foi_summary$mean),
  ciL  = exp(foi_summary$`0.025quant`),
  ciU  = exp(foi_summary$`0.975quant`)
)
est_cameroonwide_foi


# onnv only samples model
foi_summary_onnv_samples <- onnv_results_pop_grid_onnv_samples$output$summary.fixed
est_cameroonwide_foi_onnv_samples <- list(
  mean = exp(foi_summary_onnv_samples$mean),
  ciL  = exp(foi_summary_onnv_samples$`0.025quant`),
  ciU  = exp(foi_summary_onnv_samples$`0.975quant`)
)

# --- covaraites model 
foi_summary_multivariable <- onnv_results_pop_grid_multivariable$output$summary.fixed
b0 <- foi_summary_multivariable["Intercept", ]

est_cameroonwide_foi_multivariable <- list(
  mean = exp(b0$mean),
  ciL  = exp(b0$`0.025quant`),
  ciU  = exp(b0$`0.975quant`)
)
est_cameroonwide_foi_multivariable


# -- cameroon wide summary 
compute_foi_metrics <- function(foi_val, age_groups, w_age, cam_pop, total_cameroon_pop) {
  
  avg_susceptible_by_age <- numeric(nrow(age_groups))
  
  for (j in 1:nrow(age_groups)) {
    a_lower <- age_groups$age_lower[j]
    a_upper <- age_groups$age_upper[j]
    age_width <- a_upper - a_lower
    
    if (foi_val > 1e-10) {
      avg_susceptible_by_age[j] <- (1 / (foi_val * age_width)) *
        (exp(-foi_val * a_lower) -
           exp(-foi_val * a_upper))
    } else {
      avg_susceptible_by_age[j] <- 1
    }
  }
  
  avg_susceptible  <- sum(avg_susceptible_by_age * w_age)
  avg_seroprev     <- 1 - avg_susceptible
  infections       <- total_cameroon_pop * foi_val * avg_susceptible
  
  list(
    avg_susceptible = avg_susceptible,
    avg_seroprev    = avg_seroprev,
    infections      = infections
  )
}

# --- Run for mean, ciL, ciU
total_cameroon_pop <- sum(values(cam_pop), na.rm = TRUE)

metrics_mean <- compute_foi_metrics(est_cameroonwide_foi$mean, age_groups, w_age, cam_pop, total_cameroon_pop)
metrics_ciL  <- compute_foi_metrics(est_cameroonwide_foi$ciL,  age_groups, w_age, cam_pop, total_cameroon_pop)
metrics_ciU  <- compute_foi_metrics(est_cameroonwide_foi$ciU,  age_groups, w_age, cam_pop, total_cameroon_pop)

# ---summary 
cameroon_summary <- data.frame(
  metric    = c("FOI", "Avg susceptible", "Avg seroprev", "Expected infections"),
  mean      = c(est_cameroonwide_foi$mean,  metrics_mean$avg_susceptible, metrics_mean$avg_seroprev, metrics_mean$infections),
  ciL       = c(est_cameroonwide_foi$ciL,   metrics_ciL$avg_susceptible,  metrics_ciL$avg_seroprev,  metrics_ciL$infections),
  ciU       = c(est_cameroonwide_foi$ciU,   metrics_ciU$avg_susceptible,  metrics_ciU$avg_seroprev,  metrics_ciU$infections)
)

cameroon_summary_fmt <- cameroon_summary
cameroon_summary_fmt[, c("mean", "ciL", "ciU")] <- lapply(
  cameroon_summary[, c("mean", "ciL", "ciU")],
  function(x) formatC(x, format = "f", digits = 3)
)

print(cameroon_summary_fmt)

# number of individuals with a history of ONNV infection (ie seropositive)
seropos_mean <- total_cameroon_pop * metrics_mean$avg_seroprev
seropos_ciL  <- total_cameroon_pop * metrics_ciL$avg_seroprev
seropos_ciU  <- total_cameroon_pop * metrics_ciU$avg_seroprev

# print summary statement
cat(sprintf(
  "We estimated that an average of %s (95%% CI: %s-%s) 
  individuals get infected each year and that in 2020, 
  %s (95%% CI: %s-%s) individuals had a history of ONNV infection, 
  representing %.1f%% (95%% CI: %.1f%%-%.1f%%) of the population.",
  formatC(metrics_mean$infections, format = "f", digits = 0),
  formatC(metrics_ciL$infections,  format = "f", digits = 0),
  formatC(metrics_ciU$infections,  format = "f", digits = 0),
  formatC(seropos_mean, format = "f", digits = 0),
  formatC(seropos_ciL,  format = "f", digits = 0),
  formatC(seropos_ciU,  format = "f", digits = 0),
  metrics_mean$avg_seroprev * 100,
  metrics_ciL$avg_seroprev  * 100,
  metrics_ciU$avg_seroprev  * 100
))

