


# Keep the single-fit run_inla() calls above this: they are still used for the
# multivariable model, the ONNV-samples sensitivity fit, and model checking.
# ==============================================================================
 
library(posterior)
library(stringr)


# --- Source functions
source(here('R/INLAImputation_Functions.R'))

# re-run
fit_full_model <- readRDS(here("Results", "full_model_fits.rds"))
meta_data_with_labels <- read.csv(here("Results", "meta_data_with_labels.csv"))
preprocessed_data_full_model <- readRDS(here('Results/preprocessed_data_full_model.rds'))

model_data <- read.csv(here("Results", "model_data.csv"))
colnames(model_data)
# model_data must be in Stan index order so p_onnv lines up row for row
model_data <- model_data[order(model_data$stan_idx_full_model), ]
stopifnot(length(p_onnv) == nrow(model_data))
stopifnot(!is.unsorted(model_data$stan_idx_full_model))



p_onnv <- get_onnv_prob(fit_full_model, onnv_components = c(2, 4), method = "mean")

# the mean probability should be close to the hard-assigned prevalence
mean(p_onnv)
mean(model_data$ONNV_pos)



# how many individuals are genuinely uncertain?
sum(p_onnv[meta_data_with_labels$ONNV_pos == 1] < 0.7)


# 2. Time a single fit before committing to the full loop
system.time({
  test_run <- run_inla_imputed(
    p_onnv     = p_onnv,
    data       = model_data,
    cam_pop    = cam_pop,
    age_groups = age_groups,
    w_age      = w_age,
    M          = 2,
    n_draws    = 50
  )
})
# Total runtime scales as M x (one INLA fit). If a fit takes t minutes,
# budget M * t minutes. M = 20 is usually enough for a stable interval.


# ------------------------------------------------------------------------------
# 3. Full imputation run
# ------------------------------------------------------------------------------

onnv_imputed <- run_inla_imputed(
  p_onnv     = p_onnv,
  data       = model_data,
  cam_pop    = cam_pop,
  age_groups = age_groups,
  w_age      = w_age,
  M          = 20,
  n_draws    = 200,
  year_intro = 1900,
  seed       = 2024
)

saveRDS(onnv_imputed, here("Results", "ONNV_INLA_imputed.rds"))
# onnv_imputed <- readRDS(here("Results", "ONNV_INLA_imputed.rds"))

mean(onnv_imputed$cell_foi)
mean(onnv_imputed$cell_seroprev, na.rm = TRUE)
sum(onnv_imputed$cell_infections, na.rm = TRUE)

ns_imputed <- summarise_imputed(onnv_imputed)
saveRDS(ns_imputed, here("Results", "national_summary_imputed.rds"))

print(ns_imputed$foi)
print(ns_imputed$seroprev)
print(ns_imputed$infections)

# How much of the total uncertainty comes from classification rather than from
# the spatial model? Quote this in the response to the reviewer.
print(ns_imputed$variance_decomposition)

# Spread of crude prevalence across imputations -- a direct readout of how much
# the serostatus labels move.
print(ns_imputed$imputed_prevalence)


# ------------------------------------------------------------------------------
# 4. Headline numbers
# ------------------------------------------------------------------------------

total_pop <- ns_imputed$pop_covered
seropos_mean <- ns_imputed$seroprev[["mean"]] * total_pop
seropos_ciL  <- ns_imputed$seroprev[["ciL"]]  * total_pop
seropos_ciU  <- ns_imputed$seroprev[["ciU"]]  * total_pop

cat(sprintf(
  "We estimated that an average of %s (95%% CrI: %s-%s) individuals are infected
  each year and that %s (95%% CrI: %s-%s) individuals had a history of ONNV
  infection, representing %.1f%% (95%% CrI: %.1f%%-%.1f%%) of the population.\n",
  formatC(ns_imputed$infections[["mean"]], format = "f", digits = 0, big.mark = ","),
  formatC(ns_imputed$infections[["ciL"]],  format = "f", digits = 0, big.mark = ","),
  formatC(ns_imputed$infections[["ciU"]],  format = "f", digits = 0, big.mark = ","),
  formatC(seropos_mean, format = "f", digits = 0, big.mark = ","),
  formatC(seropos_ciL,  format = "f", digits = 0, big.mark = ","),
  formatC(seropos_ciU,  format = "f", digits = 0, big.mark = ","),
  ns_imputed$seroprev[["mean"]] * 100,
  ns_imputed$seroprev[["ciL"]]  * 100,
  ns_imputed$seroprev[["ciU"]]  * 100
))


# ------------------------------------------------------------------------------
# 5. Regional estimates from the SAME per-cell surface
# ------------------------------------------------------------------------------

surfaces_sf <- imputed_surfaces_sf(onnv_imputed)

cameroon_regions <- ne_states(country = "Cameroon", returnclass = "sf")
regions_sf <- cameroon_regions %>%
  dplyr::select(region = name) %>%
  st_make_valid()

infection_region <- aggregate_predictions_by_region(
  pred_sf = surfaces_sf, regions_sf = regions_sf,
  value_col = "infections", agg_type = "sum")

foi_region <- aggregate_predictions_by_region(
  pred_sf = surfaces_sf, regions_sf = regions_sf,
  value_col = "foi", agg_type = "weighted_mean")

prev_region <- aggregate_predictions_by_region(
  pred_sf = surfaces_sf, regions_sf = regions_sf,
  value_col = "seroprevalence", agg_type = "weighted_mean")

print(infection_region)
print(foi_region)
print(prev_region)

# Regional infections must now sum to the national mean exactly.
stopifnot(isTRUE(all.equal(sum(infection_region$total_value),
                           ns_imputed$infections[["mean"]],
                           tolerance = 1e-6)))

# Population-weighted regional seroprevalence must reconcile with the national
# figure (weighted means, so weight by each region's population).
stopifnot(isTRUE(all.equal(
  weighted.mean(prev_region$value_weighted, prev_region$total_population),
  ns_imputed$seroprev[["mean"]], tolerance = 1e-3)))

region_level_predictions <- list(
  foi = foi_region, prev = prev_region, infections = infection_region)
saveRDS(region_level_predictions, here("Results", "region_level_predictions.rds"))


# ------------------------------------------------------------------------------
# 6. Burden estimates -- apply the ratios to the DRAWS, not the point estimate
#
# The original code applied each probability to the mean and then to the CI
# bounds separately. Applying them to every draw propagates the uncertainty
# correctly and avoids the overwriting bug in the severe-case block.
# ------------------------------------------------------------------------------

summ <- function(x) c(mean = mean(x),
                      ciL  = unname(quantile(x, 0.025)),
                      ciU  = unname(quantile(x, 0.975)))

inf_d <- onnv_imputed$infect_draws

burden <- list(
  acute       = summ(inf_d * prob_disease),
  severe      = summ(inf_d * prob_severe),
  arthralgic  = summ(inf_d * prob_medically_attended / 2),
  deaths      = summ(inf_d * prob_ifr)
)
print(lapply(burden, round))
saveRDS(burden, here("Results", "burden_estimates_imputed.rds"))


# ------------------------------------------------------------------------------
# 7. Maps built from the same surface
# ------------------------------------------------------------------------------

p_infections <- ggplot() +
  geom_sf(data = surfaces_sf, aes(color = infections), size = 1.5, alpha = 0.5) +
  scale_color_viridis_c(option = "mako", trans = "log10",
                        name = "Annual infections",
                        labels = scales::comma_format()) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        axis.text = element_blank(), axis.ticks = element_blank())
quartz()
print(p_infections)
