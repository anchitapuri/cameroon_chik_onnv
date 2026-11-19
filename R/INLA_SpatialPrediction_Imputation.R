
 
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


# Full imputation run
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

ns_imputed <- readRDS(here('Results/national_summary_imputed.rds'))

print(ns_imputed$foi)
print(ns_imputed$seroprev)
print(ns_imputed$infections)

# How much of the total uncertainty comes from classification rather than from
# the spatial model? Quote this in the response to the reviewer.
print(ns_imputed$variance_decomposition)

# Spread of crude prevalence across imputations -- a direct readout of how much
# the serostatus labels move.
print(ns_imputed$imputed_prevalence)



# Estimates 
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

