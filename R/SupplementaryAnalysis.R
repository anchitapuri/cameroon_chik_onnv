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
library(INLA)
library(tidyverse)



# --- Source functions
source(here('R/Functions.R'))
# --- Read data files 
# original data
data_dir <- 'Data'
meta_data <- read.csv(here(data_dir, 'MFI_meta.csv'))

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
model_data <- read.csv(here('Results/model_data.csv'))



# Multivaribale INLA analysis 
# 1) An. Gambaie + Log Population 
# 2) An. Gambaie
# 3) Aedes Aegypti
# 4) Aedes albopictus
# Model comparison (AIC) between An. Gambaie, Aedes Aegypti and Albopictus 


# --- Run INLA model with covariates - An. Gam and log pop density
# extract covariates for each grid cell
cam_pop_agg <- terra::aggregate(cam_pop, fact = 10, fun = sum, na.rm = TRUE)   # ~1 km cells
dens_agg    <- cam_pop_agg / terra::cellSize(cam_pop_agg, unit = "km")         # people / km²

# Covariates at observation points for fitting
obs_pts <- terra::vect(model_data, geom = c("Longitude","Latitude"), crs = "EPSG:4326")

model_data$gam_ras      <- terra::extract(anopheles_gambiae, obs_pts)[, 2] # [,2] = values, [,1] = ID
model_data$log_dens_ras <- log(terra::extract(dens_agg,      obs_pts)[, 2] + 1)
model_data$aeg_ras      <- terra::extract(aegypti,    obs_pts)[, 2]
model_data$albo_ras     <- terra::extract(albopictus, obs_pts)[, 2]

# Covariates at prediction grid points -> for the MAP ---------------
# align with coop (prediction grid)
grid_pts <- terra::as.points(cam_pop_agg, values = TRUE, na.rm = TRUE)

covar_grid <- data.frame(
  gam_ras      = terra::extract(anopheles_gambiae, grid_pts)[, 2],
  log_dens_ras = log(terra::extract(dens_agg,      grid_pts)[, 2] + 1),
  aeg_ras      = terra::extract(aegypti,     grid_pts)[, 2],
  albo_ras     = terra::extract(albopictus,  grid_pts)[, 2]
)

# fill any grid NAs (cells where a raster has no data) so INLA give error
covar_grid$gam_ras[is.na(covar_grid$gam_ras)]  <- mean(covar_grid$gam_ras, na.rm = TRUE)
covar_grid$log_dens_ras[is.na(covar_grid$log_dens_ras)] <- mean(covar_grid$log_dens_ras, na.rm = TRUE)
covar_grid$aeg_ras[is.na(covar_grid$aeg_ras)]   <- mean(covar_grid$aeg_ras, na.rm = TRUE)
covar_grid$albo_ras[is.na(covar_grid$albo_ras)] <- mean(covar_grid$albo_ras, na.rm = TRUE)


sum(is.na(covar_grid$gam_ras))         # should be 0
sum(is.na(covar_grid$log_dens_ras))    # should be 0
sum(is.na(covar_grid$aeg_ras))         # should be 0
sum(is.na(covar_grid$albo_ras))        # should be 0


onnv_results_pop_grid_multivariable_gam_pop <- run_inla_multivariable(
  year_intro   = 1900,
  data         = model_data,
  cam_pop      = cam_pop,
  positive_col = "ONNV_pos",
  covars       = c("gam_ras", "log_dens_ras"),
  covar_grid   = covar_grid
) 


onnv_results_pop_grid_multivariable_gam <- run_inla_multivariable(
  year_intro   = 1900,
  data         = model_data,
  cam_pop      = cam_pop,
  positive_col = "ONNV_pos",
  covars       = "gam_ras",
  covar_grid   = covar_grid
) 

onnv_results_pop_grid_multivariable_aegypti <- run_inla_multivariable(
  year_intro   = 1900,
  data         = model_data,
  cam_pop      = cam_pop,
  positive_col = "ONNV_pos",
  covars       = "aeg_ras",
  covar_grid   = covar_grid
) 


onnv_results_pop_grid_multivariable_albopictus <- run_inla_multivariable(
  year_intro   = 1900,
  data         = model_data,
  cam_pop      = cam_pop,
  positive_col = "ONNV_pos",
  covars       = "albo_ras",
  covar_grid   = covar_grid
) 


# save 
saveRDS(onnv_results_pop_grid_multivariable_gam_pop, here('Results/ONNV_INLAResults_multivariable_gam_pop.rds'))
saveRDS(onnv_results_pop_grid_multivariable_gam, here('Results/ONNV_INLAResults_multivariable_gam.rds'))
saveRDS(onnv_results_pop_grid_multivariable_aegypti, here('Results/ONNV_INLAResults_multivariable_aegypti.rds'))
saveRDS(onnv_results_pop_grid_multivariable_albopictus, here('Results/ONNV_INLAResults_multivariable_albopictus.rds'))

# beta values 
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

# interpret covariate model results (cloglog link + log(age) offset => exp(beta) = FOI rate ratio)
interpret_covars_new <- function(res, digits = 3) {
  covars   <- res$covars
  covars_z <- paste0(covars, "_z")
  sf  <- res$output$summary.fixed[covars_z, , drop = FALSE]   # select by name, not position
  sds <- res$covar_sds[covars]                                # SDs used to standardise

  # posterior median and 95% CrI of beta (standardised scale)
  b_med <- sf$`0.5quant`
  b_lo  <- sf$`0.025quant`
  b_hi  <- sf$`0.975quant`

  # posterior probability that beta > 0
  P_pos <- sapply(covars_z, function(v)
    1 - INLA::inla.pmarginal(0, res$output$marginals.fixed[[v]]))

  out <- data.frame(
    covariate      = covars,
    # beta per 1 SD (the coefficient INLA actually estimated)
    beta_perSD     = round(b_med, digits),
    beta_perSD_sd  = round(sf$sd, digits),
    beta_perSD_ciL = round(b_lo, digits),
    beta_perSD_ciU = round(b_hi, digits),
    # beta per 1 natural unit (signif, since these can be very small)
    beta_perUnit     = signif(b_med / sds, digits),
    beta_perUnit_sd  = signif(sf$sd / sds, digits),
    beta_perUnit_ciL = signif(b_lo / sds, digits),
    beta_perUnit_ciU = signif(b_hi / sds, digits),
    # FOI rate ratio per 1 SD
    RR_perSD       = round(exp(b_med), digits),
    RR_perSD_ciL   = round(exp(b_lo), digits),
    RR_perSD_ciU   = round(exp(b_hi), digits),
    # FOI rate ratio per 1 natural unit
    RR_perUnit     = round(exp(b_med / sds), digits),
    RR_perUnit_ciL = round(exp(b_lo / sds), digits),
    RR_perUnit_ciU = round(exp(b_hi / sds), digits),
    P_beta_pos     = round(P_pos, digits),
    row.names = NULL
  )

  # readable sentences (computed from unrounded values)
  for (i in seq_len(nrow(out))) {
    rr    <- exp(b_med[i])
    incon <- b_lo[i] < 0 & b_hi[i] > 0
    cat(sprintf(
      paste0("%s: beta = %.3f (95%% CrI %.3f to %.3f) per SD; ",
             "a 1-SD increase is associated with a %.0f%% %s in FOI ",
             "(rate ratio %.2f, 95%% CrI %.2f-%.2f; P(beta > 0) = %.2f)%s\n"),
      covars[i], b_med[i], b_lo[i], b_hi[i],
      abs(rr - 1) * 100, ifelse(rr >= 1, "increase", "decrease"),
      rr, exp(b_lo[i]), exp(b_hi[i]), P_pos[i],
      ifelse(incon, "; CrI crosses 1, not clearly distinguishable from no effect", "")
    ))
  }
  invisible(out)
}


covar_effects_gam <- interpret_covars_new(onnv_results_pop_grid_multivariable_gam)

covar_effects_aeg <- interpret_covars_new(onnv_results_pop_grid_multivariable_aegypti)

covar_effects_albo <- interpret_covars_new(onnv_results_pop_grid_multivariable_albopictus)

# beta 
fit <- onnv_results_pop_grid_multivariable_gam$model   # adjust name
# Posterior summary for beta (mean, sd, 2.5%, median, 97.5%, mode)
beta_gam <- onnv_results_pop_grid_multivariable_gam$output$summary.fixed["gam_ras", ]
beta_gam

# Just the point estimate and 95% credible interval
beta_gam[, c("mean", "0.025quant", "0.975quant")]


# national foi prediction
national_foi_mv_gam_pop <- national_foi_covars(onnv_results_pop_grid_multivariable_gam_pop, cam_pop,
                                       covar_grid = covar_grid)

national_foi_mv_gam <- national_foi_covars(onnv_results_pop_grid_multivariable_gam, cam_pop,
                                       covar_grid = covar_grid)

national_foi_mv_gam_pop_aeg <- national_foi_covars(onnv_results_pop_grid_multivariable_aegypti, cam_pop,
                                       covar_grid = covar_grid)

national_foi_mv_gam_pop_albo <- national_foi_covars(onnv_results_pop_grid_multivariable_albopictus, cam_pop,
                                       covar_grid = covar_grid)


# save 
saveRDS(national_foi_mv_gam_pop, here('Results/national_foi_mv_gam_pop.rds'))
saveRDS(national_foi_mv_gam, here('Results/national_foi_mv_gam.rds'))
saveRDS(national_foi_mv_gam_pop_aeg, here('Results/national_foi_mv_gam_pop_aeg.rds'))
saveRDS(national_foi_mv_gam_pop_albo, here('Results/national_foi_mv_gam_pop_albo.rds'))

national_foi_mv_gam <- readRDS(here('Results/national_foi_mv_gam.rds'))
national_foi_mv_gam_pop_aeg <- readRDS(here('Results/national_foi_mv_gam_pop_aeg.rds'))
national_foi_mv_gam_pop_albo <- readRDS(here('Results/national_foi_mv_gam_pop_albo.rds'))


# print
national_foi_mv_gam_pop[c("mean","ciL","ciU")] 
national_foi_mv_gam[c("mean","ciL","ciU")]
national_foi_mv_gam_pop_aeg[c("mean","ciL","ciU")]
national_foi_mv_gam_pop_albo[c("mean","ciL","ciU")]


# AIC 
models <- list(
  "An. gambiae"   = onnv_results_pop_grid_multivariable_gam,
  "Ae. aegypti"   = onnv_results_pop_grid_multivariable_aegypti,
  "Ae. albopictus" = onnv_results_pop_grid_multivariable_albopictus
)

dic_table <- data.frame(
  model = names(models),
  DIC   = sapply(models, function(m) m$dic),
  n     = sapply(models, function(m) nrow(m$data_filtered))
)
dic_table$delta_DIC <- dic_table$DIC - min(dic_table$DIC)
dic_table <- dic_table[order(dic_table$DIC), ]
dic_table

dic_plot_data <- data.frame(
  species = c("Anopheles\ngambiae", "Aedes\naegypti", "Aedes\nalbopictus"),
  DIC     = c(onnv_results_pop_grid_multivariable_gam$dic,
              onnv_results_pop_grid_multivariable_aegypti$dic,
              onnv_results_pop_grid_multivariable_albopictus$dic)
)
dic_plot_data$species <- factor(dic_plot_data$species, levels = dic_plot_data$species)

ggplot(dic_plot_data, aes(x = species, y = DIC, group = 1)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.5, vjust = 1.5, size = 5) +
  labs(x = "Vector Species", y = "DIC") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic"))



# Covariate INLA model - FOI map 
foi_onnv_gam <- predicted_foi(onnv_results_pop_grid_multivariable_gam,
                              onnv_results_pop_grid_multivariable_gam$coop, pathogen_name = "ONNV")
foi_onnv_aegypti <- predicted_foi(onnv_results_pop_grid_multivariable_aegypti,
                                  onnv_results_pop_grid_multivariable_aegypti$coop, pathogen_name = "ONNV")
foi_onnv_albopictus <- predicted_foi(onnv_results_pop_grid_multivariable_albopictus,
                                     onnv_results_pop_grid_multivariable_albopictus$coop, pathogen_name = "ONNV")

# Italic species name under each map
label_map <- function(p, species) {
  p + labs(caption = species) +
    theme(plot.caption = element_text(face = "italic", hjust = 0.5, size = 18))
}

# --- DIC plot (same style as before, without the manual "E")
dic_plot <- ggplot(dic_plot_data, aes(x = species, y = DIC, group = 1)) +
  geom_line(linetype = "dashed") +
  geom_point() +
  labs(x = "Vector Species", y = "DIC") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, face = "italic"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(t = 10, r = 40, b = 10, l = 10, unit = "pt")
  )


# Italic species name as its own text element
species_label <- function(species, size = 22) {
  wrap_elements(full = grid::textGrob(species,
    gp = grid::gpar(fontface = "italic", fontsize = size)))
}

gap <- 0.15   # space between maps (relative to map width)

map_row <- foi_onnv_gam$plot + plot_spacer() + foi_onnv_aegypti$plot + plot_spacer() +
  foi_onnv_albopictus$plot +
  plot_layout(nrow = 1, widths = c(1, gap, 1, gap, 1))

label_row <- species_label("Anopheles gambiae") + plot_spacer() + species_label("Aedes aegypti") +
  plot_spacer() + species_label("Aedes albopictus") +
  plot_layout(nrow = 1, widths = c(1, gap, 1, gap, 1))

combined_DIC_maps_plot <- label_row / map_row  / dic_plot +
  plot_layout(heights = c(0.12, 1.3, 1))


ggsave(here("Results/onnv_maps_covariate_model_comparison.png"), combined_DIC_maps_plot, width = 12, height = 8, dpi = 300)




# Imputation INLA - all analysis + results in INLA_SpatialPrediction_Imputation.R



# All antigens + correlation 

full_data <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/base_complete_MFI_meta.csv')
colnames(full_data)

antigens <- setdiff(names(full_data)[3:26], c("ONNV_VLP", "CHIKV_sE2", "MAYV_E2"))

print(antigens)

long <- full_data %>%
  dplyr::select(all_of(antigens), ONNV_VLP) %>%
  pivot_longer(-ONNV_VLP, names_to = "antigen", values_to = "value") %>%
  filter(!is.na(value), !is.na(ONNV_VLP)) %>%
  mutate(antigen = factor(antigen, levels = antigens))

# Correlation + significance stars per antigen
cors <- long %>%
  group_by(antigen) %>%
  summarise(
    r = cor(value, ONNV_VLP, method = "pearson"),
    p = cor.test(value, ONNV_VLP, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    stars = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""),
    label = sprintf("Corr: %.2f%s", r, stars)
  )
quartz()
panel_labs <- setNames(paste0(LETTERS[seq_along(antigens)], ")   ", antigens), antigens)
corr_nonalpha_antigens <- ggplot(long, aes(x = log(value), y = log(ONNV_VLP))) +
  geom_point(colour = "#1F3A93", alpha = 0.7, size = 1.6) +
  geom_text(data = cors, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5, size = 3.2, colour = "grey30",
            inherit.aes = FALSE) +
  facet_wrap(~ antigen, scales = "free_x", ncol = 7, labeller = as_labeller(panel_labs)) +
  labs(x = "Antigen response", y = "ONNV_VLP") +
  theme_minimal(base_size = 11) +
  theme(panel.border = element_rect(fill = NA, colour = "grey85"),
    strip.text = element_text(hjust = 0)
  )

ggsave(here("Results/corr_nonalpha_antigens.png"), corr_nonalpha_antigens, width = 12, height = 8, dpi = 300)
