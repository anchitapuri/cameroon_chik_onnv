library(ggraph)
library(lhs)
library(matrixStats)
library(mvtnorm)
library(matrixcalc)
library(here)
library(emdbook)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(matrixStats)
library(stringr)
library(bayesplot)
library(posterior)
library(bayesplot)
library(data.table)
library(loo)
library(dplyr)
library(tidyr)
library(here)
library(cmdstanr)
library(patchwork)
library(mixR)
library(forcats)
library(writexl)
library(posterior)
library(brms)
library(tidybayes)


source(here('R/MultiSeroFunctions.R'))
source(here('R/Functions.R'))

# Recompile model
model_path = here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/StanModel/Full_MultiSero_Model.stan')
mod = cmdstan_model(model_path, pedantic = FALSE, force_recompile = TRUE)

# Import data file 
meta_data <- readRDS(here('Results/meta_data_clean_without_coords.rds'))
nrow(meta_data)

# This is used for the onnv samples modle 
# this is to check that removing the 920 NA chik samples doesnt change results 
meta_data_without_coords_supp_materials <- readRDS('Results/meta_data_without_coords_supp_materials.rds')
nrow(meta_data_without_coords_supp_materials)
# Remove NAs
meta_data_chik_onnv_model <- meta_data %>%
  drop_na(CHIKV_sE2, MAYV_E2, ONNV_VLP) %>%
  mutate(stan_idx_full_model = row_number()) 
nrow(meta_data_chik_onnv_model)
colnames(meta_data_chik_onnv_model)

# Number of samples per year (that the model was run on)
meta_data_chik_onnv_model %>%
  group_by(year_of_survey) %>%
  summarise(n_samples = n()) 
unique(meta_data_chik_onnv_model$year_of_survey)

# onnv samples == including 920 samples that were NA for CHIK and thus removed in the full model
meta_data_onnv_samples_model <- meta_data_without_coords_supp_materials %>%
  drop_na(MAYV_E2, ONNV_VLP) %>%
  mutate(stan_idx_full_model = row_number()) 
nrow(meta_data_onnv_samples_model) # 6155 


# Log CHIK, ONNV and MAY
cols_to_log <- c("CHIKV_sE2", "MAYV_E2", "ONNV_VLP")
new_cols_names <- paste0(cols_to_log, "_log")


meta_data_chik_onnv_model[new_cols_names] <- lapply(meta_data_chik_onnv_model[cols_to_log], log)
meta_data_onnv_samples_model[new_cols_names] <- lapply(meta_data_onnv_samples_model[cols_to_log], log)

# Extract only pathogen cols
full_model_alpha <- meta_data_chik_onnv_model %>%
  dplyr::select(CHIKV_sE2_log, MAYV_E2_log, ONNV_VLP_log)
nrow(full_model_alpha)

onnv_samples_model_alpha <- meta_data_onnv_samples_model %>%
  dplyr::select(MAYV_E2_log, ONNV_VLP_log)
nrow(onnv_samples_model_alpha)


# pathogen names for the model 
pathogens_full_model = c("ONNV_VLP_log","CHIKV_sE2_log","MAYV_E2_log")
pathogens_chik_model = c("CHIKV_sE2_log","ONNV_VLP_log", "MAYV_E2_log")
pathogens_onnv_only_model = c("ONNV_VLP_log","MAYV_E2_log")

# run with all three 
preprocessed_data_full_model <- prepare_multiplex_sero_data(
  data = full_model_alpha,
  pathogens = pathogens_full_model,
  present_pathogens = c("ONNV_VLP_log","CHIKV_sE2_log")
)

preprocessed_data_chik_model <- prepare_multiplex_sero_data(
  data = full_model_alpha,
  pathogens = pathogens_chik_model,
  present_pathogens = c("CHIKV_sE2_log")
)


preprocessed_data_onnv_model <- prepare_multiplex_sero_data(
  data = full_model_alpha,
  pathogens = pathogens_full_model,
  present_pathogens = c("ONNV_VLP_log")
)

preprocessed_data_onnv_samples_model <- prepare_multiplex_sero_data(
  data = onnv_samples_model_alpha,
  pathogens = pathogens_onnv_only_model,
  present_pathogens = c("ONNV_VLP_log")
)

saveRDS(preprocessed_data_full_model, here('Results/preprocessed_data_full_model.rds'))
saveRDS(preprocessed_data_chik_model, here('Results/preprocessed_data_chik_model.rds'))
saveRDS(preprocessed_data_onnv_model, here('Results/preprocessed_data_onnv_model.rds'))
saveRDS(preprocessed_data_onnv_samples_model, here('Results/preprocessed_data_onnv_samples_model.rds'))

#--- Fit full model 
ini_full <- init(preprocessed_data_full_model$data, nChains = 3)
ini_chik <- init(preprocessed_data_chik_model$data, nChains = 3)
ini_onnv <-  init(preprocessed_data_onnv_model$data, nChains = 3)

 
fit_full_model <- mod$sample(
  data = preprocessed_data_full_model$data, 
  chains = 3, 
  iter_sampling = 3000, 
  refresh = 100, 
  iter_warmup = 1000, 
  parallel_chains = 3,
  init = ini_full,
  save_cmdstan_config=TRUE
)


fit_chik_model <- mod_final$sample(
  data = preprocessed_data_chik_model$data, 
  chains = 3, 
  iter_sampling = 3000, 
  refresh = 100, 
  iter_warmup = 1000, 
  parallel_chains = 3,
  init = ini_chik,
  save_cmdstan_config=TRUE
)


fit_onnv_model <- mod_final$sample(
  data = preprocessed_data_onnv_model$data, 
  chains = 3, 
  iter_sampling = 3000, 
  refresh = 100, 
  iter_warmup = 1000, 
  parallel_chains = 3,
  init = ini_onnv,
  save_cmdstan_config=TRUE
)


fit_onnv_samples_model <- mod_final$sample(
  data = preprocessed_data_onnv_samples_model$data, 
  chains = 3, 
  iter_sampling = 3000, 
  refresh = 100, 
  iter_warmup = 1000, 
  parallel_chains = 3,
  init = ini_onnv,
  save_cmdstan_config=TRUE
)



#saveRDS(fit, '/Users/ap2488/Desktop/Cameroon_Analysis_2025/16thDEC_CHIK+ONNV_MultiSeroFit.rds')
fit_full_model$save_object(here('Results/full_model_fits.rds'))
fit_chik_model$save_object(here('Results/chik_model_fits.rds'))
fit_onnv_model$save_object(here('Results/onnv_model_fits.rds'))

# Read RDS files
fit_full_model <- readRDS(here('Results/full_model_fits.rds'))
fit_chik_model <- readRDS(here('Results/chik_model_fits.rds'))
fit_onnv_model <- readRDS(here('Results/onnv_model_fits.rds'))


# extract chains
chains_full <- fit_full_model$draws(format='df')
chains_df_full <- as.data.frame(chains_full)

chains_chik <- fit_chik_model$draws(format='df')
chains_df_chik <- as.data.frame(chains_chik)


chains_onnv <- fit_onnv_model$draws(format='df')
chains_df_onnv <- as.data.frame(chains_onnv)


# Plot trace plots with all chains clearly visible
color_scheme_set("mix-blue-red")
p1 <- mcmc_trace(fit_full_model$draws(c("seroAll", "lp__")))
p2 <- mcmc_trace(fit_full_model$draws(c("mu0", "mu1")))
p3 <- mcmc_trace(fit_full_model$draws(c('sd0','sd1')))
p4 <- mcmc_trace(fit_full_model$draws(c('phi','rho00')))
print(p1 + p2 + p3 + p4)



# Save posterior estimates from full model
# extract posterior estimates
phi <- extract_phi(chains_df_full, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens)
mu <- extract_mu(chains_df_full, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens)
sero <- extract_sero(chains_df_full, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens)
sds <- extract_sd(chains_df_full, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens)

print(sero)

# cross reactive titre incease 
p_CR <- plot_titer_increases_comparison(phi$phi, mu$mus1)
p_CR$plot_data
print(p_CR)

# cross reactive titre increase MAYYV
p_CR_mayv <- titer_increases_comparison_mayv(phi$phi, mu$mus1)
print(p_CR_mayv)

# --- Cluster assignment based on max probability - For INLA analysis 
N  <- preprocessed_data_full_model$data$N
nC <- preprocessed_data_full_model$data$nC
draws_post <- as_draws_df(fit_full_model$draws("post_prob"))
prob_matrix <- matrix(NA_real_, nrow = N, ncol = nC)
for (n in 1:N) {
  for (c in 1:nC) {
    prob_matrix[n, c] <- mean(draws_post[[sprintf("post_prob[%d,%d]", n, c)]])
  }
}

cluster_label <- max.col(prob_matrix, ties.method = "first") #returns leftmost if two cols have equal (and max) prob
cluster_prob <- apply(prob_matrix, 1, max)
table(cluster_label)

length(cluster_label)

cluster_df <- meta_data_chik_onnv_model |>
  dplyr::select(id, stan_idx_full_model) |>   # <-- add stan_idx here
  dplyr::mutate(
    cluster = cluster_label,
    cluster_prob = cluster_prob
  )
meta_data <- meta_data |>
  dplyr::left_join(cluster_df, by = "id")


# label 2 == ONNV pos, label 3 == CHIK pos  else 0 (Neg for both)
meta_data$ONNV_pos <- as.integer(meta_data$cluster == 2)
meta_data$CHIK_pos <- as.integer(meta_data$cluster == 3)


# --- save outputs
#save estimates
write_xlsx(
  list(
    phi  = phi$phi,
    cross_reactive_titre_increase = p_CR$plot_data,
    mu_neg   = mu$mus0,
    mu_pos = mu$mus1,
    sero = sero,
    sds  = sds
  ),
  path = here("Results/onnv_chik_model_posterior_estimates.xlsx")
)

# save meta data with cluster labels and probabilities
write.csv(meta_data, here("Results/meta_data_with_labels.csv"), row.names = FALSE)



# --- Comparison to other models
onnv_samples_fit <- readRDS(here('Results/fit_onnv_samples_model.rds'))

# --- 1) Compare estimtes of Full ONNV only model (with all samples) with ONNV+CHIK model (with 920 samples removed that were NA for CHIK)
chains_onnv_samples <- onnv_samples_fit$draws(format='df')
chains_df_onnv_samples <- as.data.frame(chains_onnv_samples)
preprocessed_data_onnv_samples_model <- readRDS(here("Results/preprocessed_data_onnv_samples_model.rds"))

# Model 2 
sero_onnv_samples <- extract_sero(chains_df_onnv_samples, preprocessed_data_onnv_samples_model$data, 
pathogens=preprocessed_data_onnv_samples_model$pathogens)
print(sero_onnv_samples)
mu_onnv_samples <- extract_mu(chains_df_onnv_samples, preprocessed_data_onnv_samples_model$data, pathogens=preprocessed_data_onnv_samples_model$pathogens)
phi_onnv_samples <- extract_phi(chains_df_onnv_samples, preprocessed_data_onnv_samples_model$data, pathogens=preprocessed_data_onnv_samples_model$pathogens)
sds_onnv_samples <- extract_sd(chains_df_onnv_samples, preprocessed_data_onnv_samples_model$data, pathogens=preprocessed_data_onnv_samples_model$pathogens)

# save estimates
write_xlsx(
  list(
    phi  = phi_onnv_samples$phi,
    mu_neg   = mu_onnv_samples$mus0,
    mu_pos = mu_onnv_samples$mus1,
    sero = sero_onnv_samples,
    sds  = sds_onnv_samples
  ),
  path = here("Results/onnv_samples_model_posterior_estimates.xlsx")
)

# cluster labels for ONNV only model (with all samples)
N_onnv_samples <- preprocessed_data_onnv_samples_model$data$N
nC_onnv_samples <- preprocessed_data_onnv_samples_model$data$nC
draws_post_onnv_samples <- as_draws_df(onnv_samples_fit$draws("post_prob"))
prob_matrix_onnv_samples <- matrix(NA_real_, nrow = N_onnv_samples, ncol = nC_onnv_samples)
for (n in 1:N_onnv_samples) {
  for (c in 1:nC_onnv_samples) {
    prob_matrix_onnv_samples[n, c] <- mean(draws_post_onnv_samples[[sprintf("post_prob[%d,%d]", n, c)]])
  }
}     

cluster_label_onnv_samples <- max.col(prob_matrix_onnv_samples, ties.method = "first") #returns leftmost if two cols have equal (and max) prob
cluster_prob_onnv_samples <- apply(prob_matrix_onnv_samples, 1, max)    
table(cluster_label_onnv_samples)

cluster_df_onnv_samples <- meta_data_onnv_samples_model |>
  dplyr::select(id, stan_idx_full_model) |>   # <-- add stan_idx here
  dplyr::mutate(
    cluster = cluster_label_onnv_samples,
    cluster_prob = cluster_prob_onnv_samples
  )

meta_data_onnv_samples <- meta_data_without_coords_supp_materials
meta_data_onnv_samples <- meta_data_onnv_samples |>
  dplyr::left_join(cluster_df_onnv_samples, by = "id")


# label 2 == ONNV pos
meta_data_onnv_samples$ONNV_pos <- as.integer(meta_data_onnv_samples$cluster == 2)
table(meta_data_onnv_samples$ONNV_pos)
nrow(meta_data_onnv_samples)

#save meta data with cluster labels and probabilities
write.csv(meta_data_onnv_samples, here("Results/meta_data_onnv_samples_with_labels.csv"), row.names = FALSE)




# --- 2) Compare ONNV+CHIK, ONNV only and CHIK only models
lli <- data.frame(model = c('Full_Model', 'CHIK_only_model', 'ONNV_only_model'),
                  par = 'LogLik',
                  med = NA, ciL = NA, ciU = NA)

lli[1, 3:5] <- quantile(chains_df_full$sumloglik,c(0.5, 0.025, 0.975))
lli[2, 3:5] <- quantile(chains_df_chik$sumloglik, c(0.5, 0.025, 0.975))
lli[3, 3:5] <- quantile(chains_df_onnv$sumloglik, c(0.5, 0.025, 0.975))


# save log-likelihood summary
write.csv(lli, here("Results/loglik_model_comparison.csv"), row.names = FALSE)



#--- 3)Compare multisero model estimates to 2D Mixture model fits 
# ---- Fit CHIK  - using 50% for comp2
fmm_normal_chik <- mixfit(full_model_alpha$CHIKV_sE2_log, ncomp = 2, family="normal")
pred.dat_normal_chik <- cbind(fmm_normal_chik$data, fmm_normal_chik$comp.prob)
chik_positive_normal <- ifelse(pred.dat_normal_chik[, 3] > 0.5, "1", "0")
table(chik_positive_normal)


# ---- Fit ONNV  - using 50% for comp2
fmm_normal_onnv <- mixfit(full_model_alpha$ONNV_VLP_log, ncomp = 2, family="normal")
pred.dat_normal_onnv <- cbind(fmm_normal_onnv$data, fmm_normal_onnv$comp.prob)
onnv_positive_normal <- ifelse(pred.dat_normal_onnv[, 3] > 0.5, "1", "0")
table(onnv_positive_normal)


# ---- Fit MAYV  - using 50% for comp2
fmm_normal_mayv <- mixfit(full_model_alpha$MAYV_E2_log, ncomp = 2, family="normal")
pred.dat_normal_mayv<-cbind(fmm_normal_mayv$data, fmm_normal_mayv$comp.prob)
mayv_positive_normal <- ifelse(pred.dat_normal_mayv[, 3] > 0.5, "1", "0")
table(mayv_positive_normal)


#all mixture model plots 
chik <- plot(fmm_normal_chik) +  theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())  
        
onnv <- plot(fmm_normal_onnv) +  theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) 

mayv <- plot(fmm_normal_mayv) +  theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())

mixture_models_fits <- chik + onnv + mayv 
print(mixture_models_fits)


# 2D mixture model fits with 50% threshold for comp2
mixture_model_data <- meta_data_full_model
mixture_model_data$CHIK_pos_mixture_model <- as.numeric(ifelse(fmm_normal_chik$comp.prob[, 2] > 0.5, "1", "0"))
mixture_model_data$ONNV_pos_mixture_model <- as.numeric(ifelse(fmm_normal_onnv$comp.prob[, 2] > 0.5, "1", "0"))
mixture_model_data$MAYV_pos_mixture_model <- as.numeric(ifelse(fmm_normal_mayv$comp.prob[, 2] > 0.5, "1", "0"))



# Proportion mosquitoes with propotion positive for each virus based on 2D mixture model fits
aegmax <- seq(0,1,0.1)
aegmin <- aegmax - 0.5
aegmin[which(aegmin<0)] <- 0

anoph_max <- seq(0, 1, 0.1)
anoph_min <- anoph_max - 0.5
anoph_min[which(anoph_min < 0)] <- 0

virus_colors <- c(
  CHIK = "#2e86ab",
  ONNV = "#b31459",
  MAYV = "#55038c"
)

# ── Run calculate_prop_by_variable for all vector × virus combinations ────────
# Aegypti
df_aeg_chik <- calculate_prop_by_variable(mixture_model_data, "aeg_pw_district", "CHIK_pos_mixture_model", aegmax, aegmin)
df_aeg_onnv <- calculate_prop_by_variable(mixture_model_data, "aeg_pw_district", "ONNV_pos_mixture_model", aegmax, aegmin)
df_aeg_mayv <- calculate_prop_by_variable(mixture_model_data, "aeg_pw_district", "MAYV_pos_mixture_model", aegmax, aegmin)

# Albopictus
df_alb_chik <- calculate_prop_by_variable(mixture_model_data, "alb_pw_district", "CHIK_pos_mixture_model", aegmax, aegmin)
df_alb_onnv <- calculate_prop_by_variable(mixture_model_data, "alb_pw_district", "ONNV_pos_mixture_model", aegmax, aegmin)
df_alb_mayv <- calculate_prop_by_variable(mixture_model_data, "alb_pw_district", "MAYV_pos_mixture_model", aegmax, aegmin)

# Funestus — update breaks if different from aeg
df_fun_chik <- calculate_prop_by_variable(mixture_model_data, "fun_pw_district", "CHIK_pos_mixture_model", anoph_max, anoph_min)
df_fun_onnv <- calculate_prop_by_variable(mixture_model_data, "fun_pw_district", "ONNV_pos_mixture_model", anoph_max, anoph_min)
df_fun_mayv <- calculate_prop_by_variable(mixture_model_data, "fun_pw_district", "MAYV_pos_mixture_model", anoph_max, anoph_min)

# Gambiae — update breaks if different from aeg
df_gam_chik <- calculate_prop_by_variable(mixture_model_data, "gam_pw_district", "CHIK_pos_mixture_model", anoph_max, anoph_min)
df_gam_onnv <- calculate_prop_by_variable(mixture_model_data, "gam_pw_district", "ONNV_pos_mixture_model", anoph_max, anoph_min)
df_gam_mayv <- calculate_prop_by_variable(mixture_model_data, "gam_pw_district", "MAYV_pos_mixture_model", anoph_max, anoph_min)


# ── Plot function: 3 viruses overlaid on one set of axes ─────────────────────
plot_vector_multi_virus <- function(obs_chik, obs_onnv, obs_mayv,
                                    xlab,
                                    colors = virus_colors) {
  
  # Combine the 3 obs dataframes, tagging each with its virus label
  combined <- dplyr::bind_rows(
    obs_chik %>% dplyr::filter(!is.nan(x)) %>% dplyr::mutate(virus = "CHIK"),
    obs_onnv %>% dplyr::filter(!is.nan(x)) %>% dplyr::mutate(virus = "ONNV"),
    obs_mayv %>% dplyr::filter(!is.nan(x)) %>% dplyr::mutate(virus = "MAYV")
  )
  combined$virus <- factor(combined$virus, levels = names(colors))

  ggplot(combined, aes(x = x, y = y, color = virus)) +
    geom_point(size = 4, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0, alpha = 0.6, linewidth = 0.6
    ) +
    scale_color_manual(
      values = colors,
      name   = NULL
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = c("0", "0.25", "0.5", "0.75", "1"),
      expand = c(0, 0)
    ) +
    coord_cartesian(ylim = c(0,0.8), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, 0.8, 0.2)) +
    labs(x = xlab, y = "Proportion positive") +
    theme(
      plot.margin    = margin(t = 10, r = 14, b = 12, l = 14),
      legend.position = "bottom", 
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )
}


# ── Build the 4 plots ─────────────────────────────────────────────────────────
prop_aeg_prev <- plot_vector_multi_virus(
  obs_chik = df_aeg_chik$obs,
  obs_onnv = df_aeg_onnv$obs,
  obs_mayv = df_aeg_mayv$obs,
  xlab     = "Proportion Aedes aegypti"
)

prop_alb_prev <- plot_vector_multi_virus(
  obs_chik = df_alb_chik$obs,
  obs_onnv = df_alb_onnv$obs,
  obs_mayv = df_alb_mayv$obs,
  xlab     = "Proportion Aedes albopictus"
)

prop_fun_prev <- plot_vector_multi_virus(
  obs_chik = df_fun_chik$obs,
  obs_onnv = df_fun_onnv$obs,
   obs_mayv = df_fun_mayv$obs,
  xlab     = "Proportion Anopheles funestus"
)

prop_gam_prev <- plot_vector_multi_virus(
  obs_chik = df_gam_chik$obs,
  obs_onnv = df_gam_onnv$obs,
  obs_mayv = df_gam_mayv$obs,
  xlab     = "Proportion Anopheles gambiae"
)


mosquito_pos_plots <- (prop_aeg_prev | prop_alb_prev) / (prop_fun_prev | prop_gam_prev) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


# save plots 
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/xStarPres/mosquito_pos_plots.png", 
       plot = mosquito_pos_plots,
       width = 10, 
       height = 8, 
       units = "in", 
       dpi = 300,
       bg = "white")



plot_vector_with_2d_model <- function(df_obs, raw_data, xlab, color, pos_col = "ONNV_pos") {

  ylab <- paste0("Proportion ", gsub("_pos", "", pos_col), "positive")

  obs_clean <- df_obs[!is.nan(df_obs$x), ]


  x_scale <- scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    expand = c(0, 0)
  )


  plot_scatter <- ggplot(obs_clean, aes(x = x, y = y)) +
    geom_point(color = color, size = 4, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0, color = color, alpha = 0.6, linewidth = 0.6
    ) +
    x_scale +
    coord_cartesian(ylim = c(0, 0.2), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, 0.2, 0.1)) +
    labs(x = xlab, y = ylab) +
    base_theme +
    theme(
      plot.margin = margin(t = 10, r = 14, b = 12, l = 14)
    )
    return(plot_scatter)
}

colnames(mixture_model_data)

prop_aeg_prev_chik <- plot_chik_2d_model(
  df_aegypti_chik_mixture_model$obs,
  mixture_model_chik_data$aeg_pw_district,
  "Proportion Aedes aegypti",
  color ="#c1518b", pos_col =  "CHIK_pos_mixture_model"
)

prop_albo_prev_chik <- plot_chik_2d_model(
  df_albopictus_chik_mixture_model$obs,
  mixture_model_chik_data$alb_pw_district,
  "Proportion Aedes albopictus",
  color = "#430726", pos_col =  "CHIK_pos_mixture_model"
)
quartz()
prop_albo_prev_chik

aedes_chik <- 
  patchwork::wrap_plots(
    prop_aeg_prev_chik, prop_albo_prev_chik,
    ncol = 2
  ) +
  patchwork::plot_layout(axes = "collect_x")


ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/xStarPres/Fig1.png", 
       plot = aedes_chik,
       width = 10, 
       height = 4, 
       units = "in", 
       dpi = 300,
       bg = "white")




# Aegypti
df_aegypti_chik_mixture_model <- calculate_prop_by_variable(
  data = mixture_model_chik_data, 
  var_col = "aeg_pw_district", 
  positive_col = "CHIK_pos_mixture_model",
  breaks_max = aegmax, 
  breaks_min = aegmin)

  # Albopictus
df_albopictus_chik_mixture_model <- calculate_prop_by_variable(
  data = mixture_model_chik_data, 
  var_col = "alb_pw_district", 
  positive_col = "CHIK_pos_mixture_model",
  breaks_max = aegmax, 
  breaks_min = aegmin)



