library(ggraph)
library(lhs)
library(matrixStats)
library(mvtnorm)
library(matrixcalc)
library(here)
#----- Functions for multivariate Gaussian mixture serology model -----#
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

source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/MultiSeroModel/MultiSeroFunctions.R'))

# Setup cmdstan
check_cmdstan_toolchain()
cmdstan_path <- "/Users/ap2488/.cmdstan/cmdstan-2.36.0"
set_cmdstan_path(cmdstan_path)

# Compile model
model_path = "/Users/ap2488/Desktop/Cameroon_Analysis_2025/Final_MultiSero.stan"
mod = cmdstan_model(model_path_final, pedantic=FALSE)

# Import data file 
meta_data <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/meta_data_without_coords.csv')
nrow(meta_data)


# Remove NAs
meta_data_full_model <- meta_data %>%
  drop_na(CHIKV_sE2, MAYV_E2, ONNV_VLP) %>%
  mutate(stan_idx_full_model = row_number()) 
nrow(meta_data_full_model)
colnames(meta_data_full_model)

# Number of samples per year (that the model was run on)
meta_data_full_model %>%
  group_by(year_of_survey) %>%
  summarise(n_samples = n()) 
unique(meta_data_full_model$year_of_survey)

meta_data_onnv_only_model <- meta_data %>%
  drop_na(MAYV_E2, ONNV_VLP) %>%
  mutate(stan_idx_full_model = row_number()) 
nrow(meta_data_onnv_only_model)


# Log CHIK, ONNV and MAY
cols_to_log <- c("CHIKV_sE2", "MAYV_E2", "ONNV_VLP")
new_cols_names <- paste0(cols_to_log, "_log")


meta_data_full_model[new_cols_names] <- lapply(meta_data_full_model[cols_to_log], log)
meta_data_onnv_only_model[new_cols_names] <- lapply(meta_data_onnv_only_model[cols_to_log], log)

# Extract only pathogen cols
full_model_alpha <- meta_data_full_model %>%
  dplyr::select(CHIKV_sE2_log, MAYV_E2_log, ONNV_VLP_log)
nrow(full_model_alpha)

onnv_only_model_alpha <- meta_data_onnv_only_model %>%
  dplyr::select(CHIKV_sE2_log, MAYV_E2_log, ONNV_VLP_log)
nrow(onnv_only_model_alpha)


# pathogen names for the model 
pathogens_full_model = c("ONNV_VLP_log","CHIKV_sE2_log","MAYV_E2_log")
pathogens_onnv_only_model = c("ONNV_VLP_log","MAYV_E2_log")

# run with all three 
preprocessed_data_full_model <- prepare_multiplex_sero_data(
  data = full_model_alpha,
  pathogens = pathogens_full_model,
  present_pathogens = c("ONNV_VLP_log","CHIKV_sE2_log")
)

# run with all three 
preprocessed_data_onnv_only_model <- prepare_multiplex_sero_data(
  data = onnv_only_model_alpha,
  pathogens = pathogens_onnv_only_model,
  present_pathogens = c("ONNV_VLP_log")
)

saveRDS(preprocessed_data_full_model, '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/preprocessed_data_full_model.rds')
saveRDS(preprocessed_data_onnv_only_model, '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/preprocessed_data_onnv_only_model.rds')

preprocessed_data_full_model <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/preprocessed_data_full_model.rds')

#--- Fit full model 
ini <- init_diffSds(preprocessed_data_full_model$data, nChains = 3)
fit_full_model <- mod$sample(
data = preprocessed_data$data, 
chains = 3, 
iter_sampling = 3000, 
refresh = 100, 
iter_warmup = 1000, 
parallel_chains = 3,
init = ini,
save_cmdstan_config=TRUE
)


#--- Fit ONNV only model 
ini <- init_diffSds(preprocessed_data_onnv_only_model$data, nChains = 3)
fit_onnv_only_model <- mod$sample(
data = preprocessed_data_onnv_only_model$data, 
chains = 3, 
iter_sampling = 3000, 
refresh = 100, 
iter_warmup = 1000, 
parallel_chains = 3,
init = ini,
save_cmdstan_config=TRUE
)


#save fits
fit_full_model$save_object('/Users/ap2488/Desktop/Cameroon_Analysis_2025/fit_full_model.rds')
fit_onnv_only_model$save_object('/Users/ap2488/Desktop/Cameroon_Analysis_2025/fit_onnv_only_model.rds')


# Read fit RDS
fit_full_model <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/adapted_full_model_fits.rds')
fit_onnv_only_model <-readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/adapted_onnv_only_model_fits.rds')

# extract chains
chains <- fit_full_model$draws(format='df')
chains_df <- as.data.frame(chains)

# Plot trace plots with all chains clearly visible
color_scheme_set("mix-blue-red")
p1 <- mcmc_trace(fit_full_model$draws(c("seroAll", "lp__")))
p2 <- mcmc_trace(fit_full_model$draws(c("mu0", "mu1")))
p3 <- mcmc_trace(fit_full_model$draws(c('sd0','sd1')))
p4 <- mcmc_trace(fit_full_model$draws(c('phi','rho00')))
quartz()
print(p1 + p2 + p3 + p4)



# Plot fits (neg component, neg-CR component, pos component)
distfits <- plot_fits(chains_df, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens, show_crossreactive_for = seq_along(preprocessed_data_full_model$pathogens))
quartz()
distfits$fitPN 

ggsave(
  filename = '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/Fig2b.png',
  plot = distfits$fitPN,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

# extract phi and mu1 posterior distributions
phi <- extract_phi(chains_df, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens)
print(phi)
mu <- extract_mu(chains_df, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens)
print(mu)
sero <- extract_sero(chains_df, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens)
print(sero)


# MAYV homologous vs cross reactive increase 
# --- Extract cross reactivity - MAYV vs ONNV and CHIK 
CR_mayv <- titer_increases_comparison_mayv(phi$phi, mu$mus1)
CR_mayv

# plot titre increease due to infection / CR for each pathogen
p_CR <- plot_titer_increases_comparison(phi$phi, mu$mus1)
print(p_CR)

ggsave(
  filename = '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/Fig2c.png',
  plot = p_CR$p,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300
)


# plot proportion pos 
p_sero <- plot_seroprevalence(chains_df)
print(p_sero)


ggsave(
  filename = '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/Fig2d.png',
  plot = p_sero,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300
)


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


cluster_assignment <- apply(prob_matrix, 1, which.max)
unique(cluster_assignment)
table(cluster_assignment)
cluster_assignment_with_uncertainty <- apply(prob_matrix, 1, max)


# add labels back to meta data 
nrow(prob_matrix)
nrow(meta_data)
nrow(meta_data_full_model)

all(meta_data_full_model$id %in% meta_data$id)

cluster_label <- max.col(prob_matrix, ties.method = "first") #returns leftmost if two cols have equal (and max) prob
cluster_prob <- apply(prob_matrix, 1, max)

cluster_df <- meta_data_full_model |>
  dplyr::select(id, stan_idx_full_model) |>   # <-- add stan_idx here
  dplyr::mutate(
    cluster = cluster_label,
    cluster_prob = cluster_prob
  )

meta_data <- meta_data |>
  dplyr::left_join(cluster_df, by = "id")


# label 2 == ONNV pos 
# label 3 == CHIK pos 
# else 0 (Neg for both)
meta_data$ONNV_pos <- as.integer(meta_data$cluster == 2)
meta_data$CHIK_pos <- as.integer(meta_data$cluster == 3)


# plot to visualise distribution, coloured by label
titres_plot <- plot_titres_coloured_by_clusters(meta_data)
print(titres_plot)


# -- Save figures and data with labels
# Figures combined 
fig2 <- distfits$fitPN / (titres_plot | p_CR$p | p_sero)  +
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1)
  )
print(fig2)
ggsave(
  filename = '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/Fig2.png',
  plot = fig2,
  width = 20,
  height = 12,
  units = "in",
  dpi = 300
)
head(meta_data)
# save file with labels 
write.csv(meta_data, "/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/final_meta_data_with_labels.csv", row.names = FALSE)


# --- Comparison to other model 


# 1) compare estimtes of ONNV only model with full model 
chains_onnv_only <- fit_onnv_only_model$draws(format='df')
chains_df_onnv_only <- as.data.frame(chains_onnv_only)

# Model 2 
sero_onnv_only <- extract_sero(chains_df_onnv_only, preprocessed_data_onnv_only_model$data, 
pathogens=preprocessed_data_onnv_only_model$pathogens)
print(sero_onnv_only)

mu_onnv_only <- extract_mu(chains_df_onnv_only, preprocessed_data_onnv_only_model$data, pathogens=preprocessed_data_onnv_only_model$pathogens)
print(mu_onnv_only)

phi_onnv_only <- extract_phi(chains_df_onnv_only, preprocessed_data_onnv_only_model$data, pathogens=preprocessed_data_onnv_only_model$pathogens)
print(phi_onnv_only)

# Model 1 
print(sero)
print(mu)
print(phi)




# 2) Compare multisero model estimates to 2D Mixture model fits 
# ---- Fit CHIK  - using 50% for comp2
fmm_normal_chik <- mixfit(full_model_alpha$CHIKV_sE2_log, ncomp = 2, family="normal")
plot(fmm_normal_chik)
#define threshold to get chik+ve and chik-ve
pred.dat_normal_chik <-cbind(fmm_normal_chik$data, fmm_normal_chik$comp.prob)
chik_positive_normal <- ifelse(pred.dat_normal_chik[, 3] > 0.5, "1", "0")
table(chik_positive_normal)


# ---- Fit ONNV  - using 50% for comp2
fmm_normal_onnv <- mixfit(full_model_alpha$ONNV_VLP_log, ncomp = 2, family="normal")
plot(fmm_normal_onnv)
#define threshold to get onnv+ve and onnv-ve
pred.dat_normal_onnv <- cbind(fmm_normal_onnv$data, fmm_normal_onnv$comp.prob)
onnv_positive_normal <- ifelse(pred.dat_normal_onnv[, 3] > 0.5, "1", "0")
table(onnv_positive_normal)


# ---- Fit MAYV  - using 50% for comp2
fmm_normal_mayv <- mixfit(full_model_alpha$MAYV_E2_log, ncomp = 2, family="normal")
plot(fmm_normal_mayv)
#define threshold to get chik+ve and chik-ve
pred.dat_normal_mayv<-cbind(fmm_normal_mayv$data, fmm_normal_mayv$comp.prob)
mayv_positive_normal <- ifelse(pred.dat_normal_mayv[, 3] > 0.5, "1", "0")
table(mayv_positive_normal)
