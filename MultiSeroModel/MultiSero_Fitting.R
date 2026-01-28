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
meta_data <- meta_data %>%
  drop_na(CHIKV_sE2, MAYV_E2, ONNV_VLP)


# Log CHIK, ONNV and MAY
cols_to_log <- c("CHIKV_sE2", "MAYV_E2", "ONNV_VLP")
new_cols_names <- paste0(cols_to_log, "_log")
meta_data[new_cols_names] <- lapply(meta_data[cols_to_log], log)

# Extract only Alpha virus cols
meta_data_alpha <- meta_data %>%
  dplyr::select(CHIKV_sE2_log, MAYV_E2_log, ONNV_VLP_log)
nrow(meta_data_alpha)


# pathogen names for the model (circulating pathogens first)
pathogens = c("ONNV_VLP_log","CHIKV_sE2_log","MAYV_E2_log")

# prepare data for stan
preprocessed_data <- prepare_multiplex_sero_data(
  data = meta_data_alpha,
  pathogens = pathogens,
  present_pathogens = c("ONNV_VLP_log","CHIKV_sE2_log")
)

#--- chain starting values
ini <- init_diffSds(preprocessed_data$data, nChains = 3)

fit <- mod$sample(
data = preprocessed_data$data, 
chains = 3, 
iter_sampling = 3000, 
refresh = 100, 
iter_warmup = 1000, 
parallel_chains = 3,
init = ini,
save_cmdstan_config=TRUE
)

#save fits
fit$save_object('/Users/ap2488/Desktop/Cameroon_Analysis_2025/final_model_fits.rds')
fit <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/redone_final_model_fits.rds')

# extract chains
chains <- fit$draws(format='df')
chains_df <- as.data.frame(chains)

# Plot trace plots with all chains clearly visible
color_scheme_set("mix-blue-red")
p1 <- mcmc_trace(fit_final_model$draws(c("seroAll", "lp__")))
p2 <- mcmc_trace(fit_final_model$draws(c("mu0", "mu1")))
p3 <- mcmc_trace(fit_final_model$draws(c('sd0','sd1')))
p4 <- mcmc_trace(fit_final_model$draws(c('phi','rho00')))

print(p1 + p2 + p3 + p4)

# Plot fits (neg component, neg-CR component, pos component)
distfits <- plot_fits(chains_df, preprocessed_data$data, pathogens=preprocessed_data$pathogens, show_crossreactive_for = seq_along(preprocessed_data$pathogens))
distfits$fitPN


# extract phi and mu1 posterior distributions
phi <- extract_phi(chains_df, preprocessed_data$data, pathogens=preprocessed_data$pathogens)
mu <- extract_mu(chains_df, preprocessed_data$data, pathogens=preprocessed_data$pathogens)


# plot titre increease due to infection / CR for each pathogen
p_CR <- plot_titer_increases_comparison(phi$phi, mu$mus1)
print(p_CR)


# plot proportion pos 
p_sero <- plot_seroprevalence(chains_df)
print(p_sero)



# plot prevelance by age group
plot_age_seroprevalence(meta_data, chains_df, component_col = 2, pathogen_name = "ONNV")





# Component labels for INLA 
N  <- preprocessed_data$data$N
nC <- preprocessed_data$data$nC
draws_post <- as_draws_df(fit_final_model$draws("post_prob"))
prob_matrix <- matrix(NA_real_, nrow = N, ncol = nC)
for (n in 1:N) {
  for (c in 1:nC) {
    prob_matrix[n, c] <- mean(draws_post[[sprintf("post_prob[%d,%d]", n, c)]])
  }
}
cluster_assignment <- apply(prob_matrix, 1, which.max)
table(cluster_assignment)

# add hard assignmnet back to data 
# 1 = Negative , 2 = ONNV pos and 3 = CHIK Pos
meta_data_alpha$label <- cluster_assignment
meta_data$label <- cluster_assignment

# label 2 == ONNV pos 
# label 3 == CHIK pos 
# else 0 (Neg for both)
meta_data$ONNV_pos <- as.integer(meta_data$label == 2)
meta_data$CHIK_pos <- as.integer(meta_data$label == 3)

# plot to visualise distribution, coloured by label
ggplot(meta_data_alpha,
       aes(x = ONNV_VLP_log,
           y = CHIKV_sE2_log,
           color = factor(label))) +
  geom_point(alpha = 0.5, size = 2) +
  labs(
    color = "Cluster",
    title = "Hard-assigned clusters"
  ) +
  theme_minimal()

# save file with labels 
write.csv(meta_data, "/Users/ap2488/Desktop/Cameroon_Analysis_2025/final_meta_data_with_labels.csv", row.names = FALSE)