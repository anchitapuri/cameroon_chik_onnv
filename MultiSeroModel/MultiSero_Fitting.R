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
source(here('/Users/ap2488/Documents/GitHub/cameroon_chik_onnv/Spatial Analysis/Functions.R'))


# Setup cmdstan
check_cmdstan_toolchain()
cmdstan_path <- "/Users/ap2488/.cmdstan/cmdstan-2.36.0"
set_cmdstan_path(cmdstan_path)

# Compile model
model_path = "/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/MultiSeroModel/Final_MultiSero.stan"
mod = cmdstan_model(model_path, pedantic=FALSE)


# Import data file 
meta_data <- read.csv('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/meta_data_without_coords.csv')
nrow(meta_data)
unique(meta_data$year_of_survey)

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
ini <- init(preprocessed_data_full_model$data, nChains = 3)


fit_full_model <- mod$sample(
data = preprocessed_data_full_model$data, 
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
fit_full_model$save_object('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/3rdMarch_full_model_fits.rds')
fit_onnv_only_model$save_object('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/adapted_onnv_only_model_fits.rds')


# Read fit RDS
fit_full_model <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/3rdMarch_full_model_fits.rds')
fit_onnv_only_model <-readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/adapted_onnv_only_model_fits.rds')
fit_naive_model <- readRDS('/Users/ap2488/Desktop/Cameroon_Analysis_2025/naive_model_fits.rds')



# extract chains
chains <- fit_full_model$draws(format='df')
chains_df <- as.data.frame(chains)


# extract chains - onnv only model
chain_onnv_only <- fit_onnv_only_model$draws(format='df')
chains_df_onnv_only <- as.data.frame(chain_onnv_only)

# extract chains - onnv only model
chains_naive_model <- fit_naive_model$draws(format='df')
chains_df_naive_model <- as.data.frame(chains_naive_model)


# Plot trace plots with all chains clearly visible
color_scheme_set("mix-blue-red")
p1 <- mcmc_trace(fit_full_model$draws(c("seroAll", "lp__")))
p2 <- mcmc_trace(fit_full_model$draws(c("mu0", "mu1")))
p3 <- mcmc_trace(fit_full_model$draws(c('sd0','sd1')))
p4 <- mcmc_trace(fit_full_model$draws(c('phi','rho00')))
quartz()
print(p1 + p2 + p3 + p4)


seq_along(preprocessed_data_full_model$pathogens)

# Plot fits (neg component, neg-CR component, pos component)
distfits <- plot_fits(chains_df, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens, show_crossreactive_for = seq_along(preprocessed_data_full_model$pathogens))
quartz()
distfits$fitPN 


ggsave(
  filename = '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/NEW_Fig2b.png',
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
sds <- extract_sd(chains_df, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens)
print(sds)

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

# After the plot_titres_coloured_by_clusters call, add this to identify and print probs for top 3 high CHIK points
high_chik_indices <- order(meta_data$CHIKV_sE2, decreasing = TRUE)[1:3]
high_chik_data <- meta_data[high_chik_indices, c("id", "CHIKV_sE2", "cluster", "cluster_prob")]
print(high_chik_data)

# To print the full probability vectors for these points
for (i in high_chik_indices) {
  cat("Sample ID:", meta_data$id[i], "\n")
  print(prob_matrix[i, ])
  cat("\n")
}


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





# CHIK vs ONNV correlation plot 
quartz()
plot(full_model_alpha$ONNV_VLP_log, full_model_alpha$CHIKV_sE2_log,col = "#003c8b", pch = 16, xlab = "Log (ONNV VLP MFI)", ylab = "Log (CHIK sE2 MFI)")


# 2) Compare multisero model estimates to 2D Mixture model fits 
# ---- Fit CHIK  - using 50% for comp2
fmm_normal_chik <- mixfit(full_model_alpha$CHIKV_sE2_log, ncomp = 2, family="normal")
plot(fmm_normal_chik)
pred.dat_normal_chik <-cbind(fmm_normal_chik$data, fmm_normal_chik$comp.prob)
chik_positive_normal <- ifelse(pred.dat_normal_chik[, 3] > 0.5, "1", "0")
table(chik_positive_normal)


# ---- Fit ONNV  - using 50% for comp2
fmm_normal_onnv <- mixfit(full_model_alpha$ONNV_VLP_log, ncomp = 2, family="normal")
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
# save plots 
ggsave("/Users/ap2488/Desktop/Cameroon_Analysis_2025/xStarPres/mixture_models_fits.png", 
       plot = mixture_models_fits,
       width = 15, 
       height = 5, 
       units = "in", 
       dpi = 300,
       bg = "white")




# 2D mixture model fits with 50% threshold for comp2
mixture_model_data <- meta_data_full_model
mixture_model_data$CHIK_pos_mixture_model <- as.numeric(ifelse(fmm_normal_chik$comp.prob[, 2] > 0.5, "1", "0"))
mixture_model_data$ONNV_pos_mixture_model <- as.numeric(ifelse(fmm_normal_onnv$comp.prob[, 2] > 0.5, "1", "0"))
mixture_model_data$MAYV_pos_mixture_model <- as.numeric(ifelse(fmm_normal_mayv$comp.prob[, 2] > 0.5, "1", "0"))


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










# compare naive model fits to multisero model fits
lli <- data.frame(model = c('Full Model', 'Naive Model'),
                  par = 'LogLik',
                  med = NA, ciL = NA, ciU = NA)

lli[1, 3:5] <- quantile(chains_df$sumloglik,c(0.5, 0.025, 0.975))
lli[2, 3:5] <- quantile(chains_df_naive_model$sumloglik, c(0.5, 0.025, 0.975))
lli



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






p <- which(preprocessed_data_full_model$pathogens == "MAYV_E2_log")
neg_idx <- which(preprocessed_data_full_model$data$infM[, p] == 0)

cross_reactive_idx <- neg_idx[
  rowSums(preprocessed_data_full_model$data$infM[neg_idx, , drop = FALSE]) > 0
]

true_neg_idx <- neg_idx[
  rowSums(preprocessed_data_full_model$data$infM[neg_idx, , drop = FALSE]) == 0
]

cat("Cross-reactive rows:\n")
print(cross_reactive_idx)

cat("Corresponding infM rows:\n")
print(preprocessed_data_full_model$data$infM[cross_reactive_idx, , drop = FALSE])

i <- 1  # pick any iteration

# simulate components exactly like inside your function
covM <- extract_covM(chains_df, preprocessed_data_full_model$data)

yy <- list()
for (c in 1:preprocessed_data_full_model$data$nC) {
  g <- paste('mu[', c, sep='')
  muu <- sapply(1:preprocessed_data_full_model$data$nP, function(pp)
    chains_df[i, paste0(g, ',', pp, ']')]
  )

  yy[[c]] <- as.data.frame(
    rmvnorm(500, mean = muu, sigma = covM[[i]][[c]])
  )
  yy[[c]]$C <- c
}

yc <- do.call(rbind, yy)

d_cr  <- density(yc[yc$C %in% cross_reactive_idx, p])
d_tn  <- density(yc[yc$C %in% true_neg_idx, p])

plot(d_tn, col="blue", main="MAYV component comparison")
lines(d_cr, col="red")
legend("topright", legend=c("True negative","Cross-reactive"),
       col=c("blue","red"), lwd=2)


infM <- preprocessed_data_full_model$data$infM
colnames(infM) <- c("ONNV","CHIKV","MAYV")
print(infM)

p <- which(colnames(infM) == "MAYV")
neg_idx <- which(infM[, p] == 0)
cross_reactive_idx <- neg_idx[rowSums(infM[neg_idx, , drop = FALSE]) > 0]
true_neg_idx <- neg_idx[rowSums(infM[neg_idx, , drop = FALSE]) == 0]

cat("neg_idx:\n"); print(neg_idx)
cat("cross_reactive_idx:\n"); print(cross_reactive_idx)
cat("true_neg_idx:\n"); print(true_neg_idx)
cat("cross-reactive infM rows:\n"); print(infM[cross_reactive_idx, , drop = FALSE])

