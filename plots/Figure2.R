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

source(here('R/MultiSeroFunctions.R'))
source(here('R/Functions.R'))


#read files 
meta_data_with_labels <- read.csv(here('Results/meta_data_with_labels.csv'))
preprocessed_data_full_model <- readRDS(here('Results/preprocessed_data_full_model.rds'))
fit_full_model <- readRDS(here('Results/full_model_fits.rds'))


# extract chains and parameters
chains_full <- fit_full_model$draws(format='df')
chains_df_full <- as.data.frame(chains_full)
mu <- extract_mu(chains_df_full, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens)
phi <- extract_phi(chains_df_full, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens)


#Figure 2a --- Plot fits (neg component, neg-CR component, pos component)
distfits <- plot_fits(chains_df_full, preprocessed_data_full_model$data, pathogens=preprocessed_data_full_model$pathogens, 
                      show_crossreactive_for = seq_along(preprocessed_data_full_model$pathogens))
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


# Fig 2b --- plot to visualise distribution, coloured by label
titres_plot <- plot_titres_coloured_by_clusters(meta_data_with_labels)
print(titres_plot)

# Fig 2c --  titre increease due to infection / CR for each pathogen
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


# Fig 2d -- plot proportion pos 
p_sero <- plot_seroprevalence(chains_df_full)
print(p_sero)


ggsave(
  filename = '/Users/ap2488/Desktop/Cameroon_Analysis_2025/FinalCode/Fig2d.png',
  plot = p_sero,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300
)



# Figures combined 
fig2 <- distfits$fitPN / (titres_plot | p_CR$p | p_sero)  +
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1)
  )
print(fig2)
ggsave(
  filename =  here('Results/Fig2.png'),
  plot = fig2,
  width = 20,
  height = 12,
  units = "in",
  dpi = 300
)


