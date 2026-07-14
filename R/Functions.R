
# ---- Function for INLA models ----
run_inla <- function(year_intro, data, cam_pop, positive_col) {
  

  # Calculate years of exposure
  data$age_intro <- data$year_of_survey - year_intro
  data$years_of_exposure <- pmin(data$age_intro, data$AgeInYears)
  
  data <- data[!is.na(data$years_of_exposure) & data$years_of_exposure > 0, ]
  
  data_points <- data %>%
  st_drop_geometry() %>%
    filter(!is.na(Easting) & !is.na(Northing)) %>%
    filter(!is.na(ONNV_pos))
  
  # Build estimation coo
  cooe <- cbind(
    Easting  = data_points$Easting,
    Northing = data_points$Northing
  )
  
  # Create mesh
  mesh <- inla.mesh.2d(
    loc      = cooe,
    max.edge = c(10, 40),
    cutoff   = 10,
    offset   = c(100, 200)
  )
  
  # Prediction Grid
  cam_pop_agg <- terra::aggregate(cam_pop, fact = 10, fun = sum, na.rm = TRUE)
  cam_pop_points <- terra::as.points(cam_pop_agg, values = TRUE, na.rm = TRUE)

  # Convert to sf and transform to UTM 33N
  cam_pop_sf <- st_as_sf(cam_pop_points)
  cam_pop_utm <- st_transform(cam_pop_sf, crs = 32633)
  coords_utm <- st_coordinates(cam_pop_utm)
  
  # For INLA - convert to km (INLA works better with smaller numbers)
  coop <- coords_utm / 1000
  colnames(coop) <- c("X", "Y")
  
  # Build SPDE model
  spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
  s.index <- inla.spde.make.index("spatial.field", spde$n.spde)
  
  # Projection matrices
  A <- inla.spde.make.A(mesh = mesh, loc = cooe)
  Ap <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coop))
  
  # Estimation stack 
  stk.e <- inla.stack(
    data = list(y = data_points[[positive_col]]),  
    A = list(1, A),
    effects = list(data.frame(Intercept = 1, age = data_points$years_of_exposure),
                   spatial.field = s.index),
    tag = "est")
  
  # Prediction stack
  stk.p <- inla.stack(
    tag = "pred", 
    data = list(y = rep(NA, nrow(coop))),
    A = list(1, Ap),
    effects = list(data.frame(Intercept = 1, age = rep(1, nrow(coop))),
                   spatial.field = s.index))
  
  # Full stack
  stk.full <- inla.stack(stk.e, stk.p)
  
  # Run INLA model
  output <- inla(y ~ -1 + Intercept + offset(log(age)) + f(spatial.field, model = spde),
                 data = inla.stack.data(stk.full),
                 family = "binomial",
                 Ntrials = 1,
                 control.family = list(link = "cloglog"),
                 control.predictor = list(A = inla.stack.A(stk.full), 
                                          compute = TRUE, 
                                          link = 1),
                 control.compute = list(dic = TRUE, config = TRUE),
                 verbose = FALSE)
  
  return(list(
    year = year_intro,
    output = output,
    dic = output$dic$dic,
    mesh = mesh,
    stk.full = stk.full,
    data_filtered = data_points,
    cooe = cooe,
    coop = coop
  ))
}

run_inla_multivariable <- function(year_intro, data, cam_pop, positive_col,
                     covars = c("gam_pw_district", "log_pop_density"),
                     covar_grid = NULL) {   

  # Years of exposure
  data$age_intro <- data$year_of_survey - year_intro
  data$years_of_exposure <- pmin(data$age_intro, data$AgeInYears)
  data <- data[!is.na(data$years_of_exposure) & data$years_of_exposure > 0, ]

  data_points <- data %>%
    st_drop_geometry() %>%
    filter(!is.na(Easting) & !is.na(Northing)) %>%
    filter(!is.na(.data[[positive_col]]))          # was hardcoded to ONNV_pos

  # Drop rows missing any covariate, then standardise (store centre/scale for prediction)
  data_points <- data_points[complete.cases(data_points[, covars]), ]
  covar_means <- sapply(covars, function(v) mean(data_points[[v]]))
  covar_sds   <- sapply(covars, function(v) sd(data_points[[v]]))
  covars_z <- paste0(covars, "_z")
  for (i in seq_along(covars)) {
    data_points[[covars_z[i]]] <- (data_points[[covars[i]]] - covar_means[i]) / covar_sds[i]
  }

  # Mesh
  cooe <- cbind(Easting = data_points$Easting, Northing = data_points$Northing)
  mesh <- inla.mesh.2d(loc = cooe, max.edge = c(10, 40), cutoff = 10, offset = c(100, 200))

  # Prediction grid
  cam_pop_agg    <- terra::aggregate(cam_pop, fact = 10, fun = sum, na.rm = TRUE)
  cam_pop_points <- terra::as.points(cam_pop_agg, values = TRUE, na.rm = TRUE)
  cam_pop_sf  <- st_as_sf(cam_pop_points)
  cam_pop_utm <- st_transform(cam_pop_sf, crs = 32633)
  coop <- st_coordinates(cam_pop_utm) / 1000
  colnames(coop) <- c("X", "Y")

  # Covariate values at grid: 0 = mean (default), else standardise supplied raw values
  if (is.null(covar_grid)) {
    Xp <- matrix(0, nrow = nrow(coop), ncol = length(covars))
  } else {
    Xp <- sapply(seq_along(covars), function(i)
      (covar_grid[[covars[i]]] - covar_means[i]) / covar_sds[i])
  }
  colnames(Xp) <- covars_z

  # SPDE
  spde    <- inla.spde2.matern(mesh = mesh, alpha = 2)
  s.index <- inla.spde.make.index("spatial.field", spde$n.spde)
  A  <- inla.spde.make.A(mesh = mesh, loc = cooe)
  Ap <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coop))

  # Estimation stack (+ covariates)
  eff.e <- cbind(data.frame(Intercept = 1, age = data_points$years_of_exposure),
                 data_points[, covars_z, drop = FALSE])
  stk.e <- inla.stack(data = list(y = data_points[[positive_col]]),
                      A = list(1, A),
                      effects = list(eff.e, spatial.field = s.index),
                      tag = "est")

  # Prediction stack (+ covariates)
  eff.p <- cbind(data.frame(Intercept = 1, age = rep(1, nrow(coop))),
                 as.data.frame(Xp))
  stk.p <- inla.stack(tag = "pred",
                      data = list(y = rep(NA, nrow(coop))),
                      A = list(1, Ap),
                      effects = list(eff.p, spatial.field = s.index))

  stk.full <- inla.stack(stk.e, stk.p)

  # Formula
  form <- as.formula(paste("y ~ -1 + Intercept +",
                           paste(covars_z, collapse = " + "),
                           "+ offset(log(age)) + f(spatial.field, model = spde)"))

  output <- inla(form,
                 data = inla.stack.data(stk.full),
                 family = "binomial", Ntrials = 1,
                 control.family    = list(link = "cloglog"),
                 control.predictor = list(A = inla.stack.A(stk.full), compute = TRUE, link = 1),
                 control.compute   = list(dic = TRUE, config = TRUE),
                 verbose = FALSE)

  list(year = year_intro, output = output, dic = output$dic$dic,
       mesh = mesh, stk.full = stk.full, data_filtered = data_points,
       cooe = cooe, coop = coop,
       covars = covars, covar_means = covar_means, covar_sds = covar_sds)
}


# --- Function to extract and plot FOI ---
predicted_foi <- function(model, coop, pathogen_name = "ONNV") { 
  # Get prediction indices
  index_pred <- inla.stack.index(model$stk.full, tag = "pred")$data
  
  # Extract the intercept
  eta_pred <- model$output$summary.linear.predictor[index_pred, "mean"]
  lambda_pred <- exp(eta_pred)
  
  # Create dataframe with UTM coordinates (km -> m)
  foi_df <- data.frame(
    X = coop[, "X"] * 1000,
    Y = coop[, "Y"] * 1000,
    foi = lambda_pred
  )
  
  # Convert to sf object for plotting
  foi_sf <- st_as_sf(
    foi_df,
    coords = c("X", "Y"),
    crs = 32633  # UTM Zone 33N
  )
  
  # Add FOI values to the sf object
  foi_sf$foi <- foi_df$foi
  
  # Plot
  p <- ggplot() +
    geom_sf(
      data = foi_sf, aes(color = foi), size = 0.5, alpha = 0.5) +
      scale_color_gradientn(
        colours = c("#1f363d",
         "#40798c",
        "#70a9a1",
        "#f46036",
        "#a4243b"
        ),
      values = scales::rescale(c(0, 0.01, 0.02, 0.03, 0.04, 0.05)),
      name = "FOI (λ)",
      limits = c(0, max(foi_sf$foi)),
      guide = guide_colorbar(
        direction = "horizontal",
        barheight = unit(0.25, "cm"),
        barwidth  = unit(6.5, "cm"),
        title.position = "top",
        label.position = "bottom",
        ticks = TRUE,
        ticks.length = unit(0.1, "cm")
      )
    ) +
    theme_minimal() + 
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 20),
      axis.text = element_blank(),      # Remove axis text (lat/long labels)
      axis.ticks = element_blank(),
      legend.text = element_text(size = 14),   # tick label size
      legend.title = element_text(size = 14)) +  # "Seroprevalence" title size) + 
    annotation_scale(
    bar_cols = c("black", "white"),  # alternating black/white like the reference
    height = unit(0.2, "cm"),
    text_family = "sans", 
    text_cex = 1.5
  )
  
  # Return foi_sf and plot invisibly
  invisible(list(
    foi_sf = foi_sf,
    foi_df = foi_df,
    plot = p
  ))
}


# ---  Function to extract and plot seroprevalence ---
predicted_seroprevalence <- function(foi_result, model, age_groups, age_weights,
                                          crs = 32633, pathogen_name = "ONNV") {
  
  # For each location and age group, calculate average prevalence within that age group
  # This requires integrating P(a) = 1 - exp(-λa) over each age interval
  
  n_locs <- nrow(foi_result$foi_sf)
  n_age_groups <- nrow(age_groups)
  
  # Matrix to store average prevalence for each location × age group
  avg_prev_mat <- matrix(0, nrow = n_locs, ncol = n_age_groups)
  
  for (i in 1:n_locs) {
    lambda <- foi_result$foi_sf$foi[i]
    
    for (j in 1:n_age_groups) {
      a_lower <- age_groups$age_lower[j]
      a_upper <- age_groups$age_upper[j]
      age_width <- a_upper - a_lower
      
      # Analytical solution for average prevalence in age interval [a_lower, a_upper]
      # Average of (1 - exp(-λa)) from a_lower to a_upper
      # = 1 - (1/λΔa) * [exp(-λa_lower) - exp(-λa_upper)]
      
      if (lambda > 1e-10) {  # Avoid division by zero
        avg_prev <- 1 - (1/(lambda * age_width)) * 
                        (exp(-lambda * a_lower) - exp(-lambda * a_upper))
      } else {
        # For very small lambda, use approximation
        avg_prev <- lambda * (a_lower + a_upper) / 2
      }
      
      avg_prev_mat[i, j] <- avg_prev
    }
  }
  
  # Now weight by age distribution
  prev_loc <- as.vector(avg_prev_mat %*% age_weights)
  
  # Create dataframe with prediction coordinates
  prev_df <- data.frame(
    X_km = model$coop[, "X"],
    Y_km = model$coop[, "Y"],
    prev = prev_loc
  )
  
  # Convert to sf object
  prev_sf <- st_as_sf(
    data.frame(X = prev_df$X_km * 1000, Y = prev_df$Y_km * 1000),
    coords = c("X", "Y"),
    crs = crs
  )
  prev_sf$prev <- prev_df$prev
  
  # Plot
  p <- ggplot() +
    geom_sf(data = prev_sf, aes(color = prev), size = 1.5, alpha = 0.5) +
    scale_color_gradientn(
        colours = c("#1f363d",
        "#40798c",
        "#70a9a1",
        "#f46036",
        "#a4243b"
        ),
      name = "Seroprevalence",
      limits = c(0, max(prev_sf$prev, na.rm = TRUE)),
      labels = scales::percent_format(accuracy = 1),
      guide = guide_colorbar(
      direction = "horizontal",
      barheight = unit(0.25, "cm"),
      barwidth  = unit(6.5, "cm"),
      title.position = "top",
      label.position = "bottom",
      ticks = TRUE,
      ticks.length = unit(0.1, "cm")
  )
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 20),
      axis.text = element_blank(),      # Remove axis text (lat/long labels)
      axis.ticks = element_blank(),
      legend.text = element_text(size = 14),   # tick label size
      legend.title = element_text(size = 14),  # "Seroprevalence" title size
      ) +
    
    annotation_scale(
      bar_cols = c("black", "white"),  # alternating black/white like the reference
      height = unit(0.2, "cm"),
      text_family = "sans", 
      text_cex = 1.5
  )
    
  
  return(list(
    plot = p,
    prev_sf = prev_sf,
    prev_range = range(prev_loc)
  ))
}





# --- Function to extract and plot annual infections ---
predicted_annual_infections <- function(foi_result, model, age_groups, age_weights, cam_pop, 
                                             agg_factor = 10, crs = 32633, pathogen_name = "ONNV") {
  
  # Aggregate population raster to manageable resolution
  cam_pop_agg <- terra::aggregate(cam_pop, fact = agg_factor, fun = sum, na.rm = TRUE)
  
  # Transform FOI sf to match cam_pop CRS before rasterizing
  foi_sf_transformed <- st_transform(foi_result$foi_sf, crs = terra::crs(cam_pop_agg))
  
  # Convert FOI predictions to raster matching the population grid
  foi_raster <- terra::rasterize(
    terra::vect(foi_sf_transformed),
    cam_pop_agg,
    field = "foi",
    fun = mean
  )
  
  # Extract values from both rasters as vectors
  cam_pop_vals <- terra::values(cam_pop_agg, mat = FALSE)
  foi_vals <- terra::values(foi_raster, mat = FALSE)
  
  # Get coordinates of raster cells for plotting later
  coords_xy <- terra::xyFromCell(cam_pop_agg, 1:terra::ncell(cam_pop_agg))
  
  n_cells <- length(cam_pop_vals)
  n_age_groups <- nrow(age_groups)
  infections_by_age <- matrix(0, nrow = n_cells, ncol = n_age_groups)
  seroprevalence_by_age <- matrix(0, nrow = n_cells, ncol = n_age_groups)
  susceptible_by_age <- matrix(0, nrow = n_cells, ncol = n_age_groups)
  
  # Calculate infections for each cell and age group
  for (i in 1:n_cells) {
    total_pop <- cam_pop_vals[i]
    lambda <- foi_vals[i]
    
    # Skip if missing data or zero population
    if (is.na(total_pop) || is.na(lambda) || total_pop == 0) next
    
    for (j in 1:n_age_groups) {
      a_lower <- age_groups$age_lower[j]
      a_upper <- age_groups$age_upper[j]
      age_width <- a_upper - a_lower
      
      # Population in this age group at this location
      N_age_loc <- total_pop * age_weights[j]
      
      # Average susceptible proportion in this age group
      # S(a) = exp(-λa)
      if (lambda > 1e-10) {
        avg_susceptible_prop <- (1/(lambda * age_width)) * 
                               (exp(-lambda * a_lower) - exp(-lambda * a_upper))
      } else {
        avg_susceptible_prop <- 1
      }
      
      # Average seroprevalence = 1 - average susceptible
      avg_seroprev_prop <- 1 - avg_susceptible_prop
      
      # Store proportions for this cell and age group
      susceptible_by_age[i, j] <- avg_susceptible_prop
      seroprevalence_by_age[i, j] <- avg_seroprev_prop
      
      # Annual infections = Population × FOI × Average susceptible
      infections_by_age[i, j] <- N_age_loc * lambda * avg_susceptible_prop
    }
  }
  
  # Sum across age groups to get total annual infections per cell
  annual_infections <- rowSums(infections_by_age)
  total_annual_infections <- sum(annual_infections)
  
  # Calculate overall susceptible and seroprevalence weighted by population and age
  # For each cell, weight by age distribution
  susceptible_prop_by_cell <- as.vector(susceptible_by_age %*% age_weights)
  seroprev_prop_by_cell <- as.vector(seroprevalence_by_age %*% age_weights)
  
  # Calculate Cameroon-wide proportions weighted by population
  total_pop_cameroon <- sum(cam_pop_vals, na.rm = TRUE)
  
  # For each cell: proportion × population = number of people
  susceptible_people <- susceptible_prop_by_cell * cam_pop_vals
  seropositive_people <- seroprev_prop_by_cell * cam_pop_vals
  
  # Sum across all cells and divide by total population
  cameroon_susceptible_prop <- sum(susceptible_people, na.rm = TRUE) / total_pop_cameroon
  cameroon_seropositive_prop <- sum(seropositive_people, na.rm = TRUE) / total_pop_cameroon
  
  cat("\n=== Cameroon-wide Summary ===\n")
  cat("Total population: ", scales::comma(round(total_pop_cameroon)), "\n")
  cat("Proportion susceptible (seronegative): ", round(cameroon_susceptible_prop * 100, 2), "%\n")
  cat("Proportion seropositive (immune): ", round(cameroon_seropositive_prop * 100, 2), "%\n")
  cat("Total annual infections: ", scales::comma(round(total_annual_infections)), "\n")
  cat("=============================\n\n")
  
  # Create sf object for plotting (only non-zero infection cells to reduce size)
  valid_cells <- which(annual_infections > 0)
  
  infections_sf <- st_as_sf(
    data.frame(
      X = coords_xy[valid_cells, 1],
      Y = coords_xy[valid_cells, 2],
      infections = annual_infections[valid_cells],
      population = cam_pop_vals[valid_cells],
      susceptible_prop = susceptible_prop_by_cell[valid_cells],
      seroprevalence = seroprev_prop_by_cell[valid_cells]
    ),
    coords = c("X", "Y"),
    crs = terra::crs(cam_pop_agg)
  )
  
  # Plot
  p <- ggplot() +
    geom_sf(data = infections_sf, aes(color = infections), 
            size = 1.5, alpha = 0.5) +
    scale_color_viridis_c(
      option = "mako",
      trans = "log10",
      name = "Annual Infections",
      labels = scales::comma_format(),
      guide = guide_colorbar(
        direction = "horizontal",
        barheight = unit(0.25, "cm"),
        barwidth  = unit(6.5, "cm"),
        title.position = "top",
        label.position = "bottom",
        ticks = TRUE,
        ticks.length = unit(0.1, "cm")
      )
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 20),
      axis.text = element_blank(),      # Remove axis text (lat/long labels)
      axis.ticks = element_blank(),
      legend.text = element_text(size = 14),   # tick label size
      legend.title = element_text(size = 14),  # "Seroprevalence" title size
    ) +
  annotation_scale(
    bar_cols = c("black", "white"),  # alternating black/white like the reference
    height = unit(0.2, "cm"),
    text_family = "sans", 
    text_cex = 1.5
  )
  
  # Calculate infections by age group (summed across all cells)
  total_infections_by_age <- colSums(infections_by_age)
  
  age_group_summary <- data.frame(
    Age_Group = cameroon_age_2025$Age,
    Age_Lower = age_groups$age_lower,
    Age_Upper = age_groups$age_upper,
    Annual_Infections = total_infections_by_age,
    Proportion = total_infections_by_age / sum(total_infections_by_age)
  )
  
  # Return plot, data, and summary statistics
  return(list(
    plot = p,
    infections_sf = infections_sf,
    infections_range = range(annual_infections[annual_infections > 0]),
    total_infections = total_annual_infections,
    infections_by_age = age_group_summary,
    infections_matrix = infections_by_age,
    population_captured = total_pop_cameroon,
    cameroon_susceptible_prop = cameroon_susceptible_prop,
    cameroon_seropositive_prop = cameroon_seropositive_prop,
    susceptible_people = sum(susceptible_people, na.rm = TRUE),
    seropositive_people = sum(seropositive_people, na.rm = TRUE)
  ))
}


# ---- Plot proportion positive by age 
plot_age_seroprevalence <- function(data, chains_df, infM, pathogen_col, pathogen_name) {
  
  # Find which components have pathogen_col = 1
  nC <- nrow(infM)
  positive_components <- which(infM[, pathogen_col] == 1)
  
  # Age groups
  age_breaks <- c(0, 5, 10, 16, 23, 31, 40, 50, 100)
  age_labels <- c("0-4", "5-9", "10-15", "16-22", "23-30", "31-39", "40-49", "50+")
  
  # Filter out NA ages and track which rows we're keeping
  data_plot <- data %>%
    mutate(original_row = row_number()) %>%  # Track original indices
    filter(!is.na(AgeInYears))  # Remove rows with NA ages
  

  kept_indices <- data_plot$original_row
  N_kept <- length(kept_indices)
  
  data_plot$age_group <- cut(
    data_plot$AgeInYears,
    breaks = age_breaks,
    labels = age_labels,
    include.lowest = TRUE,
    right = FALSE
  )
  
  # Extract probabilities ONLY for individuals we're keeping
  prob_cols_list <- lapply(positive_components, function(comp) {
    sprintf("post_prob[%d,%d]", kept_indices, comp)  # Use kept_indices instead of 1:N
  })
  
  # Sum probabilities across all positive components for each draw
  probs_all_draws <- Reduce(`+`, lapply(prob_cols_list, function(cols) {
    as.matrix(chains_df[, cols])
  }))
  
  n_draws <- nrow(probs_all_draws)
  
  # Calculate prevalence by year and age group for each draw
  prevalence_draws <- map_dfr(1:n_draws, function(draw_num) {
    probs_this_draw <- probs_all_draws[draw_num, ]
    
    data_plot %>%
      dplyr::mutate(prob_pos = probs_this_draw) %>%
      group_by(year_of_survey, age_group) %>%
      summarise(
        prevalence = mean(prob_pos, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(draw = draw_num)
  })
  
  # Summarize across draws
  obs <- prevalence_draws %>%
    group_by(year_of_survey, age_group) %>%
    summarise(
      proportion_positive = median(prevalence),
      obs_lower = quantile(prevalence, 0.025),
      obs_upper = quantile(prevalence, 0.975),
      n = first(n),
      .groups = "drop"
    )
  
  y_limits <- if (pathogen_name == "CHIK") c(0, 0.08) else c(0, 0.8)
  
  p <- ggplot(obs, aes(x = age_group)) +
    geom_errorbar(aes(ymin = obs_lower, ymax = obs_upper), 
                  width = 0.15, linewidth = 0.8,color = '#057cfc') +      
    geom_point(aes(y = proportion_positive), size = 2, color = '#057cfc') +  
    facet_wrap(~ year_of_survey, ncol = 5) +
    scale_y_continuous(limits = y_limits) +
    labs(
      x = "Age group",
      y = "Seroprevalence",
      title = paste("Model-estimated seroprevalence by age -", pathogen_name)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks = element_line(color = "black", size = 0.5),
      panel.grid = element_blank(),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
      axis.title = element_text(size = 24),
      aspect.ratio = 2,
      strip.text = element_text(size = 20),
    )
  
  print(p)
  return(list(plot = p, data = obs, draws = prevalence_draws))
}

# ---- Plot proportion positive by age stratified by gender
plot_age_seroprevalence_by_year_gender_obs_binary <- function(data_original, positive_col) {
  
  # Recreate the filtered dataset used in the model
  data_plot <- data_original
  data_plot$year_of_survey <- data_original$year_of_survey
  
  # Filter to only Sex = 1 (Male) or 2 (Female)
  data_plot <- data_plot[data_plot$Sex %in% c(1, 2), ]
  
  data_plot$sex_label <- factor(data_plot$Sex, 
                                levels = c(1, 2), 
                                labels = c("Male", "Female"))
  
  # Age groups
  age_breaks <- c(0, 5, 10, 16, 23, 31, 40, 50, 100)
  age_labels <- c("0-4", "5-9", "10-15", "16-22", "23-30", "31-39", "40-49", "50+")
  data_plot$age_group <- cut(
    data_plot$AgeInYears,
    breaks = age_breaks,
    labels = age_labels,
    include.lowest = TRUE,
    right = FALSE
  )
  
  formula_mean <- as.formula(paste(positive_col, "~ year_of_survey + age_group + sex_label"))
  formula_length <- as.formula(paste(positive_col, "~ year_of_survey + age_group + sex_label"))
  
  # Summaries by year, age group, and sex
  obs <- aggregate(formula_mean, data_plot, mean, na.rm = TRUE)
  n_by <- aggregate(formula_length, data_plot, length)
  names(n_by)[4] <- "n"
  names(obs)[4] <- "proportion_positive"
  
  obs <- merge(obs, n_by, by = c("year_of_survey", "age_group", "sex_label"))
  obs$obs_lower <- pmax(0, obs$proportion_positive - 1.96*sqrt(obs$proportion_positive*(1-obs$proportion_positive)/obs$n))
  obs$obs_upper <- pmin(1, obs$proportion_positive + 1.96*sqrt(obs$proportion_positive*(1-obs$proportion_positive)/obs$n))
  
  y_limits <- if (positive_col == "CHIK_pos") {
    c(0, 0.3)
  } else {
    c(0, 0.8)
  }
  
  # Plot with sex differentiation
  p <- ggplot() +
    geom_point(data = obs, aes(x = age_group, y = proportion_positive, color = sex_label), 
               size = 2, position = position_dodge(width = 0.5)) +
    geom_errorbar(data = obs, aes(x = age_group, ymin = obs_lower, ymax = obs_upper, color = sex_label), 
                  width = 0.15, position = position_dodge(width = 0.5)) +
    
    facet_wrap(~ year_of_survey, ncol = 5) +
    scale_y_continuous(limits = y_limits) +
    scale_color_manual(values = c("Male" = "#0f4c5c", "Female" = "#90a955")) +
    labs(
      x = "Age group",
      y = "Proportion seropositive",
      title = paste("Observed seroprevalence by age group and sex -", positive_col),
      color = "Sex"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  invisible(p)
  
  return(data_plot)
}

# --- Seroprevalence by age group by year - model fits ---
plot_age_seroprevalence_model_fits <- function(result, data, model_data, chains_df, infM, pathogen_col) {
  
  data_plot <- data
  meta_data_with_labels <- meta_data_with_labels[!is.na(meta_data_with_labels$cluster), ]
  data_plot$year_of_survey <- as.numeric(substr(data_plot$Sample, 1, 4))

  kept_ids <- data_plot$id
  
  # Age groups
  age_breaks <- c(0, 5, 10, 16, 23, 31, 40, 50, 100)
  age_labels <- c("0-4", "5-9", "10-15", "16-22", "23-30",
                  "31-39", "40-49", "50+")
  data_plot$age_group <- cut(
    data_plot$AgeInYears,
    breaks = age_breaks,
    labels = age_labels,
    include.lowest = TRUE,
    right = FALSE
  )
  
  idx_est <- inla.stack.index(result$stk.full, tag = "est")$data
  fit     <- result$output$summary.fitted.values[
    idx_est, c("mean", "0.025quant", "0.975quant")
  ]

  
  # Attach INLA predicted probabilities to each individual
  data_plot <- data_plot %>%
  dplyr::mutate(
    predicted = fit$mean,
    pred_lower = fit$`0.025quant`,
    pred_upper = fit$`0.975quant`
  )
  # ---- Observed summaries using posterior probabilities from chains ----
  # Find which components have pathogen_col = 1
  nC <- nrow(infM)
  positive_components <- which(infM[, pathogen_col] == 1)

  # Find the columns in chains_df corresponding to the kept ids
  prob_cols_list <- lapply(positive_components, function(comp) {
    sprintf("post_prob[%d,%d]", match(kept_ids, model_data$id), comp)
  })

  matched_idx <- match(kept_ids, model_data$id)
  cat("Any NA in posterior match:", any(is.na(matched_idx)), "\n")
  
  # Sum probabilities across all positive components for each draw
  probs_all_draws <- Reduce(`+`, lapply(prob_cols_list, function(cols) {
    as.matrix(chains_df[, cols])
  }))
  
  n_draws <- nrow(probs_all_draws)
  
  # Calculate prevalence by year and age group for each draw
  prevalence_draws <- purrr::map_dfr(1:n_draws, function(draw_num) {
    probs_this_draw <- probs_all_draws[draw_num, ]
    data_plot %>%
      dplyr::mutate(prob_pos = probs_this_draw) %>%
      group_by(year_of_survey, age_group) %>%
      summarise(
        prevalence = mean(prob_pos, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(draw = draw_num)
  })
  
  # Summarize observed data across draws
  obs <- prevalence_draws %>%
    group_by(year_of_survey, age_group) %>%
    summarise(
      obs_mean = median(prevalence),
      obs_lower = quantile(prevalence, 0.025),
      obs_upper = quantile(prevalence, 0.975),
      n = first(n),
      .groups = "drop"
    )
  
  # ---- Predicted summaries (mean of INLA predicted probabilities) ----
  pred <- aggregate(
    cbind(predicted, pred_lower, pred_upper) ~ year_of_survey + age_group,
    data_plot,
    mean,
    na.rm = TRUE
  )
 
  # ---- Plot ----
  p <- ggplot() +
    # observed (from posterior probabilities)
    geom_point(
      data = obs,
      aes(x = age_group, y = obs_mean, color = "Observed"),
      size = 5, 
      position = position_nudge(x = -0.1)
    ) +
    geom_errorbar(
      data = obs,
      aes(x = age_group, ymin = obs_lower, ymax = obs_upper, color = "Observed"),
      width = 0.15, 
      position = position_nudge(x = -0.1)
    ) +
    # predicted (INLA model fits)
    geom_point(
      data = pred,
      aes(x = age_group, y = predicted, color = "Estimated"),
      size = 5,
      position = position_nudge(x = 0.1)

    ) +
    geom_errorbar(
      data = pred,
      aes(x = age_group, ymin = pred_lower, ymax = pred_upper, color = "Estimated"),
      width = 0.15,
      position = position_nudge(x = 0.1)
    ) +
    scale_color_manual(
        name = "",
        values = c("Observed" = "#04989a", "Estimated" = "#c93a88")
      ) +
    facet_wrap(~ year_of_survey, ncol = 5) +
    labs(
      x = "Age group",
      y = "Proportion seropositive") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks.x = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      legend.position.inside = c(0.95, 0.5),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
      axis.title = element_text(size = 24),
      aspect.ratio = 1.2,
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 20),
      strip.text = element_text(size = 20),
      strip.background = element_rect(fill = "#ffffff"))



  print(p)
  invisible(list(plot = p, obs = obs, pred = pred, prevalence_draws = prevalence_draws))
}


# --- Seroprevalence by age group by year and gender ---
plot_age_seroprevalence_model_fits_by_gender <- function(result, data, model_data, chains_df, infM, pathogen_col) {
  
  data_plot <- data
  model_data <- model_data[!is.na(model_data$cluster), ]
  
  # Derived variables
  data_plot$year_of_survey <- as.numeric(substr(data_plot$Sample, 1, 4))
  
  # Age groups
  age_breaks <- c(0, 5, 10, 16, 23, 31, 40, 50, 100)
  age_labels <- c("0-4", "5-9", "10-15", "16-22", "23-30", "31-39", "40-49", "50+")
  
  data_plot$age_group <- cut(
    data_plot$AgeInYears,
    breaks = age_breaks,
    labels = age_labels,
    include.lowest = TRUE,
    right = FALSE
  )
  
  # Attach fitted values from INLA predictions BEFORE sex filtering
  idx_est <- inla.stack.index(result$stk.full, tag = "est")$data
  fit <- result$output$summary.fitted.values[
    idx_est, c("mean", "0.025quant", "0.975quant")
  ]
  
  data_plot <- data_plot %>%
    dplyr::mutate(
      predicted = fit$mean,
      pred_lower = fit$`0.025quant`,
      pred_upper = fit$`0.975quant`
    )
  
  # Keep only Male / Female
  data_plot <- data_plot[data_plot$Sex %in% c(1, 2), ]
  
  data_plot$sex_label <- factor(
    data_plot$Sex,
    levels = c(1, 2),
    labels = c("Male", "Female")
  )
  
  # Keep ids after filtering
  kept_ids <- data_plot$id
  
  # ---- Observed summaries from posterior probabilities ----
  positive_components <- which(infM[, pathogen_col] == 1)
  
  matched_idx <- match(kept_ids, model_data$id)
  cat("Any NA in posterior match:", any(is.na(matched_idx)), "\n")
  
  prob_cols_list <- lapply(positive_components, function(comp) {
    sprintf("post_prob[%d,%d]", matched_idx, comp)
  })
  
  all_prob_cols <- unlist(prob_cols_list)
  missing_cols <- setdiff(all_prob_cols, colnames(chains_df))
  
  cat("Number of missing columns:", length(missing_cols), "\n")
  if (length(missing_cols) > 0) {
    print(head(missing_cols, 20))
    stop("Some post_prob columns were not found in chains_df")
  }
  
  probs_all_draws <- Reduce(`+`, lapply(prob_cols_list, function(cols) {
    as.matrix(chains_df[, cols, drop = FALSE])
  }))
  
  n_draws <- nrow(probs_all_draws)
  
  prevalence_draws <- purrr::map_dfr(1:n_draws, function(draw_num) {
    probs_this_draw <- probs_all_draws[draw_num, ]
    
    data_plot %>%
      dplyr::mutate(prob_pos = probs_this_draw) %>%
      dplyr::group_by(year_of_survey, age_group, sex_label) %>%
      dplyr::summarise(
        prevalence = mean(prob_pos, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(draw = draw_num)
  })
  
  # Observed summaries
  obs <- prevalence_draws %>%
    dplyr::group_by(year_of_survey, age_group, sex_label) %>%
    dplyr::summarise(
      obs_mean = median(prevalence),
      obs_lower = quantile(prevalence, 0.025),
      obs_upper = quantile(prevalence, 0.975),
      n = dplyr::first(n),
      .groups = "drop"
    )
  
  # Predicted summaries
  pred <- aggregate(
    cbind(predicted, pred_lower, pred_upper) ~ year_of_survey + age_group + sex_label,
    data_plot,
    mean,
    na.rm = TRUE
  )
  
  # Shared dodge so bars and points align
  pd <- position_dodge(width = 0.7)
  
  # Plot
  p <- ggplot() +
    # Observed bars
    geom_col(
      data = obs,
      aes(x = age_group, y = obs_mean, fill = sex_label),
      position = pd,
      width = 0.8,
      alpha = 0.6,
      linewidth = 0.2
    ) +
    # Model fit points
    geom_point(
      data = pred,
      aes(x = age_group, y = predicted,group = sex_label),
      position = pd,
      size = 4,
      colour = "black", 
    ) +
    # Model fit CI
    geom_errorbar(
      data = pred,
      aes(x = age_group, ymin = pred_lower, ymax = pred_upper, group = sex_label),
      position = pd,
      width = 0.1,
      linewidth = 0.5,
      colour = "black",
    ) +
    facet_wrap(~ year_of_survey, ncol = 5) +
    scale_fill_manual(values = c("Male" = "#b84f74", "Female" = "#00798c")) +
    scale_colour_manual(values = c("Male" = "#b84f74", "Female" = "#00798c")) +
    labs(
      x = "Age group",
      y = "Proportion seropositive",
      fill = "Observed",
      colour = "Model fit"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
      axis.title = element_text(size = 24),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      strip.text = element_text(size = 20),
      panel.border = element_blank(),
      aspect.ratio = 1.2,
      legend.position = "bottom",
      strip.background = element_rect(fill = "#ffffff")
    )
  
  print(p)
  invisible(list(plot = p, obs = obs, pred = pred, prevalence_draws = prevalence_draws))
}

# --- Obs eroprevalence by age group by year and gender (multisero model)---
plot_age_seroprevalence_obs_only_by_gender <- function(data, model_data, chains_df, infM, pathogen_col) {
  
  data_plot <- data
  model_data <- model_data[!is.na(model_data$cluster), ]
  
  # Derived variables
  data_plot$year_of_survey <- as.numeric(substr(data_plot$Sample, 1, 4))
  
  # Age groups
  age_breaks <- c(0, 5, 10, 16, 23, 31, 40, 50, 100)
  age_labels <- c("0-4", "5-9", "10-15", "16-22", "23-30", "31-39", "40-49", "50+")
  
  data_plot$age_group <- cut(
    data_plot$AgeInYears,
    breaks = age_breaks,
    labels = age_labels,
    include.lowest = TRUE,
    right = FALSE
  )
  
  # Keep only Male / Female
  data_plot <- data_plot[data_plot$Sex %in% c(1, 2), ]
  
  data_plot$sex_label <- factor(
    data_plot$Sex,
    levels = c(1, 2),
    labels = c("Male", "Female")
  )
  
  # Keep ids after filtering
  kept_ids <- data_plot$id
  
  # ---- Observed summaries from posterior probabilities ----
  positive_components <- which(infM[, pathogen_col] == 1)
  
  matched_idx <- match(kept_ids, model_data$id)
  cat("Any NA in posterior match:", any(is.na(matched_idx)), "\n")
  
  prob_cols_list <- lapply(positive_components, function(comp) {
    sprintf("post_prob[%d,%d]", matched_idx, comp)
  })
  
  all_prob_cols <- unlist(prob_cols_list)
  missing_cols <- setdiff(all_prob_cols, colnames(chains_df))
  
  cat("Number of missing columns:", length(missing_cols), "\n")
  if (length(missing_cols) > 0) {
    print(head(missing_cols, 20))
    stop("Some post_prob columns were not found in chains_df")
  }
  
  probs_all_draws <- Reduce(`+`, lapply(prob_cols_list, function(cols) {
    as.matrix(chains_df[, cols, drop = FALSE])
  }))
  
  n_draws <- nrow(probs_all_draws)
  
  prevalence_draws <- purrr::map_dfr(1:n_draws, function(draw_num) {
    probs_this_draw <- probs_all_draws[draw_num, ]
    
    data_plot %>%
      dplyr::mutate(prob_pos = probs_this_draw) %>%
      dplyr::group_by(year_of_survey, age_group, sex_label) %>%
      dplyr::summarise(
        prevalence = mean(prob_pos, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(draw = draw_num)
  })
  
  # Observed summaries
  obs <- prevalence_draws %>%
    dplyr::group_by(year_of_survey, age_group, sex_label) %>%
    dplyr::summarise(
      obs_mean = median(prevalence),
      obs_lower = quantile(prevalence, 0.025),
      obs_upper = quantile(prevalence, 0.975),
      n = dplyr::first(n),
      .groups = "drop"
    )
  
  # Shared dodge
  pd <- position_dodge(width = 0.7)
  
  # Plot
  p <- ggplot() +
    # Observed bars
    geom_col(
      data = obs,
      aes(x = age_group, y = obs_mean, fill = sex_label),
      position = pd,
      width = 0.8,
      alpha = 0.6,
      linewidth = 0.2
    ) +
    # Observed CIs
    geom_errorbar(
      data = obs,
      aes(x = age_group, ymin = obs_lower, ymax = obs_upper, group = sex_label),
      position = pd,
      width = 0.1,
      linewidth = 0.5,
      colour = "black"
    ) +
    facet_wrap(~ year_of_survey, ncol = 5) +
    scale_fill_manual(values = c("Male" = "#b84f74", "Female" = "#00798c")) +
    labs(
      x = "Age group",
      y = "Proportion seropositive",
      fill = "Observed"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
      axis.title = element_text(size = 24),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      strip.text = element_text(size = 20),
      panel.border = element_blank(),
      aspect.ratio = 1.2,
      legend.position = "bottom",
      strip.background = element_rect(fill = "#ffffff")
    )
  
  print(p)
  invisible(list(plot = p, obs = obs, prevalence_draws = prevalence_draws))
}


# prediction FOI, seroprevelance and infections by region 
aggregate_predictions_by_region <- function(
  pred_sf,
  regions_sf,
  cam_pop = NULL,
  value_col,
  region_col = "region",
  agg_type = c("weighted_mean", "sum")
) {
  

  # --- CRS alignment ---
  if (st_crs(pred_sf) != st_crs(regions_sf)) {
    regions_sf <- st_transform(regions_sf, st_crs(pred_sf))
  }
  
  # --- spatial join ---
  joined <- st_join(pred_sf, regions_sf[, region_col], left = FALSE)
  
  if (nrow(joined) == 0) {
    stop("No prediction points fall inside regions.")
  }
  
  if (agg_type == "weighted_mean") {
    
    if (is.null(cam_pop)) {
      stop("cam_pop raster required for weighted_mean aggregation")
    }
    
    joined$population <- terra::extract(
      cam_pop,
      terra::vect(joined))[,2]
    
    # remove missing population
    joined <- joined[!is.na(joined$population) & joined$population > 0, ]
    
    # --- weighted aggregation ---
    out <- joined %>%
      dplyr::group_by(.data[[region_col]]) %>%
      dplyr::summarise(
        value_weighted = weighted.mean(
          .data[[value_col]],
          population,
          na.rm = TRUE
        ),
        value_unweighted = mean(.data[[value_col]], na.rm = TRUE),
        total_population = sum(population, na.rm = TRUE),
        n_points = dplyr::n(),
        .groups = "drop"
      )
    
  } else {
    
    # --- sum aggregation (for infections) ---
    out <- joined %>%
      dplyr::group_by(.data[[region_col]]) %>%
      dplyr::summarise(
        total_value = sum(.data[[value_col]], na.rm = TRUE),
        n_points = dplyr::n(),
        .groups = "drop"
      )
  }
  
  return(out)
}



# --- Propotion by mosquito distribution + log population (using binary serostatus)
calculate_prop_by_variable <- function(data, var_col, positive_col, breaks_max, breaks_min) {
  
  #remove NA from positive_col
  data <- data[!is.na(data[[positive_col]]), ]
  var_mid <- rep(NaN, length(breaks_max))
  prop_pos <- matrix(NaN, length(breaks_max), 3)
  
  for (i in 1:length(breaks_max)) {

      tmp <- which(
      data[[var_col]] < breaks_max[i] &
      data[[var_col]] >= breaks_min[i]
    )


    if (length(tmp) > 5) {
      prop_pos[i, 1] <- mean(data[[positive_col]][tmp], na.rm = TRUE)
      a <- prop.test(sum(data[[positive_col]][tmp]), length(tmp))
      prop_pos[i, 2:3] <- a$conf.int
      var_mid[i] <- mean(data[[var_col]][tmp], na.rm = TRUE)
    }
  }

  obs_df <- data.frame(
    x = var_mid,
    y = prop_pos[, 1],
    ymin = prop_pos[, 2],
    ymax = prop_pos[, 3]
  )

  # binomial regression
  model_df <- data[, c(var_col, positive_col)]
  model_df <- model_df[complete.cases(model_df), ]
  
  formula <- as.formula(paste(positive_col, "~", var_col))
  
  log_model <- glm(
    formula,
    family = binomial,
    data = model_df
  )

  list(
    obs = obs_df,
    log_model = log_model
  )
}


# --- Propotion by mosquito distribution + log population (using binary serostatus)
calculate_prop_by_variable_NEW <- function(data, var_col, positive_col, breaks_max, breaks_min) {

  data <- data[!is.na(data[[positive_col]]), ]
  var_mid  <- rep(NaN, length(breaks_max))
  prop_pos <- matrix(NaN, length(breaks_max), 3)

  n_bins <- length(breaks_max)

  for (i in 1:n_bins) {

    if (i == n_bins) {
      # final bin: inclusive upper edge so value == breaks_max is kept
      tmp <- which(data[[var_col]] <= breaks_max[i] &
                   data[[var_col]] >= breaks_min[i])
    } else {
      tmp <- which(data[[var_col]] <  breaks_max[i] &
                   data[[var_col]] >= breaks_min[i])
    }

    if (length(tmp) > 5) {
      prop_pos[i, 1]   <- mean(data[[positive_col]][tmp], na.rm = TRUE)
      a                <- prop.test(sum(data[[positive_col]][tmp]), length(tmp))
      prop_pos[i, 2:3] <- a$conf.int
      var_mid[i]       <- (breaks_min[i] + breaks_max[i]) / 2   # <-- midpoint, not mean
    }
  }

  obs_df <- data.frame(
    x    = var_mid,
    y    = prop_pos[, 1],
    ymin = prop_pos[, 2],
    ymax = prop_pos[, 3]
  )

  model_df <- data[, c(var_col, positive_col)]
  model_df <- model_df[complete.cases(model_df), ]

  formula   <- as.formula(paste(positive_col, "~", var_col))
  log_model <- glm(formula, family = binomial, data = model_df)

  list(obs = obs_df, log_model = log_model)
}



# --- Plot functions 
make_plot_onnv <- function(df_obs, raw_data, xlab, color, pos_col = "ONNV_pos") {

  ylab <- NULL

  obs_clean <- df_obs[!is.nan(df_obs$x), ]

  obs_clean$ymin <- pmax(obs_clean$ymin, 0)
  obs_clean$ymax <- pmin(obs_clean$ymax, 0.5)

  hist_df <- data.frame(x = as.numeric(raw_data))
  hist_df <- hist_df[!is.na(hist_df$x), , drop = FALSE]

    x_scale <- scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    expand = c(0, 0)
  )

  plot_hist <- ggplot(hist_df, aes(x = x)) +
    geom_histogram(fill = color, alpha = 0.5, bins = 30, color = NA) +
    x_scale +
    labs(x = NULL, y = NULL) +
    base_theme +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(t = 6, r = 14, b = 10, l = 14)
    )

  plot_scatter <- ggplot(obs_clean, aes(x = x, y = y)) +
    geom_point(color = color, size = 4, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0, color = color, alpha = 0.6, linewidth = 0.6
    ) +
    x_scale +
    coord_cartesian(ylim = c(0, 0.4), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, 0.4, 0.1)) +
    labs(x = xlab, y = ylab) +
    base_theme +
    theme(
      plot.margin = margin(t = 10, r = 14, b = 12, l = 14)
    )

  plot_hist / plot_scatter +
    patchwork::plot_layout(heights = c(2, 4))
  
}


make_plot_chik <- function(df_obs, raw_data, xlab, color, pos_col = "ONNV_pos") {

  ylab <- NULL

  obs_clean <- df_obs[!is.nan(df_obs$x), ]

  # truncate CIs to [0, 0.5]
  obs_clean$ymin <- pmax(obs_clean$ymin, 0)
  obs_clean$ymax <- pmin(obs_clean$ymax, 0.1)

  hist_df <- data.frame(x = as.numeric(raw_data))
  hist_df <- hist_df[!is.na(hist_df$x), , drop = FALSE]

    x_scale <- scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = c("0", "0.25", "0.5", "0.75", "1"),
    expand = c(0, 0)
  )

  plot_hist <- ggplot(hist_df, aes(x = x)) +
    geom_histogram(fill = color, alpha = 0.5, bins = 30, color = NA) +
    x_scale +
    labs(x = NULL, y = NULL) +
    base_theme +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(t = 6, r = 14, b = 10, l = 14)
    )

  plot_scatter <- ggplot(obs_clean, aes(x = x, y = y)) +
    geom_point(color = color, size = 4, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0, color = color, alpha = 0.6, linewidth = 0.6
    ) +
    x_scale +
    coord_cartesian(ylim = c(0, 0.05), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, 0.05, 0.01)) +
    labs(x = xlab, y = ylab) +
    base_theme +
    theme(
      plot.margin = margin(t = 10, r = 14, b = 12, l = 14)
    )

  plot_hist / plot_scatter +
    patchwork::plot_layout(heights = c(2, 4))
  
}



calculate_prop_by_variable_multisero_probs <- function(data, var_col, chains_df, infM, pathogen_col, breaks_max, breaks_min) {
  
  positive_components <- which(infM[, pathogen_col] == 1)
  
  data_plot <- data %>%
    filter(!is.na(.data[[var_col]]), !is.na(stan_idx_full_model))
  
  kept_indices <- data_plot$stan_idx_full_model  # <-- use this, not $id
  
  nrow(data_plot)

  prob_cols_list <- lapply(positive_components, function(comp) {
    cols <- sprintf("post_prob[%d,%d]", kept_indices, comp)
    missing_cols <- setdiff(cols, colnames(chains_df))
    if (length(missing_cols) > 0) {
      stop("Missing posterior columns: ", paste(head(missing_cols), collapse = ", "))
    }
    cols
  })
  
  probs_all_draws <- Reduce(`+`, lapply(prob_cols_list, function(cols) {
    as.matrix(chains_df[, cols])
  }))


  n_draws <- nrow(probs_all_draws)

  # --- Bin individuals by the continuous variable
  bin_indices <- lapply(seq_along(breaks_max), function(i) {
    which(data_plot[[var_col]] >= breaks_min[i] & data_plot[[var_col]] < breaks_max[i])
  })

  # --- For each draw, compute mean probability in each bin
  prevalence_draws <- map_dfr(1:n_draws, function(draw_num) {
    probs_this_draw <- probs_all_draws[draw_num, ]

    map_dfr(seq_along(breaks_max), function(i) {
      idx <- bin_indices[[i]]
      if (length(idx) <= 5) return(NULL)          # same minimum-n guard as before

      tibble(
        bin         = i,
        var_mid     = mean(data_plot[[var_col]][idx], na.rm = TRUE),
        prevalence  = mean(probs_this_draw[idx],     na.rm = TRUE),
        n           = length(idx),
        draw        = draw_num
      )
    })
  })

  # --- Summarise across draws: median + 95 % credible interval
  obs <- prevalence_draws %>%
    group_by(bin, var_mid, n) %>%
    dplyr::summarise(
      y    = median(prevalence),
      ymin = quantile(prevalence, 0.025),
      ymax = quantile(prevalence, 0.975),
      .groups = "drop"
    ) %>%
    dplyr::select(x = var_mid, y, ymin, ymax)         
 

    # posterior mean per individual
    posterior_mean_probs <- colMeans(probs_all_draws)


    #logistic regression 
    model_df <- data_plot %>%
    mutate(prob_positive = posterior_mean_probs)

    formula <- as.formula(paste("prob_positive ~", var_col))

    log_model <- glm(
      formula,
      family = quasibinomial,
      data = model_df
    )

    model_df$logit_pred <- predict(log_model, type = "response")


  return(list( obs = obs,
  model_df = model_df,
  log_model = log_model))
}
 
