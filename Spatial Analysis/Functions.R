
# ---- Function for INLA models ----
run_inla <- function(year_intro, data, cam_pop, positive_col) {
  
  # Calculate years of exposure
  data$age_intro <- data$year_of_survey - year_intro
  data$years_of_exposure <- pmin(data$age_intro, data$AgeInYears)
  
  data <- data[!is.na(data$years_of_exposure) & data$years_of_exposure > 0, ]
  
  # Check if we have enough data
  if (nrow(data) < 50) {
    warning(paste0("Year ", year_intro, " has too few valid observations (", nrow(data), "). Skipping."))
    return(NULL)
  }
  
  data_points <- data %>%
    st_drop_geometry() %>%
    filter(!is.na("Easting") & !is.na("Northing")) 
  
  
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
  
  # Estimation stack - USE THE INPUT PARAMETER
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



# --- Function to extract and plot FOI ---
extract_and_plot_foi <- function(model, coop, pathogen_name = "ONNV") { 
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
      data = foi_sf, aes(color = foi), size = 1.7, alpha = 1) +
    scale_color_viridis(
      option = "mako",
      name = "FOI (λ)",
      limits = c(0, max(foi_sf$foi))) +
    labs(
      title = paste0("Force of Infection (FOI) Predictions - ", pathogen_name),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() + 
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold")
    )
  
  print(p)
  
  # Return foi_sf and plot invisibly
  invisible(list(
    foi_sf = foi_sf,
    foi_df = foi_df,
    plot = p
  ))
}


# --- CORRECT Function to extract and plot seroprevalence ---
plot_predicted_seroprevalence <- function(foi_result, model, age_groups, age_weights,
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
    geom_sf(data = prev_sf, aes(color = prev), size = 1.7, alpha = 1, shape = 15) +
    scale_color_viridis_c(
      option = "mako",
      name = "Seroprevalence",
      limits = c(0, max(prev_sf$prev, na.rm = TRUE)),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      title = paste0("Predicted Seroprevalence - ", pathogen_name),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  
  print(p)
  
  return(list(
    plot = p,
    prev_sf = prev_sf,
    prev_range = range(prev_loc)
  ))
}


# --- Function to extract and plot annual infections ---
plot_predicted_annual_infections <- function(foi_result, model, age_groups, age_weights, cam_pop, 
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
  
  # Initialize
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
            size = 1.7, alpha = 1, shape = 15) +
    scale_color_viridis_c(
      option = "plasma",
      name = "Annual\nInfections",
      trans = "log10",
      labels = scales::comma_format()
    ) +
    labs(
      title = paste0("Predicted Annual Infections - ", pathogen_name),
      subtitle = paste0("Total: ", scales::comma(round(total_annual_infections)), 
                       #" | Seropositive: ", round(cameroon_seropositive_prop * 100, 1), 
                       "% | Susceptible: ", round(cameroon_susceptible_prop * 100, 1), "%"
       ),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  
  print(p)
  
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


# --- Prevelance Patterns 

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
  
  # Get the indices of rows we're keeping
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

# --- Seroprevalence by age group by year - model fits ---
plot_age_seroprevalence_model_fits <- function(year_intro, result, data, chains_df, infM, pathogen_col) {
  
  # Recreate the filtered dataset used in the model
  data_plot                   <- subset(data, !is.na(Latitude) & !is.na(Longitude))
  data_plot$year_of_survey    <- as.numeric(substr(data_plot$Sample, 1, 4))
  data_plot$age_intro         <- data_plot$year_of_survey - year_intro
  data_plot$years_of_exposure <- pmin(data_plot$age_intro, data_plot$AgeInYears)
  data_plot <- subset(data_plot, !is.na(years_of_exposure) & years_of_exposure > 0)
  
  # Track original row numbers before any filtering
  data_plot$original_row <- as.numeric(rownames(data_plot))
  
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
  
  # Attach fitted values from estimation stack (INLA predictions)
  idx_est <- inla.stack.index(result$stk.full, tag = "est")$data
  fit     <- result$output$summary.fitted.values[
    idx_est, c("mean", "0.025quant", "0.975quant")
  ]
  
  cat("Fit dimensions:", nrow(fit), "x", ncol(fit), "\n")

  # Basic alignment checks
  if (!is.data.frame(fit))
    stop("fit is not a data.frame. Check result$output$summary.fitted.values.")
  if (nrow(fit) != nrow(data_plot)) {
    stop(sprintf(
      "Row mismatch. fit=%d, data_plot=%d. Build data_plot in the exact order used to make the stack.",
      nrow(fit), nrow(data_plot)
    ))
  }
  
  # Attach INLA predicted probabilities to each individual
  data_plot$predicted  <- fit$mean
  data_plot$pred_lower <- fit$`0.025quant`
  data_plot$pred_upper <- fit$`0.975quant`
  
  # ---- Observed summaries using posterior probabilities from chains ----
  # Find which components have pathogen_col = 1
  nC <- nrow(infM)
  positive_components <- which(infM[, pathogen_col] == 1)
  
  # Get the indices we're keeping (after all filtering)
  kept_indices <- data_plot$original_row
  N_kept <- length(kept_indices)
  
  # Extract probabilities ONLY for individuals we're keeping
  prob_cols_list <- lapply(positive_components, function(comp) {
    sprintf("post_prob[%d,%d]", kept_indices, comp)
  })
  
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

  cat("Obs summary rows:", nrow(obs), "\n")
  print(head(obs))
  
  # ---- Predicted summaries (mean of INLA predicted probabilities) ----
  pred <- aggregate(
    cbind(predicted, pred_lower, pred_upper) ~ year_of_survey + age_group,
    data_plot,
    mean,
    na.rm = TRUE
  )
  print(head(pred))
  
  # ---- Plot ----
  p <- ggplot() +
    # observed (from posterior probabilities)
    geom_point(
      data = obs,
      aes(x = age_group, y = obs_mean, color = "Observed"),
      size = 2, 
      color = "#0d1b2a"
    ) +
    geom_errorbar(
      data = obs,
      aes(x = age_group, ymin = obs_lower, ymax = obs_upper),
      width = 0.15, 
      color = '#0d1b2a'
    ) +
    # predicted (INLA model fits)
    geom_point(
      data = pred,
       aes(x = age_group, y = predicted, color = "Estimated"),
      color = "#0a9396",
      size = 2
    ) +
    geom_errorbar(
      data = pred,
      aes(x = age_group, ymin = pred_lower, ymax = pred_upper),
      color = "#0a9396",
      width = 0.15
    ) +
    facet_wrap(~ year_of_survey, ncol = 5) +
    labs(
      x = "Age group",
      y = "Proportion seropositive",
      title = "Observed vs fitted seroprevalence by age group",
      subtitle = paste0(
        "Pathogen column: ", pathogen_col,
        " | Year of introduction: ", year_intro
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      aspect.ratio = 1
    )
  
  print(p)
  invisible(list(plot = p, obs = obs, pred = pred, prevalence_draws = prevalence_draws))
}

# --- Propotion by mosquito distribution + log population 
calculate_prop_by_variable <- function(data, var_col, positive_col, breaks_max, breaks_min) {
  var_mid <- rep(NaN, length(breaks_max))
  prop_pos <- matrix(NaN, length(breaks_max), 3)
  
  for (i in 1:length(breaks_max)) {
    tmp <- which(data[[var_col]] < breaks_max[i] & data[[var_col]] >= breaks_min[i])
    if (length(tmp) > 5) {
      prop_pos[i, 1] <- mean(data[[positive_col]][tmp], na.rm = TRUE)
      a <- prop.test(sum(data[[positive_col]][tmp]), length(tmp))
      prop_pos[i, 2:3] <- a$conf.int
      var_mid[i] <- mean(data[[var_col]][tmp], na.rm = TRUE)
    }
  }
  
  data.frame(
    x = var_mid,
    y = prop_pos[, 1],
    ymin = prop_pos[, 2],
    ymax = prop_pos[, 3]
  )
}




# OLD FUNCTIONS 

# --- Seroprevalence by age group by year - model fits ---
old_plot_age_seroprevalence_model_fits <- function(year_intro, result,data, positive_col) {
  # Recreate the filtered dataset used in the model
  data_plot                   <- subset(data, !is.na(Latitude) & !is.na(Longitude))
  data_plot$year_of_survey    <- as.numeric(substr(data_plot$Sample, 1, 4))
  data_plot$age_intro         <- data_plot$year_of_survey - year_intro
  data_plot$years_of_exposure <- pmin(data_plot$age_intro, data_plot$AgeInYears)
  data_plot <- subset(data_plot, !is.na(years_of_exposure) & years_of_exposure > 0)
  
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
  
  # Attach fitted values from estimation stack
  idx_est <- inla.stack.index(result$stk.full, tag = "est")$data
  fit     <- result$output$summary.fitted.values[
    idx_est, c("mean", "0.025quant", "0.975quant")
  ]
  
  # Basic alignment checks
  if (!is.data.frame(fit))
    stop("fit_tab is not a data.frame. Check result$output$summary.fitted.values.")
  
  if (nrow(fit) != nrow(data_plot)) {
    stop(sprintf(
      "Row mismatch. fit_tab=%d, data_plot=%d. Build data_plot in the exact order used to make the stack.",
      nrow(fit), nrow(data_plot)
    ))
  }
  
  data_plot$predicted  <- fit$mean
  data_plot$pred_lower <- fit$`0.025quant`
  data_plot$pred_upper <- fit$`0.975quant`
  
  # ---- Observed summaries (dynamic positive_col) ----
  formula_mean   <- as.formula(paste(positive_col, "~ year_of_survey + age_group"))
  formula_length <- formula_mean
  
  obs <- aggregate(formula_mean, data_plot, mean, na.rm = TRUE)
  n_by <- aggregate(formula_length, data_plot, length)
  names(n_by)[3] <- "n"
  
  obs <- merge(obs, n_by, by = c("year_of_survey", "age_group"))
  names(obs)[names(obs) == positive_col] <- "obs_mean"
  
  obs$obs_lower <- pmax(
    0,
    obs$obs_mean - 1.96 * sqrt(obs$obs_mean * (1 - obs$obs_mean) / obs$n)
  )
  obs$obs_upper <- pmin(
    1,
    obs$obs_mean + 1.96 * sqrt(obs$obs_mean * (1 - obs$obs_mean) / obs$n)
  )
  
  # ---- Predicted summaries ----
  pred <- aggregate(
    cbind(predicted, pred_lower, pred_upper) ~ year_of_survey + age_group,
    data_plot,
    mean,
    na.rm = TRUE
  )
  
  # ---- Plot ----
  p <- ggplot() +
    # observed
    geom_point(
      data = obs,
      aes(x = age_group, y = obs_mean),
      size = 2, 
      color = "#0d1b2a"
    ) +
    geom_errorbar(
      data = obs,
      aes(x = age_group, ymin = obs_lower, ymax = obs_upper),
      width = 0.15, 
      color = '#0d1b2a'
    ) +
    # predicted (line connecting dots)
    geom_line(
      data = pred,
      aes(x = age_group, y = predicted, group = 1),
      color = "#0a9396",
      linewidth = 0.8
    ) +
    geom_point(
      data = pred,
      aes(x = age_group, y = predicted),
      color = "#0a9396",
      size = 2
    ) +
    geom_errorbar(
      data = pred,
      aes(x = age_group, ymin = pred_lower, ymax = pred_upper),
      color = "#0a9396",
      width = 0.15
    ) +
    facet_wrap(~ year_of_survey, ncol = 5) +
    labs(
      x = "Age group",
      y = "Proportion seropositive",
      title = "Observed vs fitted seroprevalence by age group",
      subtitle = paste0(
        "Marker: ", positive_col,
        " | Year of introduction: ", year_intro
      )
    ) +
    theme_bw() +
    theme(axis.text.x = element_text())
  
  print(p)
  invisible(p)
}
