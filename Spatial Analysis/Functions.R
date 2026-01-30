
# ---- Function for INLA models ----
run_inla <- function(year_intro, data, cameroon, positive_col) {
  
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
  grid        <- st_make_grid(cameroon, cellsize = 0.05, what = "centers")
  grid_sf     <- st_sf(grid_id = seq_along(grid), geometry = grid)
  grid_inside <- st_intersection(grid_sf, cameroon)
  grid_inside <- grid_inside[order(grid_inside$grid_id), ]  
  grid_utm    <- st_transform(grid_inside, crs = 32633)
  coop_utm    <- st_coordinates(grid_utm)
  # For INLA - convert to km (INLA works better with smaller numbers)
  coop <- coop_utm / 1000
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
crs = 32633, pathogen_name = "ONNV") {
  # Calculate total population by age group (M + F)
  N_age <- cameroon_age_2025$M + cameroon_age_2025$F
  
  # Get lambda (FOI) at each location
  lambda_pred <- foi_result$foi_sf$foi
  
  # Initialize matrix to store infections by location and age group
  n_locations <- length(lambda_pred)
  n_age_groups <- nrow(age_groups)
  
  
  # Convert prediction locations to sf if not already
  pred_coords_sf <- st_as_sf(
    data.frame(X = model$coop[, "X"] * 1000, 
               Y = model$coop[, "Y"] * 1000),
    coords = c("X", "Y"),
    crs = crs
  )

  pop_at_locations <- terra::extract(cam_pop, pred_coords_sf, ID = FALSE)[,1]
  infections_by_age <- matrix(0, nrow = n_locations, ncol = n_age_groups)
  
  # Calculate infections for each location and age group
  # Formula: N(a) × lambda × exp(-lambda × age)
  for (i in 1:n_locations) {
    lambda <- lambda_pred[i]
    total_pop <- pop_at_locations[i]

    if (total_pop == 0 || is.na(total_pop)) next
      
      for (j in 1:n_age_groups) {
        a_lower <- age_groups$age_lower[j]
        a_upper <- age_groups$age_upper[j]
        age_width <- a_upper - a_lower

      # Population in this age group at this location
      # Distribute total population according to national age distribution
        N_age_loc <- total_pop * age_weights[j]
      
      # Average annual incidence in this age group
      # This is λ × average susceptible proportion
      # Average S(a) = (1/(λΔa)) × [exp(-λa_lower) - exp(-λa_upper)]
        if (lambda > 1e-10) {
          avg_susceptible <- (1/(lambda * age_width)) * 
                            (exp(-lambda * a_lower) - exp(-lambda * a_upper))
        } else {
          # For very small lambda, S(a) ≈ 1
          avg_susceptible <- 1
        }
        
        # Annual infections = Population × FOI × Average susceptible
        infections_by_age[i, j] <- N_age_loc * lambda * avg_susceptible
      }
    }
      

  # Sum across age groups to get total annual infections per location
  annual_infections <- rowSums(infections_by_age)
  total_annual_infections <- sum(annual_infections)

  # Create dataframe with prediction coordinates
  infections_df <- data.frame(
    X_km = model$coop[, "X"],
    Y_km = model$coop[, "Y"],
    infections = annual_infections,
    population = pop_at_locations
  )
  
  # Convert to sf object (convert km back to meters for proper CRS)
  infections_sf <- st_as_sf(
    data.frame(X = infections_df$X_km * 1000, 
               Y = infections_df$Y_km * 1000),
    coords = c("X", "Y"),
    crs = crs
  )
  infections_sf$infections <- infections_df$infections
  infections_sf$population <- infections_df$population

  
  # Plot
  p <- ggplot() +
    geom_sf(data = infections_sf, aes(color = infections), 
            size = 1.7, alpha = 1, shape = 15) +
    scale_color_viridis_c(
      option = "plasma",
      name = "Annual\nInfections",
      trans = "log10",  # Log scale often better for infections
      labels = scales::comma_format()
    ) +
    labs(
      title = paste0("Predicted Annual Infections - ", pathogen_name),
      subtitle = paste0("Total: ", scales::comma(round(total_annual_infections))),
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
  
  # Calculate infections by age group (summed across all locations)
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
    infections_matrix = infections_by_age
  ))
}

# --- Seroprevalence by age group by year - model fits ---
plot_age_seroprevalence_model_fits <- function(year_intro, result,data, positive_col) {
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
