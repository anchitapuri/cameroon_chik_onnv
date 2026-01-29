
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
  index_pred <- inla.stack.index(best_model$stk.full, tag = "pred")$data
  
  # Extract the intercept
  eta_pred <- best_model$output$summary.linear.predictor[index_pred, "mean"]
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


# --- Function to extract and plot seroprevelance ---
plot_predicted_seroprevalence <- function(foi_result,age_mid, age_weights,
                                          crs = 32633, pathogen_name = "ONNV") {
  
  # Calculate probability matrix
  prob_mat <- outer(
    foi_result$foi_sf$foi,
    age_mid,
    function(lambda, a) 1 - exp(-lambda * a)
  )
  
  # Age weighted prevalence at each location
  prev_loc <- as.vector(prob_mat %*% age_weights)
  
  # Create dataframe with prediction coordinates
  prev_df <- data.frame(
    X_km = foi_result$coop[, "X"],
    Y_km = foi_result$coop[, "Y"],
    prev = prev_loc
  )
  
  # Convert to sf object (convert km back to meters for proper CRS)
  prev_sf <- st_as_sf(
    data.frame(X = prev_df$X_km * 1000, Y = prev_df$Y_km * 1000),
    coords = c("X", "Y"),
    crs = 32633
  )
  
  # Add prevalence values to sf object
  prev_sf$prev <- prev_df$prev
  
  # Create subtitle if not provided
  if (is.null(subtitle)) {
    subtitle <- paste("Introduction year:", foi_result$year)
  }
  
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
      title = paste0("Predicted Seroprevalence -", pathogen_name),
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
  
  # Return both the plot and the prevalence data
  return(list(
    plot = p,
    prev_sf = prev_sf,
    prev_range = range(prev_loc)
  ))
}


# --- Function to extract and plot annual infections ---



# --- Seroprevalence by age group by year - model fits ---
plot_age_seroprevalence_model_fits <- function(year_intro,
                                               result,
                                               data,
                                               positive_col) {
  
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
