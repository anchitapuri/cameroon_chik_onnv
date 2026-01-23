# ---- Function for Year of Intro INLA models ----
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
  grid <- st_make_grid(cameroon, cellsize = 0.05, what = "centers")
  grid_sf   <- st_sf(grid_id = seq_along(grid), geometry = grid)
  grid_inside <- st_intersection(grid_sf, cameroon)
  grid_inside <- grid_inside[order(grid_inside$grid_id), ]   # restore original order
  
  # Project to appropriate CRS (UTM Zone 33N for Cameroon)
  grid_utm <- st_transform(grid_inside, crs = 32633)
  coop_utm <- st_coordinates(grid_utm)
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
    #data = list(y = 1 * data_points[[positive_col]]),  
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

# ---- INLA model function with covariates ----
run_inla_model_comparision <- function(year_intro, data, cameroon, anopheles_funestus, anopheles_gambiae, 
                                       positive_col, covariates = "baseline") {
  
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
  coop        <- coop_utm / 1000
  colnames(coop) <- c("X", "Y")

  
  # --- Covariates ---
  # extract mosquito data at prediction points 
  funestus_values <- extract(anopheles_funestus, vect(grid_inside), method = "bilinear")
  gambiae_values <- extract(anopheles_gambiae, vect(grid_inside), method = "bilinear")

  # create coop dataframe and add covariates 
  dp <- as.data.frame(coop)
  dp$funestus <- funestus_values$`2010_Anopheles_funestus_CMR`
  dp$gambiae  <- gambiae_values$`2010_Anopheles_gambiae_ss_CMR`
  dp$pop     <- pop_values$cmr_ppp_2020_UNadj
  
  
  # Build SPDE model
  spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
  s.index <- inla.spde.make.index("spatial.field", spde$n.spde)
  
  # Projection matrices
  A <- inla.spde.make.A(mesh = mesh, loc = cooe)
  Ap <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coop))

  
  # Prepare covariates based on input
  if (covariates == "baseline") {
    # Model 1: Baseline (intercept + offset(log(age)) + spatial)
    est_effects <- data.frame(Intercept = 1, 
                              age = data_points$years_of_exposure)
    pred_effects <- data.frame(Intercept = 1, 
                               age = rep(1, nrow(coop)))
    formula <- y ~ -1 + Intercept + offset(log(age)) + f(spatial.field, model = spde)
    
  } else if (covariates == "anopheles_funestus") {
    # Model 2: + funestus
    est_effects <- data.frame(Intercept = 1, 
                              age = data_points$years_of_exposure,
                              funestus = data_points$fun_pw_district)
    pred_effects <- data.frame(Intercept = 1, 
                               age = rep(1, nrow(coop)),
                               funestus = dp$funestus)

    formula <- y ~ -1 + Intercept + offset(log(age)) + funestus + f(spatial.field, model = spde)
    
  } else if (covariates == "anopheles_gambiae") {
    # Model 3: + gambiae
    est_effects <- data.frame(Intercept = 1, 
                              age = data_points$years_of_exposure,
                              gambiae = data_points$gam_pw_district)
    pred_effects <- data.frame(Intercept = 1, 
                               age = rep(1, nrow(coop)),
                               gambiae = dp$gambiae)
    
    formula <- y ~ -1 + Intercept + offset(log(age)) + gambiae + f(spatial.field, model = spde)
    
  } else if (covariates == "anopheles_both") {
    # Model 4: + both anopheles
    est_effects <- data.frame(Intercept = 1, 
                              age = data_points$years_of_exposure,
                              funestus = data_points$fun_pw_district,
                              gambiae = data_points$gam_pw_district)
    pred_effects <- data.frame(Intercept = 1, 
                               age = rep(1, nrow(coop)),
                               funestus = dp$funestus,
                               gambiae = dp$gambiae)
    
    formula <- y ~ -1 + Intercept + offset(log(age)) + funestus + gambiae + f(spatial.field, model = spde)
    
  } 
  
  # Estimation stack
  stk.e <- inla.stack(
    data = list(y = data_points[[positive_col]]),  
    A = list(1, A),
    effects = list(est_effects, spatial.field = s.index),
    tag = "est")
  
  # Prediction stack
  stk.p <- inla.stack(
    tag = "pred", 
    data = list(y = rep(NA, nrow(coop))),
    A = list(1, Ap),
    effects = list(pred_effects, spatial.field = s.index))
  
  # Full stack
  stk.full <- inla.stack(stk.e, stk.p)
  
  # Run INLA model
  output <- inla(formula,
                 data = inla.stack.data(stk.full),
                 family = "binomial",
                 Ntrials = 1,
                 control.family = list(link = "cloglog"),
                 control.predictor = list(A = inla.stack.A(stk.full), 
                                          compute = TRUE, 
                                          link = 1),
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                 verbose = FALSE)
  
  return(list(
    year = year_intro,
    model_type = covariates,
    output = output,
    dic = output$dic$dic,
    waic = output$waic$waic,
    mesh = mesh,
    stk.full = stk.full,
    data_filtered = data_points,
    cooe = cooe,
    coop = coop
  ))
}


# --- seroprevelance by age group by year obs data ---
plot_age_seroprevalence_by_year <- function(data, positive_col) {
  
  # Recreate the filtered dataset used in the model
  data_plot <- data
  
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
  
  # Create a formula dynamically using the positive_col parameter
  formula_mean <- as.formula(paste(positive_col, "~ year_of_survey + age_group"))
  formula_length <- as.formula(paste(positive_col, "~ year_of_survey + age_group"))
  
  # Summaries by year and age group
  obs <- aggregate(formula_mean, data_plot, mean, na.rm = TRUE)
  n_by <- aggregate(formula_length, data_plot, length)
  names(n_by)[3] <- "n"
  names(obs)[3] <- "proportion_positive"  # Rename for clarity
  
  obs <- merge(obs, n_by, by = c("year_of_survey", "age_group"))
  obs$obs_lower <- pmax(0, obs$proportion_positive - 1.96*sqrt(obs$proportion_positive*(1-obs$proportion_positive)/obs$n))
  obs$obs_upper <- pmin(1, obs$proportion_positive + 1.96*sqrt(obs$proportion_positive*(1-obs$proportion_positive)/obs$n))
  
  y_limits <- if (positive_col == "CHIK_pos") {
    c(0, 0.08)
  } else {
    c(0, 0.8)
  }
  
  # Plot
  p <- ggplot() +
    geom_point(data = obs, aes(x = age_group, y = proportion_positive), size = 2, color = '#0d1b2a') +
    geom_errorbar(data = obs, aes(x = age_group, ymin = obs_lower, ymax = obs_upper), width = 0.15, color = '#0d1b2a') +
    
    facet_wrap(~ year_of_survey, ncol = 5) +

    scale_y_continuous(limits = y_limits) +
    labs(
      x = "Age group",
      y = "Proportion seropositive",
      title = paste("Observed seroprevalence by age group -", positive_col),
    ) +
    theme_bw() +
    theme(axis.text.x = element_text())
  
  print(p)
  invisible(p)
  
  return(data_plot)
}


plot_age_seroprevalence_by_year_by_gender  <- function(data, positive_col) {
  
  # Recreate the filtered dataset used in the model
  data_plot <- data
  
  # Filter to only Sex = 1 (Male) or 2 (Female)
  data_plot <- data_plot[data_plot$Sex %in% c(1, 2), ]
  
  # Create sex labels (assuming 1 = Male, 2 = Female)
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
    c(0, 0.08)
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


# --- seroprevalence by age group by year - model fits ---
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



