
get_onnv_prob <- function(fit,
                          onnv_components = c(2, 4),
                          method = c("mean", "draws"),
                          n_keep = 100,
                          seed = 1) {
 
  method <- match.arg(method)
 
  if (method == "mean") {
 
    s <- fit$summary("post_prob", mean = mean)
 
    idx <- stringr::str_match(s$variable, "post_prob\\[(\\d+),(\\d+)\\]")
    n_idx <- as.integer(idx[, 2])
    c_idx <- as.integer(idx[, 3])
 
    keep <- c_idx %in% onnv_components
    p <- tapply(s$mean[keep], n_idx[keep], sum)
    p <- as.numeric(p[order(as.integer(names(p)))])
 
  } else {
 
    set.seed(seed)
    d <- posterior::as_draws_matrix(fit$draws("post_prob"))
 
    if (nrow(d) > n_keep) d <- d[sort(sample(nrow(d), n_keep)), , drop = FALSE]
 
    idx <- stringr::str_match(colnames(d), "post_prob\\[(\\d+),(\\d+)\\]")
    n_idx <- as.integer(idx[, 2])
    c_idx <- as.integer(idx[, 3])
 
    keep <- c_idx %in% onnv_components
    n_ind <- max(n_idx)
 
    p <- matrix(0, nrow = nrow(d), ncol = n_ind)
    for (cc in onnv_components) {
      sel <- which(c_idx == cc)
      p[, n_idx[sel]] <- p[, n_idx[sel]] + d[, sel, drop = FALSE]
    }
    rm(d); gc()
  }
 
  # sanity checks
  pm <- if (is.matrix(p)) colMeans(p) else p
  stopifnot(all(pm >= 0 - 1e-8), all(pm <= 1 + 1e-8))
  message("Mean P(ONNV+) across individuals: ", round(mean(pm), 4),
          "  (should be close to your reported ONNV prevalence)")
 
  p
}



run_inla_imputed <- function(p_onnv,
                             data,
                             cam_pop,
                             age_groups,
                             w_age,
                             M = 20,
                             n_draws = 200,
                             year_intro = 1900,
                             positive_col = "ONNV_pos",
                             agg_factor = 10,
                             seed = 2024,
                             verbose = TRUE) {
 
  stopifnot(length(w_age) == nrow(age_groups))
 
  p_is_matrix <- is.matrix(p_onnv)
  n_ind <- if (p_is_matrix) ncol(p_onnv) else length(p_onnv)
  stopifnot(n_ind == nrow(data))
 
  # --- population on the prediction grid (identical across imputations) -------
  cam_pop_agg <- terra::aggregate(cam_pop, fact = agg_factor, fun = sum, na.rm = TRUE)
  pop <- as.data.frame(terra::as.points(cam_pop_agg, values = TRUE, na.rm = TRUE))[, 1]
  ok  <- !is.na(pop)
  pop <- pop[ok]
  w   <- pop / sum(pop)
  n_cells <- length(pop)
 
  # --- age-averaged susceptible proportion for a vector of cell-level FOIs ----
stopifnot(!any(is.na(age_groups$age_lower)),
        !any(is.na(age_groups$age_upper)),
        all(age_groups$age_upper > age_groups$age_lower),
        !any(is.na(w_age)))
  L <- age_groups$age_lower; U <- age_groups$age_upper; dW <- U - L
  suscept <- function(lam) {
    s <- sapply(seq_along(L), function(j)
      ifelse(lam > 1e-10,
             (1 / (lam * dW[j])) * (exp(-lam * L[j]) - exp(-lam * U[j])),
             1))
    as.numeric(s %*% w_age)
  }
 
  # --- accumulators -----------------------------------------------------------
  foi_draws       <- numeric(M * n_draws)
  seroprev_draws  <- numeric(M * n_draws)
  infect_draws    <- numeric(M * n_draws)
  imputation_id   <- integer(M * n_draws)
 
  cell_infect_sum <- numeric(n_cells)   # running sums -> per-cell means
  cell_foi_sum    <- numeric(n_cells)
  cell_prev_sum   <- numeric(n_cells)
 
  imputed_prevalence <- numeric(M)      # crude prevalence of each imputed label
  coop <- NULL
  k <- 0L
 
  t0 <- Sys.time()
 
  for (m in seq_len(M)) {
 
    set.seed(seed + m)
 
    # ---- draw a BINARY serostatus for each individual ----
    p_m <- if (p_is_matrix) p_onnv[((m - 1) %% nrow(p_onnv)) + 1, ] else p_onnv
    md <- data
    md[[positive_col]] <- rbinom(n_ind, size = 1, prob = p_m)
    imputed_prevalence[m] <- mean(md[[positive_col]])
 
    # ---- fit INLA to the imputed binary outcome ----
    res <- run_inla(
      year_intro   = year_intro,
      data         = md,
      cam_pop      = cam_pop,
      positive_col = positive_col
    )
 
    if (is.null(coop)) {
      coop <- res$coop
      stopifnot(nrow(coop) == length(ok))
      coop <- coop[ok, , drop = FALSE]
    }
 
    # ---- posterior samples of the latent field ----
    sel  <- list(Intercept = 1, spatial.field = seq_len(res$mesh$n))
    samp <- INLA::inla.posterior.sample(n_draws, res$output, selection = sel)
    Ap   <- INLA::inla.spde.make.A(mesh = res$mesh,
                                   loc = as.matrix(res$coop))[ok, , drop = FALSE]
 
    for (d in seq_len(n_draws)) {
      lat <- samp[[d]]$latent
      b0  <- lat[grep("^Intercept",     rownames(lat)), 1]
      fld <- lat[grep("^spatial.field", rownames(lat)), 1]
 
      lam  <- exp(b0 + as.numeric(Ap %*% fld))   # FOI per cell
      S    <- suscept(lam)                        # mean susceptible per cell
      infc <- pop * lam * S                       # annual infections per cell
 
      k <- k + 1L
      foi_draws[k]      <- sum(w * lam)
      seroprev_draws[k] <- sum(w * (1 - S))
      infect_draws[k]   <- sum(infc)
      imputation_id[k]  <- m
 
      cell_infect_sum <- cell_infect_sum + infc
      cell_foi_sum    <- cell_foi_sum + lam
      cell_prev_sum   <- cell_prev_sum + (1 - S)
    }
 
    rm(samp, res); gc()
 
    if (verbose) {
      el <- difftime(Sys.time(), t0, units = "mins")
      message(sprintf("imputation %d/%d done (%.1f min elapsed, ~%.1f min remaining)",
                      m, M, as.numeric(el), as.numeric(el) / m * (M - m)))
    }
  }
 
  n_total <- M * n_draws
 
  list(
    # per-draw national quantities (pooled across imputations)
    foi_draws      = foi_draws,
    seroprev_draws = seroprev_draws,
    infect_draws   = infect_draws,
    imputation_id  = imputation_id,
 
    # per-cell posterior means, aligned to coop / pop
    cell_infections = cell_infect_sum / n_total,
    cell_foi        = cell_foi_sum / n_total,
    cell_seroprev   = cell_prev_sum / n_total,
 
    coop = coop,            # km, UTM 33N
    pop  = pop,
    pop_covered = sum(pop),
 
    M = M, n_draws = n_draws,
    imputed_prevalence = imputed_prevalence
  )
}
 
 
# --- 3. Pool the imputations into national summaries --------------------------
#
# Pooling all M * n_draws samples gives the mixture posterior directly. The
# variance decomposition is reported alongside so you can quote how much of the
# total uncertainty comes from classification.
 
summarise_imputed <- function(res) {
 
  summ <- function(x) c(mean = mean(x),
                        ciL  = unname(quantile(x, 0.025)),
                        ciU  = unname(quantile(x, 0.975)))
 
  # variance decomposition: between-imputation vs within-imputation
  decomp <- function(x, id) {
    mi <- tapply(x, id, mean)
    vi <- tapply(x, id, var)
    c(between = var(mi), within = mean(vi),
      pct_classification = 100 * var(mi) / (var(mi) + mean(vi)))
  }
 
  out <- list(
    foi        = summ(res$foi_draws),
    seroprev   = summ(res$seroprev_draws),
    infections = summ(res$infect_draws),
    pop_covered = res$pop_covered,
    variance_decomposition = list(
      foi        = decomp(res$foi_draws, res$imputation_id),
      seroprev   = decomp(res$seroprev_draws, res$imputation_id),
      infections = decomp(res$infect_draws, res$imputation_id)
    ),
    imputed_prevalence = c(mean = mean(res$imputed_prevalence),
                           min  = min(res$imputed_prevalence),
                           max  = max(res$imputed_prevalence))
  )
 
  # CONSISTENCY GUARANTEE: the per-cell means must sum to the national mean
  stopifnot(isTRUE(all.equal(sum(res$cell_infections),
                             out$infections[["mean"]],
                             tolerance = 1e-6)))
 
  out
}
 
 
# --- 4. Per-cell surfaces as sf, for mapping and regional aggregation ---------
#
# Built from the SAME per-cell means used for the national total, so regional
# sums reconcile with the national figure exactly.
 
imputed_surfaces_sf <- function(res, crs = 32623) {
 
  sf::st_as_sf(
    data.frame(
      X              = res$coop[, "X"] * 1000,   # km -> m
      Y              = res$coop[, "Y"] * 1000,
      infections     = res$cell_infections,
      foi            = res$cell_foi,
      seroprevalence = res$cell_seroprev,
      population     = res$pop
    ),
    coords = c("X", "Y"),
    crs = 32633   # UTM Zone 33N
  )
}
 
 
# --- 5. Regional aggregation that does not silently drop cells ----------------
#
# Replaces aggregate_predictions_by_region(). The original used
# st_join(..., left = FALSE), an INNER join: any prediction cell falling outside
# every region polygon was discarded without warning, so regional sums did not
# reconcile with the national total. Here unmatched cells are assigned to the
# nearest region and the number affected is reported.
 
aggregate_by_region <- function(pred_sf, regions_sf,value_col,
                                region_col = "region",agg_type = c("weighted_mean", "sum"),
                                pop_col = "population") {
 
  agg_type <- match.arg(agg_type)
 
  if (sf::st_crs(pred_sf) != sf::st_crs(regions_sf)) {
    regions_sf <- sf::st_transform(regions_sf, sf::st_crs(pred_sf))
  }
 
  joined <- sf::st_join(pred_sf, regions_sf[, region_col], left = TRUE)
 
  n_miss <- sum(is.na(joined[[region_col]]))
  if (n_miss > 0) {
    pct <- if (agg_type == "sum") {
      100 * sum(joined[[value_col]][is.na(joined[[region_col]])], na.rm = TRUE) /
        sum(joined[[value_col]], na.rm = TRUE)
    } else NA_real_
    message(n_miss, " of ", nrow(joined),
            " cells fell outside all region polygons",
            if (!is.na(pct)) sprintf(" (%.2f%% of the total)", pct) else "",
            "; assigned to nearest region.")
    nr <- sf::st_nearest_feature(joined[is.na(joined[[region_col]]), ], regions_sf)
    joined[[region_col]][is.na(joined[[region_col]])] <- regions_sf[[region_col]][nr]
  }
  stopifnot(!any(is.na(joined[[region_col]])))
 
  if (agg_type == "sum") {
    out <- joined %>%
      sf::st_drop_geometry() %>%
      dplyr::group_by(.data[[region_col]]) %>%
      dplyr::summarise(total_value = sum(.data[[value_col]], na.rm = TRUE),
                       n_cells = dplyr::n(), .groups = "drop")
  } else {
    # population-weighted mean: zero-population cells carry no weight and are
    # dropped here, so n_cells is NOT comparable with the sum aggregation above.
    out <- joined %>%
      sf::st_drop_geometry() %>%
      dplyr::filter(!is.na(.data[[pop_col]]), .data[[pop_col]] > 0) %>%
      dplyr::group_by(.data[[region_col]]) %>%
      dplyr::summarise(
        value_weighted   = weighted.mean(.data[[value_col]], .data[[pop_col]], na.rm = TRUE),
        value_unweighted = mean(.data[[value_col]], na.rm = TRUE),
        total_population = sum(.data[[pop_col]], na.rm = TRUE),
        n_cells          = dplyr::n(), .groups = "drop")
  }
 
  out
}

