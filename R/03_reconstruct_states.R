#############################################################
#############################################################
###
### Reconstruction of sea surface temperature from a fitted
### MARSS data fusion model.
###
### The fit leaves its states in the scaled, anomaly space
### the model works in: a site level anomaly and a set of
### shared factors, with the seasonal cycle and the site
### means held separately as parameters. This step puts them
### back together and returns the result to degrees C.
###
### Two views of the same fit are available.
###
### The site view is the fused reconstruction: the site's own
### latent temperature, informed by every instrument that
### observes it. It is what the model is for.
###
### An instrument view is what the model says a single
### instrument sees at each site, using that instrument's own
### intercept, its own seasonal cycle, and -- when it has
### `site_state = FALSE` -- only the shared factors rather
### than the site's own anomaly. For MUR, an interpolated
### regional product, that is the regional signal with the
### site level anomaly removed, which is the natural thing to
### compare the fused reconstruction against.
###
### Both views are returned in long format, one row per site
### and date, so they can be joined to the observations or
### handed straight to ggplot.
###
### The standard errors here are state uncertainty at fixed
### parameter values: they say how well the smoother pins
### down the states given the estimates, not how well the
### estimates themselves are pinned down. See the note on
### `reconstruct_states()` for what that leaves out.
###
### Jack H. Buckner, Oregon State University, 08/31/2026
### Generated with Claude Code
### Reviewed JHB 09/02/2026
#############################################################
#############################################################

# Present in the other pipeline files as well; repeated so this one can be
# sourced on its own.
`%||%` <- function(x, y) if (is.null(x)) y else x

# What a config's `reconstruction` block leaves unset.
.RECONSTRUCTION_DEFAULTS <- list(
  enabled               = TRUE,
  instruments           = NULL,
  formats               = c("csv", "rds"),
  parameter_uncertainty = FALSE
)

.RECONSTRUCTION_FORMATS <- c("csv", "rds")

# The four objects `write_model_outputs()` writes, minus the fit itself.
.RUN_DIR_FILES <- c(scales = "scales.rds", observations = "observations.rds",
                    design = "design.rds")


#' Recover a fitted parameter block with its names attached
#'
#' The estimates live in `fit$par$<which>`, a one column matrix, but whether it
#' carries row names depends on the optimiser: TMB and BFGS return them, and the
#' EM algorithm does not. Since `kem` is this project's default method, reading
#' `rownames(fit$par$A)` directly works on some runs and fails on others.
#'
#' The names are always available as the column names of `fit$marss$free`, one
#' per estimated parameter and in the same order as the rows of `par`. Where
#' both are present they are identical -- verified across kem, BFGS and TMB fits
#' -- so taking them from `free` is simply the form that is always there.
#'
#' @param fit A fitted `marssMLE`.
#' @param which A parameter block: "A" or "Z".
#' @return A one column matrix of estimates with row names.
.par_named <- function(fit, which) {

  p <- fit$par[[which]]
  if (is.null(p)) p <- matrix(numeric(0), ncol = 1)

  # No free parameters at all is a legitimate model -- every intercept fixed,
  # say -- and the callers that allow it handle an empty set themselves.
  if (!nrow(p)) return(matrix(numeric(0), ncol = 1,
                              dimnames = list(character(0), NULL)))

  nms <- rownames(p)
  if (is.null(nms)) nms <- dimnames(fit$marss$free[[which]])[[2]]

  if (is.null(nms) || length(nms) != nrow(p))
    stop("Cannot recover parameter names for `", which, "` from this fit: it ",
         "has ", nrow(p), " estimate(s) and ",
         if (is.null(nms)) "no names" else paste(length(nms), "name(s)"),
         ".\nThis is not a fit produced by R/run_model.R.", call. = FALSE)

  matrix(p[, 1], ncol = 1, dimnames = list(nms, NULL))
}


#' Pull named parameters out of a fitted A
#'
#' The reconstruction locates parameters by name rather than by position, the
#' same names `build_observation_intercepts()` and `build_covariate_effects()`
#' put into the model. That is what makes it independent of the parameter
#' ordering MARSS happens to produce.
#'
#' Missing names are an error rather than a silent NA. Handing a fit from one
#' site set to a design describing another is the easy mistake here -- the name
#' lookup then matches nothing and every reconstructed value comes back NA,
#' which is indistinguishable from a legitimately unobserved series.
#'
#' `optional = TRUE` allows the whole set to be absent, which happens when no
#' instrument in the model uses a given parameterisation: with no instrument
#' carrying a shared site mean there is no `mu_<i>` to find, and zero is then
#' the right value rather than a failure. A partially present set is still an
#' error, because that can only mean the names were built wrongly.
#'
#' @param A A fitted `par$A` with its names attached, from `.par_named()`.
#' @param nms The parameter names to pull, in order.
#' @param what What the names are, for the error message.
#' @param optional Allow the entire set to be absent, returning zeros.
#' @return A numeric vector as long as `nms`.
.coef_by_name <- function(A, nms, what = "parameter", optional = FALSE) {

  absent <- setdiff(nms, rownames(A))
  if (!length(absent)) return(as.numeric(A[nms, 1]))

  if (optional && length(absent) == length(nms)) return(rep(0, length(nms)))

  stop("The fit has no ", what, " named: ", paste(absent, collapse = ", "),
       ".\nEither the fit and the design describe different site sets, or the ",
       "fit was produced by a different model specification.", call. = FALSE)
}


#' Recover the seasonal covariate matrix from a fit
#'
#' Taken from the fit rather than rebuilt from the dates, so the harmonics are
#' by construction the ones the coefficients were estimated against, aligned
#' column for column with the states. `seasonal_harmonics()` in
#' 01_format_data_for_marss.R would rebuild them, but only agrees if the period
#' and the number of harmonics are both carried across correctly.
#'
#' @param fit A fitted `marssMLE`.
#' @param n_harm Number of harmonic terms, from `design$n_harmonics`.
#' @param n_time Number of time steps the fit covers.
#' @return An `n_harm` by `n_time` numeric matrix.
.harmonics_from_fit <- function(fit, n_harm, n_time) {

  d <- fit$model$fixed$d
  if (is.null(d))
    stop("The fit has no seasonal covariates (`model$fixed$d`), so the ",
         "seasonal cycle cannot be reconstructed.", call. = FALSE)

  v <- as.numeric(d)
  if (length(v) %% n_harm != 0)
    stop("The fit's seasonal covariates hold ", length(v), " values, which is ",
         "not a multiple of the ", n_harm, " harmonic term(s) the design ",
         "declares.", call. = FALSE)

  h <- matrix(v, nrow = n_harm)

  # MARSS stores a time invariant covariate once rather than T times.
  if (ncol(h) == 1L && n_time > 1L)
    h <- matrix(h[, 1], nrow = n_harm, ncol = n_time)

  if (ncol(h) != n_time)
    stop("The fit's seasonal covariates span ", ncol(h), " time step(s) but ",
         "the states span ", n_time, ".", call. = FALSE)

  h
}


#' Rebuild the factor loading matrix from the fitted Z
#'
#' `dfa_loadings()` names the free entries `lambda_<i>_<j>`, so the fitted
#' values can be folded back into a matrix by name. The dimensions are forced
#' rather than inferred: with one factor, or with a fit whose last site loads
#' on nothing, the largest index seen is smaller than the design says.
#'
#' @param fit A fitted `marssMLE`.
#' @param n_site,m_factors The design dimensions.
#' @return An `n_site` by `m_factors` numeric matrix.
.loadings <- function(fit, n_site, m_factors) {

  if (!exists("tri_from_names"))
    stop("`tri_from_names()` not found. Source R/marss_matrix_functions.R ",
         "first.", call. = FALSE)

  tri_from_names(.par_named(fit, "Z"), nrows = n_site, ncols = m_factors)
}


#' Build the observation operator for one view
#'
#' Under "site_plus_factors" a view that sees the site state gets the identity
#' in the site block; one that does not gets zeros there, which is exactly the
#' structure `build_observation_matrix()` gives an instrument with
#' `site_state = FALSE`. Under "factors_only" there is no site block at all and
#' `site_state` has nothing to act on.
#'
#' @param design The design list.
#' @param Lambda The factor loadings.
#' @param site_state Whether this view sees the site's own anomaly.
#' @return A matrix with one row per site and one column per state.
.z_block <- function(design, Lambda, site_state) {

  if (!identical(design$state_structure, "site_plus_factors")) return(Lambda)

  I <- matrix(0, nrow = design$n_site, ncol = design$n_site)
  if (isTRUE(site_state)) diag(I) <- 1
  cbind(I, Lambda)
}


#' Reconstruct one view of the fit
#'
#' @param fit A fitted `marssMLE`.
#' @param design The design list.
#' @param scales The scaling constants.
#' @param harm The harmonic matrix from `.harmonics_from_fit()`.
#' @param Lambda The loadings from `.loadings()`.
#' @param intercept,seasonality The parameterisation kinds to use.
#' @param tag The instrument tag used to build names, or NULL for the shared
#'   site terms.
#' @param site_state Whether this view sees the site's own anomaly.
#' @param optional Passed to `.coef_by_name()`.
#' @return A list with `value` and `se`, both `n_site` by `n_time`, in degrees C.
.reconstruct_view <- function(fit, design, scales, harm, Lambda,
                              intercept, seasonality, tag = NULL,
                              site_state = TRUE, optional = FALSE) {

  n_site <- design$n_site
  n_harm <- nrow(harm)
  A      <- .par_named(fit, "A")
  label  <- if (is.null(tag)) "site" else tag

  # ---- seasonal cycle ------------------------------------------------------
  # Mirrors build_covariate_effects(): "shared" uses the site's own terms,
  # "independent" a set of its own carrying the instrument tag.
  infix <- switch(seasonality,
    "shared"      = "",
    "independent" = paste0(tag, "_"),
    stop("Unknown seasonality kind for the ", label, " view: ", seasonality,
         call. = FALSE))

  C <- matrix(0, nrow = n_site, ncol = n_harm)
  for (k in seq_len(n_harm))
    C[, k] <- .coef_by_name(A, paste0("c", k, "_", infix, seq_len(n_site)),
                            paste0("seasonal coefficient(s) for the ", label,
                                   " view"), optional = optional)

  trends <- C %*% harm

  # ---- intercept -----------------------------------------------------------
  # Mirrors build_observation_intercepts(). MARSS reads the "+" of a
  # "site+instrument" name as a sum of two estimated parameters, so the fit
  # carries the site means and the instrument offset separately.
  mu <- switch(intercept,
    "zero"        = rep(0, n_site),
    "site"        = .coef_by_name(A, paste0("mu_", seq_len(n_site)),
                                  paste0("site mean(s) for the ", label,
                                         " view"), optional = optional),
    "independent" = .coef_by_name(A, paste0("mu_", tag, "_", seq_len(n_site)),
                                  paste0("site mean(s) for the ", label,
                                         " view")),
    "site+instrument" =
      .coef_by_name(A, paste0("mu_", seq_len(n_site)),
                    paste0("site mean(s) for the ", label, " view"),
                    optional = optional) +
      .coef_by_name(A, paste0("mu_", tag),
                    paste0("instrument offset for the ", label, " view")),
    stop("Unknown intercept kind for the ", label, " view: ", intercept,
         call. = FALSE))

  # ---- states --------------------------------------------------------------
  Z <- .z_block(design, Lambda, site_state)
  if (ncol(Z) != nrow(fit$states))
    stop("The design implies ", ncol(Z), " state(s) but the fit has ",
         nrow(fit$states), ". The fit and the design do not describe the same ",
         "model.", call. = FALSE)

  x <- Z %*% fit$states

  # Propagated as if the states were independent, which they are not: this
  # ignores the smoother's covariance between a site state and the shared
  # factors. See the note on reconstruct_states().
  x_se <- if (is.null(fit$states.se)) {
    warning("The fit carries no `states.se`; standard errors are returned as ",
            "NA.", call. = FALSE)
    array(NA_real_, dim = dim(x))
  } else {
    sqrt(Z^2 %*% fit$states.se^2)
  }

  # `scales$sigma`, `scales$mu` and `mu` are all length n_site and recycle down
  # the rows of an n_site by n_time matrix, which is why the row count has to
  # be exactly n_site rather than merely a multiple of it.
  list(value = scales$sigma * (x + trends + mu) + scales$mu,
       se    = scales$sigma * x_se)
}


#' Melt a reconstructed view into long format
#'
#' Site major, so the rows for one site are contiguous and in date order --
#' the layout `as.vector(t(.))` gives for free and the one a time series plot
#' wants.
#'
#' @param view Output of `.reconstruct_view()`.
#' @param sites Site keys, in state order.
#' @param dates The dates the columns span.
#' @param variable Value for the `variable` column.
#' @param file_tag Suffix `write_reconstruction()` gives this frame's file.
#' @return A data frame with date, site, variable, value and se.
.to_long <- function(view, sites, dates, variable, file_tag = "") {

  out <- data.frame(
    date     = rep(dates, times = length(sites)),
    site     = rep(sites, each = length(dates)),
    variable = variable,
    value    = as.vector(t(view$value)),
    se       = as.vector(t(view$se)),
    stringsAsFactors = FALSE
  )
  attr(out, "file_tag") <- file_tag
  out
}


#' Reconstruct sea surface temperature from a fitted model
#'
#' Returns the fused site level reconstruction, and optionally one frame per
#' named instrument giving what the model says that instrument sees at each
#' site. Every frame is long format with columns `date`, `site`, `variable`,
#' `value` and `se`, `value` and `se` in degrees C.
#'
#' Each view is assembled as
#'
#'   value = sigma * (Z x + C d + mu) + m
#'
#' where `x` are the smoothed states, `Z` selects the site state (when the view
#' sees it) and applies the factor loadings, `C d` is the seasonal cycle, `mu`
#' the intercept, and `sigma`, `m` the per site scaling constants. Which
#' parameters `Z`, `C` and `mu` are built from is what distinguishes the views,
#' and comes from the instrument specification in the design rather than from
#' anything hard coded here.
#'
#' **What the standard errors are.** They propagate the Kalman smoother's state
#' uncertainty through `Z` and the scaling, holding the parameters at their
#' estimates. Two things are therefore left out. Parameter uncertainty is not
#' included at all, so the intervals are narrower than a full posterior; that is
#' deliberate and is what `parameter_uncertainty = FALSE` in the config records.
#' And the propagation treats the states as independent, using only
#' `states.se`, so it ignores the smoother's covariance between a site state and
#' the shared factors. The full covariance is available as `fit$kf$VtT` if this
#' is ever tightened up.
#'
#' @param fit A fitted `marssMLE`, as written by `fit_marss()`.
#' @param scales The `scales.rds` contents: `mu`, `sigma` and `sites`.
#' @param observations The `observations.rds` contents; `col_dates` supplies the
#'   time axis and `site_levels` the site keys if `scales$sites` is absent.
#' @param design The `design.rds` contents: dimensions, state structure, number
#'   of harmonics and the instrument specification.
#' @param instruments Variable names to reconstruct an instrument view for, or
#'   NULL for none. Each must have an entry in `design$instruments`.
#' @param verbose Print what was reconstructed.
#'
#' @return A named list of data frames: `site` always, plus one element per
#'   entry of `instruments`, named by the variable.
#' @export
reconstruct_states <- function(fit, scales, observations, design,
                               instruments = NULL, verbose = TRUE) {

  # ---- what the design has to tell us --------------------------------------
  need <- c("n_site", "m_factors", "state_structure", "n_harmonics")
  absent <- need[vapply(need, function(k) is.null(design[[k]]), logical(1))]
  if (length(absent))
    stop("`design` is missing: ", paste(absent, collapse = ", "),
         ".\nRun directories written before the reconstruction step exists ",
         "carry a shorter design; supply the missing fields explicitly.",
         call. = FALSE)

  n_site <- design$n_site

  # ---- the time axis and the site keys -------------------------------------
  dates <- as.Date(observations$col_dates)
  if (length(dates) != ncol(fit$states))
    stop("The observation matrix spans ", length(dates), " column(s) but the ",
         "fit has ", ncol(fit$states), " state(s) in time. The fit and the ",
         "observations do not come from the same run.", call. = FALSE)

  sites <- scales$sites %||% observations$site_levels
  if (is.null(sites))
    stop("Neither `scales$sites` nor `observations$site_levels` is present, ",
         "so states cannot be labelled with a site.", call. = FALSE)
  if (length(sites) != n_site)
    stop("The design declares ", n_site, " site(s) but ", length(sites),
         " site key(s) are available.", call. = FALSE)

  if (length(scales$mu) != n_site || length(scales$sigma) != n_site)
    stop("`scales$mu` and `scales$sigma` must both be as long as the number ",
         "of sites (", n_site, ").", call. = FALSE)

  # ---- the pieces every view shares ----------------------------------------
  harm   <- .harmonics_from_fit(fit, design$n_harmonics, length(dates))
  Lambda <- .loadings(fit, n_site, design$m_factors)

  out <- list(site = .to_long(
    .reconstruct_view(fit, design, scales, harm, Lambda,
                      intercept = "site", seasonality = "shared",
                      tag = NULL, site_state = TRUE, optional = TRUE),
    sites, dates, "sst", file_tag = ""))

  # ---- one frame per requested instrument ----------------------------------
  ins <- design$instruments
  for (v in instruments) {

    spec <- ins[[v]]
    if (is.null(spec))
      stop("No instrument specification for '", v, "'. The design describes: ",
           paste(names(ins), collapse = ", "), call. = FALSE)

    out[[v]] <- .to_long(
      .reconstruct_view(fit, design, scales, harm, Lambda,
                        intercept   = spec$intercept,
                        seasonality = spec$seasonality,
                        tag         = spec$tag,
                        site_state  = isTRUE(spec$site_state)),
      sites, dates, v, file_tag = paste0("_", spec$tag, "_only"))
  }

  if (verbose)
    message("Reconstructed ", n_site, " site(s) x ", length(dates),
            " date(s): ", paste(names(out), collapse = ", "), ".")

  out
}


#' Write reconstructed states to a run directory
#'
#' The site view goes to `states.<ext>`; an instrument view goes to
#' `states_<tag>_only.<ext>`, so MUR lands as `states_mur_only.csv`.
#'
#' @param states Output of `reconstruct_states()`.
#' @param outdir Directory to write into.
#' @param formats Any of "csv" and "rds".
#' @param verbose Print what was written.
#' @return The paths written, invisibly.
#' @export
write_reconstruction <- function(states, outdir, formats = c("csv", "rds"),
                                 verbose = TRUE) {

  bad <- setdiff(formats, .RECONSTRUCTION_FORMATS)
  if (length(bad))
    stop("Unknown output format(s): ", paste(bad, collapse = ", "),
         ". Use any of: ", paste(.RECONSTRUCTION_FORMATS, collapse = ", "),
         call. = FALSE)

  # Normally the run directory the fit was written to, which already exists;
  # created here for the case where --outdir names somewhere else.
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  paths <- character()
  for (nm in names(states)) {
    df   <- states[[nm]]
    stem <- paste0("states", attr(df, "file_tag") %||% "")
    for (f in formats) {
      p <- file.path(outdir, paste0(stem, ".", f))
      if (identical(f, "csv")) utils::write.csv(df, p, row.names = FALSE)
      else                     saveRDS(df, p)
      paths <- c(paths, p)
    }
  }

  if (verbose)
    message("Wrote ", paste(basename(paths), collapse = ", "), ".")

  invisible(paths)
}


#' Resolve the reconstruction settings for a run
#'
#' Taken from the config when one is supplied, otherwise from the copy stored in
#' the design at fitting time, otherwise from the defaults. The design copy is
#' what lets a run directory be reconstructed without its config.
#'
#' @param config The `model_config` list, or NULL.
#' @param design The design list, or NULL.
#' @return A list with enabled, instruments, formats and parameter_uncertainty.
#' @export
reconstruction_settings <- function(config = NULL, design = NULL) {

  rc <- config$reconstruction %||% design$reconstruction %||% list()
  if (!is.list(rc))
    stop("`reconstruction` must be a list.", call. = FALSE)

  rc <- utils::modifyList(.RECONSTRUCTION_DEFAULTS, rc)

  # The seam for parameter uncertainty is named but not implemented: the states
  # are reconstructed at the parameter estimates. Refusing TRUE is what keeps a
  # config that asks for it from silently getting something narrower.
  if (!identical(rc$parameter_uncertainty, FALSE))
    stop("`reconstruction$parameter_uncertainty` must be FALSE: the states are ",
         "reconstructed at the parameter estimates, and no method for ",
         "propagating parameter uncertainty is implemented yet.", call. = FALSE)

  rc
}


#' Read the objects a reconstruction needs out of a run directory
#'
#' The fit defaults to the last chunk written, which is the fit the run
#' finished on; a run that stopped after its initialisation stage has only
#' `fit_init.rds`, and that is used instead.
#'
#' A run directory written before the design carried the instrument
#' specification falls back to the config copied alongside it, so earlier runs
#' can still be reconstructed.
#'
#' @param dir A model run directory.
#' @param fit_file A file name within `dir`, or NULL to choose the last chunk.
#' @param load_fit Read the fit at all. FALSE when the caller already holds it,
#'   which avoids reading a large object back off disk at the end of a run.
#' @return A list with `scales`, `observations`, `design`, and `fit` and
#'   `fit_path` when `load_fit` is TRUE.
#' @export
read_run_dir <- function(dir, fit_file = NULL, load_fit = TRUE) {

  if (!dir.exists(dir))
    stop("Run directory not found: ", dir, call. = FALSE)

  out <- list()
  for (k in names(.RUN_DIR_FILES)) {
    p <- file.path(dir, .RUN_DIR_FILES[[k]])
    if (!file.exists(p))
      stop("Run directory is missing ", .RUN_DIR_FILES[[k]], ": ", dir,
           "\nThis is not a run directory written by R/run_model.R.",
           call. = FALSE)
    out[[k]] <- readRDS(p)
  }

  if (load_fit) {
    if (is.null(fit_file)) {
      chunks <- sort(list.files(dir, pattern = "^fit_chunk_[0-9]+\\.rds$"))
      fit_file <- if (length(chunks)) chunks[length(chunks)] else "fit_init.rds"
    }
    fit_path <- file.path(dir, fit_file)
    if (!file.exists(fit_path))
      stop("Fit not found: ", fit_path, call. = FALSE)
    out$fit      <- readRDS(fit_path)
    out$fit_path <- fit_path
  }

  # ---- fill a pre-reconstruction design from the config it was fit with ----
  if (is.null(out$design$instruments)) {
    cfg_path <- file.path(dir, "config_used.R")
    if (file.exists(cfg_path) && exists("read_model_config")) {
      cfg <- read_model_config(cfg_path)
      message("The design carries no instrument specification; taking it from ",
              "config_used.R. Command line overrides used by that run are not ",
              "reflected there.")
      out$design$instruments <- out$design$instruments %||%
        cfg$structure$instruments
      out$design$state_structure <- out$design$state_structure %||%
        cfg$structure$state_structure
      out$design$m_factors <- out$design$m_factors %||% cfg$structure$m_factors
      out$design$reconstruction <- out$design$reconstruction %||%
        cfg$reconstruction
    }
  }

  out
}


#' Reconstruct a run's states and write them alongside the fit
#'
#' The single entry point both `run_model.R` and `run_reconstruct.R` call. The
#' fit may be handed in directly, which is what happens at the end of a fitting
#' run, or read from the run directory, which is what happens when a run is
#' reconstructed after the fact.
#'
#' @param dir A model run directory holding scales.rds, observations.rds and
#'   design.rds.
#' @param fit An already loaded `marssMLE`, or NULL to read one from `dir`.
#' @param fit_file A file name within `dir`, used only when `fit` is NULL.
#' @param config The `model_config` list, or NULL to use the settings stored in
#'   the design.
#' @param instruments,formats,enabled Override the corresponding resolved
#'   setting. `instruments` uses NA for "not supplied", since NULL is the
#'   meaningful value "the site view alone".
#' @param outdir Where to write; defaults to `dir`.
#' @param verbose Print progress.
#'
#' @return The paths written, invisibly, or NULL if reconstruction is disabled.
#' @export
run_reconstruction <- function(dir, fit = NULL, fit_file = NULL, config = NULL,
                               instruments = NA, formats = NULL, enabled = NULL,
                               outdir = NULL, verbose = TRUE) {

  run <- read_run_dir(dir, fit_file = fit_file, load_fit = is.null(fit))
  if (!is.null(fit)) run$fit <- fit

  rc <- reconstruction_settings(config, run$design)
  if (!identical(instruments, NA)) rc["instruments"] <- list(instruments)
  if (!is.null(formats))           rc$formats        <- formats
  if (!is.null(enabled))           rc$enabled        <- enabled

  if (!isTRUE(rc$enabled)) {
    if (verbose) message("Reconstruction is disabled; nothing written.")
    return(invisible(NULL))
  }

  states <- reconstruct_states(run$fit, run$scales, run$observations,
                               run$design, instruments = rc$instruments,
                               verbose = verbose)

  write_reconstruction(states, outdir %||% dir, formats = rc$formats,
                       verbose = verbose)
}
