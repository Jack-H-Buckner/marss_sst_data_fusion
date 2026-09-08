### ------------------------------------------------------
### UPDATED 2026-09-04, NEEDS RE-REVIEW: validation holdout
### (`data$sites_validation`). Newer than your last review of
### this file. Delete this block once you have re-read it.
### ------------------------------------------------------
#############################################################
#############################################################
###
### Specification and fitting of the MARSS data fusion model.
###
### Turns the observation matrix written by the preprocessing
### step, plus a configuration from parameters/model.R, into
### the list of fixed and estimated matrices that MARSS()
### takes, and runs the fit.
###
### Every hard coded choice in the earlier exploratory
### scripts -- which sites, which dates, which instruments
### share a mean or a seasonal cycle, how the states are
### structured -- is a configuration value here, so a fit is
### fully described by its config file.
###
### The observation equation is
###
###   y_t = Z x_t + a + D d_t + v_t,   v_t ~ N(0, R)
###
### and the state equation is
###
###   x_t = B x_{t-1} + w_t,           w_t ~ N(0, Q)
###
### where d_t holds the seasonal harmonics.
###
### Jack H. Buckner, Oregon State University, 08/30/2026
### Generated with Claude Code
### Reviewed JHB 09/02/2026
#############################################################
#############################################################

library(MARSS)

# Present in run_preprocessing.R as well; repeated so this file can be sourced
# on its own.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Internal: rank the unique values of x, so a subset of sites is renumbered
# 1..n without gaps. Defined here rather than taken from dplyr, which is not
# otherwise needed by the model step.
.dense_rank <- function(x) {
  if (all(is.na(x))) return(as.integer(x))
  match(x, sort(unique(x)))
}

# Internal: the rows a holdout withheld, as a long data frame. Keyed the same
# way as the reconstruction in 03_reconstruct_states.R -- date, site, variable
# -- so the two join on (site, date) whenever the comparison gets written.
# Values come straight from `ts_matrix`, i.e. degrees C, never scaled.
.withheld_long <- function(dat, rows, keep_cols) {

  empty <- data.frame(date = as.Date(character(0)), site = character(0),
                      variable = character(0), value = numeric(0),
                      stringsAsFactors = FALSE)
  if (!length(rows)) return(empty)

  block <- dat$ts_matrix[rows, keep_cols, drop = FALSE]
  dates <- as.Date(dat$col_dates)[keep_cols]

  out <- data.frame(
    date     = rep(dates, times = length(rows)),
    site     = rep(dat$row_site_keys[rows], each = length(dates)),
    variable = rep(dat$row_var_keys[rows],  each = length(dates)),
    value    = as.vector(t(block)),
    stringsAsFactors = FALSE)

  # Empty cells of the padded grid are gaps in the record, not withheld data.
  out <- out[!is.na(out$value), , drop = FALSE]
  out <- out[order(out$variable, out$site, out$date), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# The intercept and seasonality parameterisations a config may name.
.INTERCEPT_KINDS  <- c("site", "site+instrument", "independent", "zero")
.SEASONALITY_KINDS <- c("shared", "independent")
.STATE_STRUCTURES  <- c("site_plus_factors", "factors_only")

# What `data$sites_validation` holds out when the config does not say. In situ
# is the series a reconstruction is normally validated against, so it is the
# default; `data$validation_variables` names something else.
.VALIDATION_VARIABLES <- "insitu_sst"

# Fitting methods MARSS accepts. "kem" is the EM algorithm; everything else
# optimises the likelihood directly through optim() or nlminb(), which imposes
# a constraint on R -- see `.needs_diagonal_R()`.
.METHODS <- c("kem", "BFGS", "TMB", "BFGS_TMB", "nlminb_TMB")

# Only EM can fit a block structured R. The direct optimisers require every
# block of a variance-covariance matrix to be fixed, diagonal, or wholly
# unconstrained; a diagonal of per instrument variances with a single shared
# off-diagonal covariance is none of those.
.needs_diagonal_R <- function(method) !identical(method, "kem")

# Instruments whose day effect puts a shared covariance into the off-diagonal
# of R, which is what the direct optimisers cannot handle.
.day_effect_vars <- function(config) {
  ins <- config$structure$instruments
  if (!is.list(ins)) return(character(0))
  names(ins)[vapply(ins, function(s) isTRUE(s$day_effect), logical(1))]
}

# TMB does not estimate B: it returns whatever value B started at, without
# saying so. Holding B at a stated value is therefore a modelling choice for
# the initialisation stage, not a tuning knob, and `B_values` makes it explicit
# rather than leaving it to MARSS's default of 1 (a random walk).
#
# The free entries of B depend on the state structure, so the names are listed
# here once and validation checks `B_values` against them without having to
# build the model.
.b_param_names <- function(config) {
  if (identical(config$structure$state_structure, "site_plus_factors"))
    c("rho_chi", "rho_eta") else "rho"
}

# Turn a scalar or a named list of B values into one named numeric vector
# covering every free entry. Errors naming whatever is missing.
.resolve_b_values <- function(names, B_values) {

  if (is.numeric(B_values) && length(B_values) == 1L && is.null(names(B_values)))
    return(stats::setNames(rep(as.numeric(B_values), length(names)), names))

  if (!is.list(B_values) && !(is.numeric(B_values) && !is.null(names(B_values))))
    stop("`init$B_values` must be a single number or a named list.",
         call. = FALSE)

  vals   <- unlist(B_values)
  absent <- setdiff(names, names(vals))
  if (length(absent))
    stop("`init$B_values` has no value for: ", paste(absent, collapse = ", "),
         ". The free entries of B for this state structure are: ",
         paste(names, collapse = ", "), call. = FALSE)

  stats::setNames(as.numeric(vals[names]), names)
}

# Derive the stage 1 configuration from the full one: the initialisation method
# and its controls, a single pass, no day effects (the direct optimisers cannot
# fit the off-diagonal they add), and B held at the configured values.
#
# Dropping the day effects here is safe in a way that dropping them from the
# final fit would not be: stage 1 exists only to supply starting values, and the
# parameters it cannot estimate are simply not transferred.
.init_config <- function(config) {

  init <- config$fitting$init
  cfg  <- config

  cfg$fitting$method                <- init$method
  cfg$fitting$controls[[init$method]] <- init$controls
  cfg$fitting$chunks                <- 1L

  for (v in names(cfg$structure$instruments))
    cfg$structure$instruments[[v]]$day_effect <- FALSE

  cfg$structure$B_fixed <-
    .resolve_b_values(.b_param_names(config), init$B_values)

  cfg
}

# Locate a saved fit to start from. `path` is either a run directory or the fit
# itself; a directory resolves to its last chunk or to fit_init.rds, depending
# on which stage the caller is replacing.
#
# read_run_dir() in 03_reconstruct_states.R carries the same chunk selection.
# That duplication is deliberate: 03 has no hard dependency on this file -- its
# one use of read_model_config() is guarded by exists() -- so that it can be
# sourced on its own, and sharing four lines is not worth a source order
# dependency between them.
.resolve_fit_path <- function(path, prefer = c("chunk", "init")) {

  prefer <- match.arg(prefer)

  if (!dir.exists(path)) {
    if (!file.exists(path))
      stop("No fit found at: ", path, call. = FALSE)
    return(path)
  }

  name <- if (identical(prefer, "init")) "fit_init.rds" else {
    chunks <- sort(list.files(path, pattern = "^fit_chunk_[0-9]+\\.rds$"))
    if (length(chunks)) chunks[length(chunks)] else "fit_init.rds"
  }

  out <- file.path(path, name)
  if (!file.exists(out))
    stop("Run directory holds no ", name, ": ", path,
         "\nName the fit directly if it is somewhere else.", call. = FALSE)
  out
}

# Refuse to start from a fit of different data.
#
# transfer_params() matches on parameter names, and the names are a property of
# the model rather than of the data: a fit of a different site set, date range or
# scaling produces exactly the same names, transfers cleanly, and is silently
# wrong. Nothing downstream would notice, so it is caught here.
#
# The comparison is on values rather than row names because MARSS de-duplicates
# the site keys build_model_data() puts on y_fit into CB001-1, CB008-1,
# CB001-2 ...; undoing that needs a suffix strip that a site key ending in a
# digit would break. Comparing the matrices catches a changed site set, date
# range, instrument set or scaling in one check.
.check_same_data <- function(loaded, y, what) {

  d_old <- dim(loaded$model$data)
  d_new <- dim(y)

  if (!identical(d_old, d_new))
    stop(what, " was fit to a ", d_old[1], " x ", d_old[2],
         " observation matrix, but this run's is ", d_new[1], " x ", d_new[2],
         ".\nThe two describe different data, so its estimates cannot be ",
         "carried across. Match `data$sites` and the date range to the run ",
         "you are starting from.", call. = FALSE)

  same <- isTRUE(all.equal(unname(loaded$model$data), unname(y)))
  if (!same)
    stop(what, " was fit to observations that differ from this run's, though ",
         "they are the same shape.\nThe site or date selection matches but the ",
         "values do not -- check `scaling` and the instrument set against the ",
         "run you are starting from.", call. = FALSE)

  invisible(TRUE)
}


#' Read a model configuration
#'
#' Accepts either a path to a config R file or an already-built list. A path is
#' sourced into a private environment, so nothing lands in the global
#' environment and the file can live anywhere. Mirrors `read_config()` in
#' 00_preprocessing.R.
#'
#' @param config Path to a config file defining `model_config`, or the list
#'   itself.
#' @return The `model_config` list.
#' @export
read_model_config <- function(config) {
  if (is.list(config)) return(config)
  if (!is.character(config) || length(config) != 1)
    stop("`config` must be a file path or a configuration list.", call. = FALSE)
  if (!file.exists(config)) stop("Config file not found: ", config, call. = FALSE)

  env <- new.env(parent = globalenv())
  sys.source(normalizePath(config), envir = env)
  if (!exists("model_config", envir = env, inherits = FALSE))
    stop("Config file '", config, "' does not define `model_config`.",
         call. = FALSE)
  get("model_config", envir = env, inherits = FALSE)
}


#' Validate a model configuration against the data it will be fit to
#'
#' Every problem found is collected and reported in a single error, rather than
#' surfacing one at a time over successive runs -- most of these failures are
#' otherwise silent. A misspelled instrument name leaves those rows with a fixed
#' zero intercept and no error term, which fits happily and means nothing.
#'
#' @param config The `model_config` list.
#' @param marss_inputs Optionally the `marss_inputs.rds` contents the model will
#'   be fit to. When supplied, the site, date and variable names in the config
#'   are checked against it.
#'
#' @return `config`, invisibly, if it is valid. Otherwise an error listing every
#'   problem found.
#' @export
validate_model_config <- function(config, marss_inputs = NULL) {

  if (!is.list(config)) stop("`config` must be a list.", call. = FALSE)
  problems <- character()
  add <- function(...) problems <<- c(problems, paste0(...))

  required <- c("data", "scaling", "structure", "inits", "fitting", "output")
  absent   <- setdiff(required, names(config))
  if (length(absent))
    add("Missing required config key(s): ", paste(absent, collapse = ", "))

  # A config without a `reconstruction` block is still valid: the step then
  # runs on its defaults, which is what every config predating it does.
  optional <- "reconstruction"

  extra <- setdiff(names(config), c(required, optional))
  if (length(extra))
    warning("Unrecognised config key(s), which will be ignored: ",
            paste(extra, collapse = ", "),
            ". Check for a typo against: ",
            paste(c(required, optional), collapse = ", "),
            call. = FALSE, immediate. = TRUE)

  known_vars  <- marss_inputs$row_var_keys
  known_sites <- marss_inputs$row_site_keys

  # ---- data ----------------------------------------------------------------
  dp <- config$data
  if (!is.list(dp)) {
    add("`data` must be a list.")
  } else {
    if (!is.character(dp$input) || length(dp$input) != 1 || !nzchar(dp$input))
      add("`data$input` must be a single non-empty path.")
    if (!is.null(dp$sites)) {
      if (!is.character(dp$sites) || !length(dp$sites)) {
        add("`data$sites` must be a character vector, or NULL for all sites.")
      } else if (!is.null(known_sites)) {
        miss <- setdiff(dp$sites, known_sites)
        if (length(miss))
          add("`data$sites` names site(s) absent from the data: ",
              paste(miss, collapse = ", "), ". Available: ",
              paste(sort(unique(known_sites)), collapse = ", "))
      }
    }
    # The holdout: sites that are fitted, but with their in situ series kept out
    # of the observation matrix so it survives untouched for later validation.
    if (!is.null(dp$sites_validation)) {
      if (!is.character(dp$sites_validation) || !length(dp$sites_validation)) {
        add("`data$sites_validation` must be a character vector, or NULL for ",
            "no holdout.")
      } else {
        if (!is.null(known_sites)) {
          miss <- setdiff(dp$sites_validation, known_sites)
          if (length(miss))
            add("`data$sites_validation` names site(s) absent from the data: ",
                paste(miss, collapse = ", "), ". Available: ",
                paste(sort(unique(known_sites)), collapse = ", "))
        }
        # A validation site outside `data$sites` is not held out, it is simply
        # not fitted -- almost certainly a mistake rather than an intent.
        if (is.character(dp$sites)) {
          outside <- setdiff(dp$sites_validation, dp$sites)
          if (length(outside))
            add("`data$sites_validation` must be a subset of `data$sites`; ",
                "these are not fitted at all: ",
                paste(outside, collapse = ", "))
        }
      }
    }

    vv <- dp$validation_variables
    if (!is.null(vv)) {
      if (!is.character(vv) || !length(vv)) {
        add("`data$validation_variables` must be a character vector, or NULL ",
            "for the default (", .VALIDATION_VARIABLES, ").")
      } else if (!is.null(known_vars)) {
        miss <- setdiff(vv, known_vars)
        if (length(miss))
          add("`data$validation_variables` names variable(s) absent from the ",
              "data: ", paste(miss, collapse = ", "))
      }
    }

    for (k in c("start_date", "end_date")) {
      v <- dp[[k]]
      if (!is.null(v) && is.na(suppressWarnings(as.Date(v))))
        add("`data$", k, "` is not a date: ", v)
    }
  }

  # ---- the holdout against the rest of the config --------------------------
  if (is.list(dp) && is.character(dp$sites_validation) &&
      length(dp$sites_validation)) {

    held_vars <- dp$validation_variables %||% .VALIDATION_VARIABLES

    # Scaling is per site from one reference instrument, so holding that
    # reference out leaves the validation sites with no mean and sd to scale by.
    if (isTRUE(config$scaling$enabled) &&
        config$scaling$variable %in% held_vars)
      add("`data$validation_variables` holds out '", config$scaling$variable,
          "', which is `scaling$variable`. The validation sites would then ",
          "have no scaling reference. Hold out a different variable, or fit ",
          "with scaling disabled.")

    # Holding out every row of an instrument leaves its error and intercept
    # parameters with nothing to estimate them from, and `build_model_data()`
    # would drop it from the model with only a warning.
    if (!is.null(known_sites) && !is.null(known_vars)) {
      fitted <- if (is.null(dp$sites)) rep(TRUE, length(known_sites)) else
        known_sites %in% dp$sites
      for (v in intersect(held_vars, known_vars[fitted])) {
        left <- setdiff(unique(known_sites[fitted & known_vars == v]),
                        dp$sites_validation)
        if (!length(left))
          add("The holdout takes every row of '", v,
              "': no site would keep it. Its error and intercept parameters ",
              "would have no data. Hold it out at fewer sites, or drop it ",
              "from `structure$instruments`.")
      }
    }
  }

  # ---- scaling -------------------------------------------------------------
  sc <- config$scaling
  if (!is.list(sc)) {
    add("`scaling` must be a list.")
  } else if (isTRUE(sc$enabled)) {
    if (!is.character(sc$variable) || length(sc$variable) != 1)
      add("`scaling$variable` must name a single variable.")
    else if (!is.null(known_vars) && !sc$variable %in% known_vars)
      add("`scaling$variable` = '", sc$variable, "' is not in the data.")
    # Scaling is per site, so the reference has to reach every site rather than
    # being one nominated series.
    if (!is.null(sc$site))
      warning("`scaling$site` is no longer used: scaling is per site, taken ",
              "from `scaling$variable` at each site.",
              call. = FALSE, immediate. = TRUE)
  }

  # ---- structure -----------------------------------------------------------
  st <- config$structure
  if (!is.list(st)) {
    add("`structure` must be a list.")
  } else {
    if (!is.character(st$state_structure) ||
        length(st$state_structure) != 1 ||
        !st$state_structure %in% .STATE_STRUCTURES)
      add("`structure$state_structure` must be one of: ",
          paste(.STATE_STRUCTURES, collapse = ", "))

    m <- st$m_factors
    if (!is.numeric(m) || length(m) != 1 || is.na(m) || m < 1 ||
        m != as.integer(m)) {
      add("`structure$m_factors` must be a single positive integer.")
    } else if (!is.null(known_sites)) {
      keep <- if (is.null(config$data$sites)) known_sites else
        known_sites[known_sites %in% config$data$sites]
      n_site <- length(unique(keep))
      if (n_site && m > n_site)
        add("`structure$m_factors` (", m, ") exceeds the number of sites (",
            n_site, "): an m-factor model needs at least m sites.")
    }

    ins <- st$instruments
    if (!is.list(ins) || !length(ins) || is.null(names(ins)) ||
        !all(nzchar(names(ins)))) {
      add("`structure$instruments` must be a non-empty named list.")
    } else {
      if (!is.null(known_vars)) {
        miss <- setdiff(names(ins), known_vars)
        if (length(miss))
          add("`structure$instruments` names variable(s) absent from the data: ",
              paste(miss, collapse = ", "))
        # The dangerous direction: rows the config says nothing about get a
        # fixed zero intercept and no error term, and fit silently.
        keep_rows <- if (is.null(config$data$sites)) rep(TRUE, length(known_vars))
                     else known_sites %in% config$data$sites
        unspoken <- setdiff(unique(known_vars[keep_rows]), names(ins))
        if (length(unspoken))
          add("The data contain variable(s) with no `structure$instruments` ",
              "entry, which would enter the model with no intercept and no ",
              "error term: ", paste(unspoken, collapse = ", "))
      }
      for (v in names(ins)) {
        s <- ins[[v]]
        if (!is.list(s)) { add("`instruments$", v, "` must be a list."); next }
        if (!is.character(s$tag) || length(s$tag) != 1 || !nzchar(s$tag))
          add("`instruments$", v, "$tag` must be a single non-empty string.")
        if (!isTRUE(s$intercept %in% .INTERCEPT_KINDS))
          add("`instruments$", v, "$intercept` must be one of: ",
              paste(.INTERCEPT_KINDS, collapse = ", "))
        if (!isTRUE(s$seasonality %in% .SEASONALITY_KINDS))
          add("`instruments$", v, "$seasonality` must be one of: ",
              paste(.SEASONALITY_KINDS, collapse = ", "))
        if (!is.character(s$error) || length(s$error) != 1 || !nzchar(s$error))
          add("`instruments$", v, "$error` must name an observation variance.")
        if (!is.null(s$day_effect) &&
            (!is.logical(s$day_effect) || length(s$day_effect) != 1 ||
             is.na(s$day_effect)))
          add("`instruments$", v, "$day_effect` must be TRUE or FALSE.")
        if (!is.logical(s$site_state) || length(s$site_state) != 1 ||
            is.na(s$site_state))
          add("`instruments$", v, "$site_state` must be TRUE or FALSE.")
      }
      tags <- vapply(ins, function(s) as.character(s$tag %||% NA), character(1))
      if (anyDuplicated(tags[!is.na(tags)]))
        add("`structure$instruments` has duplicate `tag` values, which would ",
            "collide in the parameter names: ",
            paste(unique(tags[duplicated(tags)]), collapse = ", "))
    }

    v0 <- st$init_var_x0
    if (!is.numeric(v0) || length(v0) != 1 || !is.finite(v0) || v0 <= 0)
      add("`structure$init_var_x0` must be a single positive number.")
  }

  # ---- inits ---------------------------------------------------------------
  it <- config$inits
  if (!is.list(it)) {
    add("`inits` must be a list.")
  } else if (isTRUE(it$enabled)) {
    v <- it$shared_from
    if (!is.character(v) || length(v) != 1) {
      add("`inits$shared_from` must name a single variable.")
    } else {
      if (!is.null(known_vars) && !v %in% known_vars)
        add("`inits$shared_from` = '", v, "' is not in the data.")
      if (is.list(config$structure$instruments) &&
          !v %in% names(config$structure$instruments))
        add("`inits$shared_from` = '", v, "' has no `structure$instruments` ",
            "entry, so it is not part of the model.")
    }
    if (!is.logical(it$seed_independent) || length(it$seed_independent) != 1 ||
        is.na(it$seed_independent))
      add("`inits$seed_independent` must be TRUE or FALSE.")
    if (!is.numeric(it$min_obs) || length(it$min_obs) != 1 ||
        !is.finite(it$min_obs) || it$min_obs < 1)
      add("`inits$min_obs` must be a single positive number.")
  }

  # ---- fitting -------------------------------------------------------------
  ft <- config$fitting
  if (!is.list(ft)) {
    add("`fitting` must be a list.")
  } else {

    method <- ft$method
    if (!is.character(method) || length(method) != 1 ||
        !method %in% .METHODS) {
      add("`fitting$method` must be one of: ", paste(.METHODS, collapse = ", "))
      method <- NA_character_
    }

    # Each method has its own set of legal control names -- MARSS rejects a
    # control belonging to another method outright -- so the controls are keyed
    # by method rather than shared.
    if (!is.list(ft$controls) || is.null(names(ft$controls))) {
      add("`fitting$controls` must be a list keyed by method name, e.g. ",
          "list(kem = list(...), TMB = list(...)).")
    } else if (!is.na(method) && !method %in% names(ft$controls)) {
      add("`fitting$controls` has no entry for method '", method,
          "'. Present: ", paste(names(ft$controls), collapse = ", "))
    }

    if (!is.list(ft$warmup_controls))
      add("`fitting$warmup_controls` must be a list.")
    if (!is.numeric(ft$chunks) || length(ft$chunks) != 1 ||
        is.na(ft$chunks) || ft$chunks < 1)
      add("`fitting$chunks` must be a single positive integer.")

    # The check the whole method option turns on. A day effect puts a shared
    # covariance in the off-diagonal of R, and only EM can fit that; the direct
    # optimisers require each block of a variance-covariance matrix to be fixed,
    # diagonal, or wholly unconstrained. Caught here rather than several minutes
    # into a run, where MARSS reports it as a model specification problem
    # without saying which setting caused it.
    if (!is.na(method) && .needs_diagonal_R(method)) {
      with_day <- .day_effect_vars(config)
      if (length(with_day))
        add("`fitting$method` = '", method, "' cannot fit the day effect on: ",
            paste(with_day, collapse = ", "),
            ". A day effect is a shared covariance in the off-diagonal of R, ",
            "which only method 'kem' can estimate. Either set ",
            "`day_effect = FALSE` for those instruments, or fit with 'kem'.")
    }

    # ---- fitting$init ------------------------------------------------------
    # ---- starting from a fit already on disk --------------------------------
    for (k in c("resume_from", "init$from")) {
      v <- if (identical(k, "resume_from")) ft$resume_from else ft$init$from
      if (!is.null(v) && (!is.character(v) || length(v) != 1 || !nzchar(v)))
        add("`fitting$", k, "` must be a single path, or NULL.")
    }
    if (!is.null(ft$resume_from) && !is.null(ft$init$from))
      add("`fitting$init$from` and `fitting$resume_from` both name a fit to ",
          "start from. Set one: `init$from` replaces the initialisation stage, ",
          "`resume_from` continues a final fit.")

    it <- ft$init
    if (!is.null(it)) {
      # Supplying an initialisation fit and switching the stage off are
      # contradictory instructions, and quietly honouring one of them would mean
      # the run silently ignores the fit it was pointed at.
      if (is.list(it) && !is.null(it$from) && identical(it$enabled, FALSE))
        add("`fitting$init$from` names a fit to start from, but ",
            "`fitting$init$enabled` is FALSE. Drop --no-init, or drop ",
            "--init-from to fit without an initialisation stage.")

      if (!is.list(it)) {
        add("`fitting$init` must be a list, or absent.")
      } else if (isTRUE(it$enabled)) {

        if (!is.character(it$method) || length(it$method) != 1 ||
            !it$method %in% .METHODS)
          add("`fitting$init$method` must be one of: ",
              paste(.METHODS, collapse = ", "))
        if (!is.list(it$controls))
          add("`fitting$init$controls` must be a list.")

        # Checked here as well as at build time so a bad value is reported
        # alongside every other problem rather than stopping the run alone.
        ok_scalar <- is.numeric(it$B_values) && length(it$B_values) == 1L &&
                     is.null(names(it$B_values))
        if (ok_scalar) {
          if (!is.finite(it$B_values))
            add("`fitting$init$B_values` must be a finite number.")
        } else if (is.list(it$B_values) ||
                   (is.numeric(it$B_values) && !is.null(names(it$B_values)))) {
          need <- .b_param_names(config)
          have <- names(unlist(it$B_values))
          absent <- setdiff(need, have)
          if (length(absent))
            add("`fitting$init$B_values` has no value for: ",
                paste(absent, collapse = ", "),
                ". The free entries of B for state structure '",
                config$structure$state_structure %||% "?", "' are: ",
                paste(need, collapse = ", "))
        } else {
          add("`fitting$init$B_values` must be a single number or a named list.")
        }
      }
    }
  }

  # ---- output --------------------------------------------------------------
  op <- config$output
  if (!is.list(op)) {
    add("`output` must be a list.")
  } else {
    if (!is.null(op$output_root) &&
        (!is.character(op$output_root) || length(op$output_root) != 1))
      add("`output$output_root` must be a single path, or NULL.")
    if (!is.null(op$run_name) &&
        (!is.character(op$run_name) || length(op$run_name) != 1))
      add("`output$run_name` must be a single name, or NULL.")
  }

  # ---- reconstruction ------------------------------------------------------
  rc <- config$reconstruction
  if (!is.null(rc)) {
    if (!is.list(rc)) {
      add("`reconstruction` must be a list, or absent.")
    } else {
      if (!is.logical(rc$enabled) || length(rc$enabled) != 1 ||
          is.na(rc$enabled))
        add("`reconstruction$enabled` must be TRUE or FALSE.")

      if (!is.null(rc$instruments)) {
        if (!is.character(rc$instruments) || !length(rc$instruments)) {
          add("`reconstruction$instruments` must be a character vector, or ",
              "NULL for the site view alone.")
        } else if (is.list(config$structure$instruments)) {
          miss <- setdiff(rc$instruments, names(config$structure$instruments))
          if (length(miss))
            add("`reconstruction$instruments` names instrument(s) with no ",
                "`structure$instruments` entry: ", paste(miss, collapse = ", "),
                ". Available: ",
                paste(names(config$structure$instruments), collapse = ", "))
        }
      }

      bad <- setdiff(rc$formats, c("csv", "rds"))
      if (!is.character(rc$formats) || !length(rc$formats) || length(bad))
        add("`reconstruction$formats` must be a non-empty subset of: csv, rds.")

      # Named rather than implemented: the states are reconstructed at the
      # parameter estimates. Rejecting TRUE here, before a fit starts, is
      # cheaper than discovering it hours later at the reconstruction step.
      if (!identical(rc$parameter_uncertainty, FALSE))
        add("`reconstruction$parameter_uncertainty` must be FALSE: no method ",
            "for propagating parameter uncertainty is implemented yet.")
    }
  }

  if (length(problems))
    stop("Invalid model configuration (", length(problems), " problem(s)):\n",
         paste0("  - ", problems, collapse = "\n"), call. = FALSE)

  invisible(config)
}


#' Subset and scale the observations a model will be fit to
#'
#' Applies the site and date selection from the config, renumbers the surviving
#' sites 1..n, locates the rows belonging to each instrument, and scales the
#' matrix.
#'
#' `data$sites_validation` names fitted sites whose `data$validation_variables`
#' rows -- in situ by default -- are withheld from the matrix. Those sites are
#' still fitted, from their satellite series alone, so the reconstruction there
#' can be checked against data the model never saw. The withheld observations
#' come back in `validation`.
#'
#' Scaling uses one pair of constants for the whole matrix, taken from a single
#' reference series, so every series stays on a common scale and the estimated
#' biases and loadings remain comparable across instruments. Those constants are
#' saved with the fit so the reconstruction step can return predictions to
#' degrees C; with scaling disabled they are written as mu = 0, sigma = 1, so
#' that path does not have to know which was used.
#'
#' @param marss_inputs The `marss_inputs.rds` contents from the preprocessing
#'   step: `ts_matrix`, `harmonics`, `row_var_keys`, `row_site_keys`,
#'   `row_site_index`, `col_dates`.
#' @param config The `model_config` list.
#' @param verbose Print a summary of what was kept.
#'
#' @return A list with `y` (the selected observations, unscaled), `y_fit` (what
#'   is passed to MARSS), `d` (the harmonics on the same time axis), `scales`,
#'   `dates`, `row_vars`, `row_site`, `row_site_index`, `rows_by_var`,
#'   `sites_by_var`, `n_site`, `k_obs`, `n_harm`, `observations` (the
#'   subsetted input, saved for the reconstruction step) and `validation`
#'   (`sites`, `variables` and the withheld observations as a long data frame of
#'   `date`, `site`, `variable`, `value`; zero rows when nothing is held out).
#' @export
build_model_data <- function(marss_inputs, config, verbose = TRUE) {

  dat <- marss_inputs
  need <- c("ts_matrix", "harmonics", "row_var_keys", "row_site_keys",
            "row_site_index", "col_dates")
  absent <- setdiff(need, names(dat))
  if (length(absent))
    stop("`marss_inputs` is missing: ", paste(absent, collapse = ", "),
         call. = FALSE)

  # ---- rows: site selection and the validation holdout ---------------------
  # Two separate ideas. `sites` says which sites are fitted at all;
  # `sites_validation` names fitted sites whose in situ series is withheld, so
  # the reconstruction there is driven by the satellites alone and can later be
  # checked against observations the fit never saw.
  #
  # The withheld rows are dropped rather than filled with NA. Everything below
  # -- row_vars, row_site_index, rows_by_var, and through them Z, A, R and every
  # parameter name -- is derived from `rows_to_keep`, so dropping reindexes the
  # model on its own, exactly as if those sites had no in situ record. Masking
  # would instead leave an all-NA row carrying parameters nothing can estimate.
  sites     <- config$data$sites
  val_sites <- config$data$sites_validation
  val_vars  <- config$data$validation_variables %||% .VALIDATION_VARIABLES

  in_sites <- if (is.null(sites)) rep(TRUE, length(dat$row_site_keys)) else
    dat$row_site_keys %in% sites
  held_out <- if (is.null(val_sites)) rep(FALSE, length(dat$row_site_keys)) else
    dat$row_site_keys %in% val_sites & dat$row_var_keys %in% val_vars

  rows_to_keep <- in_sites & !held_out
  if (!any(rows_to_keep))
    stop("No rows left after selecting sites: ",
         paste(sites, collapse = ", "),
         if (any(held_out))
           paste0("\nThe holdout removed every remaining row: ",
                  paste(val_vars, collapse = ", "), " at ",
                  paste(val_sites, collapse = ", "), "."),
         call. = FALSE)

  # A validation site with nothing left is not being validated, it has silently
  # dropped out of the model and out of `n_site`.
  if (any(held_out)) {
    emptied <- setdiff(intersect(val_sites, dat$row_site_keys[in_sites]),
                       dat$row_site_keys[rows_to_keep])
    if (length(emptied))
      stop("The holdout leaves no observations at all at: ",
           paste(emptied, collapse = ", "),
           ".\nThose sites carry only ", paste(val_vars, collapse = ", "),
           ", so withholding it drops them from the model entirely rather ",
           "than validating them.", call. = FALSE)
  }

  # ---- columns: date selection ---------------------------------------------
  col_dates <- as.Date(dat$col_dates)
  keep_cols <- rep(TRUE, length(col_dates))
  if (!is.null(config$data$start_date))
    keep_cols <- keep_cols & col_dates > as.Date(config$data$start_date)
  if (!is.null(config$data$end_date))
    keep_cols <- keep_cols & col_dates <= as.Date(config$data$end_date)
  if (!any(keep_cols))
    stop("No columns left after applying the date range.", call. = FALSE)

  # The withheld series, over the same dates the model is fitted to and in
  # degrees C: taken from `ts_matrix` before any scaling, so it stays directly
  # comparable to the reconstruction written by 03_reconstruct_states.R.
  validation <- list(
    sites     = val_sites,
    variables = if (is.null(val_sites)) character(0) else val_vars,
    data      = .withheld_long(dat, which(in_sites & held_out), keep_cols))

  y <- dat$ts_matrix[rows_to_keep, keep_cols, drop = FALSE]
  d <- dat$harmonics[, keep_cols, drop = FALSE]

  row_vars       <- dat$row_var_keys[rows_to_keep]
  row_site       <- dat$row_site_keys[rows_to_keep]
  row_site_index <- .dense_rank(dat$row_site_index[rows_to_keep])

  n_site <- length(unique(row_site))
  k_obs  <- nrow(y)
  n_harm <- nrow(d)

  # ---- rows belonging to each instrument -----------------------------------
  # Built once here so the Z, a, D and R builders all index the same way; in the
  # exploratory scripts each rebuilt its own copy.
  instruments  <- config$structure$instruments
  rows_by_var  <- list()
  sites_by_var <- list()
  for (v in names(instruments)) {
    r <- which(row_vars == v)
    if (!length(r)) {
      warning("No rows for instrument '", v,
              "' after the site and date selection; it is dropped from the ",
              "model.", call. = FALSE, immediate. = TRUE)
      next
    }
    rows_by_var[[v]]  <- r
    sites_by_var[[v]] <- row_site_index[r]
  }
  if (!length(rows_by_var))
    stop("None of the configured instruments have any rows.", call. = FALSE)

  unspoken <- setdiff(unique(row_vars), names(rows_by_var))
  if (length(unspoken))
    stop("Rows are present for variable(s) with no `structure$instruments` ",
         "entry: ", paste(unspoken, collapse = ", "),
         "\nThey would enter the model with no intercept and no error term.",
         call. = FALSE)

  # ---- scaling -------------------------------------------------------------
  # Scaling is per site: each site's observations are put on the scale of one
  # reference instrument at that site. Every instrument at a site is scaled by
  # the same pair, so within a site the series stay comparable and the estimated
  # biases read as departures from the reference; across sites the reference
  # itself becomes mean 0, variance 1, which is what lets its own intercepts be
  # fixed at zero rather than estimated.
  site_levels <- sort(unique(row_site))
  sc <- config$scaling

  if (isTRUE(sc$enabled)) {

    ref_rows <- which(row_vars == sc$variable)
    if (!length(ref_rows))
      stop("No rows for the scaling reference variable '", sc$variable, "'.",
           call. = FALSE)

    mu <- sigma <- rep(NA_real_, n_site)
    mu[row_site_index[ref_rows]] <-
      apply(y[ref_rows, , drop = FALSE], 1, mean, na.rm = TRUE)
    sigma[row_site_index[ref_rows]] <-
      apply(y[ref_rows, , drop = FALSE], 1, stats::sd, na.rm = TRUE)

    # A site the reference does not cover, or covers too sparsely, has no scale
    # to put its other instruments on -- and would silently produce NA rows.
    bad <- which(!is.finite(mu) | !is.finite(sigma) | sigma <= 0)
    if (length(bad))
      stop("The scaling reference '", sc$variable, "' gives no usable mean and ",
           "standard deviation at: ", paste(site_levels[bad], collapse = ", "),
           ".\nName a variable that covers every selected site, or drop those ",
           "sites from `data$sites`.", call. = FALSE)

    scales <- list(mu = mu, sigma = sigma, sites = site_levels,
                   variable = sc$variable)
    y_fit  <- (y - mu[row_site_index]) / sigma[row_site_index]

  } else {
    scales <- list(mu = rep(0, n_site), sigma = rep(1, n_site),
                   sites = site_levels, variable = NA_character_)
    y_fit  <- y
  }
  rownames(y_fit) <- row_site

  dates <- list(start = min(col_dates[keep_cols]),
                end   = max(col_dates[keep_cols]),
                dt    = dat$time_step %||% 1)

  # Everything the reconstruction step needs to map rows back to sites and
  # instruments, subsetted the same way as `y`.
  observations <- list(
    ts_matrix      = y,
    row_keys       = dat$row_keys[rows_to_keep],
    row_var_keys   = row_vars,
    row_site_keys  = row_site,
    row_site_index = row_site_index,
    site_levels    = sort(unique(row_site)),
    col_dates      = col_dates[keep_cols],
    time_step      = dat$time_step %||% 1
  )

  if (verbose) {
    message("Observations: ", k_obs, " rows x ", ncol(y), " columns, ",
            n_site, " sites, ",
            format(sum(!is.na(y)), big.mark = ","), " observed cells (",
            sprintf("%.2f%%", 100 * sum(!is.na(y)) / length(y)), ").")
    message("  dates:   ", dates$start, " to ", dates$end)
    message("  rows per instrument: ",
            paste(sprintf("%s=%d", names(rows_by_var),
                          lengths(rows_by_var)), collapse = ", "))
    message("  scaling: ",
            if (isTRUE(sc$enabled))
              sprintf("per site from %s, mu %.2f-%.2f, sigma %.2f-%.2f",
                      scales$variable, min(scales$mu), max(scales$mu),
                      min(scales$sigma), max(scales$sigma))
            else "disabled")
    if (nrow(validation$data))
      message("  holdout: ", paste(validation$variables, collapse = ", "),
              " withheld at ", paste(validation$sites, collapse = ", "), " (",
              format(nrow(validation$data), big.mark = ","),
              " observations kept back for validation).")
  }

  list(y = y, y_fit = y_fit, d = d, scales = scales, dates = dates,
       row_vars = row_vars, row_site = row_site,
       row_site_index = row_site_index,
       rows_by_var = rows_by_var, sites_by_var = sites_by_var,
       n_site = n_site, k_obs = k_obs, n_harm = n_harm,
       observations = observations, validation = validation)
}


#' Build the observation matrix Z
#'
#' Under "site_plus_factors" the first `n_site` columns select each row's own
#' site state, and the remaining `m_factors` columns hold the dynamic factor
#' loadings. Under "factors_only" the site block is dropped and only the
#' loadings remain.
#'
#' An instrument with `site_state = FALSE` gets zeros in the site block: it
#' observes the seasonal terms and the shared factors, but not the site's own
#' anomaly. That is the right structure for an interpolated product such as MUR,
#' whose value at a site is a smoothed regional field rather than a measurement
#' of that site -- letting it load on the site state would make it evidence
#' about an anomaly it cannot actually see.
#'
#' The loadings come from `dfa_loadings()` in marss_matrix_functions.R, which
#' places the triangular zero restrictions that identify the factors. Every
#' instrument observing a given site shares that site's loading row, so the
#' factors describe the site rather than the instrument.
#'
#' @param md Output of `build_model_data()`.
#' @param config The `model_config` list.
#' @return A list matrix suitable for the `Z` element of a MARSS model list.
#' @export
build_observation_matrix <- function(md, config) {

  if (!exists("dfa_loadings"))
    stop("`dfa_loadings()` not found. Source R/marss_matrix_functions.R first.",
         call. = FALSE)

  m_factors <- config$structure$m_factors

  Lambda  <- dfa_loadings(md$n_site, m_factors, "lambda")
  eta_obs <- matrix(list(0), nrow = md$k_obs, ncol = m_factors)
  for (v in names(md$rows_by_var))
    eta_obs[md$rows_by_var[[v]], ] <-
      Lambda[md$sites_by_var[[v]], , drop = FALSE]

  if (identical(config$structure$state_structure, "site_plus_factors")) {
    I_sites <- matrix(0, nrow = md$n_site, ncol = md$n_site)
    diag(I_sites) <- 1
    chi_obs <- I_sites[md$row_site_index, , drop = FALSE]

    # Zeroed in place rather than assembled from row subsets, so the result
    # does not depend on the instruments happening to be in a particular order
    # down the matrix.
    for (v in names(md$rows_by_var))
      if (!isTRUE(config$structure$instruments[[v]]$site_state))
        chi_obs[md$rows_by_var[[v]], ] <- 0

    # A site state that nothing loads on is not identified: the model would
    # estimate its variance and persistence from no observations at all. This
    # happens as soon as a site's only instrument is one excluded from the site
    # block, which is easy to arrive at by narrowing `data$sites`.
    unseen <- which(colSums(chi_obs) == 0)
    if (length(unseen))
      stop("No instrument loads on the site state for: ",
           paste(sort(unique(md$row_site))[unseen], collapse = ", "),
           ".\nThose states have no observations. Either drop the site from ",
           "`data$sites`, or set `site_state = TRUE` for an instrument that ",
           "observes it.", call. = FALSE)

    return(cbind(chi_obs, eta_obs))
  }

  eta_obs
}


#' Build the observation intercepts a
#'
#' One entry per row, named according to the instrument's `intercept` kind:
#' `mu_<i>` for a shared site mean, `mu_<tag>+mu_<i>` for a site mean plus an
#' instrument offset (MARSS reads the `+` as a sum of two estimated
#' parameters), and `mu_<tag>_<i>` for an instrument with its own mean at each
#' site.
#'
#' `"zero"` fixes the intercept at 0 rather than estimating it. That is the
#' correct choice for the instrument the observations were scaled by: after
#' per-site scaling its series has mean 0 at every site by construction, so an
#' estimated intercept would be a free parameter with nothing to explain. Fixing
#' it also makes every other intercept read as a difference from that reference,
#' at the site and instrument level.
#'
#' Returned as a list rather than a character vector so that fixed and estimated
#' entries can coexist -- a character matrix would make the fixed `0` a
#' parameter named "0". MARSS treats an all-character list matrix identically to
#' a character matrix, so this costs nothing where nothing is fixed.
#'
#' @param md Output of `build_model_data()`.
#' @param config The `model_config` list.
#' @return A list of length `k_obs`: numeric entries fixed, character estimated.
#' @export
build_observation_intercepts <- function(md, config) {

  n_site <- md$n_site
  a <- rep(list(0), md$k_obs)

  for (v in names(md$rows_by_var)) {
    spec <- config$structure$instruments[[v]]
    if (identical(spec$intercept, "zero")) next   # left fixed at 0
    nms <- switch(spec$intercept,
      "site"            = paste0("mu_", seq_len(n_site)),
      "site+instrument" = paste0("mu_", spec$tag, "+mu_", seq_len(n_site)),
      "independent"     = paste0("mu_", spec$tag, "_", seq_len(n_site)),
      stop("Unknown intercept kind for '", v, "': ", spec$intercept,
           call. = FALSE))
    a[md$rows_by_var[[v]]] <- as.list(nms[md$sites_by_var[[v]]])
  }

  a
}


#' Build the seasonal covariate effects D
#'
#' One column per harmonic term, one row per observation row. An instrument with
#' "shared" seasonality uses the site's own coefficients `c<k>_<i>`; one with
#' "independent" seasonality gets its own set, `c<k>_<tag>_<i>`.
#'
#' The number of harmonic terms is taken from the covariate matrix rather than
#' assumed, so changing `harmonics` in the preprocessing config carries through
#' without touching this code.
#'
#' @param md Output of `build_model_data()`.
#' @param config The `model_config` list.
#' @return A character matrix, `k_obs` by `n_harm`.
#' @export
build_covariate_effects <- function(md, config) {

  n_site <- md$n_site
  n_harm <- md$n_harm
  D <- matrix(0, nrow = md$k_obs, ncol = n_harm)

  for (v in names(md$rows_by_var)) {
    spec   <- config$structure$instruments[[v]]
    infix  <- switch(spec$seasonality,
      "shared"      = "",
      "independent" = paste0(spec$tag, "_"),
      stop("Unknown seasonality kind for '", v, "': ", spec$seasonality,
           call. = FALSE))

    # n_site x n_harm: rows are sites, columns are harmonic terms.
    nms <- vapply(seq_len(n_harm),
                  function(k) paste0("c", k, "_", infix, seq_len(n_site)),
                  character(n_site))

    D[md$rows_by_var[[v]], ] <- nms[md$sites_by_var[[v]], , drop = FALSE]
  }

  D
}


#' Build the observation covariance R
#'
#' Each instrument gets its own variance on the diagonal. An instrument with
#' `day_effect = TRUE` additionally gets a single covariance, `cov_<tag>`,
#' between every pair of its own rows: a cloud or a wind event affects a whole
#' scene at once, so on any given day that instrument's errors are shifted
#' together across sites rather than independently.
#'
#' A day effect can only be estimated by the EM algorithm. The direct
#' optimisers require every block of a variance-covariance matrix to be fixed,
#' diagonal, or wholly unconstrained, and a diagonal of per instrument
#' variances with one shared off-diagonal covariance is none of those.
#' `validate_model_config()` rejects that combination before a run starts.
#'
#' @param md Output of `build_model_data()`.
#' @param config The `model_config` list.
#' @return A `k_obs` by `k_obs` list matrix.
#' @export
build_observation_covariance <- function(md, config) {

  R <- matrix(list(0), nrow = md$k_obs, ncol = md$k_obs)

  for (v in names(md$rows_by_var)) {
    spec <- config$structure$instruments[[v]]
    diag(R)[md$rows_by_var[[v]]] <- spec$error
  }

  for (v in names(md$rows_by_var)) {
    spec <- config$structure$instruments[[v]]
    if (!isTRUE(spec$day_effect)) next
    nm   <- paste0("cov_", spec$tag)
    rows <- md$rows_by_var[[v]]
    for (i in rows) for (j in rows) if (i != j) R[[i, j]] <- nm
  }

  R
}


#' Build the state model: B, Q, U, x0 and V0
#'
#' Under "site_plus_factors" the states are the `n_site` site level anomalies
#' followed by the `m_factors` shared factors. The site states share one
#' autocorrelation, `rho_chi`, and one innovation variance, `tau_1`; the factors
#' share `rho_eta` and have their variance fixed at 1, which together with the
#' triangular loadings identifies the factor scale.
#'
#' Under "factors_only" only the factors remain, with autocorrelation `rho` and
#' an estimated innovation variance `tau`. Freeing that variance removes the
#' scale identification, so the individual loadings from such a fit are not
#' interpretable on their own even though the states and fitted values are.
#'
#' @param md Output of `build_model_data()`.
#' @param config The `model_config` list.
#' @return A list with `B`, `Q`, `U`, `x0`, `V0` and `n_total`.
#' @export
build_state_model <- function(md, config) {

  n_eta <- config$structure$m_factors

  if (identical(config$structure$state_structure, "site_plus_factors")) {

    n_chi   <- md$n_site
    n_total <- n_chi + n_eta

    B <- matrix(list(0), nrow = n_total, ncol = n_total)
    diag(B)[seq_len(n_chi)]                <- rep("rho_chi", n_chi)
    diag(B)[n_chi + seq_len(n_eta)]        <- rep("rho_eta", n_eta)

    Q <- matrix(list(0), nrow = n_total, ncol = n_total)
    diag(Q)[seq_len(n_chi)]                <- rep("tau_1", n_chi)
    diag(Q)[n_chi + seq_len(n_eta)]        <- rep(list(1), n_eta)

  } else {

    n_total <- n_eta

    B <- matrix(list(0), nrow = n_total, ncol = n_total)
    diag(B)[seq_len(n_total)] <- rep("rho", n_total)

    Q <- matrix(list(0), nrow = n_total, ncol = n_total)
    diag(Q)[seq_len(n_total)] <- rep("tau", n_total)
  }

  # Hold B at stated values, for a method that cannot estimate it. Applied as a
  # post-pass over whatever names the structure produced, so it works for
  # rho_chi/rho_eta and for rho without branching on the structure again.
  #
  # Fixing B numerically rather than leaving it free and hoping the optimiser
  # moves it matters: TMB leaves a free B at its starting value and still
  # reports it as an estimated parameter. Fixed, it correctly disappears from
  # the parameter set, so what was and was not estimated is legible from the fit.
  B_fixed <- config$structure$B_fixed
  if (!is.null(B_fixed)) {
    for (i in seq_len(n_total)) {
      nm <- B[[i, i]]
      if (is.character(nm)) B[[i, i]] <- unname(B_fixed[[nm]])
    }
  }

  U  <- rep(list(0), n_total)
  x0 <- rep(list(0), n_total)
  V0 <- matrix(list(0), nrow = n_total, ncol = n_total)
  diag(V0) <- config$structure$init_var_x0

  list(B = B, Q = Q, U = U, x0 = x0, V0 = V0, n_total = n_total)
}


#' Assemble the full MARSS model list
#'
#' @param md Output of `build_model_data()`.
#' @param config The `model_config` list.
#' @return The list passed as the `model` argument of `MARSS()`.
#' @export
build_marss_model <- function(md, config) {

  state <- build_state_model(md, config)

  list(
    Z  = build_observation_matrix(md, config),
    A  = matrix(build_observation_intercepts(md, config), ncol = 1),
    R  = build_observation_covariance(md, config),
    B  = state$B,
    U  = matrix(state$U, ncol = 1),
    Q  = state$Q,
    D  = build_covariate_effects(md, config),
    x0 = matrix(state$x0, ncol = 1),
    V0 = state$V0,
    d  = md$d
  )
}


#' Least squares starting values for one instrument's mean and seasonality
#'
#' The EM algorithm starts every parameter at zero, leaving it to discover the
#' annual cycle -- by far the largest signal in the data -- one iteration at a
#' time. Regressing an instrument's series on the harmonics gives a starting
#' point that removes that phase of the run.
#'
#' The regression is run on the scaled matrix that is actually fitted, so the
#' coefficients are on the same scale as the parameters they seed. Sites with
#' fewer than `inits$min_obs` observations, or whose fit is rank deficient, are
#' left at the default start.
#'
#' @param md Output of `build_model_data()`.
#' @param config The `model_config` list.
#' @param variable The instrument whose rows the regression is fit to.
#' @return A data frame with one row per site and columns `site`, `mu` and
#'   `c1`..`c<n_harm>`; `NA` where a site could not be fit.
#' @export
seasonal_inits <- function(md, config, variable) {

  n_harm <- md$n_harm
  rows   <- md$rows_by_var[[variable]]
  if (is.null(rows))
    stop("'", variable, "' has no rows in the selected data.", call. = FALSE)

  out <- data.frame(site = seq_len(md$n_site), mu = NA_real_)
  for (k in seq_len(n_harm)) out[[paste0("c", k)]] <- NA_real_

  for (r in rows) {
    s  <- md$row_site_index[r]
    yy <- md$y_fit[r, ]
    ok <- !is.na(yy)
    if (sum(ok) < config$inits$min_obs) next

    X   <- t(md$d[, ok, drop = FALSE])
    cf  <- stats::coef(stats::lm(yy[ok] ~ X))
    # A rank deficient fit drops terms and returns NA coefficients; seeding from
    # it would put NA into the parameter vector and stop the run.
    if (length(cf) == n_harm + 1L && all(is.finite(cf)))
      out[s, 2:(2 + n_harm)] <- cf
  }

  out
}


#' Collect starting values for every parameter the warm start can seed
#'
#' Two things are seeded, and they are kept separate because they are different
#' parameters.
#'
#' The shared site terms, `mu_<i>` and `c<k>_<i>`, come from a single
#' instrument named by `inits$shared_from`. Any instrument may be named; the
#' question is only which series is the best available description of the
#' site's own seasonal cycle. A gap free product covering every site is usually
#' a better starting point than an exact record covering two of them, even
#' though it is a proxy, since this only sets where the optimiser begins.
#'
#' Instruments with their own mean or seasonality -- `mu_<tag>_<i>`,
#' `c<k>_<tag>_<i>` -- are seeded from their own rows when
#' `inits$seed_independent` is TRUE. Those parameters otherwise start at zero,
#' which for an instrument with an independent seasonal cycle means starting
#' with no seasonal cycle at all.
#'
#' An instrument whose intercept is "site+instrument" contributes nothing here:
#' its offset `mu_<tag>` is a separate parameter, and the site part of its mean
#' is already covered by the shared terms.
#'
#' @param md Output of `build_model_data()`.
#' @param config The `model_config` list.
#' @param verbose Report where each group of starting values came from.
#'
#' @return A named numeric vector of parameter names to starting values, using
#'   the model's own names without a MARSS matrix prefix.
#' @export
build_inits_seeds <- function(md, config, verbose = TRUE) {

  seeds  <- numeric(0)
  n_harm <- md$n_harm
  terms  <- c("mu", paste0("c", seq_len(n_harm)))

  collect <- function(pars, name_for) {
    got <- numeric(0)
    for (i in seq_len(md$n_site)) {
      for (nm in terms) {
        value <- pars[[nm]][i]
        key   <- name_for(nm, i)
        if (is.null(key) || is.null(value) || !is.finite(value)) next
        got[key] <- value
      }
    }
    got
  }

  # ---- shared site terms ---------------------------------------------------
  shared <- config$inits$shared_from
  pars   <- seasonal_inits(md, config, shared)
  got    <- collect(pars, function(nm, i) paste0(nm, "_", i))
  seeds  <- c(seeds, got)
  if (verbose)
    message("  shared site terms: ", length(got), " from ", shared,
            " (", sum(!is.na(pars$mu)), " of ", md$n_site, " sites).")

  # ---- instrument specific terms -------------------------------------------
  if (isTRUE(config$inits$seed_independent)) {
    for (v in names(md$rows_by_var)) {
      spec <- config$structure$instruments[[v]]
      own_mu   <- identical(spec$intercept, "independent")
      own_seas <- identical(spec$seasonality, "independent")
      if (!own_mu && !own_seas) next

      pars <- seasonal_inits(md, config, v)
      got  <- collect(pars, function(nm, i) {
        if (nm == "mu" && !own_mu)  return(NULL)
        if (nm != "mu" && !own_seas) return(NULL)
        paste0(nm, "_", spec$tag, "_", i)
      })
      seeds <- c(seeds, got)
      if (verbose)
        message("  ", v, " own terms: ", length(got), " from its own rows (",
                sum(!is.na(pars$mu)), " of ", md$n_site, " sites).")
    }
  }

  seeds
}


#' Seed a fitted MARSS object's parameter vector with starting values
#'
#' `MARSSvectorizeparam()` is used to recover the parameter names and their
#' ordering from a throwaway fit, the named entries are overwritten, and the
#' vector is converted back into an inits list.
#'
#' MARSS folds the covariate term `D d` into the observation intercept when it
#' converts to its internal form, so a seasonal coefficient may be named either
#' `A.c1_1` or `D.c1_1` depending on the version and the model. Both are tried,
#' rather than assuming one.
#'
#' @param fit A fitted `marssMLE`, typically from a two iteration warm up run.
#' @param seeds Named numeric vector from `build_inits_seeds()`.
#' @param verbose Report how many parameters were seeded.
#'
#' @return A `marssMLE` whose parameters carry the starting values.
#' @export
apply_inits <- function(fit, seeds, verbose = TRUE) {

  init_vec <- MARSSvectorizeparam(fit)
  seeded   <- 0L

  for (nm in names(seeds)) {
    # Whichever form MARSS gave the parameter, seed it there.
    keys <- paste0(c("A.", "D."), nm)
    key  <- keys[keys %in% names(init_vec)]
    if (!length(key)) next
    init_vec[key[1]] <- seeds[[nm]]
    seeded <- seeded + 1L
  }

  if (verbose)
    message("Warm start: seeded ", seeded, " of ", length(seeds),
            " starting values into the parameter vector.")
  if (length(seeds) > 0 && seeded == 0)
    warning("No starting values matched a parameter name, so the fit begins ",
            "from the MARSS defaults.", call. = FALSE, immediate. = TRUE)

  MARSSvectorizeparam(fit, init_vec)
}


#' Carry parameter estimates from one fit to another by name
#'
#' The initialisation stage fits a different model from the one that follows --
#' B is held fixed and the day effects are dropped -- so its estimates cannot be
#' copied across by position. They are matched on the parameter names MARSS
#' itself uses (`A.mu_1`, `R.sigma_2_lst`, `Z.lambda_1_1`, ...). Anything the
#' first stage did not estimate simply has no name to match and keeps the
#' second stage's own starting value.
#'
#' This is deliberately distinct from `apply_inits()`, which takes unprefixed
#' names from the least squares warm start and has to guess between the `A.` and
#' `D.` forms. Here both sides are full MARSS names and must match exactly.
#'
#' @param from_fit The fit supplying estimates.
#' @param to_fit A fit of the target model, used for its parameter ordering.
#' @param verbose Report how many parameters carried over.
#' @param from What to call the source in that report. The default suits the
#'   two stage run; a resumed fit is not an initialisation fit and says so.
#'
#' @return A `marssMLE` for the target model carrying the transferred values.
#' @export
transfer_params <- function(from_fit, to_fit, verbose = TRUE,
                            from = "the initialisation fit") {

  src <- MARSSvectorizeparam(from_fit)
  dst <- MARSSvectorizeparam(to_fit)

  common <- intersect(names(src), names(dst))
  kept   <- setdiff(names(dst), names(src))

  dst[common] <- src[common]

  if (verbose) {
    message("Transferred ", length(common), " of ", length(dst),
            " parameters from ", from, ".")
    if (length(kept))
      message("  left at defaults (not estimable in that stage): ",
              paste(kept, collapse = ", "))
  }
  if (!length(common))
    warning("No parameter names matched between the two stages, so the ",
            "initialisation fit contributed nothing.",
            call. = FALSE, immediate. = TRUE)

  MARSSvectorizeparam(to_fit, dst)
}


#' Fit the model, optionally through a fast initialisation stage first
#'
#' The final fit is run as repeated short calls rather than one long one. Each
#' chunk is saved as it completes, so a run that is interrupted -- or that is
#' still climbing after days of EM iterations -- leaves a usable fit behind, and
#' the log likelihood trace shows whether it is still improving. It stops as
#' soon as MARSS reports convergence.
#'
#' Chunking is an EM idiom. The direct optimisers reach their own stopping rule
#' inside one call, so with `method` other than "kem" a run normally finishes at
#' the first chunk and `chunks` is never reached.
#'
#' When `fitting$init$enabled` is TRUE the run has two stages. The first fits a
#' reduced model -- B held at `init$B_values`, day effects dropped -- with a fast
#' method, typically TMB. Neither restriction is a compromise on the result: the
#' stage exists only to supply starting values, and the parameters it cannot
#' estimate are the ones not transferred. The second stage then fits the model
#' the config actually describes, starting from those values.
#'
#' The throwaway fit used to recover a parameter ordering always runs under EM,
#' whichever method the real fit uses: the ordering is a property of the model,
#' not of the optimiser, and EM will produce it for any model the config can
#' express.
#'
#' @param md Output of `build_model_data()`.
#' @param model_list Output of `build_marss_model()` for the full model.
#' @param config The `model_config` list.
#' @param outdir Directory to write the fits and `loglik.csv` into.
#' @param verbose Print progress.
#'
#' @return A list with the final `fit`, the `trace` data frame, and `init_fit`
#'   (NULL when no initialisation stage ran).
#' @export
fit_marss <- function(md, model_list, config, outdir, verbose = TRUE) {

  y        <- md$y_fit
  method   <- config$fitting$method %||% "kem"
  controls <- config$fitting$controls[[method]]
  init_cfg <- config$fitting$init

  trace <- data.frame(stage = character(), method = character(),
                      chunk = integer(), iterations = integer(),
                      logLik = numeric(), convergence = integer())

  # Rewritten after every fit, so an interrupted run still has its trace.
  record <- function(trace, stage, method, chunk, fit) {
    # The direct optimisers do not report an iteration count the way EM does.
    iters <- if (is.null(fit$numIter)) NA_integer_ else as.integer(fit$numIter)
    trace <- rbind(trace, data.frame(
      stage = stage, method = method, chunk = chunk, iterations = iters,
      logLik = fit$logLik, convergence = fit$convergence))
    utils::write.csv(trace, file.path(outdir, "loglik.csv"), row.names = FALSE)
    message(sprintf("%-5s | chunk %2s | iters %5s | logLik %.3f",
                    stage, chunk, if (is.na(iters)) "-" else iters, fit$logLik))
    trace
  }

  # ---- a fit to start from -------------------------------------------------
  # Either replaces stage 1 with one already on disk, or picks up where a stage
  # 2 left off. Both are the same mechanism -- load a fit, transfer its
  # parameters by name into this model, start there -- and differ only in which
  # file is loaded and what the trace calls it.
  loaded      <- NULL
  start_stage <- NULL
  if (!is.null(init_cfg$from) || !is.null(config$fitting$resume_from)) {

    resuming    <- !is.null(config$fitting$resume_from)
    start_stage <- if (resuming) "resume" else "init"
    src <- .resolve_fit_path(
      if (resuming) config$fitting$resume_from else init_cfg$from,
      prefer = if (resuming) "chunk" else "init")

    if (verbose) message("Starting from ", src, " ...")
    loaded <- readRDS(src)
    if (!inherits(loaded, "marssMLE"))
      stop("Not a fitted MARSS model: ", src, call. = FALSE)
    .check_same_data(loaded, y, paste0("The fit at ", src))

    # Copied in so the new run directory is self-contained, the way every other
    # run directory in this project is.
    file.copy(src, file.path(outdir, if (resuming) "fit_resumed.rds"
                                     else "fit_init.rds"), overwrite = TRUE)
    trace <- record(trace, start_stage, loaded$method %||% "?", NA_integer_,
                    loaded)
  }

  # ---- least squares warm start --------------------------------------------
  # Feeds whichever fit runs first. Starting the seasonal terms at zero is
  # expensive for every method, not just EM.
  #
  # A loaded fit supersedes it: `seeds` is consumed only by stage 1 and by the
  # branch of stage 2 that runs when there is no fit to start from, so with one
  # loaded these regressions would be fit and then discarded.
  seeds <- NULL
  if (isTRUE(config$inits$enabled) && is.null(loaded)) {
    if (verbose) message("Fitting warm start regressions ...")
    seeds <- build_inits_seeds(md, config, verbose = verbose)
  }

  ordering_fit <- function(mod) {
    MARSS(y, mod, method = "kem", control = config$fitting$warmup_controls)
  }

  # ---- stage 1: initialisation ---------------------------------------------
  # A loaded fit stands in for whichever stage it came from, so stage 2 needs no
  # branch of its own: it already knows how to start from an `init_fit`.
  init_fit <- loaded
  if (isTRUE(init_cfg$enabled) && is.null(loaded)) {

    icfg  <- .init_config(config)
    imod  <- build_marss_model(md, icfg)
    imeth <- icfg$fitting$method

    if (verbose) {
      dropped <- .day_effect_vars(config)
      message("Stage 1: ", imeth, ", B fixed at ",
              paste(sprintf("%s=%s", names(icfg$structure$B_fixed),
                            icfg$structure$B_fixed), collapse = ", "),
              if (length(dropped))
                paste0(", day effects dropped (",
                       paste(dropped, collapse = ", "), ")") else "", ".")
    }

    iinits <- NULL
    if (!is.null(seeds)) {
      if (verbose) message("Recovering the parameter ordering ...")
      iinits <- stats::coef(apply_inits(ordering_fit(imod), seeds,
                                        verbose = verbose), type = "list")
    }

    iargs <- list(y = y, model = imod, method = imeth,
                  control = icfg$fitting$controls[[imeth]])
    if (!is.null(iinits)) iargs$inits <- iinits

    init_fit <- do.call(MARSS, iargs)
    saveRDS(init_fit, file.path(outdir, "fit_init.rds"))
    trace <- record(trace, "init", imeth, NA_integer_, init_fit)
  }

  # ---- stage 2: the model the config describes -----------------------------
  if (verbose) message("Stage 2: ", method, ".")

  inits <- NULL
  if (!is.null(init_fit)) {
    # The ordering fit gives the full model's parameter layout; the values then
    # come from stage 1 wherever the names line up.
    inits <- stats::coef(
      transfer_params(init_fit, ordering_fit(model_list), verbose = verbose,
                      from = if (identical(start_stage, "resume"))
                        "the resumed fit" else "the initialisation fit"),
      type = "list")
  } else if (!is.null(seeds)) {
    if (verbose) message("Recovering the parameter ordering ...")
    inits <- stats::coef(apply_inits(ordering_fit(model_list), seeds,
                                     verbose = verbose), type = "list")
  }

  # ---- where B starts ------------------------------------------------------
  # B is free here, and nothing above sets it: stage 1 holds it at a fixed value
  # so it is not a parameter there and has no name to transfer, and the warm
  # start seeds only the seasonal terms. Left alone it would inherit whatever
  # the throwaway ordering fit reached in two EM iterations from MARSS's own
  # default of 1 -- a random walk, and an arbitrary place to begin that moves
  # whenever `warmup_controls` is touched. Under a method that cannot estimate B
  # that arbitrary value is also the answer, reported as though it had been
  # estimated.
  #
  # So every method starts B where the config says, at the same value stage 1
  # held it at. Only the first chunk is seeded; later chunks continue from the
  # one before.
  #
  # A resumed fit is the exception. Unlike an initialisation fit it does carry
  # an estimated B, which transfers by name like everything else, and resetting
  # it to the configured start would throw away the very progress the resume
  # exists to keep.
  if (!is.null(init_cfg$B_values) && is.null(config$structure$B_fixed) &&
      !identical(start_stage, "resume")) {
    b_names <- .b_param_names(config)
    b_vals  <- .resolve_b_values(b_names, init_cfg$B_values)
    if (is.null(inits)) inits <- list()
    inits$B <- matrix(b_vals, ncol = 1, dimnames = list(b_names, NULL))
    if (verbose)
      message("B starts at ",
              paste(sprintf("%s=%s", b_names, b_vals), collapse = ", "), ".")
  }

  fit <- NULL
  for (chunk in seq_len(config$fitting$chunks)) {

    args <- list(y = y, model = model_list, method = method,
                 control = controls)
    start <- if (is.null(fit)) inits else stats::coef(fit, type = "list")
    if (!is.null(start)) args$inits <- start

    fit <- do.call(MARSS, args)

    saveRDS(fit, file.path(outdir, sprintf("fit_chunk_%02d.rds", chunk)))
    trace <- record(trace, "fit", method, chunk, fit)

    if (isTRUE(fit$convergence == 0)) {
      message("Converged after ", chunk, " chunk(s).")
      break
    }
  }

  if (!isTRUE(fit$convergence == 0))
    message("Stopped at the chunk limit (", config$fitting$chunks,
            ") without convergence; raise `fitting$chunks` to continue.")

  list(fit = fit, trace = trace, init_fit = init_fit)
}


#' Describe the model a fit was produced by
#'
#' The dimensions and parameterisation choices the reconstruction step needs to
#' interpret a fit, and the record of how it was fitted.
#'
#' `instruments` and `reconstruction` are the *resolved* config values rather
#' than a pointer at config_used.R, which matters because run_model.R copies the
#' config file before applying the command line overrides: a run launched with
#' --structure or --m-factors is not described by the file sitting next to it,
#' but is described by this.
#'
#' @param md Output of `build_model_data()`.
#' @param config The `model_config` list.
#' @return The design list, as written to design.rds.
#' @export
build_design <- function(md, config) {

  init <- config$fitting$init

  list(n_site          = md$n_site,
       m_factors       = config$structure$m_factors,
       state_structure = config$structure$state_structure,
       n_harmonics     = md$n_harm,
       instruments     = config$structure$instruments,
       reconstruction  = config$reconstruction,
       # The holdout, so a finished run says which sites were fitted without
       # their in situ series without anyone having to read config_used.R --
       # which does not describe a run launched with --sites-validation.
       sites_validation     = md$validation$sites,
       validation_variables = md$validation$variables,
       method          = config$fitting$method %||% "kem",
       day_effects     = .day_effect_vars(config),
       init_method     = if (isTRUE(init$enabled)) init$method else NA_character_,
       init_B_values   = if (isTRUE(init$enabled))
         .resolve_b_values(.b_param_names(config), init$B_values)
         else NULL)
}


#' Write the fixed record of a fit
#'
#' These four objects are the contract with the reconstruction step in
#' 03_reconstruct_states.R: the observation matrix and its row keys, the
#' scaling constants needed to return predictions to degrees C, the date range
#' the columns span, and the model dimensions.
#'
#' A run that held data out also writes `validation.csv` and `validation.rds`:
#' the withheld observations in degrees C, keyed the same way as `states.csv`
#' so the two join on `(site, date)`. Comparing them is not yet automated.
#'
#' @param md Output of `build_model_data()`.
#' @param config The `model_config` list.
#' @param outdir Directory to write into.
#' @param verbose Print what was written.
#' @return The paths written, invisibly.
#' @export
write_model_outputs <- function(md, config, outdir, verbose = TRUE) {

  design <- build_design(md, config)

  paths <- c(observations = file.path(outdir, "observations.rds"),
             scales       = file.path(outdir, "scales.rds"),
             dates        = file.path(outdir, "dates.rds"),
             design       = file.path(outdir, "design.rds"))

  saveRDS(md$observations, paths[["observations"]])
  saveRDS(md$scales,       paths[["scales"]])
  saveRDS(md$dates,        paths[["dates"]])
  saveRDS(design,          paths[["design"]])

  # Only when something was actually withheld: an empty validation.csv in a run
  # directory would read as "the holdout found nothing", not "there was none".
  if (!is.null(md$validation) && nrow(md$validation$data)) {
    paths <- c(paths,
               validation_csv = file.path(outdir, "validation.csv"),
               validation_rds = file.path(outdir, "validation.rds"))
    utils::write.csv(md$validation$data, paths[["validation_csv"]],
                     row.names = FALSE)
    saveRDS(md$validation$data, paths[["validation_rds"]])
  }

  if (verbose)
    message("Wrote ", paste(basename(paths), collapse = ", "), ".")

  invisible(paths)
}
