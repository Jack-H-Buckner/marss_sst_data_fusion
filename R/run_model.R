#!/usr/bin/env Rscript
#############################################################
#############################################################
###
### Command line entry point for the model fitting step.
###
###   Rscript R/run_model.R --config parameters/model.R
###
### Reads the observation matrix written by the preprocessing
### step, builds the MARSS model described by the config, and
### fits it as a series of restarted chunks, writing each one
### alongside the data and constants the reconstruction step
### needs to interpret it.
###
### Run with --help for the full option list.
###
### Jack H. Buckner, Oregon State University, 08/30/2026
### Generated with Claude Code
###
#############################################################
#############################################################

# ---------------------------------------------------------------------------
# Locate the project root from this script's own path, so `source()` and the
# default paths work regardless of the working directory the user invoked from.
# ---------------------------------------------------------------------------
.script_path <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(normalizePath(f[1], mustWork = FALSE))
  NA_character_   # sourced interactively; fall back to the working directory
}

SCRIPT <- .script_path()
ROOT   <- normalizePath(if (is.na(SCRIPT)) getwd() else dirname(dirname(SCRIPT)))

# A relative path in the config is resolved against the project root, not the
# caller's working directory -- otherwise the same config would read and write
# different places depending on where the run was launched from.
.is_absolute <- function(path) grepl("^(/|~|[A-Za-z]:)", path)
.from_root   <- function(path) if (.is_absolute(path)) path else file.path(ROOT, path)

DEFAULTS <- list(
  config  = file.path(ROOT, "parameters", "example_model.R"),
  outroot = file.path(ROOT, "models")
)

USAGE <- "
Usage: Rscript R/run_model.R [options]

Fits the configured MARSS data fusion model to the observation matrix written
by the preprocessing step, and writes every chunk of the fit together with the
data, scaling constants and design it was fit to.

Options:
  --config PATH      Config file defining `model_config`.
                     Default: parameters/model.R
  --input PATH       marss_inputs.rds, or the preprocessing run directory that
                     holds it. Overrides the config's `data$input`.
  --outdir PATH      Output directory, used verbatim. Overrides both
                     --run-name and the config's `output` block.
                     Default: <output_root>/<run-name>
  --run-name NAME    Name of the run subdirectory within the config's
                     `output_root`. Overrides the config's `run_name`.
                     Default: the config file stem plus today's date.
  --sites LIST       Comma separated site keys to fit, overriding the config.
                     Pass 'all' to fit every site present.
  --start-date DATE  Earliest date to fit, YYYY-MM-DD. Pass 'none' to remove
                     the config's lower bound.
  --end-date DATE    Latest date to fit, YYYY-MM-DD. Pass 'none' to remove it.
  --structure NAME   State structure: 'site_plus_factors' or 'factors_only'.
  --m-factors N      Number of shared dynamic factors.
  --method NAME      Method for the final fit: 'kem' (the EM algorithm, the
                     default and the only one that can fit a day effect),
                     'BFGS', 'TMB', 'BFGS_TMB' or 'nlminb_TMB'.
  --no-day-effects   Drop the per instrument day effects, i.e. the shared
                     off-diagonal covariances in R. Required by every method
                     other than 'kem'.
  --init-method NAME Method for the fast initialisation stage, which fits a
                     reduced model (B fixed, no day effects) and hands its
                     estimates to the final fit. Default 'TMB'.
  --init-B VALUE     Value B is held at during initialisation, and the value
                     the final fit starts B from. Default 0.975.
  --no-init          Skip the initialisation stage and fit directly.
  --init-from PATH   Start from a saved initialisation fit instead of running
                     the stage again: a run directory, or a fit_init.rds. Its
                     estimates are transferred by name, so the final --method
                     may differ from the one that run used.
  --resume-from PATH Start the final fit from a saved fit: a run directory, or
                     a fit_chunk_NN.rds. Continues a run that hit the chunk
                     limit, or hands a converged fit to another optimiser.
                     Mutually exclusive with --init-from.
  --chunks N         Maximum number of restarted fitting chunks.
  --maxit N          Iterations per chunk, for the selected method.
  --no-scaling       Fit the observations in degrees C rather than scaling them.
  --no-warm-start    Skip the least squares warm start for the seasonal terms.
                     Distinct from --no-init: this is the regression that seeds
                     the seasonal parameters, which feeds whichever fit is first.
  --no-figures       Skip the observation coverage figure.
  --overwrite        Allow writing into a non-empty output directory.
  --quiet            Suppress progress messages.
  --help             Show this message and exit.

Outputs, under the output directory:
  fit_init.rds              the initialisation stage fit, if one ran, or a copy
                            of the one --init-from named
  fit_resumed.rds           a copy of the fit --resume-from started from
  fit_chunk_NN.rds          the fitted model after each chunk
  loglik.csv                stage, method, log likelihood and iteration count
  observations.rds          the observation matrix and its row keys
  scales.rds                the constants that return predictions to degrees C
  dates.rds                 the date range the columns span
  design.rds                the model dimensions and parameterisation
  states.csv, states.rds    reconstructed SST in degrees C, long format
  states_<tag>_only.*       the same as seen by one instrument alone, for each
                            instrument named in the config's `reconstruction`
  coverage.png              which cells of the observation matrix are observed
  config_used.R             a copy of the configuration this run used
  run_log.txt               the console log and sessionInfo()

The reconstruction runs at the end of a fit, controlled by the `reconstruction`
block of the config. R/run_reconstruct.R re-runs it against a finished run
directory without refitting.

Where a run is written is set by the `output` block of the config
(`output_root`, `run_name`); the flags above override it. A relative
`output_root` is taken from the project root, not the working directory.
"

# ---------------------------------------------------------------------------
# Argument parsing. Base R only -- the project has no dependency manager, and a
# handful of flags does not justify adding one.
# ---------------------------------------------------------------------------
parse_args <- function(argv) {

  # `sites`, `start_date` and `end_date` use NA to mean "not supplied, use the
  # config". A real NULL means "no restriction", which is why those are stored
  # through `[` below: `opts$sites <- NULL` would drop the element instead of
  # setting it.
  opts <- list(config = NULL, input = NULL, outdir = NULL, run_name = NULL,
               sites = NA, start_date = NA, end_date = NA,
               structure = NULL, m_factors = NULL, method = NULL,
               chunks = NULL, maxit = NULL,
               init = NA, init_method = NULL, init_B = NULL,
               init_from = NULL, resume_from = NULL,
               scaling = NA, inits = NA, day_effects = NA, figures = TRUE,
               overwrite = FALSE, verbose = TRUE)

  takes_value <- c("--config", "--input", "--outdir", "--run-name", "--sites",
                   "--start-date", "--end-date", "--structure", "--m-factors",
                   "--method", "--chunks", "--maxit", "--init-method",
                   "--init-B", "--init-from", "--resume-from")

  as_count <- function(flag, val) {
    n <- suppressWarnings(as.integer(val))
    if (is.na(n) || n < 1)
      stop(flag, " must be a positive integer, got: ", val, call. = FALSE)
    n
  }
  as_number <- function(flag, val) {
    x <- suppressWarnings(as.numeric(val))
    if (is.na(x) || !is.finite(x))
      stop(flag, " must be a number, got: ", val, call. = FALSE)
    x
  }
  as_date_or_null <- function(flag, val) {
    if (tolower(val) %in% c("none", "null", "")) return(NULL)
    if (is.na(suppressWarnings(as.Date(val))))
      stop(flag, " must be a date as YYYY-MM-DD, or 'none', got: ", val,
           call. = FALSE)
    val
  }

  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]

    # Accept --flag=value as well as --flag value.
    val <- NULL
    if (grepl("^--[^=]+=", a)) {
      val <- sub("^--[^=]+=", "", a)
      a   <- sub("=.*$", "", a)
    } else if (a %in% takes_value) {
      if (i == length(argv))
        stop(a, " requires a value.", call. = FALSE)
      i   <- i + 1L
      val <- argv[i]
    }

    switch(a,
      "--config"     = opts$config    <- val,
      "--input"      = opts$input     <- val,
      "--outdir"     = opts$outdir    <- val,
      "--run-name"   = opts$run_name  <- val,
      "--sites"      = opts["sites"]  <- list(
                          if (tolower(val) %in% c("all", "none", "null", ""))
                            NULL
                          else trimws(strsplit(val, ",", fixed = TRUE)[[1]])),
      "--start-date" = opts["start_date"] <- list(as_date_or_null(a, val)),
      "--end-date"   = opts["end_date"]   <- list(as_date_or_null(a, val)),
      "--structure"  = opts$structure <- val,
      "--m-factors"  = opts$m_factors <- as_count(a, val),
      "--method"     = opts$method    <- val,
      "--chunks"     = opts$chunks    <- as_count(a, val),
      "--maxit"      = opts$maxit     <- as_count(a, val),
      "--init-method" = opts$init_method <- val,
      "--init-B"     = opts$init_B    <- as_number(a, val),
      "--init-from"  = opts$init_from <- val,
      "--resume-from" = opts$resume_from <- val,
      "--no-init"    = opts$init      <- FALSE,
      "--no-scaling" = opts$scaling   <- FALSE,
      "--no-warm-start" = opts$inits  <- FALSE,
      "--no-day-effects" = opts$day_effects <- FALSE,
      "--no-figures" = opts$figures   <- FALSE,
      "--overwrite"  = opts$overwrite <- TRUE,
      "--quiet"      = opts$verbose   <- FALSE,
      "--help"       = { cat(USAGE); quit(save = "no", status = 0) },
      "-h"           = { cat(USAGE); quit(save = "no", status = 0) },
      stop("Unknown option: ", a, "\nRun with --help for the option list.",
           call. = FALSE)
    )
    i <- i + 1L
  }

  # Caught here rather than in validate_model_config() so the message names the
  # flags the user actually typed.
  if (!is.null(opts$init_from) && !is.null(opts$resume_from))
    stop("--init-from and --resume-from both name a fit to start from. Use ",
         "--init-from to replace the initialisation stage, or --resume-from ",
         "to continue a final fit.", call. = FALSE)

  opts
}


# Apply the command line overrides to the config, so that everything downstream
# reads one settled object and the copy saved as config_used.R is not the only
# incomplete record of the run.
apply_overrides <- function(cfg, opts) {

  if (!is.null(opts$input))                cfg$data$input       <- opts$input
  if (!identical(opts$sites, NA))          cfg$data["sites"]      <- list(opts$sites)
  if (!identical(opts$start_date, NA))     cfg$data["start_date"] <- list(opts$start_date)
  if (!identical(opts$end_date, NA))       cfg$data["end_date"]   <- list(opts$end_date)

  if (identical(opts$scaling, FALSE))      cfg$scaling$enabled  <- FALSE
  if (identical(opts$inits, FALSE))        cfg$inits$enabled    <- FALSE

  if (!is.null(opts$structure))            cfg$structure$state_structure <- opts$structure
  if (!is.null(opts$m_factors))            cfg$structure$m_factors       <- opts$m_factors
  if (!is.null(opts$method))               cfg$fitting$method            <- opts$method

  if (identical(opts$day_effects, FALSE))
    for (v in names(cfg$structure$instruments))
      cfg$structure$instruments[[v]]$day_effect <- FALSE

  if (identical(opts$init, FALSE))          cfg$fitting$init$enabled  <- FALSE
  if (!is.null(opts$init_method))           cfg$fitting$init$method   <- opts$init_method
  if (!is.null(opts$init_B))                cfg$fitting$init$B_values <- opts$init_B

  # Resolved against the project root like every other path, so the same command
  # names the same fit whatever directory it is run from.
  if (!is.null(opts$init_from))    cfg$fitting$init$from   <- .from_root(opts$init_from)
  if (!is.null(opts$resume_from))  cfg$fitting$resume_from <- .from_root(opts$resume_from)

  if (!is.null(opts$chunks))               cfg$fitting$chunks   <- opts$chunks
  if (!is.null(opts$maxit)) {
    # Applied only to the selected method's controls; the others are left as
    # the config has them, so config_used.R stays a faithful record.
    m <- cfg$fitting$method %||% "kem"
    cfg$fitting$controls[[m]]$maxit <- opts$maxit
    # minit above maxit is rejected by MARSS, and a --maxit below the config's
    # minit is a request for a short run, not a contradiction.
    if (!is.null(cfg$fitting$controls[[m]]$minit))
      cfg$fitting$controls[[m]]$minit <-
        min(cfg$fitting$controls[[m]]$minit, opts$maxit)
  }

  cfg
}


# `input` may name the marss_inputs.rds itself or the run directory holding it.
resolve_input <- function(input) {
  path <- .from_root(input)
  if (dir.exists(path)) path <- file.path(path, "marss_inputs.rds")
  if (!file.exists(path))
    stop("Observation matrix not found: ", path,
         "\nPass --input with the marss_inputs.rds written by ",
         "R/run_preprocessing.R, or the run directory holding it.",
         call. = FALSE)
  normalizePath(path)
}


# The coverage figure replaces the interactive image(y) calls the exploratory
# scripts relied on. It is a convenience, not part of the fit, so a missing
# plotting package downgrades it to a message rather than failing the run.
write_coverage_plot <- function(md, outdir) {
  ok <- tryCatch({
    # 00_diagnostics.R takes its palette from 00_preprocessing.R, as its header
    # notes, so the two are sourced in that order.
    suppressPackageStartupMessages({
      source(file.path(ROOT, "R", "00_preprocessing.R"))
      source(file.path(ROOT, "R", "00_diagnostics.R"))
    })
    p <- plot_marss_coverage(md$observations)
    if (is.null(p)) return(invisible(FALSE))
    ggplot2::ggsave(file.path(outdir, "coverage.png"), p,
                    width = 9, height = 6, dpi = 150)
    TRUE
  }, error = function(e) {
    message("Skipping the coverage figure: ", conditionMessage(e))
    FALSE
  })
  if (isTRUE(ok)) message("Wrote coverage.png.")
  invisible(ok)
}


main <- function(argv) {

  opts   <- parse_args(argv)
  config <- opts$config %||% DEFAULTS$config
  if (!file.exists(config)) stop("Config file not found: ", config, call. = FALSE)

  # The model library must be loaded before the output location can be resolved,
  # because the location comes from the config and `read_model_config()` lives
  # in 02_model_specification.R. Nothing is logged yet -- there is no run
  # directory to log into until this is settled.
  suppressPackageStartupMessages({
    source(file.path(ROOT, "R", "marss_matrix_functions.R"))
    source(file.path(ROOT, "R", "02_model_specification.R"))
    source(file.path(ROOT, "R", "03_reconstruct_states.R"))
  })

  cfg <- apply_overrides(read_model_config(config), opts)

  # Precedence: --outdir wins outright, then --run-name, then the config's
  # `output` block, then a name derived from the config file and the date.
  run_name <- opts$run_name %||% cfg$output$run_name %||%
    paste0(sub("\\.[^.]*$", "", basename(config)), "_", format(Sys.Date()))

  out_root <- cfg$output$output_root
  out_root <- if (is.null(out_root)) DEFAULTS$outroot else .from_root(out_root)

  outdir <- opts$outdir %||% file.path(out_root, run_name)

  # Refuse to scatter a new run through the artifacts of an old one. The
  # exploratory scripts wrote every structure to the same fixed path, so a
  # second run silently overwrote the first.
  if (dir.exists(outdir) && length(list.files(outdir)) && !opts$overwrite)
    stop("Output directory is not empty: ", outdir,
         "\nPass --overwrite to write into it anyway.", call. = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Copy every message and warning into a log, so a run can be audited from its
  # own outputs, while still letting them reach the console. A `sink()` would
  # capture them but hide the progress of a fit that takes hours.
  log_file <- file.path(outdir, "run_log.txt")
  log_con  <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  to_log <- function(...) cat(..., file = log_con, sep = "")

  tee <- function(expr) withCallingHandlers(
    expr,
    message = function(m) {
      to_log(conditionMessage(m))
      if (!opts$verbose) invokeRestart("muffleMessage")
    },
    warning = function(w) {
      to_log("WARNING: ", conditionMessage(w), "\n")
    })

  tee({
    message("Model run: ", run_name)
    message("  config:  ", config)
    message("  outdir:  ", outdir)
    message("  started: ", format(Sys.time()))

    input <- resolve_input(cfg$data$input)
    message("Reading ", input, " ...")
    marss_inputs <- readRDS(input)

    validate_model_config(cfg, marss_inputs)
    day <- .day_effect_vars(cfg)
    message("Structure: ", cfg$structure$state_structure, ", ",
            cfg$structure$m_factors, " factor(s), day effects: ",
            if (length(day)) paste(day, collapse = ", ") else "none", ".")
    message("Fitting: ",
            if (isTRUE(cfg$fitting$init$enabled))
              paste0(cfg$fitting$init$method, " initialisation, then ") else "",
            cfg$fitting$method %||% "kem", ".")

    md         <- build_model_data(marss_inputs, cfg, verbose = TRUE)
    model_list <- build_marss_model(md, cfg)
    message("Model: ", ncol(model_list$Z), " states, ",
            nrow(model_list$Z), " observation rows.")

    # Written before the fit starts, so an interrupted run still leaves the
    # record needed to interpret whatever chunks it managed to produce.
    file.copy(config, file.path(outdir, "config_used.R"), overwrite = TRUE)
    write_model_outputs(md, cfg, outdir, verbose = TRUE)
    if (opts$figures) write_coverage_plot(md, outdir)

    res <- fit_marss(md, model_list, cfg, outdir, verbose = TRUE)

    # Reconstruction is a post-processing step on a fit that is already safely
    # on disk, so a failure here must not take the run down with it -- a fit can
    # be hours of work, and run_reconstruct.R can retry against the run
    # directory once the cause is fixed.
    tryCatch(
      run_reconstruction(outdir, fit = res$fit, config = cfg, verbose = TRUE),
      error = function(e)
        message("WARNING: reconstruction failed: ", conditionMessage(e),
                "\nThe fit is written; retry with ",
                "Rscript R/run_reconstruct.R --run ", outdir))

    message("Finished: ", format(Sys.time()))
  })

  to_log("\n--- sessionInfo ---\n",
         paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"),
         "\n")

  invisible(NULL)
}


`%||%` <- function(x, y) if (is.null(x)) y else x

if (!interactive()) {
  ok <- tryCatch({ main(commandArgs(trailingOnly = TRUE)); TRUE },
                 error = function(e) {
                   message("ERROR: ", conditionMessage(e))
                   FALSE
                 })
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}
