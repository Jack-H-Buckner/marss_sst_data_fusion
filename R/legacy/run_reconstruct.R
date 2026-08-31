#!/usr/bin/env Rscript
#############################################################
#############################################################
###
### Command line entry point for the reconstruction step.
###
###   Rscript R/run_reconstruct.R --run models/v2_bfgs
###
### Reads a finished model run directory and writes the
### reconstructed temperatures it implies, in degrees C and
### in long format. This is the same step run_model.R runs at
### the end of a fit; it exists on its own so a run can be
### reconstructed without being refitted, which matters when
### the fit took hours and only the reconstruction settings
### have changed.
###
### Run with --help for the full option list.
###
### Jack H. Buckner, Oregon State University, 08/31/2026
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

.is_absolute <- function(path) grepl("^(/|~|[A-Za-z]:)", path)
.from_root   <- function(path) if (.is_absolute(path)) path else file.path(ROOT, path)

USAGE <- "
Usage: Rscript R/run_reconstruct.R --run PATH [options]

Reconstructs sea surface temperature from a fitted model and writes it in long
format, one row per site and date, in degrees C.

The site view -- the fused reconstruction, informed by every instrument that
observes a site -- is always written. An instrument view is what the model says
that instrument alone sees at each site, using its own intercept and seasonal
cycle and, when it has `site_state = FALSE`, only the shared factors.

Options:
  --run PATH         The model run directory to reconstruct. Required. It must
                     hold scales.rds, observations.rds, design.rds and a fit.
  --fit NAME         Which fit within the run directory to use.
                     Default: the highest numbered fit_chunk_NN.rds, or
                     fit_init.rds if the run has no chunk.
  --instruments LIST Comma separated variables to write an instrument view for,
                     overriding the run's configuration. Pass 'none' for the
                     site view alone.
  --formats LIST     Comma separated output formats: csv, rds, or both.
  --outdir PATH      Where to write. Default: the run directory itself.
  --config PATH      Configuration to take the `reconstruction` block from,
                     instead of the copy stored in the run's design.rds.
  --quiet            Suppress progress messages.
  --help             Show this message and exit.

Outputs, under the output directory:
  states.csv, states.rds    the fused site level reconstruction
  states_<tag>_only.*       one instrument's view, for each instrument
                            requested; MUR lands as states_mur_only.csv

Columns: date, site, variable, value, se. `se` is state uncertainty at the
parameter estimates -- it excludes parameter uncertainty, and propagates the
states as if they were independent. See R/03_reconstruct_states.R.
"

# ---------------------------------------------------------------------------
# Argument parsing. Base R only, mirroring run_model.R.
# ---------------------------------------------------------------------------
parse_args <- function(argv) {

  # `instruments` uses NA to mean "not supplied, use the run's configuration".
  # A real NULL means "the site view alone", which is why it is stored through
  # `[` below: `opts$instruments <- NULL` would drop the element instead.
  opts <- list(run = NULL, fit = NULL, instruments = NA, formats = NULL,
               outdir = NULL, config = NULL, verbose = TRUE)

  takes_value <- c("--run", "--fit", "--instruments", "--formats", "--outdir",
                   "--config")

  as_list <- function(val) trimws(strsplit(val, ",", fixed = TRUE)[[1]])

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
      "--run"    = opts$run    <- val,
      "--fit"    = opts$fit    <- val,
      "--instruments" = opts["instruments"] <- list(
                          if (tolower(val) %in% c("none", "null", ""))
                            NULL
                          else as_list(val)),
      "--formats" = opts$formats <- as_list(val),
      "--outdir" = opts$outdir <- val,
      "--config" = opts$config <- val,
      "--quiet"  = opts$verbose <- FALSE,
      "--help"   = { cat(USAGE); quit(save = "no", status = 0) },
      "-h"       = { cat(USAGE); quit(save = "no", status = 0) },
      stop("Unknown option: ", a, "\nRun with --help for the option list.",
           call. = FALSE)
    )
    i <- i + 1L
  }

  if (is.null(opts$run))
    stop("--run is required: the model run directory to reconstruct.\n",
         "Run with --help for the option list.", call. = FALSE)

  opts
}


main <- function(argv) {

  opts <- parse_args(argv)

  suppressPackageStartupMessages({
    source(file.path(ROOT, "R", "marss_matrix_functions.R"))
    source(file.path(ROOT, "R", "02_model_specification.R"))
    source(file.path(ROOT, "R", "03_reconstruct_states.R"))
  })

  run_dir <- .from_root(opts$run)

  # The settings come from the run's own design.rds unless a config is named;
  # the flags then override whichever of those applies, so a run can be
  # reconstructed differently without editing anything on disk.
  cfg <- if (is.null(opts$config)) NULL else
    read_model_config(.from_root(opts$config))

  if (opts$verbose) message("Reconstructing ", run_dir)

  run_reconstruction(run_dir, fit_file = opts$fit, config = cfg,
                     instruments = opts$instruments,
                     formats     = opts$formats,
                     # Not a flag: asking for a reconstruction by running this
                     # script is a clearer statement of intent than whatever
                     # the run's config defaulted to.
                     enabled     = TRUE,
                     outdir = if (is.null(opts$outdir)) NULL else
                       .from_root(opts$outdir),
                     verbose = opts$verbose)

  invisible(NULL)
}


if (!interactive()) {
  ok <- tryCatch({ main(commandArgs(trailingOnly = TRUE)); TRUE },
                 error = function(e) {
                   message("ERROR: ", conditionMessage(e))
                   FALSE
                 })
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}
