#!/usr/bin/env Rscript
#############################################################
#############################################################
###
### Command line entry point for the preprocessing step.
###
###   Rscript R/run_preprocessing.R --input data/raw/x.csv
###
### Reads the raw long format extract, applies the quality
### control configured in parameters/preprocessing.R, and
### writes the cleaned data alongside a complete visual and
### tabular record of what was removed and why.
###
### Run with --help for the full option list.
###
### Jack H. Buckner, Oregon State University
### Generated with Claude Code
### Reviewed JHB 09/02/2026
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

# A relative `output_root` is resolved against the project root, not the caller's
# working directory -- otherwise the same config would write to a different place
# depending on where the run was launched from.
.is_absolute <- function(path) grepl("^(/|~|[A-Za-z]:)", path)

# Used only when the config carries no `output` block.
DEFAULTS <- list(
  input   = file.path(ROOT, "data", "raw", "tasi_salmon_20_yrs.csv"),
  config  = file.path(ROOT, "parameters", "preprocessing.R"),
  outroot = file.path(ROOT, "outputs")
)

USAGE <- "
Usage: Rscript R/run_preprocessing.R [options]

Applies the configured preprocessing QC to a raw long format extract and writes
the cleaned data, the removed observations, diagnostic tables and figures.

Options:
  --input PATH       Input CSV.
                     Default: data/raw/tasi_salmon_20_yrs.csv
  --config PATH      Config file defining `preprocess_config`.
                     Default: parameters/preprocessing.R
  --outdir PATH      Output directory, used verbatim. Overrides both
                     --run-name and the config's `output` block.
                     Default: <output_root>/<run-name>
  --run-name NAME    Name of the run subdirectory within the config's
                     `output_root`. Overrides the config's `run_name`.
                     Default: the input file stem plus today's date.
  --group-id NAME    Override the config `group_id`. Pass 'none' for pooled
                     statistics only. Also applies to the match-up variance
                     breakdown, unless the config set its `by` independently.
  --no-diagnostics   Skip the diagnostic tables and figures.
  --no-variance      Skip the match-up variance summary.
  --no-marss         Skip building the MARSS observation and harmonics
                     matrices.
  --overwrite        Allow writing into a non-empty output directory.
  --quiet            Suppress progress messages.
  --help             Show this message and exit.

Outputs, under the output directory:
  clean_data.rds            the cleaned long format data
  removed_all.csv           every removed observation, tagged with its stage
  tables/*.csv              per-variable, per-group and per-stage summaries,
                            including match_up_variance.csv
  diagnostics/*.png         one figure per screening step and per variable
  diagnostics_report.pdf    every figure in one document
  marss_inputs.rds          the MARSS observation matrix, the seasonal
                            harmonics matrix and the row metadata
  config_used.R             a copy of the configuration this run used
  run_log.txt               the console log and sessionInfo()

Where a run is written is set by the `output` block of the config
(`output_root`, `run_name`); the flags above override it. A relative
`output_root` is taken from the project root, not the working directory.
"

# ---------------------------------------------------------------------------
# Argument parsing. Base R only -- the project has no dependency manager, and a
# handful of flags does not justify adding one.
# ---------------------------------------------------------------------------
parse_args <- function(argv) {

  # `group_id = NA` means "not supplied, use the config". A real NULL means
  # "pooled diagnostics only", which is why it is stored through `[` below:
  # `opts$group_id <- NULL` would drop the element instead of setting it.
  opts <- list(input = NULL, config = NULL, outdir = NULL, run_name = NULL,
               group_id = NA, diagnostics = TRUE, variance = TRUE,
               marss = TRUE, overwrite = FALSE, verbose = TRUE)

  takes_value <- c("--input", "--config", "--outdir", "--run-name", "--group-id")
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
      "--input"          = opts$input    <- val,
      "--config"         = opts$config   <- val,
      "--outdir"         = opts$outdir   <- val,
      "--run-name"       = opts$run_name <- val,
      "--group-id"       = opts["group_id"] <- list(
                             if (tolower(val) %in% c("none", "null", ""))
                               NULL else val),
      "--no-diagnostics" = opts$diagnostics <- FALSE,
      "--no-variance"    = opts$variance    <- FALSE,
      "--no-marss"       = opts$marss       <- FALSE,
      "--overwrite"      = opts$overwrite   <- TRUE,
      "--quiet"          = opts$verbose     <- FALSE,
      "--help"           = { cat(USAGE); quit(save = "no", status = 0) },
      "-h"               = { cat(USAGE); quit(save = "no", status = 0) },
      stop("Unknown option: ", a, "\nRun with --help for the option list.",
           call. = FALSE)
    )
    i <- i + 1L
  }
  opts
}


main <- function(argv) {

  opts <- parse_args(argv)

  input  <- opts$input  %||% DEFAULTS$input
  config <- opts$config %||% DEFAULTS$config

  if (!file.exists(input))  stop("Input file not found: ", input, call. = FALSE)
  if (!file.exists(config)) stop("Config file not found: ", config, call. = FALSE)

  # The pipeline must be loaded before the output location can be resolved,
  # because the location comes from the config and `read_config()` lives in
  # 00_preprocessing.R. Nothing is logged yet -- there is no run directory to
  # log into until this is settled.
  suppressPackageStartupMessages({
    source(file.path(ROOT, "R", "00_preprocessing.R"))
    source(file.path(ROOT, "R", "00_diagnostics.R"))
    source(file.path(ROOT, "R", "01_format_data_for_marss.R"))
  })
  cfg <- read_config(config)

  # Precedence: --outdir wins outright, then --run-name, then the config's
  # `output` block, then a name derived from the input file and the date.
  run_name <- opts$run_name %||% cfg$output$run_name %||%
    paste0(sub("\\.[^.]*$", "", basename(input)), "_", format(Sys.Date()))

  out_root <- cfg$output$output_root
  out_root <- if (is.null(out_root)) DEFAULTS$outroot else
                if (.is_absolute(out_root)) out_root else
                  file.path(ROOT, out_root)

  outdir <- opts$outdir %||% file.path(out_root, run_name)

  # Refuse to scatter a new run through the artifacts of an old one.
  if (dir.exists(outdir) && length(list.files(outdir)) && !opts$overwrite)
    stop("Output directory is not empty: ", outdir,
         "\nPass --overwrite to write into it anyway.", call. = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Copy every message and warning into a log, so a run can be audited from its
  # own outputs, while still letting them reach the console. A `sink()` would
  # capture them but hide the progress of a run that takes minutes.
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
    message("Preprocessing run: ", run_name)
    message("  input:   ", input)
    message("  config:  ", config)
    message("  outdir:  ", outdir)
    message("  started: ", format(Sys.time()))

    message("Reading ", input, " ...")
    dat <- utils::read.csv(input, stringsAsFactors = FALSE)
    message("  ", format(nrow(dat), big.mark = ","), " rows, ",
            ncol(dat), " columns.")

    # `cfg` is already read; passing the list avoids sourcing the config twice.
    result <- run_pipeline(dat, cfg,
                           group_id    = opts$group_id,
                           diagnostics = opts$diagnostics,
                           variance    = opts$variance,
                           marss       = opts$marss,
                           verbose     = TRUE)

    file.copy(config, file.path(outdir, "config_used.R"), overwrite = TRUE)
    write_pipeline_outputs(result, outdir, result$config, verbose = TRUE)

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
