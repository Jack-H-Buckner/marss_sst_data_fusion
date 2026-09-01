#!/usr/bin/env Rscript
#############################################################
#############################################################
###
### Knits the reconstructed states from the three example
### runs -- east Storm Bay, west Storm Bay and Macquarie
### Harbour -- into one long format table.
###
###   Rscript examples/combine_states.R
###
### Each run writes its own states.csv, already in long
### format with one row per site and date. This script stacks
### them, labels every row with the run it came from, and
### joins the site metadata in data/raw/site_locations_combined.csv
### so each row also carries the site's name and position.
###
### The run label matters as more than provenance: CB001 is
### fit by both Storm Bay runs, so site alone does not
### identify a series in the combined table. The pair
### (model, site) does.
###
### The runs need not share a time step -- as things stand
### east Storm Bay is daily and the other two are 5-daily,
### having been preprocessed with multi_day_preprocessing.R.
### Stacking them is still the right move, since every row
### carries its own date, but anything that counts rows or
### assumes an even spacing has to do so within a model. The
### step of each run is reported as it is read.
###
### Jack H. Buckner, Oregon State University, 09/01/2026
### Generated with Claude Code
###
#############################################################
#############################################################

# ---------------------------------------------------------------------------
# Locate the project root from this script's own path, so the default paths
# resolve the same way regardless of the working directory it is invoked from.
# Mirrors R/run_reconstruct.R.
# ---------------------------------------------------------------------------
.script_path <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(normalizePath(f[1], mustWork = FALSE))
  NA_character_   # sourced interactively; fall back to the working directory
}

SCRIPT <- .script_path()
ROOT   <- normalizePath(if (is.na(SCRIPT)) getwd() else dirname(dirname(SCRIPT)))


############################################################
### What to combine.
###
### `MODELS` names the subdirectories of examples/models to
### read, in the order they should appear in the output. Each
### holds one directory per run, named <config>_<date>; the
### most recent by name is used, so re-running a model and
### re-running this script is enough to refresh the table.
###
### `STATES_FILE` selects which reconstruction to read.
### "states.csv" is the fused site view -- the site's own
### latent temperature, informed by every instrument that
### observes it. Switch to "states_mur_only.csv" for MUR's
### view of the same sites; the two differ by exactly what
### the higher resolution instruments add.
###
### `SITES_FILE` supplies the names and coordinates. Its `id`
### column is the site key used throughout the model.
############################################################
MODELS      <- c("east_storm_bay", "west_storm_bay", "macquarie_harbor")
MODELS_ROOT <- file.path(ROOT, "examples", "models")
STATES_FILE <- "states.csv"
SITES_FILE  <- file.path(ROOT, "data", "raw", "site_locations_combined.csv")
OUTPUT      <- file.path(MODELS_ROOT, "states_combined.csv")


############################################################
### The most recent run within a model directory.
###
### Run directories are named <config>_<YYYY-MM-DD>, so the
### last in sorted order is the newest. Only directories that
### actually hold the states file are considered, so a run
### that was interrupted before reconstruction is skipped
### rather than silently returning nothing.
############################################################
latest_run <- function(model) {
  model_dir <- file.path(MODELS_ROOT, model)
  if (!dir.exists(model_dir))
    stop("No such model directory: ", model_dir, call. = FALSE)

  runs <- list.dirs(model_dir, recursive = FALSE, full.names = TRUE)
  runs <- runs[file.exists(file.path(runs, STATES_FILE))]
  if (!length(runs))
    stop("No run under ", model_dir, " holds a ", STATES_FILE, call. = FALSE)

  sort(runs)[length(runs)]
}


############################################################
### Read one run's states, tagged with where it came from.
############################################################
read_states <- function(model) {
  run <- latest_run(model)
  x   <- utils::read.csv(file.path(run, STATES_FILE), stringsAsFactors = FALSE)
  x$date <- as.Date(x$date)

  # The modal gap between consecutive dates at one site; a run written on an
  # uneven grid reports the step it mostly uses, which is all this line is for.
  gaps <- diff(sort(unique(x$date)))
  step <- if (length(gaps)) as.integer(names(which.max(table(gaps)))) else NA

  message(sprintf("  %-17s %-33s %6d rows, %d sites, %d-day step, %s to %s",
                  model, basename(run), nrow(x), length(unique(x$site)),
                  step, min(x$date), max(x$date)))

  data.frame(model = model, run = basename(run), x, stringsAsFactors = FALSE)
}

message("Reading reconstructions:")
states <- do.call(rbind, lapply(MODELS, read_states))


############################################################
### Site names and positions.
###
### Kept to the columns that describe where a site is and
### what it is called; the rest of site_locations_combined.csv
### is provenance for the site list itself. `notes` is carried
### because two of the fixed stations -- CB001 and CB002 --
### have a documented coordinate conflict, and a combined
### table that quietly drops that warning invites the reader
### to trust those two positions more than they should.
############################################################
sites <- utils::read.csv(SITES_FILE, stringsAsFactors = FALSE)
sites <- sites[, c("id", "site_name", "site_type", "lat", "lon",
                   "region", "notes")]
names(sites)[names(sites) == "id"] <- "site"

missing <- setdiff(unique(states$site), sites$site)
if (length(missing))
  warning("No location on record for: ", paste(missing, collapse = ", "),
          " -- these rows get NA coordinates.", call. = FALSE)

combined <- merge(states, sites, by = "site", all.x = TRUE, sort = FALSE)


############################################################
### Order the output so it reads as a table rather than a
### merge result: models in the order given above, then site,
### then time.
############################################################
combined$model <- factor(combined$model, levels = MODELS)
combined <- combined[order(combined$model, combined$site, combined$date), ]
combined$model <- as.character(combined$model)

combined <- combined[, c("model", "run", "site", "site_name", "site_type",
                         "region", "lat", "lon", "date", "variable",
                         "value", "se", "notes")]
rownames(combined) <- NULL

utils::write.csv(combined, OUTPUT, row.names = FALSE, na = "")

message(sprintf("\nWrote %s\n  %d rows, %d site-model series, %s to %s",
                sub(paste0("^", ROOT, "/"), "", OUTPUT),
                nrow(combined),
                nrow(unique(combined[, c("model", "site")])),
                min(combined$date), max(combined$date)))
