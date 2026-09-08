#############################################################
#############################################################
###
### This file defines constants used in the data processing
### steps before model fitting. The data input into this 
### program have already been pre processed by the 
### coastal_sst_data package, which has it own parameters.
###
### These parameters are used to filter outliers and other pre-
### processing steps they are documented here as a one-stop-shop
### for understand the assumption that go into the model. 
###
### Jack H. Buckner, Oregon State University, 08/27/2026
###
#############################################################
#############################################################


############################################################
### Grouping variable used for the cross-group diagnostics
### and for the match-up variance summary.
###
### Must name a column present in the input data (e.g.
### "point_id" for a per-site breakdown, or "aoi" for a
### coarser regional one), or be NULL to report pooled
### statistics only. This field is mandatory.
############################################################
group_id <- "point_id"

############################################################
### Define the variable that are needed for the model and
### subsequent processing steps.
###
### These values depend on the data set used to generate the 
### and must match the extraction process in coastal_sst_data
### to work as intended.
############################################################
variables <- list(eco_sst_v002_clean = list(variable = "eco_sst_v002_clean", stat = "nanmean", radius_m = 300),
                  eco_wind_speed_era5 = list(variable = "eco_wind_speed_era5", stat = "nanmean", radius_m = 300),
                  eco_hour_v002 = list(variable = "eco_hour_v002", stat = "nearest", radius_m = 0),
                  lst_sst_clean = list(variable = "lst_sst_clean", stat = "nanmean", radius_m = 300),
                  lst_wind_speed_era5 = list(variable = "lst_wind_speed_era5", stat = "nanmean", radius_m = 300),
                  lst_hour = list(variable = "lst_hour", stat = "nearest", radius_m = 0),
                  modis_sst = list(variable = "modis_sst_aqua", stat = "nanmean", radius_m = 1000),
                  mur_sst = list(variable = "mur_sst", stat = "nanmean", radius_m = 1000),
                  insitu_sst = list(variable = "insitu_sst", stat = "nanmean", radius_m = 300))
value_var <- "value"

value_var_params <- list(variables=variables,value_var=value_var)
############################################################
### A list of variables that need to be converted from
### Kelvin to centigrade
############################################################
variables       <- c("eco_sst_v002_clean","lst_sst_clean","modis_sst","mur_sst")
from        <- "K"     # convert from kelvin 
to          <- "C"     # convert to celsius 
value_var   <-"value"  # column with values
check_range <- TRUE    # check if values look plausible to validate units
strict      <- TRUE    # fail loudly if range is implausible 
quiet       <- FALSE   # don't print summaries for this process

convert_units_params <- list(variables=variables,from=from,to=to,value_var=value_var,
                             check_range=check_range,strict=strict,quiet=quiet)

############################################################
### Variables to filter for physical realism
############################################################
variables       <- c("eco_sst_v002_clean","lst_sst_clean")
lower_threshold <- 0  # lower bound for plausible values (freezing)
upper_threshold <- 40 # upper bound 

physics_filter_params <- list(
  variables        = variables,
  lower_threshold  = lower_threshold,
  upper_threshold  = upper_threshold
)



############################################################
### During low wind conditions the ocean can form and 
### insulated layer at the surface causing TIR estimates of
### ocean temperatures to differ from the bulk ocean temperature.
### to avoid these errors we eliminate low wind observations 
### particularly those made during the day time. 
### min_wind_speed_ms is the low wind speed threshold for 
### keeping these data
### The large skin temperature errors occur at high solar 
### zenith angles and late in the afternoon. 
############################################################
variables <- list(
  list(target = "eco_sst_v002_clean",  wind = "eco_wind_speed_era5", time = "eco_hour_v002"),
  list(target = "lst_sst_clean",  wind = "lst_wind_speed_era5", time = "lst_hour"))

early      = 9 # Beginning of high solar exposure
late       = 18 # End of high solar exposure
threshold  = 2.0 # Low winds speed threshold to filter
local_time = FALSE # if false assume UTC and convert 
lon_var    = "lon" # longitude variable to for time conversion
keep_na    = TRUE # keep missing values 

wind_filter_params <- list(
  variables = variables,
  early = early,
  late = late,
  threshold = threshold,
  local_time = local_time,
  lon_var = lon_var,
  keep_na = keep_na
)

############################################################
### Filtering ECOSTRESS outliers 
### ECOSTRESS has a small number of very extreme high 
### temperature readings that are almost surely spurious. 
### eco_outlier_nsd  defines how large these outliers must be
### in terms of standard deviations to be removed. 
############################################################
variables <- c("eco_sst_v002_clean", "lst_sst_clean")
nsd <- c(3.75,3.3) # standard deviations to flag
harmonics  <-  1 # number of harmonics to capture seasonality 
time_var <- "time" # time variable for trends 
group_var <- "point_id" # add site level variable to trend model 
robust = TRUE # use robust mthod for estimating the standard deviation

outlier_params <- list(
  variables = variables,
  nsd = nsd, 
  harmonics = harmonics,
  time_var = time_var,
  group_var = group_var,
  robust = robust
)



############################################################
### Cold outliers are very likely to be cloud contaminated 
### pixels. We apply a second filter to these by comparing 
### the ecostress and landast images to modis images 
### matchups that fall below the cold outlier threshold are
### removed as likely cloud contaminated pixels.
###
### High temperature outliers might also survive the wind 
### filter, these can be filtered with the optional 
### hot_outlier_nsd
############################################################

variables <- list(eco = list(target = "eco_sst_v002_clean", reference = "modis_sst"),
                  lst = list(target = "lst_sst_clean", reference = "modis_sst"))
lower   = 10#2.58 # remove outliers in the 1% of the distribution below the mean
upper   = 10 #2.58 # remove only extreme hot outliers since there is not clear cause
slope     = 1 # assume a 1:1 relationship between instruments, see ?match_up_filter
intercept = NA  # estimate intercepts for each instrument
passes  = 1 # run one iteration
scale   = "mad" # estimate of the sd that is robust to outlines
center  = "median"  # use median rather tan mean for center to improve cold tails
min_n   = 10 # minimum match-ups required to fit a group

matchup_filter_params <- list(
  variables = variables, lower = lower, upper = upper, slope = slope,
  intercept = intercept, passes = passes, scale = scale, center = center, min_n = min_n
)


############################################################
### Diagnostics
###
### Controls the visual and tabular record of the QC. The
### distribution of every variable is plotted for the full
### sample and broken down by `group_id`, with the values
### removed by each filter shown in red.
###
### When there are more than `max_facets` groups a panel of
### histograms becomes unreadable, so the per-group view
### switches to ordered boxplots plus a companion bar chart
### of the percent removed per group. The tables always carry
### the full per-group breakdown regardless.
############################################################
enabled    <- TRUE   # produce diagnostics at all
bins       <- 100    # histogram bins, matches the per-process plots
max_facets <- 12     # group count above which facets become boxplots
max_points <- 20000  # thinning for the per-group point / box panels
width      <- 9      # png width, inches
height     <- 6      # png height, inches
dpi        <- 150    # png resolution

diagnostics_params <- list(
  enabled = enabled, bins = bins, max_facets = max_facets,
  max_points = max_points, width = width, height = height, dpi = dpi
)


############################################################
### Match-up variance
###
### Summarises the variability between each target and each
### reference around the 1:1 line, on the cleaned data. This
### is the source of the observation error estimates used to
### set the MARSS R matrix, so it is recorded as a table.
###
### `by` controls the breakdown; it defaults to `group_id`.
### Set it to NULL for pooled statistics only.
############################################################
enabled     <- TRUE
targets     <- c("eco_sst_v002_clean","lst_sst_clean")
references  <- c("modis_sst","mur_sst","insitu_sst")
pairs       <- "cross"    # every target against every reference
by          <- group_id   # breakdown column, NULL for pooled only
id_cols     <- NULL       # NULL: every column other than variable / value
valid_range <- NULL       # optional c(lower, upper) gate applied before stats
trim        <- NULL       # optional robust SDs for the extra sd_trim column
min_n       <- 30         # pairs below this return NA rather than erroring

match_up_variance_params <- list(
  enabled = enabled, targets = targets, references = references,
  pairs = pairs, by = by, id_cols = id_cols, valid_range = valid_range,
  trim = trim, min_n = min_n
)


############################################################
### Formatting for the MARSS model.
###
### The cleaned observations become a matrix with time in the
### columns and one row per site and instrument. Rows are
### blocked by variable in the order given below, with sites
### sorted inside each block. The model specification indexes
### rows by that order, so it is part of the interface rather
### than an incidental detail.
###
### A second matrix carries the seasonal basis on the same
### time axis, as a MARSS covariate with 2 * harmonics rows.
### `harmonics` and `period` use the same convention as the
### seasonal outlier screen above, so the season is described
### the same way when screening and when fitting.
############################################################
enabled   <- TRUE
variables <- c("insitu_sst", "lst_sst_clean", "eco_sst_v002_clean",
               "modis_sst", "mur_sst")  # row block order, deliberate
time_step <- 5        # width of each column in days
harmonics <- 2        # sin/cos pairs; the matrix has 2 * this many rows
period    <- 365.25   # length of the seasonal cycle in days
site_var  <- "point_id"
time_var  <- "time"

marss_format_params <- list(
  enabled = enabled, variables = variables, time_step = time_step,
  harmonics = harmonics, period = period,
  site_var = site_var, time_var = time_var
)


############################################################
### Where a run is written.
###
### A run directory holds the cleaned data, every removed
### observation, the diagnostic tables and figures, and a
### copy of this file, so a result can always be traced back
### to the assumptions that produced it.
###
### `output_root` is taken relative to the project root
### unless it is an absolute path, so a run lands in the same
### place regardless of the working directory it was launched
### from. `run_name` names the subdirectory within that root;
### NULL derives it from the input file name and the date.
###
### The command line flags --outdir and --run-name override
### these values.
############################################################
output_root <- "outputs"  # run directories are created here
run_name    <- NULL       # NULL: <input file name>_<date>

output_params <- list(output_root = output_root, run_name = run_name)


preprocess_config <- list(
  group_id          = group_id,
  value_var         = value_var_params,
  convert_units     = convert_units_params,
  physics_filter    = physics_filter_params,
  wind_filter       = wind_filter_params,
  outliers          = outlier_params,
  matchup_filter    = matchup_filter_params,
  diagnostics       = diagnostics_params,
  match_up_variance = match_up_variance_params,
  marss_format      = marss_format_params,
  output            = output_params
)





