#' Convert a long-format observation table into a MARSS-ready time series matrix
#'
#' Rows are variable-by-site series, grouped into contiguous blocks by variable
#' (in the order given by `variables`), with sites sorted within each block.
#' Columns are consecutive bins of `time_step` days spanning the full observed
#' range, so the time axis is regular and gaps appear as NA (which is what MARSS
#' expects). A variable-site series with no non-missing value anywhere in the
#' record is dropped entirely.
#'
#' @param data          long-format data frame.
#' @param variables     character vector of variables to include, in the order
#'                      the row blocks should appear. NULL = all present, sorted.
#' @param time_step     width of each column in days. 1 = one column per day.
#' @param value.var,variable.var,site.var,time.var
#'                      column names holding the value, variable name, site id
#'                      and observation date.
#' @param statistic.var column naming the statistic, when the same variable is
#'                      present under several summaries (mean, sd, count, ...).
#' @param statistic     which statistic to pull. Required if `statistic.var` set.
#' @param agg.fun       used only when several observations fall in one bin
#'                      (i.e. `time_step` > 1). Defaults to the mean; use `sum`
#'                      for count-type variables.
#'
#' @return list with ts_matrix, row_keys, row_var_keys, row_site_keys,
#'         row_site_index, plus site_levels, col_dates, time_step and
#'         dropped (the variable-site series that were omitted).
convert_long_format_marss <- function(data, variables = NULL, time_step,
                                      value.var = "value",
                                      variable.var = "variable",
                                      site.var = "site_id",
                                      time.var = "date",
                                      statistic.var = NULL,
                                      statistic = NULL,
                                      agg.fun = function(x) mean(x, na.rm = TRUE)) {
  
  ## ---- 1. validate arguments ---------------------------------------------
  data <- as.data.frame(data)
  
  need <- c(value.var, variable.var, site.var, time.var)
  absent <- need[!need %in% names(data)]
  if (length(absent))
    stop("column(s) not found in `data`: ", paste(absent, collapse = ", "))
  
  if (missing(time_step) || is.null(time_step))
    stop("`time_step` (column width, in days) must be supplied.")
  if (!is.numeric(time_step) || length(time_step) != 1L ||
      is.na(time_step) || time_step < 1)
    stop("`time_step` must be a single number of days >= 1.")
  
  if (!is.null(statistic.var)) {
    if (!statistic.var %in% names(data))
      stop("`statistic.var` column not found in `data`: ", statistic.var)
    if (is.null(statistic))
      stop("`statistic` must be supplied when `statistic.var` is set.")
  } else if (!is.null(statistic)) {
    stop("`statistic` was supplied without `statistic.var`.")
  }
  
  ## ---- 2. subset to the requested statistic and variables ----------------
  if (!is.null(statistic.var)) {
    keep <- as.character(data[[statistic.var]]) %in% as.character(statistic)
    if (!any(keep))
      stop("no rows have ", statistic.var, " equal to ",
           paste(statistic, collapse = "/"), ".")
    data <- data[keep, , drop = FALSE]
  }
  
  vv <- as.character(data[[variable.var]])
  present <- unique(vv)
  if (is.null(variables)) {
    variables <- sort(present)
  } else {
    variables <- as.character(variables)
    gone <- setdiff(variables, present)
    if (length(gone))
      warning("variable(s) requested but not present in `data`: ",
              paste(gone, collapse = ", "), call. = FALSE)
    variables <- intersect(variables, present)   # keeps the caller's order
  }
  if (!length(variables))
    stop("none of the requested variables are present in `data`.")
  
  keep <- vv %in% variables
  data <- data[keep, , drop = FALSE]
  vv   <- vv[keep]
  site <- as.character(data[[site.var]])
  
  ## ---- 3. build the regular time axis ------------------------------------
  tt <- data[[time.var]]
  if (inherits(tt, "POSIXt")) tt <- as.Date(tt)
  if (!inherits(tt, "Date"))  tt <- as.Date(as.character(tt))
  if (anyNA(tt))
    stop(sum(is.na(tt)), " value(s) in `", time.var,
         "` could not be parsed as dates.")
  
  t0    <- min(tt)
  bin   <- as.integer(floor(as.numeric(tt - t0) / time_step))   # 0-based
  n_col <- max(bin) + 1L
  col_dates <- t0 + (seq_len(n_col) - 1L) * time_step           # bin start date
  
  ## ---- 4. warn if several records share one variable/site/date -----------
  if (is.null(statistic.var) &&
      anyDuplicated(data.frame(vv, site, tt, stringsAsFactors = FALSE)) > 0L)
    warning("`data` holds more than one record per variable/site/date; these ",
            "will be combined with `agg.fun`. If the table stores several ",
            "statistics per variable, set `statistic.var` and `statistic`.",
            call. = FALSE)
  
  ## ---- 5. drop missing values --------------------------------------------
  ## Dropping NA here is what implements the omission rule: a variable-site
  ## series with nothing but NA leaves no rows behind and never gets a row.
  val <- suppressWarnings(as.numeric(data[[value.var]]))
  ok  <- !is.na(val)
  val <- val[ok]; vv <- vv[ok]; site <- site[ok]; bin <- bin[ok]
  if (!length(val))
    stop("no non-missing values remain after filtering.")
  
  ## ---- 6. lay out the rows: blocks by variable, sites sorted within ------
  row_var_keys <- row_site_keys <- character(0)
  for (v in variables) {
    s <- sort(unique(site[vv == v]))
    row_var_keys  <- c(row_var_keys,  rep(v, length(s)))
    row_site_keys <- c(row_site_keys, s)
  }
  row_keys <- paste(row_var_keys, row_site_keys, sep = "_")
  n_row    <- length(row_keys)
  
  site_levels    <- sort(unique(site))
  row_site_index <- match(row_site_keys, site_levels)
  
  ## which variable-site pairs were asked for but had no data at all
  all_sites <- sort(unique(as.character(data[[site.var]])))
  wanted    <- expand.grid(variable = variables, site = all_sites,
                           KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  dropped   <- wanted[!paste(wanted$variable, wanted$site, sep = "_") %in% row_keys, ,
                      drop = FALSE]
  rownames(dropped) <- NULL
  
  ## ---- 7. fill the matrix -------------------------------------------------
  mat <- matrix(NA_real_, nrow = n_row, ncol = n_col,
                dimnames = list(row_keys, format(col_dates, "%Y-%m-%d")))
  
  r      <- match(paste(vv, site, sep = "_"), row_keys)
  cellid <- bin * n_row + r          # column-major linear index into `mat`
  
  if (anyDuplicated(cellid) == 0L) {
    mat[cellid] <- val               # one observation per cell, no aggregation
  } else {
    cells <- sort(unique(cellid))
    mat[cells] <- as.numeric(tapply(val, cellid, agg.fun))
  }
  mat[is.nan(mat)] <- NA_real_
  
  ## ---- 8. return ----------------------------------------------------------
  out <- list(ts_matrix      = mat,
              row_keys       = row_keys,
              row_var_keys   = row_var_keys,
              row_site_keys  = row_site_keys,
              row_site_index = row_site_index,
              site_levels    = site_levels,
              col_dates      = col_dates,
              time_step      = time_step,
              dropped        = dropped)
  return(out)
}



#' Seasonal harmonic basis on a set of dates
#'
#' Builds the sin/cos pairs used to represent an annual cycle, as a matrix laid
#' out the way MARSS wants a covariate: one row per harmonic term, one column
#' per time step.
#'
#' The basis is deliberately identical to the one `.seasonal_resid()` uses to
#' screen seasonal outliers in `R/00_preprocessing.R` -- same `2 * pi * k * doy
#' / period` terms, same `sin1, cos1, sin2, ...` naming. The screening step and
#' the fitting step then describe the season the same way, rather than one using
#' `6.28 / 365` and the other `2 * pi / 365.25`.
#'
#' @param dates Vector of dates, or anything `as.Date()` accepts. Normally the
#'   `col_dates` of the observation matrix, so the two share a time axis.
#' @param harmonics Number of sin/cos pairs. The result has `2 * harmonics` rows.
#' @param period Length of the cycle in days. Defaults to 365.25.
#'
#' @return A `2 * harmonics` by `length(dates)` numeric matrix, with rownames
#'   `sin1, cos1, sin2, cos2, ...` and the dates as colnames.
#' @importFrom lubridate yday
#' @export
seasonal_harmonics <- function(dates, harmonics = 2, period = 365.25) {

  if (!is.numeric(harmonics) || length(harmonics) != 1 || !is.finite(harmonics) ||
      harmonics < 1 || harmonics != round(harmonics))
    stop("`harmonics` must be a single whole number >= 1.")
  if (!is.numeric(period) || length(period) != 1 || !is.finite(period) ||
      period <= 0)
    stop("`period` must be a single positive number of days.")

  d <- dates
  if (inherits(d, "POSIXt")) d <- as.Date(d)
  if (!inherits(d, "Date"))  d <- as.Date(as.character(d))
  if (!length(d)) stop("`dates` is empty.")
  if (anyNA(d)) stop(sum(is.na(d)), " date(s) could not be parsed.")

  if (2 * harmonics >= length(d))
    warning("`harmonics` gives ", 2 * harmonics, " terms for only ", length(d),
            " time step(s); the basis cannot be identified.", call. = FALSE)

  doy <- lubridate::yday(d)

  X <- do.call(rbind, lapply(seq_len(harmonics), function(k) {
    rbind(sin(2 * pi * k * doy / period),
          cos(2 * pi * k * doy / period))
  }))

  dimnames(X) <- list(paste0(c("sin", "cos"), rep(seq_len(harmonics), each = 2)),
                      format(d, "%Y-%m-%d"))

  # MARSS rejects missing values in a covariate matrix. Dates are complete by
  # construction here, so this only guards against a pathological `period`.
  if (any(!is.finite(X))) {
    warning(sum(!is.finite(X)), " non-finite harmonic value(s) set to 0.",
            call. = FALSE)
    X[!is.finite(X)] <- 0
  }
  X
}


#' Build the matrices a MARSS fit needs from cleaned long format data
#'
#' Composes the two pieces the model step consumes: the observation matrix from
#' `convert_long_format_marss()`, and the seasonal covariate matrix from
#' `seasonal_harmonics()` on that matrix's own time axis, so the two are aligned
#' column for column.
#'
#' Rows of the observation matrix are blocked by variable in the order given by
#' `variables`, with sites sorted inside each block. That order is what the model
#' specification indexes by, so it is part of the interface rather than an
#' incidental detail.
#'
#' @param data Cleaned long format data, e.g. `run_pipeline()$data`.
#' @param variables Variables to include, in the order the row blocks should
#'   appear.
#' @param time_step Width of each column in days. 1 gives one column per day.
#' @param harmonics,period Passed to `seasonal_harmonics()`.
#' @param site_var,time_var,value_var,variable_var Column names in `data`. The
#'   defaults match what the preprocessing pipeline produces, which is *not*
#'   what `convert_long_format_marss()` defaults to.
#' @param agg_fun Used only when several observations land in one bin.
#'
#' @return The full `convert_long_format_marss()` list -- `ts_matrix`,
#'   `row_keys`, `row_var_keys`, `row_site_keys`, `row_site_index`,
#'   `site_levels`, `col_dates`, `time_step`, `dropped` -- with `harmonics` (the
#'   covariate matrix), `row_metadata`, `n_harmonics` and `period` added.
#' @export
build_marss_inputs <- function(data,
                               variables,
                               time_step    = 1,
                               harmonics    = 2,
                               period       = 365.25,
                               site_var     = "point_id",
                               time_var     = "time",
                               value_var    = "value",
                               variable_var = "variable",
                               agg_fun      = function(x) mean(x, na.rm = TRUE)) {

  stopifnot(is.data.frame(data))
  if (!is.character(variables) || !length(variables))
    stop("`variables` must be a non-empty character vector.")

  out <- convert_long_format_marss(
    data, variables = variables, time_step = time_step,
    value.var = value_var, variable.var = variable_var,
    site.var = site_var, time.var = time_var, agg.fun = agg_fun)

  out$harmonics   <- seasonal_harmonics(out$col_dates, harmonics, period)
  out$n_harmonics <- harmonics
  out$period      <- period

  out$row_metadata <- data.frame(
    row        = seq_along(out$row_keys),
    row_key    = out$row_keys,
    variable   = out$row_var_keys,
    site       = out$row_site_keys,
    site_index = out$row_site_index,
    stringsAsFactors = FALSE)

  # The model specification indexes rows and columns by position, so a mismatch
  # here would be silent and wrong rather than an error.
  stopifnot(
    ncol(out$ts_matrix) == length(out$col_dates),
    ncol(out$harmonics) == ncol(out$ts_matrix),
    nrow(out$harmonics) == 2 * harmonics,
    nrow(out$ts_matrix) == nrow(out$row_metadata),
    identical(rownames(out$ts_matrix), out$row_metadata$row_key))

  out
}

