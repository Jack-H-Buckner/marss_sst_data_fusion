library(dplyr)
library(reshape2)
library(ggplot2)

# The configuration is not sourced here. It arrives through the `config`
# argument of `run_pipeline()`, so this file can be sourced from any working
# directory and leaves nothing in the global environment.

#' Fill colours for retained and removed observations
#'
#' Every screening diagnostic in this project uses the same encoding: retained
#' observations in grey, removed observations in firebrick. Kept as one constant
#' so the per-process plots and the whole-dataset plots in `00_diagnostics.R`
#' cannot drift apart.
#' @export
RETAINED_PAL <- c(`TRUE` = "grey40", `FALSE` = "firebrick")

#' Collapse variable/stat column pairs into a single variable column
#'
#' Selects rows of `data` matching each specification in `vars`, labels them with
#' the corresponding name from `vars`, drops the columns used for matching, and
#' renames `value_var` to `value`.
#'
#' @param data A long format data frame.
#' @param vars A named list of named lists. Each inner list gives the column
#'   values identifying one variable, e.g.
#'   `list(var_1 = list(variable = "lst_sst", stat = "nanmean"),
#'         var_2 = list(variable = "lst_wind", stat = "nearest"))`.
#' @param value_var Character scalar naming the column holding the values.
#'
#' @return A data frame with `variable` and `value` columns replacing the
#'   matching columns and `value_var`.
#' @export
get_variable_value_format <- function(data, vars, value_var) {
  stopifnot(
    is.data.frame(data),
    is.list(vars), length(vars) > 0,
    !is.null(names(vars)), all(nzchar(names(vars))),
    is.character(value_var), length(value_var) == 1
  )
  
  key_cols <- unique(unlist(lapply(vars, names)))
  missing_cols <- setdiff(c(key_cols, value_var), names(data))
  if (length(missing_cols)) {
    stop("Columns not found in `data`: ", paste(missing_cols, collapse = ", "))
  }
  
  keep     <- rep(FALSE, nrow(data))
  var_name <- rep(NA_character_, nrow(data))
  
  for (i in seq_along(vars)) {
    spec <- vars[[i]]
    keep_var <- rep(TRUE, nrow(data))
    for (nm in names(spec)) {
      keep_var <- keep_var & (data[[nm]] %in% spec[[nm]])
    }
    if (any(keep & keep_var)) {
      warning("Rows match more than one entry of `vars`; later entries win.")
    }
    var_name[keep_var] <- names(vars)[i]
    keep <- keep | keep_var
  }
  
  if (!any(keep)) warning("No rows matched any entry of `vars`.")
  
  data_vars <- data[keep, setdiff(names(data), key_cols), drop = FALSE]
  data_vars$variable <- var_name[keep]
  data_vars$value    <- data[[value_var]][keep]
  if (!identical(value_var, "value")) data_vars[[value_var]] <- NULL
  
  rownames(data_vars) <- NULL
  data_vars
}
  

#' Convert temperature variables between units in a long format data frame
#'
#' Converts the `value` entries of named temperature variables in place. Rows are
#' modified where they sit rather than being split out and re-bound, so no rows
#' can be lost or reordered and columns other than the value column are untouched.
#'
#' Expects data in variable/value long format, e.g. the output of
#' `get_variable_value_format()`. Variables are selected by name alone.
#'
#' A plausibility check on the input range guards against the most damaging
#' mistake with unit conversion, which is running it twice. Sea surface
#' temperatures near 288 K become 15 C on the first pass and -258 C on the
#' second; the check catches this before the data is corrupted.
#'
#' @param data Long format data frame with a `variable` column and a value column.
#' @param variables Character vector of variable names to convert. Names not
#'   present in the data raise a warning and are skipped.
#' @param from,to Units, one of "K", "C" or "F". Default "K" to "C".
#' @param value_var Name of the value column. Defaults to "value".
#' @param check_range If TRUE (default), verify the values look plausible in
#'   `from` units before converting, and in `to` units afterwards.
#' @param strict If TRUE (default), a failed range check is an error. If FALSE it
#'   is a warning and the conversion proceeds.
#' @param quiet If TRUE, suppress the per-variable summary message.
#'
#' @return `data` with the selected values converted, in the original row order
#'   and with all other columns unchanged.
#' @export
convert_temperature <- function(data,
                                variables,
                                from        = "K",
                                to          = "C",
                                value_var   = "value",
                                check_range = TRUE,
                                strict      = TRUE,
                                quiet       = FALSE) {
  
  stopifnot(is.data.frame(data), is.character(variables), length(variables) > 0,
            is.logical(check_range), is.logical(strict), is.logical(quiet))
  
  from <- toupper(as.character(from)); to <- toupper(as.character(to))
  ok_u <- c("K", "C", "F")
  if (length(from) != 1 || !from %in% ok_u) stop("`from` must be one of K, C, F.")
  if (length(to)   != 1 || !to   %in% ok_u) stop("`to` must be one of K, C, F.")
  if (identical(from, to)) {
    if (!quiet) message("`from` and `to` are both ", from, "; returning data unchanged.")
    return(data)
  }
  
  if (!"variable" %in% names(data)) stop("`data` has no `variable` column.")
  if (!value_var %in% names(data)) stop("Value column not found: ", value_var)
  if (!is.numeric(data[[value_var]]))
    stop("`", value_var, "` is not numeric; cannot convert.")
  
  dup <- unique(variables[duplicated(variables)])
  if (length(dup)) {
    warning("Duplicated entries in `variables`: ", paste(dup, collapse = ", "))
    variables <- unique(variables)
  }
  absent <- setdiff(variables, unique(data$variable))
  if (length(absent)) {
    warning("Not present in `data$variable`, skipped: ", paste(absent, collapse = ", "))
    variables <- setdiff(variables, absent)
  }
  if (!length(variables)) {
    if (!quiet) message("Nothing to convert; returning data unchanged.")
    return(data)
  }
  
  # `%in%` rather than `==` so a missing `variable` excludes the row from
  # conversion instead of propagating NA into the mask.
  sel   <- data$variable %in% variables
  vals  <- data[[value_var]][sel]
  n_val <- sum(!is.na(vals))
  if (n_val == 0) {
    if (!quiet) message("All selected values are NA; returning data unchanged.")
    return(data)
  }
  
  # ---- plausibility --------------------------------------------------------
  bands <- list(K = c(150, 400), C = c(-125, 125), F = c(-190, 260))
  
  if (check_range) {
    b  <- bands[[from]]
    md <- stats::median(vals, na.rm = TRUE)
    if (md < b[1] || md > b[2]) {
      guess <- names(bands)[vapply(bands, function(z) md >= z[1] && md <= z[2],
                                   logical(1))]
      msg <- paste0("Values do not look like ", from, ": median is ",
                    signif(md, 6), ", expected roughly ", b[1], " to ", b[2], ".",
                    if (length(guess)) paste0(" They look like ", guess[1],
                                              " already - has this already been converted?"))
      if (strict) stop(msg, "\n  Set strict = FALSE to convert anyway.") else warning(msg)
    }
  }
  
  # ---- convert -------------------------------------------------------------
  cv <- function(x, from, to) {
    k <- switch(from, K = x, C = x + 273.15, F = (x - 32) * 5/9 + 273.15)
    switch(to,       K = k, C = k - 273.15, F = (k - 273.15) * 9/5 + 32)
  }
  
  before <- vals
  data[[value_var]][sel] <- cv(vals, from, to)
  after <- data[[value_var]][sel]
  
  if (check_range) {
    b  <- bands[[to]]
    md <- stats::median(after, na.rm = TRUE)
    if (md < b[1] || md > b[2]) {
      msg <- paste0("Converted values are implausible as ", to, ": median ",
                    signif(md, 6), ".")
      if (strict) stop(msg) else warning(msg)
    }
  }
  
  if (!quiet) {
    sel_var <- data$variable[sel]
    per <- do.call(rbind, lapply(variables, function(v) {
      x <- before[sel_var == v]; y <- after[sel_var == v]
      if (!sum(!is.na(x))) return(NULL)
      data.frame(variable = v, n = sum(!is.na(x)),
                 from_min = round(min(x, na.rm = TRUE), 2),
                 from_max = round(max(x, na.rm = TRUE), 2),
                 to_min   = round(min(y, na.rm = TRUE), 2),
                 to_max   = round(max(y, na.rm = TRUE), 2))
    }))
    message("Converted ", from, " -> ", to, " (", format(n_val, big.mark = ","),
            " values):")
    message(paste(utils::capture.output(print(per, row.names = FALSE)),
                  collapse = "\n"))
  }
  
  data
}


#' Convert temperature variables from kelvin to celsius
#'
#' Thin wrapper on [convert_temperature()].
#'
#' @inheritParams convert_temperature
#' @return `data` with the selected values converted to celsius.
#' @export
convert_from_K_to_C <- function(data, variables, ...) {
  convert_temperature(data, variables, from = "K", to = "C", ...)
}



#' Filter extreme outliers from an ECOSTRESS time series
#'
#' Remove physically implausible value from the emperature data set. 
#'
#' @param data Long format data frame, e.g. output of `get_variable_value_format()`.
#' @param var String naming the variable (in the `variable` column) to filter.
#' @param lower_threshold lower bound for admissible data
#' @param upper_threshold upper bound for admissible data
#' @param keep_na If TRUE (default), rows of `var` with a missing value are passed
#'   through rather than dropped. A missing value is untestable, not implausible,
#'   and dropping it here would silently remove rows that no other stage removes.
#' @return A list with `data` (full data set, implausible values removed),
#'   `removed` (the dropped rows, with the bounds that rejected them and a
#'   `reason`), `plot` (histogram of the screened values, removed in red), and
#'   `summary` (counts by outcome).
trim_implausibe_values <- function(data, var, lower_threshold = 0,
                                   upper_threshold = 40, keep_na = TRUE) {

  stopifnot(is.data.frame(data),
            is.character(var), length(var) == 1,
            is.numeric(lower_threshold), length(lower_threshold) == 1,
            is.numeric(upper_threshold), length(upper_threshold) == 1,
            lower_threshold < upper_threshold,
            is.logical(keep_na), length(keep_na) == 1)

  need <- c("variable", "value")
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "))


  data <- as.data.frame(data)
  ord  <- seq_len(nrow(data))
  
  # Split target variable off, without NSE
  is_target <- !is.na(data$variable) & data$variable == var
  if (!any(is_target)) stop("No rows in `data` with variable == '", var, "'.")
  
  dat_var   <- data[is_target, , drop = FALSE]
  ord_var   <- ord[is_target]
  data_rest <- data[!is_target, , drop = FALSE]
  ord_rest  <- ord[!is_target]
  
  # Rows the model cannot use at all
  unusable <- is.na(dat_var$value) 
  
  dat_pass <- dat_var[!unusable, , drop = FALSE]
  ord_pass <- ord_var[!unusable]
  dat_skip <- dat_var[unusable, , drop = FALSE]
  ord_skip <- ord_var[unusable]
  
  below   <- dat_pass$value <= lower_threshold
  above   <- dat_pass$value >= upper_threshold
  keep    <- !below & !above

  removed <- dat_pass[!keep, , drop = FALSE]
  if (nrow(removed)) {
    removed$lower_threshold <- lower_threshold
    removed$upper_threshold <- upper_threshold
    removed$reason <- ifelse(below[!keep],
                             "below lower threshold", "above upper threshold")
  }

  plot_dat <- data.frame(value = dat_pass$value, retained = keep)

  kept_idx <- c(ord_rest, ord_pass[keep], if (keep_na) ord_skip)
  if (!keep_na && nrow(dat_skip)) {
    message("Dropped ", nrow(dat_skip), " row(s) of '", var,
            "' with missing values.")
  }
  data_out <- data[sort(kept_idx), , drop = FALSE]
  rownames(data_out) <- rownames(removed) <- NULL

  smry <- data.frame(
    reason = c("retained", "missing value (untestable)",
               "below lower threshold", "above upper threshold"),
    n = c(sum(keep), nrow(dat_skip), sum(below), sum(above)),
    row.names = NULL)
  smry$pct_of_var <- round(100 * smry$n / max(nrow(dat_var), 1), 2)
  attr(smry, "n_var")     <- nrow(dat_var)
  attr(smry, "n_tested")  <- nrow(dat_pass)
  attr(smry, "n_removed") <- sum(!keep) + if (keep_na) 0L else nrow(dat_skip)

  plt <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = value, fill = retained)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::scale_fill_manual(values = RETAINED_PAL, name = "Retained") +
    ggplot2::geom_vline(xintercept = c(lower_threshold, upper_threshold),
                        linetype = 2, colour = "grey40") +
    ggplot2::labs(x = var, y = "Count",
                  title = paste0("Physical range screening: ", var),
                  subtitle = sprintf("Bounds: %g to %g",
                                     lower_threshold, upper_threshold)) +
    ggplot2::theme_classic()

  list(data = data_out, removed = removed, plot = plt, summary = smry)
}



#' Filter extreme outliers from an ECOSTRESS time series
#'
#' Fits a harmonic seasonal trend with a per-location intercept, standardises the
#' residuals, and removes observations whose scaled residual exceeds a threshold.
#' Runs one pass per element of `thresholds`, refitting the trend on the surviving
#' data each time. The default two-pass schedule uses a loose first threshold to
#' strip gross outliers that distort the fit, then the target threshold.
#'
#' @param data Long format data frame, e.g. output of `get_variable_value_format()`.
#' @param var String naming the variable (in the `variable` column) to filter.
#' @param time_var Name of the time column. Defaults to "time".
#' @param group_var Name of the location column. Defaults to "point_id".
#' @param thresholds Numeric vector of thresholds, in standard deviations, applied
#'   in order. Defaults to `c(5, 2.58)`.
#' @param harmonics Number of sin/cos pairs in the seasonal fit. Defaults to 1.
#' @param robust If TRUE, standardise residuals by median/MAD instead of mean/SD.
#' @param keep_na If TRUE, rows of `var` with missing values are passed through
#'   rather than dropped. Defaults to TRUE.
#' @param period Length of the seasonal cycle in days. Defaults to 365.25.
#'
#' @return A list with `data` (full data set, outliers removed), `removed` (the
#'   dropped rows with their scaled residual and the pass that dropped them),
#'   `plot` (per-pass residual histograms), and `summary` (rows in/out per pass).
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap scale_fill_manual labs theme_classic
#' @importFrom lubridate yday
#' @importFrom stats lm residuals sd mad median na.exclude as.formula
#' @export
trim_ecostress_outliers <- function(data,
                                    var,
                                    time_var   = "time",
                                    group_var  = "point_id",
                                    thresholds = c(5, 2.58),
                                    harmonics  = 1,
                                    robust     = FALSE,
                                    keep_na    = TRUE,
                                    period     = 365.25) {
  
  stopifnot(is.data.frame(data),
            is.character(var), length(var) == 1,
            is.numeric(thresholds), length(thresholds) > 0,
            all(is.finite(thresholds)), all(thresholds > 0),
            harmonics >= 1)
  
  need <- c("variable", "value", time_var, group_var)
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "))
  
  if (is.unsorted(rev(thresholds))) {
    warning("`thresholds` is not decreasing; passes run in the order given.")
  }
  
  data <- as.data.frame(data)
  ord  <- seq_len(nrow(data))
  
  # Split target variable off, without NSE
  is_target <- !is.na(data$variable) & data$variable == var
  if (!any(is_target)) stop("No rows in `data` with variable == '", var, "'.")
  
  dat_var   <- data[is_target, , drop = FALSE]
  ord_var   <- ord[is_target]
  data_rest <- data[!is_target, , drop = FALSE]
  ord_rest  <- ord[!is_target]
  
  # Rows the model cannot use at all
  doy_all  <- lubridate::yday(as.Date(dat_var[[time_var]]))
  unusable <- is.na(dat_var$value) | is.na(doy_all) | is.na(dat_var[[group_var]])
  if (any(is.na(doy_all))) {
    warning(sum(is.na(doy_all)), " row(s) have an unparseable `", time_var, "`.")
  }
  
  dat_pass <- dat_var[!unusable, , drop = FALSE]
  ord_pass <- ord_var[!unusable]
  dat_skip <- dat_var[unusable, , drop = FALSE]
  ord_skip <- ord_var[unusable]
  
  keep    <- rep(TRUE, nrow(dat_pass))
  final_r <- rep(NA_real_, nrow(dat_pass))
  cut_by  <- rep(NA_integer_, nrow(dat_pass))
  diag    <- vector("list", length(thresholds))
  smry    <- data.frame(pass      = seq_along(thresholds),
                        threshold = thresholds,
                        n_in      = NA_integer_,
                        n_removed = NA_integer_)
  
  for (p in seq_along(thresholds)) {
    idx <- which(keep)
    cur <- dat_pass[idx, , drop = FALSE]
    smry$n_in[p] <- nrow(cur)
    
    r <- .seasonal_resid(cur, time_var, group_var, harmonics, period, robust)
    
    if (all(is.na(r))) {
      warning("Pass ", p, ": seasonal fit failed; no rows removed.")
      ok <- rep(TRUE, nrow(cur))
    } else {
      # NA residuals are kept, not dropped -- don't discard what you couldn't test
      ok <- is.na(r) | abs(r) < thresholds[p]
    }
    
    final_r[idx]      <- r
    cut_by[idx[!ok]]  <- p
    keep[idx[!ok]]    <- FALSE
    smry$n_removed[p] <- sum(!ok)
    
    diag[[p]] <- data.frame(pass = p, threshold = thresholds[p],
                            resid = r, retained = ok)
  }
  
  removed <- dat_pass[!keep, , drop = FALSE]
  removed$resid        <- final_r[!keep]
  removed$removed_pass <- cut_by[!keep]
  
  kept_idx <- c(ord_rest, ord_pass[keep], if (keep_na) ord_skip)
  if (!keep_na && nrow(dat_skip)) {
    message("Dropped ", nrow(dat_skip), " row(s) of '", var, "' with missing values.")
  }
  data_out <- data[sort(kept_idx), , drop = FALSE]
  rownames(data_out) <- NULL
  
  plot_dat <- do.call(rbind, diag)
  plot_dat <- plot_dat[!is.na(plot_dat$resid), , drop = FALSE]
  plot_dat$pass_lab <- factor(
    sprintf("Pass %d  (threshold %.2f)", plot_dat$pass, plot_dat$threshold),
    levels = sprintf("Pass %d  (threshold %.2f)", seq_along(thresholds), thresholds)
  )
  
  plt <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = resid, fill = retained)) +
    ggplot2::geom_histogram(bins = 100) +
    ggplot2::facet_wrap(~ pass_lab, ncol = 1, scales = "free") +
    ggplot2::scale_fill_manual(values = RETAINED_PAL, name = "Retained") +
    ggplot2::labs(x = "Scaled residual", y = "Count",
                  title = paste0("Outlier screening: ", var)) +
    ggplot2::theme_classic()
  
  list(data = data_out, removed = removed, plot = plt, summary = smry)
}


# Internal: scaled residuals from a harmonic seasonal fit with per-group intercepts
.seasonal_resid <- function(dat, time_var, group_var, harmonics, period, robust) {
  
  doy <- lubridate::yday(as.Date(dat[[time_var]]))
  
  X <- do.call(cbind, lapply(seq_len(harmonics), function(k) {
    cbind(sin(2 * pi * k * doy / period),
          cos(2 * pi * k * doy / period))
  }))
  colnames(X) <- paste0(c("sin", "cos"), rep(seq_len(harmonics), each = 2))
  
  md      <- data.frame(value = dat$value, X)
  grp     <- factor(dat[[group_var]])
  use_grp <- nlevels(grp) > 1
  if (use_grp) md$grp <- grp
  
  rhs   <- c(colnames(X), if (use_grp) "grp")
  n_par <- ncol(X) + if (use_grp) nlevels(grp) else 1
  
  if (nrow(md) < n_par + 2) {
    warning("Too few observations (", nrow(md), ") for ", n_par, " parameters.")
    return(rep(NA_real_, nrow(dat)))
  }
  
  fit <- tryCatch(
    stats::lm(stats::as.formula(paste("value ~", paste(rhs, collapse = " + "))),
              data = md, na.action = stats::na.exclude),
    error = function(e) { warning("Seasonal fit failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(dat)))
  
  r <- as.vector(stats::residuals(fit))
  
  if (robust) {
    ctr <- stats::median(r, na.rm = TRUE)
    scl <- stats::mad(r, na.rm = TRUE)
    if (!is.finite(scl) || scl == 0) scl <- stats::sd(r, na.rm = TRUE)
  } else {
    ctr <- mean(r, na.rm = TRUE)
    scl <- stats::sd(r, na.rm = TRUE)
  }
  if (!is.finite(scl) || scl == 0) {
    warning("Residual scale is zero or undefined; returning NA residuals.")
    return(rep(NA_real_, nrow(dat)))
  }
  
  (r - ctr) / scl
}



#' Filter skin temperatures likely to have diverged from bulk SST
#'
#' Removes observations of a target temperature variable collected under low wind
#' speed within a daytime solar-time window, where a diurnal warm layer is likely
#' to have formed and skin temperature diverged from bulk SST. An observation is
#' rejected only if its local solar time falls inside `(early, late)` *and* its
#' wind speed is at or below `threshold`; calm nights and windy days are retained.
#'
#' @param data Long format data frame with `variable` and `value` columns.
#' @param vars Named list giving the target, wind and hour variable names.
#' @param early,late Bounds of the solar-time window, in decimal hours.
#' @param threshold Wind speed in m/s at or below which an observation inside the
#'   window is rejected.
#' @param local_time If TRUE, `vars$time` already holds local solar time.
#' @param lon_var Longitude column, degrees east; -180..180 or 0..360 both work.
#' @param id_cols Columns identifying one observation. Defaults to every column
#'   other than `variable` and `value`. Exclude any column that varies between
#'   the target and its covariates (e.g. an extraction radius), or the pivot will
#'   place them on separate rows and nothing will be testable.
#' @param drop_covariates If TRUE, wind and hour are removed alongside the target.
#' @param keep_na If TRUE (default), observations with a missing wind speed or
#'   hour are retained rather than removed.
#'
#' @return A list with `data`, `removed`, `summary`, `plot` and `plot_values`.
#'   `summary` carries attributes `n_obs`, `n_keys`, `n_rows_removed`,
#'   `n_target_removed` and `n_target_total`.
#' @importFrom dplyr filter semi_join left_join all_of
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
#' @export
wind_filter <- function(data, vars, early = 10, late = 18, threshold = 2.0,
                        local_time = FALSE, lon_var = "lon", id_cols = NULL,
                        drop_covariates = FALSE, keep_na = TRUE) {
  
  stopifnot(is.data.frame(data), is.list(vars),
            is.numeric(early), is.numeric(late), is.numeric(threshold),
            length(early) == 1, length(late) == 1, length(threshold) == 1,
            is.finite(early), is.finite(late), is.finite(threshold),
            is.logical(local_time), is.logical(drop_covariates), is.logical(keep_na))
  
  miss_spec <- setdiff(c("target", "wind", "time"), names(vars))
  if (length(miss_spec)) stop("`vars` is missing: ", paste(miss_spec, collapse = ", "))
  var_names <- unlist(vars[c("target", "wind", "time")], use.names = FALSE)
  if (anyDuplicated(var_names)) stop("`vars` entries must name distinct variables.")
  
  if (early >= late) stop("`early` (", early, ") must be less than `late` (", late, ").")
  if (early < 0 || late > 24) stop("`early` and `late` must lie within 0-24 hours.")
  
  need <- c("variable", "value", if (!local_time) lon_var)
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "))
  
  if (is.null(id_cols)) id_cols <- setdiff(names(data), c("variable", "value"))
  if (!length(id_cols)) stop("No id columns available; supply `id_cols`.")
  miss_id <- setdiff(id_cols, names(data))
  if (length(miss_id)) stop("`id_cols` not in `data`: ", paste(miss_id, collapse = ", "))
  
  absent <- setdiff(var_names, unique(data$variable))
  if (length(absent)) stop("Variables not present in `data$variable`: ",
                           paste(absent, collapse = ", "))
  
  row_key <- ".wf_row_id"
  if (row_key %in% names(data)) stop("`data` already has a column named ", row_key, ".")
  
  wide <- tidyr::pivot_wider(
    dplyr::filter(data, .data$variable %in% var_names),
    id_cols = dplyr::all_of(id_cols), names_from = "variable", values_from = "value")
  
  listy <- var_names[vapply(wide[var_names], is.list, logical(1))]
  if (length(listy)) stop("Duplicate observations per key for: ",
                          paste(listy, collapse = ", "),
                          ". Deduplicate `data`, or pass a narrower `id_cols`.")
  
  hour <- wide[[vars$time]]; wind <- wide[[vars$wind]]
  if (!is.numeric(hour)) stop("`", vars$time, "` is not numeric after pivoting.")
  if (!is.numeric(wind)) stop("`", vars$wind, "` is not numeric after pivoting.")
  
  has_any  <- !is.na(hour) | !is.na(wind) | !is.na(wide[[vars$target]])
  n_any    <- sum(has_any)
  n_usable <- sum(!is.na(hour) & !is.na(wind))
  if (n_usable == 0) {
    stop("No observation has both `", vars$time, "` and `", vars$wind,
         "` after pivoting - the id columns are not aligning the variables.\n",
         "  id_cols: ", paste(id_cols, collapse = ", "), "\n",
         "  Pass a narrower `id_cols` (e.g. the spatial and time keys only).")
  }
  if (n_any && n_usable < 0.5 * n_any)
    warning(round(100 * (1 - n_usable / n_any)), "% of observations that have ",
            "any data lack a matching wind speed or hour; check `id_cols`.")
  if (n_any < nrow(wide))
    message("Note: ", format(nrow(wide) - n_any, big.mark = ","), " of ",
            format(nrow(wide), big.mark = ","), " key combinations are empty ",
            "(no target, wind or hour); `data` is a padded grid.")
  
  h_rng <- range(hour, na.rm = TRUE)
  if (h_rng[1] < 0 || h_rng[2] > 24)
    warning("`", vars$time, "` spans ", signif(h_rng[1], 4), " to ", signif(h_rng[2], 4),
            "; is it really decimal hours?")
  if (max(wind, na.rm = TRUE) > 60)
    warning("`", vars$wind, "` reaches ", signif(max(wind, na.rm = TRUE), 4),
            "; is `threshold` in the same units (m/s)?")
  
  solar_time <- if (local_time) hour %% 24 else (hour + wide[[lon_var]] / 15) %% 24
  
  in_window <- solar_time > early & solar_time < late
  calm      <- wind <= threshold
  risky     <- in_window & calm
  retained  <- !risky
  retained[is.na(retained)] <- keep_na
  bad <- !retained
  
  reason <- rep(NA_character_, nrow(wide))
  reason[bad & !is.na(risky)]                       <- "calm within window"
  reason[bad & is.na(risky) & in_window %in% TRUE]  <- "within window, wind unknown"
  reason[bad & is.na(risky) & calm      %in% TRUE]  <- "calm, time unknown"
  reason[bad & is.na(reason)]                       <- "missing covariates"
  
  diag_df <- wide[, id_cols, drop = FALSE]
  diag_df$solar_time <- solar_time; diag_df$wind_speed <- wind
  diag_df$target_value <- wide[[vars$target]]
  diag_df$retained <- retained; diag_df$reason <- reason
  
  target_set <- if (drop_covariates) var_names else vars$target
  in_target  <- data$variable %in% target_set
  is_drop    <- rep(FALSE, nrow(data))
  if (any(bad) && any(in_target)) {
    rows <- which(in_target)
    cand <- data[rows, id_cols, drop = FALSE]; cand[[row_key]] <- rows
    hit <- dplyr::semi_join(cand, wide[bad, id_cols, drop = FALSE], by = id_cols)
    is_drop[hit[[row_key]]] <- TRUE
  }
  
  removed <- data[is_drop, , drop = FALSE]
  if (nrow(removed)) removed <- dplyr::left_join(
    removed, diag_df[, c(id_cols, "solar_time", "wind_speed", "reason")], by = id_cols)
  data_out <- data[!is_drop, , drop = FALSE]
  rownames(data_out) <- rownames(removed) <- NULL
  
  lev <- c("calm within window", "within window, wind unknown",
           "calm, time unknown", "missing covariates")
  smry <- data.frame(reason = c("retained", lev),
                     n = c(sum(retained & has_any),
                           vapply(lev, function(l) sum(reason %in% l & has_any), integer(1))),
                     row.names = NULL)
  smry$pct_of_obs <- round(100 * smry$n / max(n_any, 1), 2)
  attr(smry, "n_obs")            <- n_any
  attr(smry, "n_keys")           <- nrow(wide)
  attr(smry, "n_rows_removed")   <- sum(is_drop)
  attr(smry, "n_target_removed") <- sum(!is.na(removed$value))
  attr(smry, "n_target_total")   <- sum(!is.na(wide[[vars$target]]))
  
  if (!any(bad)) message("No observations rejected. In-window: ",
                         sum(in_window, na.rm = TRUE), "; calm: ", sum(calm, na.rm = TRUE),
                         "; both: ", sum(risky, na.rm = TRUE), " (of ", n_usable, " testable).")
  
  x_lab <- if (local_time) "Local solar time (h)" else "Local solar time (h, from longitude)"
  pal <- RETAINED_PAL
  plt <- ggplot2::ggplot(diag_df[!is.na(diag_df$solar_time) & !is.na(diag_df$wind_speed), ],
                         ggplot2::aes(x = solar_time, y = wind_speed)) +
    ggplot2::annotate("rect", xmin = early, xmax = late, ymin = -Inf, ymax = threshold,
                      fill = "firebrick", alpha = 0.12) +
    ggplot2::geom_point(ggplot2::aes(colour = retained), alpha = 0.4, size = 0.8) +
    ggplot2::geom_vline(xintercept = c(early, late), linetype = 2, colour = "grey40") +
    ggplot2::geom_hline(yintercept = threshold, linetype = 2, colour = "grey40") +
    ggplot2::scale_colour_manual(values = pal, name = "Retained") +
    ggplot2::labs(x = x_lab, y = "Wind speed (m/s)",
                  title = paste0("Wind / solar-time screening: ", vars$target),
                  subtitle = "Shaded area is the rejection region") + ggplot2::theme_classic()
  plt_val <- ggplot2::ggplot(diag_df[!is.na(diag_df$target_value), ],
                             ggplot2::aes(x = target_value, fill = retained)) +
    ggplot2::geom_histogram(bins = 60, position = "identity", alpha = 0.55) +
    ggplot2::scale_fill_manual(values = pal, name = "Retained") +
    ggplot2::labs(x = vars$target, y = "Count",
                  title = "Target distribution, retained vs removed") + ggplot2::theme_classic()
  
  list(data = data_out, removed = removed, summary = smry,
       plot = plt, plot_values = plt_val)
}





#' Filter a temperature variable against a co-located reference observation
#'
#' Regresses a target variable on a reference variable across match-ups, then
#' removes target observations whose residual falls outside an asymmetric window
#' around the fitted relationship. Bounds are asymmetric by default because the
#' two tails mean different things: residuals far below the fit are usually cloud
#' contamination (spuriously cold), while residuals above it may be genuine skin
#' warming, so the warm side is given more room.
#'
#' Target observations with no reference match-up are retained by default; a
#' missing reference makes the residual unknown, not out of range.
#'
#' @param data Long format data frame with `variable` and `value` columns.
#' @param vars Named list with `target` (filtered) and `reference` (comparison),
#'   e.g. `list(target = "eco_sst_valid", reference = "modis_sst_aqua")`.
#' @param slope Slope of the target-reference relationship. Fixed at 1 by default,
#'   which is the physically correct value when both variables measure the same
#'   quantity, and which makes the fit far more resistant to contamination: only a
#'   location parameter is estimated, so outliers can shift the line but cannot
#'   tilt it. Estimating the slope from contaminated match-ups biases it toward
#'   zero (regression dilution, worsened by the very outliers being screened for),
#'   which widens residuals at the ends of the range and hides outliers there.
#'   Set `slope = NA` to estimate it instead.
#' @param intercept Offset between the two variables. `NA` (default) estimates it
#'   using `center`, from `target - slope * reference`.
#'   Set to `0` to screen on the raw difference with no bias correction.
#' @param lower,upper Residual bounds in standard deviations. An observation is
#'   removed if its residual is below `-lower * s` or above `upper * s`, where `s`
#'   is the residual scale. Defaults 2.58 and 3.0.
#' @param id_cols Columns identifying one observation. Defaults to every column
#'   other than `variable` and `value`. Exclude any column that varies between
#'   the target and the reference (e.g. an extraction radius), or the pivot will
#'   place them on separate rows and no match-ups will be found.
#' @param by Optional column(s) to fit separately within, e.g. `"aoi"`. Groups
#'   with fewer than `min_n` match-ups fall back to the pooled fit.
#' @param passes Number of fit-and-filter passes. Each pass refits on the
#'   surviving match-ups. Defaults to 1.
#' @param scale Scale estimator used to set the thresholds. "mad" (default) is
#'   `stats::mad()`, the median absolute deviation rescaled by 1.4826 so that it
#'   estimates the same quantity as `sd` for normal residuals. It has a breakdown
#'   point of 50%, so the very outliers being screened for cannot inflate the
#'   bounds meant to catch them. "sd" is the ordinary standard deviation, which a
#'   single bad retrieval can move without limit. Note that `mad` gives a smaller
#'   scale on contaminated data, so the same `lower`/`upper` values remove more:
#'   thresholds tuned against `sd` should be widened when switching.
#' @param center Location estimator. "median" (default) estimates the offset as
#'   the median of `target - slope * reference`, and fits with `MASS::rlm()` when
#'   the slope is estimated. "mean" uses the mean and `stats::lm()`.
#' @param robust Deprecated. `TRUE` sets scale = "mad", center = "median";
#'   `FALSE` sets scale = "sd", center = "mean".
#' @param keep_na If TRUE (default), target observations whose reference is
#'   missing are retained.
#' @param min_n Minimum match-ups required to fit. Defaults to 30.
#'
#' @return A list with `data` (filtered), `removed` (dropped rows with residual,
#'   reference value and reason), `fit` (per-group coefficients, scale, n, and
#'   which coefficients were fixed),
#'   `summary` (counts), `plot` (target vs reference) and `plot_resid`.
#' @importFrom dplyr filter semi_join left_join all_of
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
#' @export
match_up_filter <- function(data,
                            vars,
                            lower   = 2.58,
                            upper   = 3.0,
                            slope     = 1,
                            intercept = NA,
                            id_cols = NULL,
                            by      = NULL,
                            passes  = 1,
                            scale   = c("mad", "sd"),
                            center  = c("median", "mean"),
                            robust  = NULL,
                            keep_na = TRUE,
                            min_n   = 30) {
  
  stopifnot(is.data.frame(data), is.list(vars),
            is.numeric(lower), is.numeric(upper),
            length(lower) == 1, length(upper) == 1,
            is.finite(lower), is.finite(upper), lower > 0, upper > 0,
            is.numeric(passes), length(passes) == 1, passes >= 1,
            is.logical(keep_na))
  scale <- match.arg(scale); center <- match.arg(center)
  if (!is.null(robust)) {
    if (!is.logical(robust) || length(robust) != 1 || is.na(robust))
      stop("`robust` must be TRUE or FALSE, or NULL to use `scale`/`center`.")
    scale  <- if (robust) "mad"    else "sd"
    center <- if (robust) "median" else "mean"
    message("`robust` is deprecated; using scale = \"", scale,
            "\", center = \"", center, "\".")
  }
  fix_b <- !is.null(slope)     && !all(is.na(slope))
  fix_a <- !is.null(intercept) && !all(is.na(intercept))
  if (fix_b && (!is.numeric(slope)     || length(slope) != 1     || !is.finite(slope)))
    stop("`slope` must be a single finite number, or NA to estimate it.")
  if (fix_a && (!is.numeric(intercept) || length(intercept) != 1 || !is.finite(intercept)))
    stop("`intercept` must be a single finite number, or NA to estimate it.")
  if (fix_b && fix_a && passes > 1) {
    message("Both coefficients are fixed; extra passes only shrink the scale.")
  }
  
  miss_spec <- setdiff(c("target", "reference"), names(vars))
  if (length(miss_spec)) stop("`vars` is missing: ", paste(miss_spec, collapse = ", "))
  tgt <- vars$target; ref <- vars$reference
  if (identical(tgt, ref)) stop("`target` and `reference` must be different variables.")
  
  miss <- setdiff(c("variable", "value"), names(data))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "))
  
  if (is.null(id_cols)) id_cols <- setdiff(names(data), c("variable", "value"))
  miss_id <- setdiff(c(id_cols, by), names(data))
  if (length(miss_id)) stop("Columns not in `data`: ", paste(miss_id, collapse = ", "))
  if (length(by) && !all(by %in% id_cols))
    stop("`by` columns must also be in `id_cols`: ", paste(setdiff(by, id_cols), collapse = ", "))
  
  absent <- setdiff(c(tgt, ref), unique(data$variable))
  if (length(absent)) stop("Variables not present in `data$variable`: ",
                           paste(absent, collapse = ", "))
  
  row_key <- ".mu_row_id"
  if (row_key %in% names(data)) stop("`data` already has a column named ", row_key, ".")
  
  fix_b_pre <- !is.null(slope) && !all(is.na(slope))
  if (center == "median" && !fix_b_pre && !requireNamespace("MASS", quietly = TRUE)) {
    warning("MASS not installed; the slope will be fitted by lm() rather than rlm().")
  }
  
  # ---- wide view -----------------------------------------------------------
  wide <- tidyr::pivot_wider(
    dplyr::filter(data, .data$variable %in% c(tgt, ref)),
    id_cols = dplyr::all_of(id_cols), names_from = "variable", values_from = "value")
  
  listy <- c(tgt, ref)[vapply(wide[c(tgt, ref)], is.list, logical(1))]
  if (length(listy)) stop("Duplicate observations per key for: ",
                          paste(listy, collapse = ", "),
                          ". Deduplicate `data`, or pass a narrower `id_cols`.")
  
  y <- wide[[tgt]]; x <- wide[[ref]]
  if (!is.numeric(y)) stop("`", tgt, "` is not numeric after pivoting.")
  if (!is.numeric(x)) stop("`", ref, "` is not numeric after pivoting.")
  
  n_tgt   <- sum(!is.na(y))
  n_match <- sum(!is.na(y) & !is.na(x))
  if (n_tgt == 0) stop("No non-missing `", tgt, "` values.")
  if (n_match < min_n) {
    stop("Only ", n_match, " match-ups between `", tgt, "` and `", ref,
         "` (min_n = ", min_n, ").\n",
         "  If this is unexpected, the id columns may not be aligning them.\n",
         "  id_cols: ", paste(id_cols, collapse = ", "))
  }
  if (n_match < 0.2 * n_tgt)
    message("Only ", n_match, " of ", n_tgt, " `", tgt,
            "` observations have a reference match-up (",
            round(100 * n_match / n_tgt, 1), "%); the rest are ",
            if (keep_na) "retained untested." else "removed (keep_na = FALSE).")
  
  grp <- if (length(by)) interaction(wide[by], drop = TRUE, sep = " | ")
  else factor(rep("(all)", nrow(wide)))
  
  # ---- iterate -------------------------------------------------------------
  eligible <- !is.na(y) & !is.na(x)   # usable as a match-up
  keep     <- rep(TRUE, nrow(wide))   # match-up survives so far
  resid    <- rep(NA_real_, nrow(wide))
  fit_log  <- list()
  
  for (p in seq_len(passes)) {
    use <- eligible & keep
    fp  <- .mu_fit(y, x, grp, use, scale, center, min_n, slope, intercept)
    resid <- fp$resid
    fit_log[[p]] <- cbind(pass = p, fp$stats)
    
    lo <- -lower * fp$scale
    hi <-  upper * fp$scale
    ok <- is.na(resid) | (resid > lo & resid < hi)
    keep <- keep & ok
  }
  
  fit <- do.call(rbind, fit_log)
  rownames(fit) <- NULL
  
  # ---- classify ------------------------------------------------------------
  scale_v <- .mu_fit(y, x, grp, eligible & keep, scale, center, min_n,
                     slope, intercept)$scale
  bad_tested <- eligible & !keep
  reason <- rep(NA_character_, nrow(wide))
  reason[bad_tested & resid <= -lower * scale_v] <- "below fit (cold / cloud)"
  reason[bad_tested & resid >=  upper * scale_v] <- "above fit (warm)"
  reason[bad_tested & is.na(reason)]             <- "outside bounds"
  
  no_match <- !is.na(y) & is.na(x)
  drop_obs <- bad_tested | (no_match & !keep_na)
  reason[no_match & !keep_na] <- "no reference match-up"
  
  diag_df <- wide[, id_cols, drop = FALSE]
  diag_df$target_value <- y
  diag_df$ref_value    <- x
  diag_df$resid        <- resid
  diag_df$group        <- as.character(grp)
  diag_df$retained     <- !drop_obs
  diag_df$reason       <- reason
  
  # ---- remove from the long frame -----------------------------------------
  in_target <- data$variable %in% tgt
  is_drop   <- rep(FALSE, nrow(data))
  if (any(drop_obs) && any(in_target)) {
    rows <- which(in_target)
    cand <- data[rows, id_cols, drop = FALSE]; cand[[row_key]] <- rows
    hit  <- dplyr::semi_join(cand, wide[drop_obs, id_cols, drop = FALSE], by = id_cols)
    is_drop[hit[[row_key]]] <- TRUE
  }
  
  removed <- data[is_drop, , drop = FALSE]
  if (nrow(removed)) removed <- dplyr::left_join(
    removed, diag_df[, c(id_cols, "ref_value", "resid", "reason")], by = id_cols)
  data_out <- data[!is_drop, , drop = FALSE]
  rownames(data_out) <- rownames(removed) <- NULL
  
  # ---- summary -------------------------------------------------------------
  lev <- c("below fit (cold / cloud)", "above fit (warm)",
           "outside bounds", "no reference match-up")
  smry <- data.frame(
    reason = c("retained (tested)", "retained (no match-up)", lev),
    n = c(sum(eligible & keep), sum(no_match & keep_na),
          vapply(lev, function(l) sum(reason %in% l), integer(1))),
    row.names = NULL)
  smry$pct_of_target <- round(100 * smry$n / n_tgt, 2)
  attr(smry, "n_target")   <- n_tgt
  attr(smry, "n_match_up") <- n_match
  attr(smry, "n_removed")  <- sum(is_drop)
  
  # ---- plots ---------------------------------------------------------------
  pal <- RETAINED_PAL
  pd  <- diag_df[eligible, , drop = FALSE]
  last <- fit[fit$pass == max(fit$pass), , drop = FALSE]
  
  bands <- do.call(rbind, lapply(seq_len(nrow(last)), function(i) {
    rng <- range(pd$ref_value[pd$group == last$group[i]], na.rm = TRUE)
    if (!all(is.finite(rng))) return(NULL)
    data.frame(group = last$group[i], x = rng,
               y   = last$intercept[i] + last$slope[i] * rng,
               ylo = last$intercept[i] + last$slope[i] * rng - lower * last$scale[i],
               yhi = last$intercept[i] + last$slope[i] * rng + upper * last$scale[i])
  }))
  
  plt <- ggplot2::ggplot(pd, ggplot2::aes(x = ref_value, y = target_value)) +
    ggplot2::geom_point(ggplot2::aes(colour = retained), alpha = 0.4, size = 0.8) +
    ggplot2::geom_line(data = bands, ggplot2::aes(x = x, y = y), colour = "grey20") +
    ggplot2::geom_line(data = bands, ggplot2::aes(x = x, y = ylo),
                       colour = "firebrick", linetype = 2) +
    ggplot2::geom_line(data = bands, ggplot2::aes(x = x, y = yhi),
                       colour = "firebrick", linetype = 2) +
    ggplot2::scale_colour_manual(values = pal, name = "Retained") +
    ggplot2::labs(x = ref, y = tgt,
                  title = paste0("Match-up screening: ", tgt, " vs ", ref),
                  subtitle = sprintf("Bounds: -%.2f / +%.2f %s%s", lower, upper,
                                     if (scale == "mad") "robust SD (MAD)" else "SD",
                                     if (fix_b) sprintf("  |  slope fixed at %g", slope) else "")) +
    ggplot2::theme_classic()
  if (length(by)) plt <- plt + ggplot2::facet_wrap(~ group)
  
  plt_resid <- ggplot2::ggplot(pd[!is.na(pd$resid), ],
                               ggplot2::aes(x = resid, fill = retained)) +
    ggplot2::geom_histogram(bins = 80, position = "identity", alpha = 0.6) +
    ggplot2::scale_fill_manual(values = pal, name = "Retained") +
    ggplot2::labs(x = "Residual from fit", y = "Count",
                  title = "Match-up residuals") + ggplot2::theme_classic()
  if (length(by)) plt_resid <- plt_resid + ggplot2::facet_wrap(~ group, scales = "free_y")
  
  list(data = data_out, removed = removed, fit = fit, summary = smry,
       plot = plt, plot_resid = plt_resid)
}


# Internal: fit target ~ reference within groups, return residuals and scale
.mu_fit <- function(y, x, grp, use, scale, center, min_n,
                    slope = NA, intercept = NA) {
  
  resid <- rep(NA_real_, length(y))
  scl   <- rep(NA_real_, length(y))   # not `scale`: that name is an argument here
  stats <- list()
  
  pooled <- .mu_fit1(y[use], x[use], scale, center, slope, intercept)
  if (is.null(pooled)) stop("Pooled fit of target on reference failed.")
  
  for (g in levels(grp)) {
    idx <- which(grp == g & use)
    f   <- if (length(idx) >= min_n) .mu_fit1(y[idx], x[idx], scale, center, slope, intercept) else NULL
    if (is.null(f)) {
      if (length(idx)) message("Group '", g, "': ", length(idx),
                               " match-ups (< min_n); using the pooled fit.")
      f <- pooled
    }
    gi <- which(grp == g)
    resid[gi] <- y[gi] - (f$intercept + f$slope * x[gi])
    scl[gi]   <- f$scale
    stats[[g]] <- data.frame(group = g, n = length(idx), intercept = f$intercept,
                             slope = f$slope, scale = f$scale, fixed = f$fixed,
                             scale_est = scale, center_est = center)
  }
  list(resid = resid, scale = scl, stats = do.call(rbind, stats))
}

.mu_fit1 <- function(y, x, scale, center, slope = NA, intercept = NA) {
  
  ok <- !is.na(y) & !is.na(x)
  y <- y[ok]; x <- x[ok]
  if (length(y) < 3) return(NULL)
  
  fix_b <- !is.null(slope)     && !all(is.na(slope))
  fix_a <- !is.null(intercept) && !all(is.na(intercept))
  use_rlm <- center == "median" && requireNamespace("MASS", quietly = TRUE)
  
  if (fix_b && fix_a) {
    # Nothing to estimate: residuals are the raw offset from the stated line.
    b <- slope; a <- intercept
    r <- y - a - b * x
    lab <- "slope + intercept"
    
  } else if (fix_b) {
    # Location-only fit. The offset is a median/mean of the differences, which
    # a contaminated tail can shift but cannot tilt.
    b <- slope
    d <- y - b * x
    a <- if (center == "median") stats::median(d) else mean(d)
    r <- d - a
    lab <- "slope"
    
  } else if (fix_a) {
    a  <- intercept
    md <- data.frame(yy = y - a, x = x)
    if (length(unique(x)) < 2) return(NULL)
    fit <- tryCatch({
      if (use_rlm) MASS::rlm(yy ~ 0 + x, data = md, maxit = 50)
      else         stats::lm(yy ~ 0 + x, data = md)
    }, error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) return(NULL)
    b <- unname(stats::coef(fit)[1]); r <- as.vector(stats::residuals(fit))
    lab <- "intercept"
    
  } else {
    if (length(unique(x)) < 2) return(NULL)
    md <- data.frame(y = y, x = x)
    fit <- tryCatch({
      if (use_rlm) MASS::rlm(y ~ x, data = md, maxit = 50)
      else         stats::lm(y ~ x, data = md)
    }, error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) return(NULL)
    cf <- stats::coef(fit)
    a <- unname(cf[1]); b <- unname(cf[2]); r <- as.vector(stats::residuals(fit))
    lab <- "none"
  }
  
  # mad() takes its own median centre, so the scale does not depend on `a`.
  s <- if (scale == "mad") stats::mad(r, na.rm = TRUE) else stats::sd(r, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    s <- stats::sd(r, na.rm = TRUE)   # ties or heavy discretisation give mad == 0
  }
  if (!is.finite(s) || s == 0) return(NULL)
  
  list(intercept = a, slope = b, scale = s, fixed = lab)
}




#' Summarise variability between target and reference variables around 1:1
#'
#' For each target-reference pair, computes match-up statistics on the difference
#' `d = target - reference`. Because the comparison is against a fixed 1:1 line,
#' no slope is estimated and no regression dilution arises. Both parametric and
#' robust statistics are returned, since a small number of failed retrievals can
#' dominate a variance while barely moving a median.
#'
#' Statistics returned per pair (and per `by` group):
#' \describe{
#'   \item{n}{match-ups with both values present}
#'   \item{bias}{mean of `d` -- systematic offset from 1:1}
#'   \item{sd}{standard deviation of `d` -- scatter about the bias-adjusted 1:1 line}
#'   \item{rmsd}{root mean square of `d` -- total departure from the exact 1:1 line}
#'   \item{median_bias}{median of `d`, resistant to outliers}
#'   \item{rsd}{robust SD, `mad(d)`, already scaled to match `sd` for normal data}
#'   \item{p025, p975}{2.5th and 97.5th percentiles of `d`}
#'   \item{r}{Pearson correlation between target and reference}
#'   \item{slope_ols}{OLS slope of target on reference, as a diagnostic only. It is
#'     attenuated by measurement error in the reference and by outliers, so a value
#'     below 1 is expected and is not evidence against the 1:1 assumption.}
#'   \item{n_trim, sd_trim}{if `trim` is set, the count and SD after dropping
#'     match-ups more than `trim` robust SDs from the median difference}
#' }
#'
#' `sd` and `rsd` should agree closely when the differences are well behaved. A
#' large gap between them means the variance is being driven by a few points, and
#' `rsd` is the number to trust.
#'
#' @param data Long format data frame with `variable` and `value` columns.
#' @param targets Character vector of target variable names.
#' @param references Character vector of reference variable names.
#' @param pairs If "cross" (default), every target is compared with every
#'   reference. If "elementwise", `targets` and `references` are paired in order
#'   and must have the same length. Self-comparisons are always dropped.
#' @param id_cols Columns identifying one observation. Defaults to every column
#'   other than `variable` and `value`. Exclude any column that varies between the
#'   variables being compared (e.g. an extraction radius), or the pivot will place
#'   them on separate rows and no match-ups will be found.
#' @param by Optional column(s) to summarise within, e.g. `"aoi"`. A pooled
#'   `(all)` row is always included.
#' @param valid_range Optional length-2 numeric. Values outside it are treated as
#'   missing before any statistic is computed. Use this to exclude failed
#'   retrievals, e.g. `c(271, 313)` for sea surface temperature in kelvin.
#' @param trim Optional number of robust SDs beyond which match-ups are excluded
#'   from the extra `sd_trim` column. The main statistics are never trimmed.
#' @param min_n Pairs with fewer match-ups than this return `NA` statistics rather
#'   than erroring. Defaults to 10.
#'
#' @return A list with `table` (one row per pair and group), `differences` (the
#'   per-match-up differences, for plotting or further work), `plot` (target vs
#'   reference with the 1:1 line) and `plot_resid` (distribution of differences).
#' @importFrom dplyr filter all_of
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
#' @export
match_up_variance <- function(data,
                              targets,
                              references,
                              pairs       = c("cross", "elementwise"),
                              id_cols     = NULL,
                              by          = NULL,
                              valid_range = NULL,
                              trim        = NULL,
                              min_n       = 10) {
  
  pairs <- match.arg(pairs)
  stopifnot(is.data.frame(data),
            is.character(targets), length(targets) > 0,
            is.character(references), length(references) > 0)
  
  miss <- setdiff(c("variable", "value"), names(data))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "))
  
  if (is.null(id_cols)) id_cols <- setdiff(names(data), c("variable", "value"))
  miss_id <- setdiff(c(id_cols, by), names(data))
  if (length(miss_id)) stop("Columns not in `data`: ", paste(miss_id, collapse = ", "))
  if (length(by) && !all(by %in% id_cols))
    stop("`by` columns must also be in `id_cols`: ",
         paste(setdiff(by, id_cols), collapse = ", "))
  
  if (!is.null(valid_range)) {
    if (!is.numeric(valid_range) || length(valid_range) != 2 ||
        any(!is.finite(valid_range)) || valid_range[1] >= valid_range[2])
      stop("`valid_range` must be two increasing finite numbers.")
  }
  if (!is.null(trim) && (!is.numeric(trim) || length(trim) != 1 ||
                         !is.finite(trim) || trim <= 0))
    stop("`trim` must be a single positive number, or NULL.")
  
  all_vars <- unique(c(targets, references))
  absent   <- setdiff(all_vars, unique(data$variable))
  if (length(absent)) stop("Variables not present in `data$variable`: ",
                           paste(absent, collapse = ", "))
  
  # ---- pair list -----------------------------------------------------------
  if (pairs == "elementwise") {
    if (length(targets) != length(references))
      stop("`targets` and `references` must be the same length when ",
           "pairs = \"elementwise\".")
    pr <- data.frame(target = targets, reference = references,
                     stringsAsFactors = FALSE)
  } else {
    pr <- expand.grid(target = targets, reference = references,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  }
  pr <- pr[pr$target != pr$reference, , drop = FALSE]
  if (!nrow(pr)) stop("No comparable pairs (every pair was a self-comparison).")
  rownames(pr) <- NULL
  
  # ---- one pivot for all variables ----------------------------------------
  wide <- tidyr::pivot_wider(
    dplyr::filter(data, .data$variable %in% all_vars),
    id_cols = dplyr::all_of(id_cols), names_from = "variable", values_from = "value")
  
  listy <- all_vars[vapply(wide[all_vars], is.list, logical(1))]
  if (length(listy)) stop("Duplicate observations per key for: ",
                          paste(listy, collapse = ", "),
                          ". Deduplicate `data`, or pass a narrower `id_cols`.")
  
  not_num <- all_vars[!vapply(wide[all_vars], is.numeric, logical(1))]
  if (length(not_num)) stop("Not numeric after pivoting: ",
                            paste(not_num, collapse = ", "))
  
  if (!is.null(valid_range)) {
    n_gated <- 0L
    for (v in all_vars) {
      bad <- !is.na(wide[[v]]) & (wide[[v]] < valid_range[1] | wide[[v]] > valid_range[2])
      n_gated <- n_gated + sum(bad)
      wide[[v]][bad] <- NA_real_
    }
    if (n_gated) message(n_gated, " value(s) outside `valid_range` set to NA.")
  }
  
  grp <- if (length(by)) as.character(interaction(wide[by], drop = TRUE, sep = " | "))
  else rep("(all)", nrow(wide))
  groups <- c("(all)", if (length(by)) sort(unique(grp)))
  
  # ---- accumulate ----------------------------------------------------------
  rows  <- list()
  diffs <- list()
  thin  <- getOption("match_up_variance.max_points", 20000L)
  
  for (i in seq_len(nrow(pr))) {
    y <- wide[[pr$target[i]]]; x <- wide[[pr$reference[i]]]
    d <- y - x
    for (g in groups) {
      sel <- (g == "(all)" | grp == g) & !is.na(d)
      rows[[length(rows) + 1L]] <- .muv_stats(
        d[sel], y[sel], x[sel], pr$target[i], pr$reference[i], g, trim, min_n)
    }
    keep <- which(!is.na(d))
    if (length(keep)) {
      if (length(keep) > thin) keep <- sort(sample(keep, thin))
      diffs[[length(diffs) + 1L]] <- data.frame(
        target = pr$target[i], reference = pr$reference[i],
        pair = paste(pr$target[i], "vs", pr$reference[i]),
        group = grp[keep], target_value = y[keep], ref_value = x[keep],
        diff = d[keep], stringsAsFactors = FALSE)
    }
  }
  
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  dif <- if (length(diffs)) do.call(rbind, diffs) else
    data.frame(target = character(), reference = character(), pair = character(),
               group = character(), target_value = numeric(),
               ref_value = numeric(), diff = numeric())
  rownames(dif) <- NULL
  
  thin_hit <- any(tab$n[tab$group == "(all)"] > thin, na.rm = TRUE)
  if (thin_hit) message("Plots subsample to ", format(thin, big.mark = ","),
                        " points per pair; `table` uses all match-ups.")
  if (any(tab$n < min_n, na.rm = TRUE))
    message(sum(tab$n < min_n, na.rm = TRUE),
            " pair/group combination(s) had fewer than min_n = ", min_n,
            " match-ups; statistics returned as NA.")
  
  # ---- plots ---------------------------------------------------------------
  if (nrow(dif)) {
    plt <- ggplot2::ggplot(dif, ggplot2::aes(x = ref_value, y = target_value)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "firebrick") +
      ggplot2::geom_point(alpha = 0.25, size = 0.7, colour = "grey30") +
      ggplot2::facet_wrap(~ pair, scales = "free") +
      ggplot2::labs(x = "Reference", y = "Target",
                    title = "Match-ups against the 1:1 line") +
      ggplot2::theme_classic()
    plt_r <- ggplot2::ggplot(dif, ggplot2::aes(x = diff)) +
      ggplot2::geom_histogram(bins = 80, fill = "grey40") +
      ggplot2::geom_vline(xintercept = 0, colour = "firebrick") +
      ggplot2::facet_wrap(~ pair, scales = "free") +
      ggplot2::labs(x = "Target - reference", y = "Count",
                    title = "Distribution of differences") +
      ggplot2::theme_classic()
    if (length(by)) {
      plt   <- plt   + ggplot2::aes(colour = group)
      plt_r <- plt_r + ggplot2::aes(fill = group)
    }
  } else {
    plt <- plt_r <- NULL
    warning("No match-ups found for any pair; check `id_cols`.")
  }
  
  list(table = tab, differences = dif, plot = plt, plot_resid = plt_r)
}


# Internal: statistics for one pair within one group
.muv_stats <- function(d, y, x, target, reference, group, trim, min_n) {
  
  n   <- length(d)
  out <- data.frame(target = target, reference = reference, group = group, n = n,
                    bias = NA_real_, sd = NA_real_, rmsd = NA_real_,
                    median_bias = NA_real_, rsd = NA_real_,
                    p025 = NA_real_, p975 = NA_real_,
                    r = NA_real_, slope_ols = NA_real_,
                    stringsAsFactors = FALSE)
  if (!is.null(trim)) { out$n_trim <- NA_integer_; out$sd_trim <- NA_real_ }
  if (n < min_n) return(out)
  
  out$bias        <- mean(d)
  out$sd          <- stats::sd(d)
  out$rmsd        <- sqrt(mean(d^2))
  out$median_bias <- stats::median(d)
  out$rsd         <- stats::mad(d)
  qq              <- stats::quantile(d, c(0.025, 0.975), names = FALSE)
  out$p025        <- qq[1]; out$p975 <- qq[2]
  
  if (stats::sd(x) > 0 && stats::sd(y) > 0) {
    out$r         <- stats::cor(y, x)
    out$slope_ols <- unname(stats::coef(stats::lm(y ~ x))[2])
  }
  
  if (!is.null(trim)) {
    s  <- if (is.finite(out$rsd) && out$rsd > 0) out$rsd else out$sd
    ok <- abs(d - out$median_bias) <= trim * s
    out$n_trim  <- sum(ok)
    out$sd_trim <- if (sum(ok) >= 2) stats::sd(d[ok]) else NA_real_
  }
  out
}

#' Validate a preprocessing configuration against the data it will be run on
#'
#' Every problem found is collected and reported in a single error, rather than
#' surfacing one at a time over successive runs. Most of the failures this
#' catches are silent ones: a misspelled stage key makes every argument resolve
#' to `NULL`, and a variable named with its raw extraction name rather than its
#' post-`get_variable_value_format()` label simply matches nothing.
#'
#' @param config The `preprocess_config` list.
#' @param data The raw long format data frame the pipeline will be run on. Used
#'   to check that `group_id` and the key columns actually exist.
#'
#' @return `config`, invisibly, if it is valid. Otherwise an error listing every
#'   problem found.
#' @export
validate_config <- function(config, data) {

  if (!is.list(config)) stop("`config` must be a list.")
  problems <- character()
  add <- function(...) problems <<- c(problems, paste0(...))

  required <- c("group_id", "value_var", "convert_units", "physics_filter",
                "wind_filter", "outliers", "matchup_filter")
  optional <- c("diagnostics", "match_up_variance", "marss_format", "output")

  absent <- setdiff(required, names(config))
  if (length(absent))
    add("Missing required config key(s): ", paste(absent, collapse = ", "))

  # A misspelled stage key does not error, it silently makes every argument of
  # that stage resolve to NULL, so this warning is issued immediately rather
  # than being deferred to the end of the run where it is easy to miss.
  extra <- setdiff(names(config), c(required, optional))
  if (length(extra))
    warning("Unrecognised config key(s), which will be ignored: ",
            paste(extra, collapse = ", "),
            ". Check for a typo against: ",
            paste(c(required, optional), collapse = ", "),
            call. = FALSE, immediate. = TRUE)

  # ---- group_id ------------------------------------------------------------
  # NULL is a legitimate value (pooled diagnostics only), but the key must be
  # present, so test membership rather than NULL-ness.
  if ("group_id" %in% names(config)) {
    g <- config$group_id
    if (!is.null(g)) {
      if (!is.character(g) || length(g) != 1 || is.na(g) || !nzchar(g)) {
        add("`group_id` must be a single non-empty string, or NULL.")
      } else if (is.data.frame(data) && !g %in% names(data)) {
        add("`group_id` = '", g, "' is not a column of the data. Available: ",
            paste(names(data), collapse = ", "))
      }
    }
  }

  # ---- value_var -----------------------------------------------------------
  vars <- config$value_var$variables
  known <- character()
  if (!is.list(vars) || !length(vars) || is.null(names(vars)) ||
      !all(nzchar(names(vars)))) {
    add("`value_var$variables` must be a non-empty named list.")
  } else {
    known <- names(vars)
    if (anyDuplicated(known))
      add("`value_var$variables` has duplicate names: ",
          paste(unique(known[duplicated(known)]), collapse = ", "))
    bad <- known[!vapply(vars, function(s) is.list(s) && length(s) > 0 &&
                           !is.null(names(s)), logical(1))]
    if (length(bad))
      add("`value_var$variables` entries must be named lists of column values: ",
          paste(bad, collapse = ", "))
    if (is.data.frame(data)) {
      key_cols <- unique(unlist(lapply(vars, names)))
      miss <- setdiff(c(key_cols, config$value_var$value_var), names(data))
      if (length(miss))
        add("Columns required by `value_var` are not in the data: ",
            paste(miss, collapse = ", "))
    }
  }

  # Variables are referenced downstream by their post-format label, i.e. the
  # names of value_var$variables, not the raw `variable` column values.
  check_known <- function(x, where) {
    if (!length(known)) return(invisible(NULL))
    unknown <- setdiff(x[!is.na(x)], known)
    if (length(unknown))
      add("`", where, "` names variable(s) not defined in `value_var$variables`: ",
          paste(unknown, collapse = ", "))
  }

  # ---- convert_units -------------------------------------------------------
  cu <- config$convert_units
  if (!is.list(cu)) {
    add("`convert_units` must be a list.")
  } else {
    check_known(cu$variables, "convert_units$variables")
    for (k in c("from", "to")) {
      v <- cu[[k]]
      if (!is.character(v) || length(v) != 1 || !v %in% c("K", "C", "F"))
        add("`convert_units$", k, "` must be one of \"K\", \"C\", \"F\".")
    }
  }

  # ---- physics_filter ------------------------------------------------------
  pf <- config$physics_filter
  if (!is.list(pf)) {
    add("`physics_filter` must be a list.")
  } else {
    check_known(pf$variables, "physics_filter$variables")
    lo <- pf$lower_threshold; hi <- pf$upper_threshold
    if (!is.numeric(lo) || length(lo) != 1 || !is.finite(lo) ||
        !is.numeric(hi) || length(hi) != 1 || !is.finite(hi)) {
      add("`physics_filter$lower_threshold` and `upper_threshold` must be ",
          "single finite numbers.")
    } else if (lo >= hi) {
      add("`physics_filter$lower_threshold` (", lo, ") must be below ",
          "`upper_threshold` (", hi, ").")
    }
  }

  # ---- wind_filter ---------------------------------------------------------
  wf <- config$wind_filter
  if (!is.list(wf)) {
    add("`wind_filter` must be a list.")
  } else {
    if (!is.list(wf$variables) || !length(wf$variables)) {
      add("`wind_filter$variables` must be a non-empty list of ",
          "target / wind / time triples.")
    } else {
      for (i in seq_along(wf$variables)) {
        spec <- wf$variables[[i]]
        need <- c("target", "wind", "time")
        if (!is.list(spec) || length(setdiff(need, names(spec)))) {
          add("`wind_filter$variables[[", i, "]]` needs ",
              paste(need, collapse = ", "), ".")
        } else {
          check_known(unlist(spec[need]),
                      paste0("wind_filter$variables[[", i, "]]"))
        }
      }
    }
    if (!is.numeric(wf$early) || !is.numeric(wf$late) || wf$early >= wf$late)
      add("`wind_filter$early` must be a number below `late`.")
    if (!is.numeric(wf$threshold) || wf$threshold <= 0)
      add("`wind_filter$threshold` must be a positive wind speed.")
    if (isFALSE(wf$local_time) && is.data.frame(data) &&
        !is.null(wf$lon_var) && !wf$lon_var %in% names(data))
      add("`wind_filter$lon_var` = '", wf$lon_var,
          "' is not a column of the data; it is required when ",
          "`local_time = FALSE`.")
  }

  # ---- outliers ------------------------------------------------------------
  ol <- config$outliers
  if (!is.list(ol)) {
    add("`outliers` must be a list.")
  } else {
    check_known(ol$variables, "outliers$variables")
    if (!is.numeric(ol$nsd) || !length(ol$nsd) || any(!is.finite(ol$nsd)) ||
        any(ol$nsd <= 0)) {
      add("`outliers$nsd` must be positive finite thresholds.")
    } else if (is.unsorted(rev(ol$nsd))) {
      add("`outliers$nsd` should be decreasing, so a loose pass strips the ",
          "gross outliers that would distort the tighter pass's fit. Got: ",
          paste(ol$nsd, collapse = ", "))
    }
    if (is.data.frame(data)) {
      miss <- setdiff(c(ol$time_var, ol$group_var), names(data))
      if (length(miss))
        add("`outliers` time_var / group_var not in the data: ",
            paste(miss, collapse = ", "))
    }
  }

  # ---- matchup_filter ------------------------------------------------------
  mf <- config$matchup_filter
  if (!is.list(mf)) {
    add("`matchup_filter` must be a list.")
  } else {
    if (!is.list(mf$variables) || !length(mf$variables)) {
      add("`matchup_filter$variables` must be a non-empty list of ",
          "target / reference pairs.")
    } else {
      for (i in seq_along(mf$variables)) {
        spec <- mf$variables[[i]]
        need <- c("target", "reference")
        if (!is.list(spec) || length(setdiff(need, names(spec)))) {
          add("`matchup_filter$variables[[", i, "]]` needs ",
              paste(need, collapse = ", "), ".")
        } else {
          check_known(unlist(spec[need]),
                      paste0("matchup_filter$variables[[", i, "]]"))
          if (identical(spec$target, spec$reference))
            add("`matchup_filter$variables[[", i,
                "]]` compares a variable with itself.")
        }
      }
    }
    for (k in c("lower", "upper")) {
      v <- mf[[k]]
      if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v <= 0)
        add("`matchup_filter$", k, "` must be a single positive number.")
    }
    if (!is.null(mf$scale) && !mf$scale %in% c("mad", "sd"))
      add("`matchup_filter$scale` must be \"mad\" or \"sd\".")
    if (!is.null(mf$center) && !mf$center %in% c("median", "mean"))
      add("`matchup_filter$center` must be \"median\" or \"mean\".")
    if (!is.null(mf$passes) && (!is.numeric(mf$passes) || mf$passes < 1))
      add("`matchup_filter$passes` must be at least 1.")
  }

  # ---- diagnostics ---------------------------------------------------------
  dg <- config$diagnostics
  if (!is.null(dg)) {
    if (!is.list(dg)) {
      add("`diagnostics` must be a list.")
    } else {
      for (k in c("bins", "max_facets", "max_points", "width", "height", "dpi")) {
        v <- dg[[k]]
        if (!is.null(v) && (!is.numeric(v) || length(v) != 1 ||
                            !is.finite(v) || v <= 0))
          add("`diagnostics$", k, "` must be a single positive number.")
      }
    }
  }

  # ---- marss_format --------------------------------------------------------
  mf <- config$marss_format
  if (!is.null(mf)) {
    if (!is.list(mf)) {
      add("`marss_format` must be a list.")
    } else if (isTRUE(mf$enabled)) {
      if (!is.character(mf$variables) || !length(mf$variables)) {
        add("`marss_format$variables` must be a non-empty character vector.")
      } else {
        check_known(mf$variables, "marss_format$variables")
        if (anyDuplicated(mf$variables))
          add("`marss_format$variables` has duplicate entries: ",
              paste(unique(mf$variables[duplicated(mf$variables)]),
                    collapse = ", "))
      }
      ts <- mf$time_step
      if (!is.numeric(ts) || length(ts) != 1 || !is.finite(ts) || ts < 1)
        add("`marss_format$time_step` must be a single number of days >= 1.")
      hh <- mf$harmonics
      if (!is.numeric(hh) || length(hh) != 1 || !is.finite(hh) || hh < 1 ||
          hh != round(hh))
        add("`marss_format$harmonics` must be a single whole number >= 1.")
      pd <- mf$period
      if (!is.numeric(pd) || length(pd) != 1 || !is.finite(pd) || pd <= 0)
        add("`marss_format$period` must be a single positive number of days.")
      if (is.data.frame(data)) {
        miss <- setdiff(c(mf$site_var, mf$time_var), names(data))
        if (length(miss))
          add("`marss_format` site_var / time_var not in the data: ",
              paste(miss, collapse = ", "))
      }
    }
  }

  # ---- output --------------------------------------------------------------
  op <- config$output
  if (!is.null(op)) {
    if (!is.list(op)) {
      add("`output` must be a list.")
    } else {
      r <- op$output_root
      if (!is.null(r) && (!is.character(r) || length(r) != 1 || is.na(r) ||
                          !nzchar(r)))
        add("`output$output_root` must be a single non-empty string.")
      n <- op$run_name
      if (!is.null(n)) {
        if (!is.character(n) || length(n) != 1 || is.na(n) || !nzchar(n)) {
          add("`output$run_name` must be a single non-empty string, or NULL.")
        } else if (grepl("[/\\\\]", n)) {
          # A separator here would silently nest the run somewhere other than
          # `output_root`, which is exactly what a reader of the config would
          # not expect. Use `output_root` to change where runs live.
          add("`output$run_name` must not contain a path separator: '", n,
              "'. Set `output$output_root` to change where runs are written.")
        }
      }
    }
  }

  # ---- match_up_variance ---------------------------------------------------
  mv <- config$match_up_variance
  if (!is.null(mv)) {
    if (!is.list(mv)) {
      add("`match_up_variance` must be a list.")
    } else if (isTRUE(mv$enabled)) {
      if (!is.character(mv$targets) || !length(mv$targets))
        add("`match_up_variance$targets` must be a non-empty character vector.")
      else check_known(mv$targets, "match_up_variance$targets")
      if (!is.character(mv$references) || !length(mv$references))
        add("`match_up_variance$references` must be a non-empty character vector.")
      else check_known(mv$references, "match_up_variance$references")
      if (!is.null(mv$pairs) && !mv$pairs %in% c("cross", "elementwise"))
        add("`match_up_variance$pairs` must be \"cross\" or \"elementwise\".")
      if (identical(mv$pairs, "elementwise") &&
          length(mv$targets) != length(mv$references))
        add("`match_up_variance$targets` and `references` must be the same ",
            "length when pairs = \"elementwise\".")
      # `by` must survive into the pivot, so it has to be an id column. With the
      # default id_cols (everything but variable/value) that means a data column.
      if (!is.null(mv$by)) {
        pool <- if (is.null(mv$id_cols)) {
          if (is.data.frame(data)) setdiff(names(data),
                                           c("variable", "value", "stat", "radius_m"))
          else NULL
        } else mv$id_cols
        if (!is.null(pool) && !all(mv$by %in% pool))
          add("`match_up_variance$by` must be among the id columns. ",
              "Not available: ", paste(setdiff(mv$by, pool), collapse = ", "))
      }
      if (!is.null(mv$valid_range) &&
          (!is.numeric(mv$valid_range) || length(mv$valid_range) != 2 ||
           any(!is.finite(mv$valid_range)) ||
           mv$valid_range[1] >= mv$valid_range[2]))
        add("`match_up_variance$valid_range` must be two increasing finite ",
            "numbers, or NULL.")
    }
  }

  if (length(problems)) {
    stop("Invalid preprocessing configuration (", length(problems),
         " problem(s)):\n", paste0("  - ", problems, collapse = "\n"),
         call. = FALSE)
  }

  invisible(config)
}


#' Read a preprocessing configuration
#'
#' Accepts either a path to a config R file or an already-built list. A path is
#' sourced into a private environment, so nothing lands in the global
#' environment and the file can live anywhere.
#'
#' @param config Path to a config file defining `preprocess_config`, or the list
#'   itself.
#' @return The `preprocess_config` list.
#' @export
read_config <- function(config) {
  if (is.list(config)) return(config)
  if (!is.character(config) || length(config) != 1)
    stop("`config` must be a file path or a configuration list.")
  if (!file.exists(config)) stop("Config file not found: ", config)

  env <- new.env(parent = globalenv())
  sys.source(normalizePath(config), envir = env)
  if (!exists("preprocess_config", envir = env, inherits = FALSE))
    stop("Config file '", config, "' does not define `preprocess_config`.")
  get("preprocess_config", envir = env, inherits = FALSE)
}


# Internal: normalise one stage's heterogeneous return list into a common record.
# The stages return different elements (`fit` only from the match-up filter,
# `plot_values` only from the wind filter), so anything absent stays NULL.
.stage_record <- function(stage, variable, res, n_in) {

  removed <- res$removed
  if (is.null(removed)) removed <- res$data[0, , drop = FALSE]

  plots <- list(res$plot, res$plot_values, res$plot_resid)
  names(plots) <- c("plot", "plot_values", "plot_resid")
  plots <- plots[!vapply(plots, is.null, logical(1))]

  list(stage     = stage,
       variable  = variable,
       n_in      = n_in,
       n_out     = nrow(res$data),
       n_removed = n_in - nrow(res$data),
       removed   = removed,
       summary   = res$summary,
       fit       = res$fit,
       plots     = plots)
}


#' Run the full preprocessing pipeline
#'
#' Applies the configured quality control stages in order -- variable formatting,
#' unit conversion, physical range, wind / solar time, seasonal outliers, and
#' match-up screening -- then optionally summarises the variability between
#' targets and references on the cleaned data.
#'
#' Unlike calling the stages by hand, this keeps every diagnostic each stage
#' produces: the removed rows (tagged with the stage that removed them), the
#' per-stage summaries and fits, and the per-stage plots. That record is what
#' `build_diagnostics()` in `00_diagnostics.R` turns into the tables and figures.
#'
#' @param data Raw long format data frame, as extracted by `coastal_sst_data`.
#' @param config Path to a config file defining `preprocess_config`, or the list.
#' @param group_id Optional override of `config$group_id`. Use `NA` to mean "not
#'   supplied"; pass `NULL` explicitly to force pooled-only diagnostics.
#' @param diagnostics If TRUE (default) and `config$diagnostics$enabled` is TRUE,
#'   build the whole-dataset diagnostic tables and plots. Requires
#'   `00_diagnostics.R` to have been sourced.
#' @param variance If TRUE (default) and `config$match_up_variance$enabled` is
#'   TRUE, run `match_up_variance()` on the cleaned data.
#' @param marss If TRUE (default) and `config$marss_format$enabled` is TRUE,
#'   build the MARSS observation and harmonics matrices. Requires
#'   `01_format_data_for_marss.R` to have been sourced.
#' @param verbose If TRUE (default), report progress as each stage runs.
#'
#' @return A list with `data` (the cleaned long frame), `removed` (every removed
#'   row, tagged with `stage` and `stage_variable`), `stages` (the per-stage
#'   records), `removal_summary`, `diagnostics`, `match_up_variance`, `marss`
#'   and the resolved `config`.
#' @export
run_pipeline <- function(data,
                         config,
                         group_id    = NA,
                         diagnostics = TRUE,
                         variance    = TRUE,
                         marss       = TRUE,
                         verbose     = TRUE) {

  stopifnot(is.data.frame(data))
  say <- function(...) if (verbose) message(...)

  cfg <- read_config(config)

  # NA means "not supplied"; NULL is a real value meaning pooled-only. Assign
  # through `[` -- `cfg$group_id <- NULL` would delete the key rather than set
  # it, and validate_config would then report it as missing.
  if (!identical(group_id, NA)) {
    # `match_up_variance$by` is written as `by <- group_id` in the config, so it
    # binds at source time. If it is still tracking `group_id`, carry the
    # override through to it; a `by` deliberately set to something else is left
    # alone.
    if (identical(cfg$match_up_variance$by, cfg$group_id))
      cfg$match_up_variance["by"] <- list(group_id)
    cfg["group_id"] <- list(group_id)
  }

  validate_config(cfg, data)
  gid <- cfg$group_id

  stages  <- list()
  record  <- function(stage, variable, res, n_in) {
    stages[[length(stages) + 1L]] <<- .stage_record(stage, variable, res, n_in)
    res$data
  }

  # ---- 1. variable / value formatting --------------------------------------
  say("Formatting variables (", length(cfg$value_var$variables), " defined) ...")
  dat <- get_variable_value_format(data,
                                   cfg$value_var$variables,
                                   cfg$value_var$value_var)
  n_formatted <- nrow(dat)
  say("  ", format(n_formatted, big.mark = ","), " rows retained of ",
      format(nrow(data), big.mark = ","), ".")

  # ---- 2. unit conversion ---------------------------------------------------
  say("Converting units (", cfg$convert_units$from, " -> ",
      cfg$convert_units$to, ") ...")
  dat <- convert_temperature(
    dat, cfg$convert_units$variables,
    from        = cfg$convert_units$from,
    to          = cfg$convert_units$to,
    value_var   = cfg$convert_units$value_var,
    check_range = cfg$convert_units$check_range,
    strict      = cfg$convert_units$strict,
    quiet       = cfg$convert_units$quiet %||% !verbose)

  # ---- 3. physical range ----------------------------------------------------
  for (v in cfg$physics_filter$variables) {
    say("Physical range screening: ", v, " ...")
    n_in <- nrow(dat)
    res  <- trim_implausibe_values(
      dat, v,
      lower_threshold = cfg$physics_filter$lower_threshold,
      upper_threshold = cfg$physics_filter$upper_threshold,
      keep_na         = cfg$physics_filter$keep_na %||% TRUE)
    dat <- record("physics_filter", v, res, n_in)
    say("  removed ", format(n_in - nrow(dat), big.mark = ","), " row(s).")
  }

  # ---- 4. wind / solar time -------------------------------------------------
  for (spec in cfg$wind_filter$variables) {
    say("Wind / solar-time screening: ", spec$target, " ...")
    n_in <- nrow(dat)
    res  <- wind_filter(
      dat, spec,
      early           = cfg$wind_filter$early,
      late            = cfg$wind_filter$late,
      threshold       = cfg$wind_filter$threshold,
      local_time      = cfg$wind_filter$local_time,
      lon_var         = cfg$wind_filter$lon_var,
      id_cols         = cfg$wind_filter$id_cols,
      drop_covariates = cfg$wind_filter$drop_covariates %||% FALSE,
      keep_na         = cfg$wind_filter$keep_na %||% TRUE)
    dat <- record("wind_filter", spec$target, res, n_in)
    say("  removed ", format(n_in - nrow(dat), big.mark = ","), " row(s).")
  }

  # ---- 5. seasonal outliers -------------------------------------------------
  for (v in cfg$outliers$variables) {
    say("Seasonal outlier screening: ", v, " ...")
    n_in <- nrow(dat)
    res  <- trim_ecostress_outliers(
      dat, v,
      thresholds = cfg$outliers$nsd,
      time_var   = cfg$outliers$time_var,
      group_var  = cfg$outliers$group_var,
      harmonics  = cfg$outliers$harmonics,
      robust     = cfg$outliers$robust,
      keep_na    = cfg$outliers$keep_na %||% TRUE)
    dat <- record("outliers", v, res, n_in)
    say("  removed ", format(n_in - nrow(dat), big.mark = ","), " row(s).")
  }

  # ---- 6. match-up screening ------------------------------------------------
  for (spec in cfg$matchup_filter$variables) {
    say("Match-up screening: ", spec$target, " vs ", spec$reference, " ...")
    n_in <- nrow(dat)
    res  <- match_up_filter(
      dat, spec,
      lower     = cfg$matchup_filter$lower,
      upper     = cfg$matchup_filter$upper,
      slope     = cfg$matchup_filter$slope,
      intercept = cfg$matchup_filter$intercept,
      id_cols   = cfg$matchup_filter$id_cols,
      by        = cfg$matchup_filter$by,
      passes    = cfg$matchup_filter$passes,
      scale     = cfg$matchup_filter$scale,
      center    = cfg$matchup_filter$center,
      keep_na   = cfg$matchup_filter$keep_na %||% TRUE,
      min_n     = cfg$matchup_filter$min_n)
    dat <- record("matchup_filter", spec$target, res, n_in)
    say("  removed ", format(n_in - nrow(dat), big.mark = ","), " row(s).")
  }

  # ---- collect the removals ------------------------------------------------
  removed <- .bind_removed(stages)

  removal_summary <- data.frame(
    step          = seq_along(stages),
    stage         = vapply(stages, `[[`, character(1), "stage"),
    variable      = vapply(stages, `[[`, character(1), "variable"),
    n_in          = vapply(stages, `[[`, numeric(1), "n_in"),
    n_out         = vapply(stages, `[[`, numeric(1), "n_out"),
    n_removed     = vapply(stages, `[[`, numeric(1), "n_removed"),
    stringsAsFactors = FALSE)
  removal_summary$pct_removed <-
    round(100 * removal_summary$n_removed / pmax(removal_summary$n_in, 1), 3)

  say("Pipeline complete: ", format(nrow(dat), big.mark = ","), " of ",
      format(n_formatted, big.mark = ","), " formatted rows retained (",
      format(sum(removal_summary$n_removed), big.mark = ","), " removed).")

  out <- list(data              = dat,
              removed           = removed,
              stages            = stages,
              removal_summary   = removal_summary,
              n_formatted       = n_formatted,
              diagnostics       = NULL,
              match_up_variance = NULL,
              marss             = NULL,
              config            = cfg)

  # ---- 7. match-up variance on the cleaned data -----------------------------
  mv <- cfg$match_up_variance
  if (variance && isTRUE(mv$enabled)) {
    say("Match-up variance on the cleaned data ...")
    present <- unique(dat$variable)
    tg <- intersect(mv$targets, present)
    rf <- intersect(mv$references, present)
    if (!length(tg) || !length(rf)) {
      warning("Skipping match_up_variance: no configured ",
              if (!length(tg)) "targets" else "references",
              " survive in the cleaned data.")
    } else {
      if (length(tg) < length(mv$targets) || length(rf) < length(mv$references))
        warning("match_up_variance: dropping variable(s) absent from the ",
                "cleaned data: ",
                paste(setdiff(c(mv$targets, mv$references), present),
                      collapse = ", "))
      out$match_up_variance <- match_up_variance(
        dat, targets = tg, references = rf,
        pairs       = mv$pairs %||% "cross",
        id_cols     = mv$id_cols,
        by          = mv$by,
        valid_range = mv$valid_range,
        trim        = mv$trim,
        min_n       = mv$min_n %||% 10)
    }
  }

  # ---- 8. MARSS formatting --------------------------------------------------
  mf <- cfg$marss_format
  if (marss && isTRUE(mf$enabled)) {
    if (!exists("build_marss_inputs", mode = "function")) {
      warning("Skipping MARSS formatting: `build_marss_inputs()` not found. ",
              "Source R/01_format_data_for_marss.R first.")
    } else {
      present <- unique(dat$variable)
      vars    <- intersect(mf$variables, present)   # keeps the config's order
      if (!length(vars)) {
        warning("Skipping MARSS formatting: none of the configured variables ",
                "survive in the cleaned data.")
      } else {
        if (length(vars) < length(mf$variables))
          warning("MARSS formatting: dropping variable(s) absent from the ",
                  "cleaned data: ",
                  paste(setdiff(mf$variables, present), collapse = ", "))
        say("Formatting for MARSS (", length(vars), " variable(s), ",
            mf$harmonics %||% 2, " harmonic(s), ",
            mf$time_step %||% 1, "-day columns) ...")
        out$marss <- build_marss_inputs(
          dat, variables = vars,
          time_step = mf$time_step %||% 1,
          harmonics = mf$harmonics %||% 2,
          period    = mf$period    %||% 365.25,
          site_var  = mf$site_var  %||% "point_id",
          time_var  = mf$time_var  %||% "time")
        say("  ", nrow(out$marss$ts_matrix), " series x ",
            ncol(out$marss$ts_matrix), " time steps; harmonics ",
            nrow(out$marss$harmonics), " x ", ncol(out$marss$harmonics), ".")
        if (nrow(out$marss$dropped))
          say("  ", nrow(out$marss$dropped),
              " site/variable series had no data and were dropped.")
      }
    }
  }

  # ---- 9. diagnostics -------------------------------------------------------
  if (diagnostics && !isFALSE(cfg$diagnostics$enabled)) {
    if (!exists("build_diagnostics", mode = "function")) {
      warning("Skipping diagnostics: `build_diagnostics()` not found. ",
              "Source R/00_diagnostics.R first.")
    } else {
      say("Building diagnostics",
          if (is.null(gid)) " (pooled only)" else paste0(" by ", gid), " ...")
      out$diagnostics <- build_diagnostics(out, cfg)
    }
  }

  out
}


# Internal: stack every stage's removed rows into one frame.
# The stages report different diagnostic columns (resid, reason, ref_value,
# solar_time, wind_speed); bind_rows fills the gaps with NA so a single frame can
# answer "what was removed, from where, and why".
.bind_removed <- function(stages) {

  parts <- lapply(stages, function(s) {
    r <- s$removed
    if (is.null(r) || !nrow(r)) return(NULL)
    r$stage          <- s$stage
    r$stage_variable <- s$variable
    r
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]

  if (!length(parts)) {
    return(data.frame(variable = character(), value = numeric(),
                      stage = character(), stage_variable = character(),
                      stringsAsFactors = FALSE))
  }

  out <- as.data.frame(dplyr::bind_rows(parts))
  # Lead with the provenance columns; the rest keep their original order.
  lead <- intersect(c("stage", "stage_variable", "variable", "value"), names(out))
  out  <- out[, c(lead, setdiff(names(out), lead)), drop = FALSE]
  rownames(out) <- NULL
  out
}


# Internal: default for a NULL config value.
`%||%` <- function(x, y) if (is.null(x)) y else x
