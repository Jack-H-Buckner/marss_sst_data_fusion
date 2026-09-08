#############################################################
#############################################################
###
### Scoring a fit against the in situ data it was not shown.
###
### A run made with `data$sites_validation` withheld the in
### situ series at some of its sites and wrote those
### observations to validation.csv. This step reads them back
### and asks how well the fit predicted them, against two
### baselines: the reference instrument as the model sees it,
### and that same instrument's raw observations.
###
### Everything is read from the run directory. observations.rds
### carries the unscaled observation matrix and its row keys,
### so which instruments reached which site -- and the raw
### reference series itself -- come out of the run rather than
### out of a data file that may have moved on since the fit.
###
### The sites scored are the ones design.rds records as held
### out. A site fitted *with* its in situ record cannot be
### scored here: its in situ trained the model, and its error
### would be an in sample number wearing an out of sample
### label.
###
### MSE and bias both carry a bootstrap interval throughout,
### off one set of replicates, so a replicate's bias and its
### MSE describe the same resampled holdout. Replicates are
### paired across models -- one resampled index set scores
### all three -- so the model to model differences shed the
### shared sampling noise and are far tighter than differencing
### two independent intervals. Where more than one site is
### involved the resampling is two stage, over sites and then
### over observations within each drawn site, so the interval
### carries site to site variation rather than treating every
### observation as exchangeable.
###
### Every interval is reported at two levels, 95% and 90% by
### default (`alpha` and `alpha2`), from the same replicates:
### the columns without a suffix are the outer level, those
### ending in `2` the inner one, and the figures draw the inner
### interval as a heavy bar inside the outer whisker. With this
### few sites the outer level is wide enough that a comparison
### can look inconclusive while the bulk of the replicates sit
### well clear of zero; the second level says so.
###
### Sites can be split by which instruments reached them --
### `group_var = "has_modis_sst"` scores sites with and without
### MODIS separately and puts an interval on the gap between
### them. `min_obs` sets how much of an instrument a site needs
### before it counts as having it, and defaults to 1: a single
### scene. Where every held out site has the instrument that
### split has nothing to contrast, and raising `min_obs` to a
### count that falls between the sites turns it into the
### comparison that is actually available -- dense coverage
### against thin.
###
#############################################################
#############################################################

`%||%` <- function(x, y) if (is.null(x)) y else x

# The three predictions, in the order they are reported.
.HOLDOUT_MODELS <- c("full", "ref_model", "ref_raw")


#' The withheld observations and what each prediction said
#'
#' Joins validation.csv to the fused reconstruction, to the reference
#' instrument's own reconstruction, and to the reference instrument's raw
#' observations, on `(site, date)`.
#'
#' The raw reference comes from `observations.rds`, which stores the unscaled
#' matrix the model was given. That keeps the comparison self contained: no
#' preprocessing output, no raw CSV, no unit conversion to get wrong.
#'
#' @param run_dir A finished run directory.
#' @param reference Variable to use as the baseline, or NULL for the run's
#'   `reconstruction$instruments`.
#' @return A data frame of `site`, `date`, `obs`, `pred_full`, `se_full`,
#'   `pred_ref_model`, `pred_ref_raw`, with a `reference` attribute naming the
#'   baseline variable and its tag.
#' @export
holdout_predictions <- function(run_dir, reference = NULL) {

  design <- readRDS(file.path(run_dir, "design.rds"))
  val_sites <- design$sites_validation
  if (is.null(val_sites) || !length(val_sites))
    stop("Nothing was held out of this run: design.rds has no ",
         "`sites_validation`. There is nothing to score.\n  ", run_dir,
         call. = FALSE)

  ref <- reference %||% design$reconstruction$instruments[1]
  if (is.na(ref) || is.null(ref))
    stop("No reference instrument: the run reconstructed no instrument view, ",
         "and none was named.", call. = FALSE)
  tag <- design$instruments[[ref]]$tag %||% ref

  ref_file <- file.path(run_dir, paste0("states_", tag, "_only.csv"))
  if (!file.exists(ref_file))
    stop("No reconstruction for the reference instrument '", ref, "': ",
         basename(ref_file), " is not in the run directory. Re-run the ",
         "reconstruction with '", ref, "' in `reconstruction$instruments`.",
         call. = FALSE)

  truth <- utils::read.csv(file.path(run_dir, "validation.csv"),
                           stringsAsFactors = FALSE)
  truth <- truth[truth$variable %in% design$validation_variables,
                 c("site", "date", "value")]
  names(truth)[3] <- "obs"

  states <- utils::read.csv(file.path(run_dir, "states.csv"),
                            stringsAsFactors = FALSE)
  states <- states[states$variable == "sst", c("site", "date", "value", "se")]
  names(states)[3:4] <- c("pred_full", "se_full")

  ref_states <- utils::read.csv(ref_file, stringsAsFactors = FALSE)
  ref_states <- ref_states[ref_states$variable == ref, c("site", "date", "value")]
  names(ref_states)[3] <- "pred_ref_model"

  raw <- .observed_long(run_dir, ref)
  names(raw)[3] <- "pred_ref_raw"

  out <- merge(truth, states,     by = c("site", "date"))
  out <- merge(out,   ref_states, by = c("site", "date"))
  out <- merge(out,   raw,        by = c("site", "date"))
  out <- out[order(out$site, out$date), , drop = FALSE]
  rownames(out) <- NULL

  # A join that quietly drops rows shrinks a denominator instead of failing,
  # and the MSE it produces still looks entirely reasonable.
  cols <- c("obs", "pred_full", "pred_ref_model", "pred_ref_raw")
  if (nrow(out) != nrow(truth))
    stop("Joining the predictions to the withheld observations lost ",
         nrow(truth) - nrow(out), " of ", nrow(truth), " row(s). The run ",
         "directory is inconsistent.", call. = FALSE)
  if (anyNA(out[cols]))
    stop("Missing values among the joined predictions; MSE would be computed ",
         "on a silently smaller sample.", call. = FALSE)

  # A held out site with no withheld observation is not an inconsistent run
  # directory, it is a site whose in situ record does not reach the fitted
  # dates: holding it out removed nothing and left nothing to score. That is
  # reported and dropped. A *scored* site that was never held out is the
  # genuine inconsistency, and still stops the run.
  scored <- sort(unique(out$site))
  extra  <- setdiff(scored, val_sites)
  if (length(extra))
    stop("validation.csv carries observations at site(s) design.rds did not ",
         "hold out: ", paste(extra, collapse = ", "), ". Scoring them would ",
         "put in sample error under an out of sample label. The run ",
         "directory is inconsistent.\n  ", run_dir, call. = FALSE)
  if (!length(scored))
    stop("None of the held out site(s) -- ", paste(val_sites, collapse = ", "),
         " -- has a ", paste(design$validation_variables, collapse = "/"),
         " observation inside the run's dates, so there is nothing to score.",
         "\n  ", run_dir, call. = FALSE)

  attr(out, "reference")     <- ref
  attr(out, "reference_tag") <- tag
  attr(out, "run_dir")       <- run_dir
  attr(out, "scored_sites")  <- scored
  attr(out, "empty_sites")   <- setdiff(val_sites, scored)
  out
}


# Internal: one variable's observed series out of observations.rds, long. This
# is the unscaled matrix the model was handed, so the values are degrees C.
.observed_long <- function(run_dir, variable) {

  obs <- readRDS(file.path(run_dir, "observations.rds"))
  r   <- which(obs$row_var_keys == variable)
  if (!length(r))
    stop("No '", variable, "' rows in observations.rds: that instrument ",
         "reached none of the fitted sites.", call. = FALSE)

  out <- data.frame(
    site  = rep(obs$row_site_keys[r], each = ncol(obs$ts_matrix)),
    date  = rep(as.character(obs$col_dates), times = length(r)),
    value = as.vector(t(obs$ts_matrix[r, , drop = FALSE])),
    stringsAsFactors = FALSE)
  out[!is.na(out$value), , drop = FALSE]
}


#' Which instruments reached which site
#'
#' Counts the observed cells each instrument contributed at each site, from the
#' matrix the model was actually given. `has_<variable>` is TRUE where that
#' instrument contributed at least `min_obs` of them -- which is what makes
#' `group_var = "has_modis_sst"` a meaningful split rather than a guess from
#' the raw data.
#'
#' `min_obs` is the lever on what "has" means. At the default of 1 a single
#' scene counts, so `has_` asks only whether the instrument reached the site at
#' all. That lumps a site with a handful of scenes in with one that has
#' hundreds, and those are not the same claim about coverage: raise `min_obs`
#' to ask instead which sites the instrument covered densely enough to be worth
#' fusing. It applies to every instrument's `has_` column at once, since a
#' threshold that means one thing for MODIS and another for ECOSTRESS would not
#' be readable off the table.
#'
#' @param run_dir A finished run directory.
#' @param sites Restrict to these sites, or NULL for every fitted site.
#' @param min_obs Observations an instrument needs at a site for that site's
#'   `has_<variable>` to be TRUE. 1, the default, is "anything at all".
#' @return A data frame with one row per site: `site`, then `n_<variable>` and
#'   `has_<variable>` for every instrument in the run, with a `min_obs`
#'   attribute recording the threshold used.
#' @export
site_coverage <- function(run_dir, sites = NULL, min_obs = 1) {

  if (length(min_obs) != 1 || !is.numeric(min_obs) || is.na(min_obs) ||
      min_obs < 1)
    stop("`min_obs` must be a single number of observations, at least 1. ",
         "It is the count at which an instrument counts as having reached a ",
         "site.", call. = FALSE)

  obs <- readRDS(file.path(run_dir, "observations.rds"))
  n   <- rowSums(!is.na(obs$ts_matrix))

  keep  <- if (is.null(sites)) rep(TRUE, length(obs$row_site_keys)) else
    obs$row_site_keys %in% sites
  site_levels <- sort(unique(obs$row_site_keys[keep]))
  variables   <- unique(obs$row_var_keys)

  out <- data.frame(site = site_levels, stringsAsFactors = FALSE)
  for (v in variables) {
    counts <- setNames(rep(0L, length(site_levels)), site_levels)
    r <- which(keep & obs$row_var_keys == v)
    if (length(r)) counts[obs$row_site_keys[r]] <- as.integer(n[r])
    out[[paste0("n_", v)]]   <- unname(counts[site_levels])
    out[[paste0("has_", v)]] <- unname(counts[site_levels]) >= min_obs
  }
  attr(out, "min_obs") <- min_obs
  out
}


#' Skill of one prediction against the withheld observations
#'
#' `bias` is signed, prediction minus observation, so positive means the
#' prediction runs warm. `cover95` is the share of observations inside the
#' prediction's 95% interval and is only meaningful where `se` is a calibrated
#' predictive standard error -- a reconstruction run with
#' `parameter_uncertainty = FALSE` does not produce one.
#'
#' @param obs Observed values.
#' @param pred Predicted values.
#' @param se Predictive standard errors, or NULL to skip `cover95`.
#' @return A one row data frame: `n`, `bias`, `MSE`, `RMSE`, `MAE`, `cor`,
#'   `cover95`.
#' @export
skill_metrics <- function(obs, pred, se = NULL) {
  e <- pred - obs
  data.frame(n = length(e), bias = mean(e), MSE = mean(e^2),
             RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)),
             cor = stats::cor(obs, pred),
             cover95 = if (is.null(se)) NA_real_ else mean(abs(e) <= 1.96 * se))
}


# Internal: signed errors of the three predictions, one column each. Squared,
# these are what MSE averages; unsquared, what bias averages.
.errors <- function(preds) {
  cbind(full      = preds$pred_full      - preds$obs,
        ref_model = preds$pred_ref_model - preds$obs,
        ref_raw   = preds$pred_ref_raw   - preds$obs)
}


# Internal: bootstrap MSE and bias for each column of `err`.
#
# `units` is a list of row-index vectors, one per resampling unit -- a site.
# With several units the resampling is two stage: draw units with replacement,
# then draw observations with replacement within each drawn unit, so the
# spread across units enters the interval. With a single unit it reduces to an
# ordinary bootstrap over that unit's observations.
#
# One index set per replicate serves every column and both statistics, which
# is what pairs the models and makes their differences comparable, and what
# keeps a replicate's bias and MSE describing the same resampled holdout.
.boot_stats <- function(err, units, B) {

  dn  <- list(NULL, colnames(err))
  out <- list(MSE  = matrix(NA_real_, B, ncol(err), dimnames = dn),
              bias = matrix(NA_real_, B, ncol(err), dimnames = dn))
  sq  <- err^2
  nu  <- length(units)
  for (b in seq_len(B)) {
    drawn <- if (nu == 1L) units else units[sample.int(nu, nu, replace = TRUE)]
    idx   <- unlist(lapply(drawn, function(r) sample(r, length(r),
                                                     replace = TRUE)),
                    use.names = FALSE)
    out$MSE[b, ]  <- colMeans(sq[idx, , drop = FALSE])
    out$bias[b, ] <- colMeans(err[idx, , drop = FALSE])
  }
  out
}


# Internal: percentile interval. Not BCa -- the resampling is hierarchical and
# paired, so boot::boot.ci()'s acceleration term does not apply to it cleanly,
# and a quantile of the replicates is the thing that is easy to state and check.
.ci <- function(x, alpha) stats::quantile(x, c(alpha / 2, 1 - alpha / 2),
                                          names = FALSE, na.rm = TRUE)


# Internal: point estimates plus intervals for one set of rows.
#
# The resampling units are the sites present in `preds`, derived here rather
# than passed in so the indices cannot fall out of step with the rows: one
# site gives an ordinary bootstrap over its observations, several give the two
# stage scheme over sites and then observations.
#
# Both levels come off the same replicates, so they are two quantiles of one
# bootstrap rather than two bootstraps that could disagree with each other.
.score_rows <- function(preds, B, alpha, alpha2) {

  units <- unname(split(seq_len(nrow(preds)), preds$site))

  est <- rbind(
    cbind(model = "full",      skill_metrics(preds$obs, preds$pred_full,
                                             preds$se_full)),
    cbind(model = "ref_model", skill_metrics(preds$obs, preds$pred_ref_model)),
    cbind(model = "ref_raw",   skill_metrics(preds$obs, preds$pred_ref_raw)))

  draws <- .boot_stats(.errors(preds), units, B)

  est$dMSE <- est$MSE - est$MSE[est$model == "full"]
  est <- .add_ci(est, draws, alpha,  "")
  est <- .add_ci(est, draws, alpha2, "2")

  list(estimates = est, draws = draws)
}


# Internal: one interval level's bounds, appended with `suffix` on each name.
.add_ci <- function(est, draws, alpha, suffix) {

  mse  <- t(apply(draws$MSE,  2, .ci, alpha = alpha))
  bias <- t(apply(draws$bias, 2, .ci, alpha = alpha))
  # Paired, replicate by replicate: the shared sampling noise cancels.
  d    <- t(apply(draws$MSE - draws$MSE[, "full"], 2, .ci, alpha = alpha))
  nm   <- function(x) paste0(x, suffix)

  est[[nm("bias_lo")]] <- bias[est$model, 1]
  est[[nm("bias_hi")]] <- bias[est$model, 2]
  est[[nm("mse_lo")]]  <- mse[est$model, 1]
  est[[nm("mse_hi")]]  <- mse[est$model, 2]
  # RMSE bounds are the square roots of the MSE bounds; the transform is
  # monotone, so percentile intervals map through it exactly.
  est[[nm("rmse_lo")]] <- sqrt(mse[est$model, 1])
  est[[nm("rmse_hi")]] <- sqrt(mse[est$model, 2])
  est[[nm("dmse_lo")]] <- ifelse(est$model == "full", NA_real_, d[est$model, 1])
  est[[nm("dmse_hi")]] <- ifelse(est$model == "full", NA_real_, d[est$model, 2])
  est
}


# Columns are ordered so an interval reads outwards from the estimate: the
# inner level sits between the point estimate and the outer bounds.
.SKILL_COLS <- c("n",
                 "bias_lo", "bias_lo2", "bias", "bias_hi2", "bias_hi",
                 "mse_lo", "mse_lo2", "MSE", "mse_hi2", "mse_hi",
                 "rmse_lo", "rmse_lo2", "RMSE", "rmse_hi2", "rmse_hi",
                 "MAE", "cor", "cover95",
                 "dmse_lo", "dmse_lo2", "dMSE", "dmse_hi2", "dmse_hi")


#' Score a run against its withheld in situ data
#'
#' Pooled across the held out sites, per site, and -- when `group_var` names a
#' site attribute -- per group of sites, with the gap between two groups given
#' its own interval.
#'
#' A held out site whose in situ record does not reach the run's dates has no
#' withheld observation to score and is reported and dropped, not scored as an
#' empty set; `settings$empty_sites` names any such site.
#'
#' @param run_dir A finished run directory.
#' @param group_var A column of `site_coverage()`, typically
#'   `"has_modis_sst"`, or NULL for no grouping. A split with only one
#'   non-empty group is reported and then skipped: there is nothing to contrast.
#' @param min_obs Observations an instrument needs at a site before that site
#'   counts as having it; passed to `site_coverage()`. 1, the default, means a
#'   single scene counts. Raising it is what turns a `has_` split that every
#'   held out site falls on the same side of into a contrast between the sites
#'   the instrument covered densely and the rest.
#' @param reference Baseline instrument, or NULL for the run's reconstruction.
#' @param B Bootstrap replicates.
#' @param seed Passed to `set.seed()`, so a run is reproducible.
#' @param alpha 1 - alpha is the outer interval level, reported in the columns
#'   without a suffix.
#' @param alpha2 1 - alpha2 is the second, inner interval level, reported in
#'   the columns ending in `2` and drawn as the heavy bar in the figures. Both
#'   levels are quantiles of the same replicates.
#' @param verbose Print the tables and what was and was not scored.
#' @return A list of `overall`, `by_site`, `by_group`, `group_contrast`,
#'   `preds`, `coverage`, `draws` and `settings`.
#' @export
validate_holdout <- function(run_dir, group_var = NULL, min_obs = 1,
                             reference = NULL,
                             B = 10000, seed = 1, alpha = 0.05, alpha2 = 0.10,
                             verbose = TRUE) {

  preds  <- holdout_predictions(run_dir, reference)
  design <- readRDS(file.path(run_dir, "design.rds"))
  ref    <- attr(preds, "reference")
  tag    <- attr(preds, "reference_tag")

  # The sites actually scored: those held out that have a withheld observation
  # inside the run's dates. `empty` are the rest -- held out, but with nothing
  # of theirs in the window to score against.
  held_out   <- design$sites_validation
  val_sites  <- attr(preds, "scored_sites")
  empty      <- attr(preds, "empty_sites")
  coverage   <- site_coverage(run_dir, val_sites, min_obs = min_obs)
  all_fitted <- readRDS(file.path(run_dir, "observations.rds"))$site_levels
  not_scored <- setdiff(all_fitted, held_out)

  if (verbose) {
    message("Run: ", run_dir)
    message("  scored (held out): ", length(val_sites), " site(s) -- ",
            paste(val_sites, collapse = ", "))
    if (length(empty))
      message("  held out but not scorable, no ",
              paste(design$validation_variables, collapse = "/"),
              " inside the run's dates: ", paste(empty, collapse = ", "))
    if (length(not_scored))
      message("  NOT scored (fitted with their ",
              paste(design$validation_variables, collapse = "/"),
              " record): ", paste(not_scored, collapse = ", "))
    message("  ", nrow(preds), " withheld observation(s); baseline '", ref,
            "'.")
  }

  set.seed(seed)

  overall <- .score_rows(preds, B, alpha, alpha2)
  ov <- overall$estimates
  ov$site <- "(all)"
  ov <- ov[c("site", "model", .SKILL_COLS)]

  by_site <- list()
  site_draws <- list()
  for (s in val_sites) {
    r <- preds$site == s
    one <- .score_rows(preds[r, , drop = FALSE], B, alpha, alpha2)
    e <- one$estimates
    e$site <- s
    by_site[[s]] <- e[c("site", "model", .SKILL_COLS)]
    site_draws[[s]] <- one$draws
  }
  by_site <- do.call(rbind, by_site)
  rownames(by_site) <- NULL

  # ---- grouped by a site attribute ----------------------------------------
  by_group <- NULL; group_draws <- list(); contrast <- NULL
  if (!is.null(group_var)) {
    if (!group_var %in% names(coverage))
      stop("`group_var` = '", group_var, "' is not a site attribute. ",
           "Available: ", paste(setdiff(names(coverage), "site"),
                                collapse = ", "), call. = FALSE)

    g <- setNames(coverage[[group_var]], coverage$site)
    preds$group <- .group_label(g[preds$site], group_var, min_obs)
    levels_g <- unique(preds$group)

    if (length(levels_g) < 2) {
      if (verbose) {
        message("  '", group_var, "' takes one value across the held out ",
                "sites (", levels_g, "), so there is nothing to contrast; ",
                "the grouped comparison is skipped.")
        # The counts are the way out of this: a threshold between them splits
        # the sites even where every one of them has the instrument.
        n_col <- sub("^has_", "n_", group_var)
        if (n_col %in% names(coverage))
          message("    ", n_col, " at those sites: ",
                  paste(sprintf("%s %d", coverage$site, coverage[[n_col]]),
                        collapse = ", "),
                  " (threshold: min_obs = ", min_obs,
                  "). Raise `min_obs` to a count between them to split them.")
      }
    } else {
      out <- list()
      for (lv in levels_g) {
        r <- preds$group == lv
        one <- .score_rows(preds[r, , drop = FALSE], B, alpha, alpha2)
        e <- one$estimates
        e$group  <- lv
        e$sites  <- length(unique(preds$site[r]))
        out[[lv]] <- e[c("group", "sites", "model", .SKILL_COLS)]
        group_draws[[lv]] <- one$draws
      }
      by_group <- do.call(rbind, out)
      rownames(by_group) <- NULL
      contrast <- .group_contrast(group_draws, levels_g, alpha, alpha2,
                                  by_group)
    }
  }

  res <- list(overall = ov, by_site = by_site, by_group = by_group,
              group_contrast = contrast, preds = preds, coverage = coverage,
              draws = list(overall = overall$draws, by_site = site_draws,
                           by_group = group_draws),
              settings = list(run_dir = run_dir, reference = ref,
                              reference_tag = tag, group_var = group_var,
                              min_obs = min_obs,
                              B = B, seed = seed, alpha = alpha,
                              alpha2 = alpha2, held_out = held_out,
                              validation_variables =
                                design$validation_variables,
                              val_sites = val_sites, empty_sites = empty,
                              not_scored = not_scored,
                              labels = holdout_labels(tag)))

  if (verbose) print_holdout(res)
  res
}


# Internal: readable group names. A logical `has_modis_sst` becomes
# "with modis_sst" / "without modis_sst"; anything else prints as itself.
#
# Above a threshold of 1 the split is no longer presence and absence -- every
# site in the lean group may still have the instrument -- so the labels name
# the count instead, and a figure built from them says what it actually
# compared rather than claiming sites have no MODIS when they have eleven
# scenes of it.
.group_label <- function(x, group_var, min_obs = 1) {
  if (is.logical(x)) {
    what <- sub("^has_", "", group_var)
    if (min_obs <= 1)
      ifelse(x, paste("with", what), paste("without", what))
    else
      ifelse(x, sprintf("%s, %g+ obs", what, min_obs),
             sprintf("%s, under %g obs", what, min_obs))
  } else as.character(x)
}


# Internal: the between-group gap in dMSE, with an interval.
#
# The groups are disjoint sets of sites, so their replicates are independent
# and differencing them per replicate is legitimate. Within a replicate the
# models are still paired, so what is being differenced is two clean estimates
# of "how much worse than the fused model".
#
# The reported dMSE and contrast are the point estimates from `by_group`, not
# the means of the draws -- a bootstrap mean is not the estimate, and reporting
# it here would put a different number in this table than in that one. The
# draws supply the interval and nothing else.
.group_contrast <- function(group_draws, levels_g, alpha, alpha2, by_group) {

  if (length(levels_g) != 2) return(NULL)
  a <- group_draws[[levels_g[1]]]$MSE
  b <- group_draws[[levels_g[2]]]$MSE
  dmse <- function(g, m) by_group$dMSE[by_group$group == g & by_group$model == m]

  out <- list()
  for (m in c("ref_model", "ref_raw")) {
    da <- a[, m] - a[, "full"]
    db <- b[, m] - b[, "full"]
    ci  <- .ci(db - da, alpha)
    ci2 <- .ci(db - da, alpha2)
    est_a <- dmse(levels_g[1], m)
    est_b <- dmse(levels_g[2], m)
    out[[m]] <- data.frame(
      model = m,
      group_a = levels_g[1], dMSE_a = est_a,
      group_b = levels_g[2], dMSE_b = est_b,
      contrast = est_b - est_a,
      lo = ci[1], lo2 = ci2[1], hi2 = ci2[2], hi = ci[2],
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}


#' Display labels for the three predictions
#' @param tag The reference instrument's tag, e.g. "mur".
#' @return A named character vector over `full`, `ref_model`, `ref_raw`.
#' @export
holdout_labels <- function(tag) {
  up <- toupper(tag)
  c(full = "Full model", ref_model = paste0(up, " only (model)"),
    ref_raw = paste0(up, " only (raw)"))
}


#' Print the tables from `validate_holdout()`
#' @param res Output of `validate_holdout()`.
#' @return `res`, invisibly.
#' @export
print_holdout <- function(res) {

  s <- res$settings
  show <- function(d, keys) {
    d <- d[c(keys, "model", "n",
             "bias_lo", "bias_lo2", "bias", "bias_hi2", "bias_hi",
             "mse_lo", "mse_lo2", "MSE", "mse_hi2", "mse_hi",
             "RMSE", "MAE", "cor", "cover95",
             "dmse_lo", "dmse_lo2", "dMSE", "dmse_hi2", "dmse_hi")]
    d$model <- s$labels[d$model]
    num <- vapply(d, is.numeric, logical(1))
    d[num] <- lapply(d[num], round, 4)
    print(d, row.names = FALSE)
  }

  cat("\n=== Pooled over the ", length(s$val_sites),
      " held out site(s) ==================\n", sep = "")
  show(res$overall, "site")

  cat("\n=== By site ==============================================\n")
  show(res$by_site, "site")

  if (!is.null(res$by_group)) {
    cat("\n=== By ", s$group_var,
        if ((s$min_obs %||% 1) > 1)
          paste0(", at ", s$min_obs, " observation(s) or more") else "",
        " ===========================================\n", sep = "")
    show(res$by_group, c("group", "sites"))

    if (!is.null(res$group_contrast)) {
      cat("\nGap between the groups in how far the baseline sits above the ",
          "full model.\nA positive contrast means the fused model gains more ",
          "in '", res$group_contrast$group_b[1], "'.\n", sep = "")
      d <- res$group_contrast
      d$model <- s$labels[d$model]
      num <- vapply(d, is.numeric, logical(1))
      d[num] <- lapply(d[num], round, 4)
      print(d, row.names = FALSE)
    }
  }

  cat("\n", 100 * (1 - s$alpha), "% and ", 100 * (1 - s$alpha2),
      "% percentile intervals from ", format(s$B, big.mark = ","),
      " bootstrap replicates:\n",
      "the columns ending in 2 are the ", 100 * (1 - s$alpha2),
      "% level, the rest the ", 100 * (1 - s$alpha), "%.\n",
      "bias is signed, prediction minus observation, so an interval clear of\n",
      "zero is a prediction that runs systematically warm or cold.\n",
      "dMSE is a baseline's MSE minus the full model's, resampled in step with\n",
      "it, so a wholly positive interval means the full model is better.\n",
      sep = "")

  if (length(s$empty_sites))
    cat("\n", length(s$empty_sites), " of the ", length(s$held_out),
        " held out site(s) had no ", paste(s$validation_variables,
                                           collapse = "/"),
        " observation inside the run's dates and\nare absent from every table ",
        "above: ", paste(s$empty_sites, collapse = ", "), ".\n", sep = "")

  n_per <- table(res$preds$site)
  cat("\nThe holdout is ", nrow(res$preds), " observation(s) across ",
      length(s$val_sites), " site(s) (",
      paste(sort(unique(as.integer(n_per))), collapse = "/"),
      " each). Intervals built on this few\nsites and observations are ",
      "indicative rather than tight.\n", sep = "")

  invisible(res)
}


#' Figures for `validate_holdout()`
#'
#' @param res Output of `validate_holdout()`.
#' @return A named list of ggplot objects: `series`, `scatter`, `bias`, `mse`,
#'   `dmse` and `group_dmse`. The last is NULL when no grouping was requested
#'   or the split had only one group.
#' @export
plot_holdout <- function(res) {

  s   <- res$settings
  lab <- s$labels
  as_model <- function(d) { d$model <- factor(lab[d$model], levels = lab); d }

  long <- res$preds
  long <- do.call(rbind, lapply(.HOLDOUT_MODELS, function(m) {
    data.frame(site = long$site, date = long$date, obs = long$obs,
               model = m, pred = long[[paste0("pred_", m)]],
               stringsAsFactors = FALSE)
  }))
  long <- as_model(long)

  run <- s$run_dir
  series <- rbind(
    .series(file.path(run, "states.csv"), "sst", s$val_sites, "full"),
    .series(file.path(run, paste0("states_", s$reference_tag, "_only.csv")),
            s$reference, s$val_sites, "ref_model"))
  series <- as_model(series)

  p <- list()

  p$series <- ggplot2::ggplot(series,
      ggplot2::aes(as.Date(date), value, colour = model)) +
    ggplot2::geom_line(alpha = 0.8) +
    ggplot2::geom_point(data = res$preds,
                        ggplot2::aes(as.Date(date), obs),
                        inherit.aes = FALSE, size = 1.1) +
    ggplot2::facet_wrap(~site) +
    ggplot2::labs(x = NULL, y = "SST (deg C)", colour = NULL,
                  title = "Reconstructions and the withheld in situ data",
                  subtitle = "Points are observations the model never saw") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "bottom")

  p$scatter <- ggplot2::ggplot(long, ggplot2::aes(pred, obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::facet_grid(model ~ site) +
    ggplot2::labs(x = "Predicted (deg C)", y = "Withheld observation (deg C)",
                  title = "Predicted against withheld") +
    ggplot2::theme_classic()

  mse <- as_model(rbind(res$overall, res$by_site))
  mse$site <- factor(mse$site, levels = c("(all)", s$val_sites))

  p$bias <- ggplot2::ggplot(mse, ggplot2::aes(model, bias)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = bias_lo, ymax = bias_hi),
                           width = 0.15) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = bias_lo2, ymax = bias_hi2),
                            linewidth = 1.4) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~site, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Bias (deg C, prediction minus observation)",
                  title = "Out of sample bias against the withheld in situ",
                  subtitle = paste0(.interval_note(s),
                                    "; an interval clear of zero is a ",
                                    "systematic offset")) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  p$mse <- ggplot2::ggplot(mse, ggplot2::aes(model, MSE)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mse_lo, ymax = mse_hi),
                           width = 0.15) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = mse_lo2, ymax = mse_hi2),
                            linewidth = 1.4) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~site, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "MSE (deg C squared)",
                  title = "Out of sample MSE against the withheld in situ data",
                  subtitle = paste0(.interval_note(s),
                                    "; '(all)' resamples sites and then ",
                                    "observations")) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  p$dmse <- ggplot2::ggplot(mse[mse$model != lab[["full"]], , drop = FALSE],
                            ggplot2::aes(model, dMSE)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = dmse_lo, ymax = dmse_hi),
                           width = 0.15) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = dmse_lo2, ymax = dmse_hi2),
                            linewidth = 1.4) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~site, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "MSE minus the full model's",
                  title = "How much worse than the fused model",
                  subtitle = paste0("Above zero favours the fused model. ",
                                    .interval_note(s))) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  p$group_dmse <- if (is.null(res$by_group)) NULL else
    .plot_group_dmse(res, lab)

  p
}


# Internal: the grouped comparison -- how much the baseline loses to the fused
# model in each group of sites, with the per site values shown behind the
# group estimate so a group of three sites cannot pass itself off as a point.
.plot_group_dmse <- function(res, lab) {

  s  <- res$settings
  g  <- setNames(res$coverage[[s$group_var]], res$coverage$site)
  bs <- res$by_site
  bs$group <- .group_label(g[bs$site], s$group_var, s$min_obs %||% 1)
  bs <- bs[bs$model != "full", , drop = FALSE]
  bs$model <- factor(lab[bs$model], levels = lab)

  bg <- res$by_group[res$by_group$model != "full", , drop = FALSE]
  bg$model <- factor(lab[bg$model], levels = lab)

  n_sites <- tapply(res$by_group$sites, res$by_group$group, function(x) x[1])
  bg$label <- sprintf("%s (%d sites)", bg$group, n_sites[bg$group])
  bs$label <- sprintf("%s (%d sites)", bs$group, n_sites[bs$group])

  ggplot2::ggplot(bg, ggplot2::aes(label, dMSE)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_point(data = bs, ggplot2::aes(label, dMSE),
                        colour = "grey55", size = 1.6,
                        position = ggplot2::position_jitter(width = 0.08,
                                                            height = 0)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = dmse_lo, ymax = dmse_hi),
                           width = 0.12) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = dmse_lo2, ymax = dmse_hi2),
                            linewidth = 1.5) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::facet_wrap(~model) +
    ggplot2::labs(x = NULL, y = "MSE minus the full model's",
                  title = paste0("How much the fused model gains, split by ",
                                 sub("^has_", "", s$group_var),
                                 if ((s$min_obs %||% 1) > 1)
                                   sprintf(" (%g obs or more)", s$min_obs)
                                 else ""),
                  subtitle = paste0("Black: the group estimate, heavy bar ",
                                    100 * (1 - s$alpha2), "% and whisker ",
                                    100 * (1 - s$alpha),
                                    "%. Grey: the individual sites.",
                                    "\nAbove zero favours the fused model.")) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))
}


.interval_note <- function(s)
  paste0("Whisker ", 100 * (1 - s$alpha), "%, heavy bar ",
         100 * (1 - s$alpha2), "% bootstrap interval")


# Internal: one reconstruction file as a plotting series.
.series <- function(path, variable, sites, model) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$variable == variable & d$site %in% sites, ]
  data.frame(site = d$site, date = d$date, value = d$value, model = model,
             stringsAsFactors = FALSE)
}


#' Write the validation tables, draws and figures into a run directory
#'
#' @param res Output of `validate_holdout()`.
#' @param outdir Where to write, defaulting to the run directory.
#' @param figures Write the PNGs as well as the tables.
#' @param verbose Say what was written.
#' @return The paths written, invisibly.
#' @export
write_holdout_validation <- function(res, outdir = NULL, figures = TRUE,
                                     verbose = TRUE) {

  outdir <- outdir %||% res$settings$run_dir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  paths <- character()

  put <- function(name, obj) {
    p <- file.path(outdir, name)
    if (grepl("\\.csv$", name)) utils::write.csv(obj, p, row.names = FALSE)
    else saveRDS(obj, p)
    paths <<- c(paths, p)
  }

  put("validation_skill_overall.csv", res$overall)
  put("validation_skill_by_site.csv", res$by_site)
  put("validation_site_coverage.csv", res$coverage)
  if (!is.null(res$by_group))       put("validation_skill_by_group.csv",
                                        res$by_group)
  if (!is.null(res$group_contrast)) put("validation_group_contrast.csv",
                                        res$group_contrast)
  put("validation_bootstrap.rds",
      res[c("draws", "preds", "coverage", "settings")])

  if (figures) {
    p <- plot_holdout(res)
    sizes <- list(series = c(10, 7), scatter = c(10, 7), bias = c(10, 7),
                  mse = c(10, 7), dmse = c(10, 7), group_dmse = c(9, 5))
    files <- c(series = "validation_timeseries.png",
               scatter = "validation_scatter.png",
               bias = "validation_bias_ci.png",
               mse = "validation_mse_ci.png",
               dmse = "validation_dmse_ci.png",
               group_dmse = "validation_group_dmse.png")
    for (k in names(files)) {
      if (is.null(p[[k]])) next
      f <- file.path(outdir, files[[k]])
      ggplot2::ggsave(f, p[[k]], width = sizes[[k]][1], height = sizes[[k]][2],
                      dpi = 150)
      paths <- c(paths, f)
    }
  }

  if (verbose)
    message("Wrote ", paste(basename(paths), collapse = ", "), "\n  into ",
            outdir, ".")

  invisible(paths)
}
