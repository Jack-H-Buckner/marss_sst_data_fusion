#############################################################
#############################################################
###
### Diagnostics for the preprocessing pipeline.
###
### The stage functions in 00_preprocessing.R each produce a
### plot of what they screened. Those are per-process views:
### they show the residual or the covariate the filter acted
### on, not the variable itself.
###
### This file adds the complementary whole-dataset view. For
### every variable it reports the distribution of the full
### sample and the breakdown across groups (`group_id`), in
### both figures and tables, always showing the values that
### were removed alongside those that were kept.
###
### Source after R/00_preprocessing.R -- RETAINED_PAL and the
### grey / firebrick encoding come from there, so the two sets
### of figures cannot drift apart.
###
### Jack H. Buckner, August 2026, Oregon State University 
### Generated with Claude Code
### reviewed JHB 9/1/26
###
#############################################################
#############################################################

library(dplyr)
library(ggplot2)


# ---------------------------------------------------------------------------
# Small statistics helpers
#
# Every summary is computed over a subset that can legitimately be empty (a
# variable no filter touched has no removed values; a site can have no valid
# observations at all). Base R is inconsistent there -- mean() gives NaN, min()
# gives -Inf with a warning -- so everything goes through these, which return NA.
# ---------------------------------------------------------------------------

.st <- function(x, f) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  v <- suppressWarnings(f(x))
  if (length(v) != 1 || !is.finite(v)) NA_real_ else as.numeric(v)
}

.qt <- function(x, p) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, p, na.rm = TRUE))
}


#' Kept and removed observations in one frame
#'
#' Reduces the cleaned data and the stage-tagged removals to the four columns
#' every diagnostic needs. Doing this once, rather than per variable, keeps a
#' large input from being copied repeatedly.
#'
#' @param kept The cleaned long data frame.
#' @param removed The stacked removals, from `run_pipeline()$removed`.
#' @param group_id Column naming the group, or NULL for pooled only.
#' @return A data frame with `variable`, `group`, `value` and `retained`.
#' @keywords internal
kept_removed_frame <- function(kept, removed, group_id = NULL) {

  stopifnot(is.data.frame(kept))
  grp <- function(d) {
    if (is.null(group_id)) return(rep("(all)", nrow(d)))
    if (!group_id %in% names(d))
      stop("`group_id` = '", group_id, "' is not a column of the data.")
    as.character(d[[group_id]])
  }

  out <- data.frame(variable = as.character(kept$variable),
                    group    = grp(kept),
                    value    = as.numeric(kept$value),
                    retained = TRUE,
                    stringsAsFactors = FALSE)

  if (!is.null(removed) && nrow(removed)) {
    out <- rbind(out, data.frame(
      variable = as.character(removed$variable),
      group    = grp(removed),
      value    = as.numeric(removed$value),
      retained = FALSE,
      stringsAsFactors = FALSE))
  }
  out
}


# Internal: the full statistic block for one (variable, group) cell.
#
# Two denominators matter here and they differ by an order of magnitude. The
# input is a padded grid: every variable has a row at every key, most of them
# empty. `n_input` counts those rows; `n_present` counts the rows that actually
# carry a value. A removal rate against `n_input` is diluted by the padding and
# understates the screening, so `pct_removed` uses `n_present` -- the same
# denominator the distribution plots show -- and the row-based rate is kept
# alongside it as `pct_of_rows`.
.value_stats <- function(value, retained) {

  kv <- value[retained]
  rv <- value[!retained]
  n  <- length(value)

  n_kept_present <- sum(!is.na(kv))
  n_rm_present   <- sum(!is.na(rv))
  n_present      <- n_kept_present + n_rm_present

  data.frame(
    n_input     = n,
    n_present   = n_present,
    n_kept      = length(kv),
    n_removed   = length(rv),
    pct_removed = round(100 * n_rm_present / max(n_present, 1), 3),
    pct_of_rows = round(100 * length(rv) / max(n, 1), 3),
    n_missing   = sum(is.na(kv)),
    # distribution of what survived
    mean   = .st(kv, mean),
    sd     = .st(kv, stats::sd),
    median = .st(kv, stats::median),
    mad    = .st(kv, stats::mad),
    min    = .st(kv, min),
    p01    = .qt(kv, 0.01),
    p05    = .qt(kv, 0.05),
    p25    = .qt(kv, 0.25),
    p75    = .qt(kv, 0.75),
    p95    = .qt(kv, 0.95),
    p99    = .qt(kv, 0.99),
    max    = .st(kv, max),
    # distribution of what was removed
    rm_n      = sum(!is.na(rv)),
    rm_mean   = .st(rv, mean),
    rm_sd     = .st(rv, stats::sd),
    rm_median = .st(rv, stats::median),
    rm_min    = .st(rv, min),
    rm_max    = .st(rv, max),
    stringsAsFactors = FALSE)
}


# Internal: apply .value_stats over the cells defined by `keys`.
# The key values are read back from the data rather than parsed out of the
# interaction label, so a group name containing the separator cannot corrupt
# them.
.summarise_cells <- function(kr, keys) {

  f     <- interaction(kr[keys], drop = TRUE)
  parts <- split(seq_len(nrow(kr)), f)
  parts <- parts[lengths(parts) > 0]

  rows <- lapply(parts, function(idx) {
    cbind(kr[idx[1], keys, drop = FALSE],
          .value_stats(kr$value[idx], kr$retained[idx]))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Distribution of every variable, in the full sample and across groups
#'
#' One row per variable for the pooled sample, plus one row per variable and
#' group when `group_id` is set. The pooled rows carry the group label `(all)`,
#' the same convention `match_up_variance()` uses.
#'
#' Statistics of the surviving values (`mean` ... `max`) describe what the model
#' will actually see. The `rm_*` block describes what was taken out, which is
#' what tells you whether a filter caught a distinct population or simply
#' trimmed the tail of a healthy one.
#'
#' @param kept The cleaned long data frame.
#' @param removed The stacked removals, from `run_pipeline()$removed`.
#' @param group_id Column naming the group, or NULL for pooled statistics only.
#' @return A data frame with `variable`, `group` and the statistic columns.
#' @export
variable_summary <- function(kept, removed, group_id = NULL) {

  kr <- kept_removed_frame(kept, removed, group_id)
  if (!nrow(kr)) return(data.frame())

  pooled <- .summarise_cells(kr, "variable")
  pooled <- cbind(pooled["variable"], group = "(all)",
                  pooled[setdiff(names(pooled), "variable")],
                  stringsAsFactors = FALSE)

  if (is.null(group_id)) {
    out <- pooled
  } else {
    by_grp <- .summarise_cells(kr, c("variable", "group"))
    out    <- rbind(pooled, by_grp[names(pooled)])
  }

  # Pooled row first within each variable, then groups in order.
  out <- out[order(out$variable, out$group != "(all)", out$group), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' What each stage removed, per variable and group
#'
#' The removal cascade: which filter took out what, expressed both as a count
#' and as a share of every observation of that variable that entered the
#' pipeline. Reading it top to bottom shows whether one stage is doing all the
#' work, or whether a later stage is removing what an earlier one should have.
#'
#' `n_removed` counts rows, which is what the pipeline acts on. `n_present`
#' counts only the removed rows that actually carried a value, and
#' `pct_of_present` expresses those against the observations of that variable
#' that exist -- not against the empty cells of the padded grid, which would
#' dilute every rate by the padding.
#'
#' @param removed The stacked removals, from `run_pipeline()$removed`.
#' @param kept The cleaned long data frame, used for the denominators.
#' @param group_id Column naming the group, or NULL for pooled counts only.
#' @return A data frame with `stage`, `variable`, `group`, `n_removed`,
#'   `n_present` and `pct_of_present`.
#' @export
removal_summary <- function(removed, kept, group_id = NULL) {

  cols <- c("stage", "variable", "group", "n_removed", "n_present",
            "pct_of_present")
  if (is.null(removed) || !nrow(removed))
    return(setNames(data.frame(character(), character(), character(),
                               integer(), integer(), numeric(),
                               stringsAsFactors = FALSE), cols))

  kr <- kept_removed_frame(kept, removed, group_id)
  kr <- kr[!is.na(kr$value), , drop = FALSE]

  # Denominator: every observation of the variable that entered the pipeline and
  # carried a value, i.e. what survived plus everything any stage removed.
  tot_pooled <- table(kr$variable)
  tot_group  <- table(kr$variable, kr$group)

  grp <- if (is.null(group_id)) rep("(all)", nrow(removed))
         else as.character(removed[[group_id]])

  base <- data.frame(stage    = as.character(removed$stage),
                     variable = as.character(removed$variable),
                     group    = grp,
                     present  = !is.na(removed$value),
                     stringsAsFactors = FALSE)

  # `table` over the same cells twice: all removed rows, then only those that
  # carried a value.
  tally <- function(keys) {
    a <- as.data.frame(table(base[keys]), stringsAsFactors = FALSE)
    names(a)[names(a) == "Freq"] <- "n_removed"
    b <- as.data.frame(table(base[base$present, keys, drop = FALSE]),
                       stringsAsFactors = FALSE)
    names(b)[names(b) == "Freq"] <- "n_present"
    a <- merge(a, b, by = keys, all.x = TRUE)
    a$n_present[is.na(a$n_present)] <- 0L
    a[a$n_removed > 0, , drop = FALSE]
  }

  pooled <- tally(c("stage", "variable"))
  pooled$group <- "(all)"
  pooled$pct_of_present <-
    round(100 * pooled$n_present / as.numeric(tot_pooled[pooled$variable]), 3)
  out <- pooled[cols]

  if (!is.null(group_id)) {
    bg <- tally(c("stage", "variable", "group"))
    if (nrow(bg)) {
      denom <- tot_group[cbind(bg$variable, bg$group)]
      bg$pct_of_present <- round(100 * bg$n_present / as.numeric(denom), 3)
      out <- rbind(out, bg[cols])
    }
  }

  out <- out[order(out$variable, out$stage,
                   out$group != "(all)", out$group), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Full-sample distribution of one variable, removed values in red
#'
#' @param kr Output of `kept_removed_frame()`, already subset to one variable.
#' @param var The variable name, for the title.
#' @param bins Histogram bins. Defaults to 100, matching the per-process plots.
#' @return A ggplot object, or NULL if there is nothing to plot.
#' @export
plot_variable_distribution <- function(kr, var, bins = 100) {

  d <- kr[!is.na(kr$value), , drop = FALSE]
  if (!nrow(d)) return(NULL)

  n_rm <- sum(!d$retained)
  sub  <- sprintf("%s of %s observations removed (%.2f%%)",
                  format(n_rm, big.mark = ","),
                  format(nrow(d), big.mark = ","),
                  100 * n_rm / nrow(d))

  ggplot2::ggplot(d, ggplot2::aes(x = value, fill = retained)) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::scale_fill_manual(values = RETAINED_PAL, name = "Retained") +
    ggplot2::labs(x = var, y = "Count",
                  title = paste0("Distribution: ", var),
                  subtitle = sub) +
    ggplot2::theme_classic()
}


#' Distribution of one variable across groups, removed values in red
#'
#' The form adapts to how many groups there are. Up to `max_facets` the panel is
#' a grid of histograms, which reads the same way as the full-sample plot. Above
#' that a grid becomes unreadable, so it switches to boxplots of the retained
#' values ordered by median, with a companion bar chart of the percent removed
#' per group in the same order. Together those answer the two questions a facet
#' grid would: where does this group sit, and how much of it survived.
#'
#' @param kr Output of `kept_removed_frame()`, already subset to one variable.
#' @param var The variable name, for the title.
#' @param bins Histogram bins.
#' @param max_facets Group count above which the panel becomes boxplots.
#' @param max_points Maximum points drawn per group panel.
#' @return A named list of ggplot objects, empty if there is nothing to plot.
#' @export
plot_variable_distribution_by_group <- function(kr, var, bins = 100,
                                                max_facets = 12,
                                                max_points = 20000) {

  d <- kr[!is.na(kr$value), , drop = FALSE]
  if (!nrow(d)) return(list())

  groups   <- sort(unique(d$group))
  n_groups <- length(groups)

  # ---- few groups: a real facet grid of histograms -------------------------
  if (n_groups <= max_facets) {
    plt <- ggplot2::ggplot(d, ggplot2::aes(x = value, fill = retained)) +
      ggplot2::geom_histogram(bins = bins) +
      ggplot2::facet_wrap(~ group, scales = "free_y") +
      ggplot2::scale_fill_manual(values = RETAINED_PAL, name = "Retained") +
      ggplot2::labs(x = var, y = "Count",
                    title = paste0("Distribution by group: ", var),
                    subtitle = paste0(n_groups, " group(s)")) +
      ggplot2::theme_classic()
    return(list(by_group = plt))
  }

  # ---- many groups: boxplots ordered by median, plus a removal bar ---------
  pct <- do.call(rbind, lapply(groups, function(g) {
    s <- d$group == g
    data.frame(group = g, n = sum(s), n_removed = sum(!d$retained[s]),
               pct_removed = 100 * sum(!d$retained[s]) / sum(s),
               med = .st(d$value[s & d$retained], stats::median),
               stringsAsFactors = FALSE)
  }))
  # Groups with no retained values have no median to sort on; park them last.
  ord <- pct$group[order(pct$med, na.last = TRUE)]

  # Cap points *per group*: one busy site can dominate the total, so the cheap
  # whole-frame check only decides whether the per-group split is worth doing.
  box_dat <- d[d$retained, , drop = FALSE]
  if (nrow(box_dat) > max_points) {
    n_before <- nrow(box_dat)
    keep <- unlist(lapply(split(seq_len(nrow(box_dat)), box_dat$group),
                          function(i) if (length(i) > max_points)
                            sample(i, max_points) else i), use.names = FALSE)
    box_dat <- box_dat[sort(keep), , drop = FALSE]
    if (nrow(box_dat) < n_before)
      message("  thinned '", var, "' to ", format(max_points, big.mark = ","),
              " retained value(s) per group for the boxplot panel (",
              format(n_before - nrow(box_dat), big.mark = ","), " dropped).")
  }
  box_dat$group <- factor(box_dat$group, levels = ord)
  pct$group     <- factor(pct$group,     levels = ord)

  p_box <- ggplot2::ggplot(box_dat, ggplot2::aes(x = group, y = value)) +
    ggplot2::geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3,
                          fill = RETAINED_PAL[["TRUE"]], colour = "grey20",
                          linewidth = 0.3) +
    ggplot2::labs(x = NULL, y = var,
                  title = paste0("Retained values by group: ", var),
                  subtitle = paste0(n_groups,
                                    " groups, ordered by median")) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1, size = 6))

  p_pct <- ggplot2::ggplot(pct, ggplot2::aes(x = group, y = pct_removed)) +
    ggplot2::geom_col(fill = RETAINED_PAL[["FALSE"]]) +
    ggplot2::geom_hline(yintercept = 100 * sum(!d$retained) / nrow(d),
                        linetype = 2, colour = "grey40") +
    ggplot2::labs(x = NULL, y = "Removed (%)",
                  title = paste0("Removed by group: ", var),
                  subtitle = "Same group order; dashed line is the overall rate") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1, size = 6))

  list(by_group = p_box, by_group_removed = p_pct)
}


#' The removal cascade: how much each stage took out of each variable
#'
#' Bars count removed observations that carried a value (`n_present`). A stage
#' can also drop empty rows of the padded grid -- the wind filter rejects a whole
#' key, whether or not the target was populated there -- and counting those would
#' make a stage look like it removed far more data than it did.
#'
#' @param tbl Pooled rows of `removal_summary()`.
#' @return A ggplot object, or NULL if nothing was removed.
#' @export
plot_removal_cascade <- function(tbl) {

  d <- tbl[tbl$group == "(all)" & tbl$n_present > 0, , drop = FALSE]
  if (!nrow(d)) return(NULL)

  ggplot2::ggplot(d, ggplot2::aes(x = variable, y = n_present, fill = stage)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Observations removed", fill = "Stage",
                  title = "Removals by stage",
                  subtitle = paste("Removed values only; empty cells of the",
                                   "padded grid are not counted")) +
    ggplot2::theme_classic()
}


#' Observation coverage of the MARSS matrix
#'
#' Which cells of the observation matrix carry a value and which are gaps, with
#' the series ordered as they appear in the matrix so the variable blocks read as
#' bands. This is the figure that shows at a glance which series are too sparse
#' to support the model, before a fit spends an hour discovering it.
#'
#' @param marss The list from `build_marss_inputs()`.
#' @param max_cols Time steps to draw. Wider matrices are thinned by taking a
#'   regular subsample of columns; the aggregate pattern survives, and drawing
#'   7,000+ columns as tiles would not be legible anyway.
#' @return A ggplot object, or NULL if there is nothing to draw.
#' @export
plot_marss_coverage <- function(marss, max_cols = 1500) {

  m <- marss$ts_matrix
  if (is.null(m) || !length(m)) return(NULL)

  cols <- seq_len(ncol(m))
  thinned <- FALSE
  if (ncol(m) > max_cols) {
    cols <- unique(round(seq(1, ncol(m), length.out = max_cols)))
    thinned <- TRUE
  }

  d <- expand.grid(row = seq_len(nrow(m)), ci = seq_along(cols),
                   KEEP.OUT.ATTRS = FALSE)
  d$observed <- !is.na(m[cbind(d$row, cols[d$ci])])
  d$date     <- marss$col_dates[cols][d$ci]
  d          <- d[d$observed, , drop = FALSE]

  n_obs <- sum(!is.na(m))
  sub   <- sprintf("%s of %s cells observed (%.2f%%)%s",
                   format(n_obs, big.mark = ","),
                   format(length(m), big.mark = ","),
                   100 * n_obs / length(m),
                   if (thinned) "; time axis subsampled for display" else "")

  # Variable block boundaries, so the bands are readable as instruments.
  vb <- marss$row_var_keys
  brk <- which(c(FALSE, vb[-1] != vb[-length(vb)]))
  lab <- tapply(seq_along(vb), factor(vb, levels = unique(vb)), mean)

  ggplot2::ggplot(d, ggplot2::aes(x = date, y = row)) +
    ggplot2::geom_raster(fill = RETAINED_PAL[["TRUE"]]) +
    ggplot2::geom_hline(yintercept = brk - 0.5, colour = "firebrick",
                        linewidth = 0.3) +
    ggplot2::scale_y_reverse(breaks = round(lab), labels = names(lab)) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "MARSS observation coverage",
                  subtitle = sub) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
}


#' Build the whole-dataset diagnostic tables and figures
#'
#' @param result The list returned by `run_pipeline()`.
#' @param config The resolved configuration, for `group_id` and
#'   `diagnostics` settings.
#' @return A list with `tables` and `plots`.
#' @export
build_diagnostics <- function(result, config) {

  gid <- config$group_id
  par <- config$diagnostics
  bins       <- par$bins       %||% 100
  max_facets <- par$max_facets %||% 12
  max_points <- par$max_points %||% 20000

  kr <- kept_removed_frame(result$data, result$removed, gid)

  # `result$removal_summary` is written separately as stage_counts.csv, so it is
  # deliberately not repeated here.
  tables <- list(
    variable_summary = variable_summary(result$data, result$removed, NULL),
    removal_by_stage = removal_summary(result$removed, result$data, gid))

  if (!is.null(gid)) {
    tables$variable_summary_by_group <-
      variable_summary(result$data, result$removed, gid)
  }

  plots <- list()
  for (v in sort(unique(kr$variable))) {
    sub <- kr[kr$variable == v, , drop = FALSE]

    p <- plot_variable_distribution(sub, v, bins = bins)
    if (!is.null(p)) plots[[paste0("dist_", v)]] <- p

    if (!is.null(gid)) {
      ps <- plot_variable_distribution_by_group(sub, v, bins = bins,
                                                max_facets = max_facets,
                                                max_points = max_points)
      for (nm in names(ps)) plots[[paste0(nm, "_", v)]] <- ps[[nm]]
    }
  }

  casc <- plot_removal_cascade(tables$removal_by_stage)
  if (!is.null(casc)) plots$removal_cascade <- casc

  list(tables = tables, plots = plots)
}


# Internal: make a string safe to use as a file name.
.slug <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)


#' Write every pipeline artifact to disk
#'
#' Lays out the run as tables, individual figures, and one combined report:
#'
#' ```
#' <outdir>/clean_data.rds
#'          removed_all.csv
#'          tables/*.csv
#'          diagnostics/*.png
#'          diagnostics_report.pdf
#' ```
#'
#' The per-process plots produced by the individual stages are written alongside
#' the whole-dataset ones, so the report holds the complete record of the run.
#'
#' @param result The list returned by `run_pipeline()`.
#' @param outdir Directory to write into. Created if absent.
#' @param config The resolved configuration.
#' @param verbose If TRUE (default), report what is written.
#' @return `outdir`, invisibly.
#' @export
write_pipeline_outputs <- function(result, outdir, config, verbose = TRUE) {

  say <- function(...) if (verbose) message(...)
  par <- config$diagnostics
  w   <- par$width  %||% 9
  h   <- par$height %||% 6
  dpi <- par$dpi    %||% 150

  tdir <- file.path(outdir, "tables")
  ddir <- file.path(outdir, "diagnostics")
  for (p in c(outdir, tdir)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

  # ---- data ----------------------------------------------------------------
  saveRDS(result$data, file.path(outdir, "clean_data.rds"))
  say("Wrote clean_data.rds (", format(nrow(result$data), big.mark = ","),
      " rows).")

  utils::write.csv(result$removed, file.path(outdir, "removed_all.csv"),
                   row.names = FALSE)
  say("Wrote removed_all.csv (", format(nrow(result$removed), big.mark = ","),
      " rows).")

  # ---- tables --------------------------------------------------------------
  wtab <- function(x, name) {
    if (is.null(x) || !NROW(x)) return(invisible(NULL))
    utils::write.csv(as.data.frame(x), file.path(tdir, name), row.names = FALSE)
    say("  tables/", name)
  }

  wtab(result$removal_summary, "stage_counts.csv")
  for (nm in names(result$diagnostics$tables))
    wtab(result$diagnostics$tables[[nm]], paste0(.slug(nm), ".csv"))

  # Per-stage summaries and fits, one file per stage occurrence.
  for (i in seq_along(result$stages)) {
    s   <- result$stages[[i]]
    tag <- sprintf("stage_%02d_%s_%s", i, .slug(s$stage), .slug(s$variable))
    wtab(s$summary, paste0(tag, "_summary.csv"))
    wtab(s$fit,     paste0(tag, "_fit.csv"))
  }

  if (!is.null(result$match_up_variance)) {
    wtab(result$match_up_variance$table, "match_up_variance.csv")
  }

  # ---- MARSS inputs --------------------------------------------------------
  if (!is.null(result$marss)) {
    saveRDS(result$marss, file.path(outdir, "marss_inputs.rds"))
    say("Wrote marss_inputs.rds (", nrow(result$marss$ts_matrix), " series x ",
        ncol(result$marss$ts_matrix), " time steps, ",
        nrow(result$marss$harmonics), " harmonic rows).")
    wtab(result$marss$row_metadata, "marss_row_metadata.csv")
    wtab(result$marss$dropped,      "marss_dropped_series.csv")
  }

  # ---- figures -------------------------------------------------------------
  plots <- list()
  for (i in seq_along(result$stages)) {
    s <- result$stages[[i]]
    for (nm in names(s$plots)) {
      key <- sprintf("stage_%02d_%s_%s_%s", i, .slug(s$stage),
                     .slug(s$variable), nm)
      plots[[key]] <- s$plots[[nm]]
    }
  }
  plots <- c(plots, result$diagnostics$plots)
  if (!is.null(result$match_up_variance)) {
    plots$match_up_variance_scatter <- result$match_up_variance$plot
    plots$match_up_variance_diff    <- result$match_up_variance$plot_resid
  }
  if (!is.null(result$marss)) {
    plots$marss_coverage <- tryCatch(plot_marss_coverage(result$marss),
      error = function(e) {
        warning("Could not build the MARSS coverage figure: ",
                conditionMessage(e))
        NULL
      })
  }
  plots <- plots[!vapply(plots, is.null, logical(1))]

  if (!length(plots)) {
    say("No figures to write.")
    return(invisible(outdir))
  }

  dir.create(ddir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(plots)) {
    f <- file.path(ddir, paste0(.slug(nm), ".png"))
    ok <- tryCatch({
      suppressMessages(ggplot2::ggsave(f, plots[[nm]], width = w, height = h,
                                       dpi = dpi))
      TRUE
    }, error = function(e) {
      warning("Could not write ", basename(f), ": ", conditionMessage(e))
      FALSE
    })
    if (ok) say("  diagnostics/", basename(f))
  }

  pdf_file <- file.path(outdir, "diagnostics_report.pdf")
  grDevices::pdf(pdf_file, width = w, height = h, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (nm in names(plots)) {
    tryCatch(suppressMessages(print(plots[[nm]])),
             error = function(e)
               warning("Could not render ", nm, " into the report: ",
                       conditionMessage(e)))
  }
  say("Wrote diagnostics_report.pdf (", length(plots), " pages).")

  invisible(outdir)
}
