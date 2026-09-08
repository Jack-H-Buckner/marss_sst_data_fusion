
#!/usr/bin/env Rscript
# Fill missing modis_sst_aqua for sites that have none, from a nearby site that
# does.
#
# This is the MODIS counterpart of imputs_mur_to_nearest.R, and the difference
# between the two products drives the difference between the scripts. MUR is a
# gapless L4 field, so any site that has it has it on every date and the only
# question is which donor is closest. MODIS is a swath product thinned by cloud,
# so donors differ in *how much* they have as well as how far away they are: the
# nearest site may carry 53 usable days where one a little further out carries
# 114. Taking the nearest donor blindly can therefore hand a site a much emptier
# series than it could have had.
#
# So the donor is chosen on two controls:
#
#   MAX_DONOR_KM        no site further than this is considered at all.
#   MIN_DONOR_FRACTION  if the nearest candidate's observation count falls
#                       below this fraction of the best covered candidate
#                       within MAX_DONOR_KM, the better covered one is used
#                       instead. 1 always takes the best covered candidate,
#                       0 always takes the nearest.
#
# The three nearest donors are printed for every site with their distance and
# their observation count, so the choice can be checked against what was on
# offer rather than taken on trust.
#
# Every filled value is labelled with its donor, the distance it travelled and
# that donor's coverage, so imputed rows can always be excluded again.

suppressMessages(library(data.table))

# ---------------- CONFIG ----------------
CSV_PATH           <- "figures/EDA/coastal_sst_nc_mur_imputed_clouds_R.csv"
VARIABLE           <- "modis_sst_aqua"
STAT               <- "nanmean"
MAX_DONOR_KM       <- 25      # NA for no limit
MIN_DONOR_FRACTION <- 0.5     # 0 to always take the nearest candidate
N_REPORT           <- 3       # nearest donors listed per site
OUT_CSV            <- "figures/EDA/coastal_sst_nc_modis_and_modis_imputed_R.csv"
# ----------------------------------------
# Point CSV_PATH at the output of imputs_mur_to_nearest.R to stack the two
# fills; the provenance columns of each are kept separately.

haversine_km <- function(lat, lon, lats, lons) {
  R <- 6371
  p <- function(x) x * pi / 180
  dlat <- p(lats - lat); dlon <- p(lons - lon)
  a <- sin(dlat / 2)^2 + cos(p(lat)) * cos(p(lats)) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

d <- fread(CSV_PATH)
d[, time := as.IDate(time)]

sites <- d[, .(lat = first(lat), lon = first(lon)), by = point_id]
setkey(sites, point_id)

# A site "has" the variable if at least one value is non-NA. Unlike MUR, the
# count matters as well as the fact, so it is carried through.
cover <- d[variable == VARIABLE & stat == STAT,
           .(n_valid = sum(!is.na(value))), by = point_id]
have <- sort(cover[n_valid > 0, point_id])
lack <- sort(setdiff(sites$point_id, have))
n_obs <- setNames(cover$n_valid, cover$point_id)

cat(sprintf("%d sites have %s, %d do not\n", length(have), VARIABLE, length(lack)))
cat(sprintf("Donor coverage: %d to %d observation(s); limit %s km, switch below %.2f x the best\n\n",
            min(n_obs[have]), max(n_obs[have]),
            ifelse(is.na(MAX_DONOR_KM), "no", format(MAX_DONOR_KM)),
            MIN_DONOR_FRACTION))

# Donor series, one row per (donor, date).
donors <- d[variable == VARIABLE & stat == STAT & point_id %in% have,
            .(point_id, time, donor_value = value)]
if (anyDuplicated(donors, by = c("point_id", "time")))
  stop("Duplicate (point_id, time) rows among the donors: the join below would ",
       "expand rather than align.", call. = FALSE)
setkey(donors, point_id, time)

don_sites <- sites[point_id %in% have]

# Provenance columns, so imputed values are never mistaken for observed ones.
d[, `:=`(modis_imputed = FALSE, modis_donor = NA_character_,
         modis_donor_km = NA_real_, modis_donor_n = NA_integer_)]

filled <- 0L; n_switched <- 0L; n_skipped <- 0L
for (site in lack) {
  s    <- sites[point_id == site]
  dist <- haversine_km(s$lat, s$lon, don_sites$lat, don_sites$lon)
  ord  <- order(dist)                      # nearest first
  cand <- data.table(point_id = don_sites$point_id[ord],
                     km       = dist[ord],
                     n        = as.integer(n_obs[don_sites$point_id[ord]]))

  cat(sprintf("%s\n", site))
  for (k in seq_len(min(N_REPORT, nrow(cand))))
    cat(sprintf("   %d. %-34s %7.2f km  %4d obs%s\n", k, cand$point_id[k],
                cand$km[k], cand$n[k],
                if (!is.na(MAX_DONOR_KM) && cand$km[k] > MAX_DONOR_KM)
                  "   (beyond the limit)" else ""))

  within <- if (is.na(MAX_DONOR_KM)) cand else cand[km <= MAX_DONOR_KM]
  if (!nrow(within)) {
    cat(sprintf("   -> SKIPPED, nearest donor is %.1f km (limit %g km)\n\n",
                cand$km[1], MAX_DONOR_KM))
    n_skipped <- n_skipped + 1L
    next
  }

  # The nearest candidate unless it is thin: `best` is the best covered site
  # inside the limit, and ties there fall to the nearer one because `within` is
  # already in distance order.
  near <- within[1]
  best <- within[which.max(n)]
  switched <- near$n < MIN_DONOR_FRACTION * best$n
  pick <- if (switched) best else near

  if (switched) {
    cat(sprintf("   -> %s, %.2f km, %d obs (nearest has %d, below %.2f x %d)\n",
                pick$point_id, pick$km, pick$n, near$n, MIN_DONOR_FRACTION,
                best$n))
    n_switched <- n_switched + 1L
  } else {
    cat(sprintf("   -> %s, %.2f km, %d obs (nearest)\n",
                pick$point_id, pick$km, pick$n))
  }

  idx <- d[, which(variable == VARIABLE & stat == STAT & point_id == site)]
  if (!length(idx)) {
    cat("   -> SKIPPED, this site has no ", VARIABLE, " rows to fill\n\n",
        sep = "")
    n_skipped <- n_skipped + 1L
    next
  }
  # Align the donor onto this site's own dates. NAs come across as NAs: the
  # recipient inherits the donor's cloud gaps, it does not inherit coverage it
  # never had.
  vals <- donors[.(pick$point_id, d$time[idx]), on = .(point_id, time),
                 donor_value]
  stopifnot(length(vals) == length(idx))

  set(d, i = idx, j = "value",          value = vals)
  set(d, i = idx, j = "modis_imputed",  value = TRUE)
  set(d, i = idx, j = "modis_donor",    value = pick$point_id)
  set(d, i = idx, j = "modis_donor_km", value = round(pick$km, 2))
  set(d, i = idx, j = "modis_donor_n",  value = pick$n)

  got <- sum(!is.na(vals))
  filled <- filled + got
  cat(sprintf("      filled %d of %d date(s)\n\n", got, length(idx)))
}

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
fwrite(d, OUT_CSV)
cat(sprintf("Filled %s site-days across %d site(s) (%d switched off the nearest donor, %d skipped)\nSaved %s\n",
            format(filled, big.mark = ","), length(lack) - n_skipped,
            n_switched, n_skipped, OUT_CSV))
