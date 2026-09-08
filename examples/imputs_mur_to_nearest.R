
#!/usr/bin/env Rscript
# Fill missing mur_sst for sites that have none, using the nearest site that does.
#
# MUR is a 1 km gridded L4 product, so a nearby site's value is a defensible
# stand-in over short distances. Every filled value is labelled with its donor
# and the distance it travelled, so imputed rows can always be excluded again.

suppressMessages(library(data.table))

# ---------------- CONFIG ----------------
CSV_PATH     <- "data/raw/coastal_sst_data_nc_clouds.csv"
VARIABLE     <- "mur_sst"
STAT         <- "nanmean"
MAX_DONOR_KM <- 25            # NA for no limit
OUT_CSV      <- "figures/EDA/coastal_sst_nc_mur_imputed_clouds_R.csv"
# ----------------------------------------

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

# A site "has" the variable if at least one value is non-NA.
cover <- d[variable == VARIABLE & stat == STAT,
           .(n_valid = sum(!is.na(value))), by = point_id]
have <- sort(cover[n_valid > 0, point_id])
lack <- sort(setdiff(sites$point_id, have))
cat(sprintf("%d sites have %s, %d do not\n", length(have), VARIABLE, length(lack)))

# Donor series, one row per (donor, date).
donors <- d[variable == VARIABLE & stat == STAT & point_id %in% have,
            .(point_id, time, donor_value = value)]
setkey(donors, point_id, time)

# Provenance columns, so imputed values are never mistaken for observed ones.
d[, `:=`(mur_imputed = FALSE, mur_donor = NA_character_, mur_donor_km = NA_real_)]

filled <- 0L
for (site in lack) {
  s <- sites[point_id == site]
  dist <- haversine_km(s$lat, s$lon, sites[point_id %in% have]$lat,
                       sites[point_id %in% have]$lon)
  i    <- which.min(dist)
  best <- sites[point_id %in% have][i]$point_id
  km   <- dist[i]
  
  if (!is.na(MAX_DONOR_KM) && km > MAX_DONOR_KM) {
    cat(sprintf("  %-34s SKIPPED, nearest donor %s is %.1f km (limit %g km)\n",
                site, best, km, MAX_DONOR_KM))
    next
  }
  
  idx  <- d[, which(variable == VARIABLE & stat == STAT & point_id == site)]
  # Align the donor onto this site's own dates.
  vals <- donors[.(best, d$time[idx]), on = .(point_id, time), donor_value]
  
  set(d, i = idx, j = "value",        value = vals)
  set(d, i = idx, j = "mur_imputed",  value = TRUE)
  set(d, i = idx, j = "mur_donor",    value = best)
  set(d, i = idx, j = "mur_donor_km", value = round(km, 2))
  
  filled <- filled + sum(!is.na(vals))
  cat(sprintf("  %-34s <- %-34s %6.2f km  %d days\n",
              site, best, km, sum(!is.na(vals))))
}

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
fwrite(d, OUT_CSV)
cat(sprintf("\nFilled %s site-days\nSaved %s\n",
            format(filled, big.mark = ","), OUT_CSV))