#############################################################
#############################################################
###
### Formal validation of the North Carolina runs against the
### in situ data held out of the fit.
###
### All the work is in R/04_validate_holdout.R. This file only
### says which runs to score and how to split their sites.
###
### Two runs are scored:
###
###   north_carolina          4 held out sites, all of which
###                           have MODIS
###   north_carolina_imputed  7 held out sites, 4 with MODIS
###                           and 3 without
###
### The second is the one that can answer whether fusing helps
### more where MODIS is available: `group_var = "has_modis_sst"`
### scores the two sets of sites separately and puts an
### interval on the gap between them. In the first run every
### held out site has MODIS, so that split has nothing to
### contrast and is skipped with a message rather than an
### error -- which is what makes the same call work on both.
###
### MIN_OBS is how much MODIS a site needs before it counts as
### having it. At 1, the default, one scene is enough and the
### split is presence against absence. Raising it splits on
### coverage instead: the Tasmanian runs are the case that
### needs this, since every held out site there has some MODIS
### -- 11 to 607 scenes -- so at MIN_OBS = 1 the groups collapse
### to one and the grouped figure is skipped. A threshold that
### falls between the sites (100 splits east_coast 3 against 2)
### gives back the comparison the data can actually support:
### densely covered sites against thinly covered ones, not with
### against without. The skip message prints the per site counts
### so the threshold can be chosen from them.
###
### MIN_OBS applies to every run in RUNS, so score the North
### Carolina runs at 1 -- where the split is genuinely presence
### and absence -- and the Tasmanian ones separately.
###
### Every interval is reported at two levels, ALPHA and
### ALPHA2: the columns without a suffix are the 95% level and
### those ending in 2 the 90%, and the figures draw the 90% as
### a heavy bar inside the 95% whisker.
###
### To add more sites, hold them out in the model config
### (`data$sites_validation`), refit, and point RUNS at the new
### directory. Nothing here needs to change: the held out
### sites and each site's instrument coverage are read from
### design.rds and observations.rds.
###
### north_carolina_validation.R is the exploratory version of
### this and is left alone.
###
#############################################################
#############################################################

library(ggplot2)
source("R/04_validate_holdout.R")


RUNS <- c(
  # north_carolina =
  #   "examples/models/north_carolina/north_carolina_2026-09-04",
  # north_carolina_imputed =
  #   "examples/models/north_carolina_imputed/north_carolina_2026-09-04",
  # north_carolina_clouds =
  #   "examples/models/north_carolina_clouds/north_carolina_2026-09-04",
  # north_carolina_diff_calib =
  #   "examples/models/north_carolina_diff_calib/north_carolina_diff_calib_2026-09-04",
  # north_carolina_modis_imputed =
  #   "examples/models/north_carolina_modis_imputed/north_carolina_2026-09-05",
  # tasi_modis_sites =
  #   "examples/models/tasi_modis_sites/tasi_modis_sites_2026-09-05",
  # tasi_modis =
  #   "examples/models/tasi_modis/tasi_modis_2026-09-05",
  tasi_east_coast =
    "examples/models/east_coast/east_coast_2026-09-05",
  tasi_storm_bay =
    "examples/models/storm_bay/storm_bay_models_2026-09-05")

GROUP_VAR <- "has_modis_sst"   # any has_<instrument> from site_coverage()
MIN_OBS   <- 100                 # observations that count as "has"; see above
B         <- 10000
SEED      <- 1
ALPHA     <- 0.05              # outer level: columns without a suffix
ALPHA2    <- 0.10              # inner level: columns ending in 2


results <- list()
for (nm in names(RUNS)) {
  cat("\n\n##########################################################\n")
  cat("### ", nm, "\n", sep = "")
  cat("##########################################################\n")

  res <- validate_holdout(RUNS[[nm]], group_var = GROUP_VAR,
                          min_obs = MIN_OBS,
                          B = B, seed = SEED, alpha = ALPHA, alpha2 = ALPHA2)
  write_holdout_validation(res)
  results[[nm]] <- res
}


############################################################
### Which sites carry MODIS, for the record.
############################################################
cat("\n\n=== MODIS coverage at the held out sites =================\n")
for (nm in names(results)) {
  cv <- results[[nm]]$coverage
  cat("\n", nm, ":\n", sep = "")
  print(cv[c("site", "n_modis_sst", "has_modis_sst", "n_mur_sst",
             "n_eco_sst_v002_clean", "n_lst_sst_clean")], row.names = FALSE)
}
