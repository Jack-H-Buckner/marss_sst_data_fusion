

source("R/00_preprocessing.R")
source("R/00_diagnostics.R")

dat <- read.csv("data/raw/tasi_salmon_20_yrs.csv")
res <- run_pipeline(dat, "parameters/preprocessing.R")

res$data             # cleaned long format data
res$removed          # every removed row, tagged with the stage that removed it
res$removal_summary  # rows in / out / removed per stage
res$diagnostics$tables$variable_summary

write_pipeline_outputs(res, "outputs/macquarie_harbor", res$config)   # optional

cfg <- read_model_config("parameters/macquarie_harbor_sites.R")
dat <- readRDS(cfg$data$input)
validate_model_config(cfg, dat)

md    <- build_model_data(dat, cfg)     # subset, scale, index rows by instrument
model <- build_marss_model(md, cfg)     # Z, A, R, B, U, Q, D, x0, V0, d