#############################################################
#############################################################
###
### Seasonal differences variant of the MARSS data fusion
### model.
###
### This is a complete configuration in its own right, not a
### patch on model.R: a config file is the archived record of
### a run, copied into the run directory as config_used.R, so
### it has to be readable on its own.
###
### It differs from model.R in exactly one value:
###
###   state_structure <- "factors_only"
###
### which drops the per site latent states and fits the
### shared dynamic factors alone. What is left of the site
### level signal then lives entirely in the seasonal terms,
### so this run answers a narrower question: how much of the
### between instrument disagreement is a difference in
### seasonality rather than a difference in day to day
### anomalies.
###
### Jack H. Buckner, Oregon State University, 08/30/2026
### Generated with Claude Code
#############################################################
#############################################################


############################################################
### The observations to fit. Identical to model.R, so the
### two structures are compared on the same data. See the
### note in model.R about the missing "SF0017".
############################################################
input      <- "outputs/tasi_salmon_20_yrs_2026-08-29/marss_inputs.rds"
sites      <- c("CB001", "CB008", "SF001", "SF025", "SF050")
start_date <- "2013-01-01"
end_date   <- NULL

data_params <- list(
  input = input, sites = sites,
  start_date = start_date, end_date = end_date
)


############################################################
### Scaling of the observations. Per site, from MUR at each
### site, which is why MUR's own intercepts are fixed at zero
### below. See model.R for the full rationale.
############################################################
enabled  <- TRUE
variable <- "mur_sst"

scaling_params <- list(enabled = enabled, variable = variable)


############################################################
### Structure of the state model.
###
### "factors_only": no per site states, only `m_factors`
### shared dynamic factors. Because the factors are now the
### only states, their innovation variance is estimated
### rather than fixed at 1.
###
### Note that this weakens the identification described in
### R/marss_matrix_functions.R: the triangular loading
### constraints alone do not pin down the factor scale once
### Q is free, so the loadings and the innovation variance
### trade off against each other. The states and the fitted
### values are still identified; the individual loadings are
### not, and should not be interpreted on their own.
############################################################
state_structure <- "factors_only"
m_factors       <- 2


############################################################
### How each instrument enters the observation model.
### Identical to model.R -- see that file for the meaning of
### each field and why MUR is given independent terms.
############################################################
instruments <- list(
  insitu_sst = list(
    tag = "insitu", intercept = "site", seasonality = "shared",
    error = "sigma_2_insitu", day_effect = FALSE, site_state = TRUE),
  lst_sst_clean = list(
    tag = "lst", intercept = "site+instrument", seasonality = "shared",
    error = "sigma_2_lst", day_effect = TRUE, site_state = TRUE),
  eco_sst_v002_clean = list(
    tag = "eco", intercept = "site+instrument", seasonality = "shared",
    error = "sigma_2_eco", day_effect = TRUE, site_state = TRUE),
  modis_sst = list(
    tag = "modis", intercept = "site+instrument", seasonality = "shared",
    error = "sigma_2_modis", day_effect = FALSE, site_state = TRUE),
  mur_sst = list(
    tag = "mur", intercept = "zero", seasonality = "independent",
    error = "sigma_2_mur", day_effect = FALSE, site_state = FALSE)
)

init_var_x0 <- 1.0

structure_params <- list(
  state_structure = state_structure, m_factors = m_factors,
  instruments = instruments, init_var_x0 = init_var_x0
)


############################################################
### Warm start for the seasonal parameters. The shared site
### terms are seeded from MUR, which is gap free at every
### site, and MUR's own terms from its own rows. See model.R.
############################################################
enabled          <- TRUE
shared_from      <- "mur_sst"
seed_independent <- TRUE
min_obs          <- 30

inits_params <- list(
  enabled = enabled, shared_from = shared_from,
  seed_independent = seed_independent, min_obs = min_obs
)


############################################################
### Fitting. Identical to model.R -- see that file for what
### `method` selects and why the controls are keyed by it.
### This structure has far fewer states, so the chunks run
### faster, but the same cap is kept so the two runs are
### directly comparable.
############################################################
method <- "kem"

controls <- list(

  kem = list(
    trace               = 1,
    maxit               = 20,
    minit               = 5,
    abstol              = 0.001,
    conv.test.slope.tol = 0.5,
    safe                = TRUE
  ),

  BFGS = list(
    trace  = 0,
    maxit  = 5000,
    REPORT = 100,
    reltol = 1e-8
  ),

  TMB = list(
    trace      = 0,
    maxit      = 5000,
    tmb.silent = TRUE
  )
)

warmup_controls <- list(minit = 1, maxit = 2, safe = TRUE)
chunks          <- 75


############################################################
### Initialisation stage. A fast TMB pass on a reduced model
### supplies starting values for the fit above. See model.R
### for why B has to be held fixed and the day effects
### dropped for that stage.
###
### Under this structure the only free entry of B is `rho`,
### so a scalar covers it; list(rho = 0.9) would be the
### equivalent named form.
############################################################
enabled       <- TRUE
init_method   <- "TMB"      # distinct name: `method` above is the final fit
B_values      <- 0.9
init_controls <- list(
  trace      = 0,
  maxit      = 5000,
  tmb.silent = TRUE
)

init_params <- list(
  enabled = enabled, method = init_method,
  B_values = B_values, controls = init_controls
)

fitting_params <- list(
  method = method, controls = controls,
  warmup_controls = warmup_controls, chunks = chunks,
  init = init_params
)


############################################################
### Where a run is written. `run_name` is left NULL so the
### directory is named after this file and the date, keeping
### it distinct from the model.R run.
############################################################
output_root <- "models"
run_name    <- NULL

output_params <- list(output_root = output_root, run_name = run_name)


model_config <- list(
  data      = data_params,
  scaling   = scaling_params,
  structure = structure_params,
  inits     = inits_params,
  fitting   = fitting_params,
  output    = output_params
)
