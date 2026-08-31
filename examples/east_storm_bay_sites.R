#############################################################
#############################################################
###
### This file defines the constants used to specify and fit
### the MARSS data fusion model. The observations it consumes
### have already been screened by the preprocessing step,
### which has its own parameters in preprocessing.R.
###
### The model is a state space model in which each site has a
### latent temperature, observed with instrument specific
### bias, seasonality and error. These parameters are
### documented here as a one-stop-shop for understanding the
### assumptions that go into the fit.
###
### Jack H. Buckner, Oregon State University, 08/30/2026
### 
#############################################################
#############################################################


############################################################
### The observations to fit.
###
### `input` is the marss_inputs.rds written by the
### preprocessing step, or the run directory that holds it.
### It is taken relative to the project root unless it is an
### absolute path, so the same config resolves to the same
### data regardless of the working directory.
###
### `sites` selects the subset of sites to fit, by their
### site key; NULL fits every site present. Restricting the
### site set is the main lever on run time, since the number
### of free parameters grows with it.
###
### `start_date` and `end_date` trim the time axis. The early
### record is in situ and MUR only, so a later start buys a
### denser observation matrix at the cost of years of data.
### NULL leaves that end of the record untrimmed.
###
### NOTE: the exploratory scripts this config replaces listed
### a sixth site, "SF0017", which is not a site key in the
### data -- it matched no rows, so the fits they produced
### used the five sites below. The data do contain "SF017";
### add it here if that is what was meant, bearing in mind
### that it changes the fit and makes it non-comparable with
### the runs already in models/.
############################################################
input      <- "outputs/tasi_salmon_20_yrs_2026-08-29/marss_inputs.rds"
sites      <- c("CB001", "CB004", "SF001", "SF013", "SF007","SF025","SF008")
start_date <- "2006-01-01"
end_date   <- NULL

data_params <- list(
  input = input, sites = sites,
  start_date = start_date, end_date = end_date
)


############################################################
### Scaling of the observations.
###
### The optimizers are better behaved when the observations
### are close to unit scale. Scaling is done per site: every
### instrument at a site is centred and scaled by the mean
### and standard deviation of one reference instrument at
### that same site.
###
### Using one pair per site rather than one pair overall is
### what makes the rest of the model readable. Within a site
### every series is on a common scale, so the estimated
### biases are departures from the reference; and the
### reference itself becomes mean 0, variance 1 at every
### site, which is why its own intercepts are fixed at zero
### below rather than estimated. Every other parameter then
### describes a difference from the reference, at the site
### and instrument level.
###
### `variable` must therefore be an instrument that covers
### every selected site with enough observations to give a
### mean and a standard deviation -- MUR, being a gap free
### interpolated product, is the only one here that does.
### The run refuses to start if it does not.
###
### The constants are written to scales.rds, as vectors in
### the order of `scales$sites`, so the reconstruction step
### can return predictions to degrees C. With scaling
### disabled they are written as mu = 0, sigma = 1, so that
### path is unconditional.
############################################################
enabled  <- TRUE
variable <- "mur_sst"

scaling_params <- list(enabled = enabled, variable = variable)


############################################################
### Structure of the state model.
###
### "site_plus_factors" gives every site its own AR(1) state
### and adds `m_factors` shared dynamic factors, whose
### variance is fixed at 1 for identification. This is the
### full model: site level anomalies plus correlated,
### basin scale departures.
###
### "factors_only" drops the per site states and keeps the
### shared factors alone, with a free innovation variance.
### The site level signal then has to come through the
### seasonal terms, so the two structures differ in what a
### site specific anomaly is allowed to look like.
###
### `m_factors` must not exceed the number of sites; see the
### identification note in R/marss_matrix_functions.R.
############################################################
state_structure <- "site_plus_factors"
m_factors       <- 3


############################################################
### How each instrument enters the observation model.
###
### One entry per variable, keyed by the variable name used
### in the preprocessing config. The order does not matter;
### rows are located by name.
###
###   tag          short name used to build parameter names
###   intercept    "site"            shared site mean, mu_<i>
###                "site+instrument" site mean plus an
###                                  instrument offset,
###                                  mu_<tag>+mu_<i>
###                "independent"     its own mean per site,
###                                  mu_<tag>_<i>
###                "zero"            fixed at 0, not
###                                  estimated -- for the
###                                  instrument the data were
###                                  scaled by
###   seasonality  "shared"      site seasonality, c<k>_<i>
###                "independent" its own, c<k>_<tag>_<i>
###   error        name of its observation variance
###   day_effect   TRUE to give the instrument a day effect:
###                a single covariance cov_<tag> shared
###                between every pair of its own rows
###   site_state   TRUE if the instrument observes the site's
###                own latent anomaly. FALSE puts zeros in
###                the site block of Z, so the instrument
###                sees only the seasonal terms and the
###                shared factors. Ignored under the
###                "factors_only" structure, which has no
###                site states.
###
### MUR is the reference the observations are scaled by, so
### its intercept is fixed at zero: after per-site scaling
### its series has mean 0 at every site by construction, and
### an estimated intercept would be a free parameter with
### nothing left to explain. Every other intercept is then a
### difference from MUR at that site.
###
### It keeps an independent seasonality, and gets no site
### state, because it is an interpolated product: its value
### at a site is a smoothed regional field rather than a
### measurement of that site. Its seasonal cycle is regional,
### so forcing it to share the site terms would bias them;
### and it cannot see a site level anomaly, so loading it on
### the site state would make it evidence about something it
### does not observe.
###
### The high resolution TIR instruments are given a day
### effect because a cloud or a wind event affects every site
### in a scene at once: on any given day their errors shift
### together across sites rather than independently, and
### without the covariance that common shift is read as
### evidence about the state.
###
### A day effect can only be estimated by the EM algorithm.
### It puts a shared covariance in the off-diagonal of R, and
### the direct optimisers require each block of a
### variance-covariance matrix to be fixed, diagonal, or
### wholly unconstrained. Set these to FALSE to fit with
### `method = "TMB"`; the run refuses to start otherwise.
############################################################
instruments <- list(
  insitu_sst = list(
    tag = "insitu", intercept = "site", seasonality = "shared",
    error = "sigma_2_insitu", day_effect = FALSE, site_state = TRUE),
  lst_sst_clean = list(
    tag = "lst", intercept = "site+instrument", seasonality = "shared",
    error = "sigma_2_lst", day_effect = FALSE, site_state = TRUE),
  eco_sst_v002_clean = list(
    tag = "eco", intercept = "site+instrument", seasonality = "shared",
    error = "sigma_2_eco", day_effect = FALSE, site_state = TRUE),
  modis_sst = list(
    tag = "modis", intercept = "site", seasonality = "shared",
    error = "sigma_2_modis", day_effect = FALSE, site_state = TRUE),
  mur_sst = list(
    tag = "mur", intercept = "zero", seasonality = "independent",
    error = "sigma_2_mur", day_effect = FALSE, site_state = FALSE)
)

############################################################
### Prior variance on the initial state. The states are
### scaled, so 1.0 is a diffuse but proper prior.
############################################################
init_var_x0 <- 1.0

structure_params <- list(
  state_structure = state_structure, m_factors = m_factors,
  instruments = instruments, init_var_x0 = init_var_x0
)


############################################################
### Warm start for the seasonal parameters.
###
### The EM algorithm starts every parameter at zero, which
### leaves it to discover the annual cycle -- by far the
### largest signal in the data -- one iteration at a time.
### Seeding the site means and harmonic coefficients with an
### ordinary least squares fit removes that phase of the run.
###
### `shared_from` names the rows used to seed the shared site
### terms, mu_<i> and c<k>_<i>. Any instrument may be named.
### MUR is used here because it is gap free at every site: an
### interpolated product is only a proxy for the site's own
### seasonal cycle, but a proxy at all five sites is a better
### starting point than an exact answer at the two sites with
### in situ records. This is only where the optimiser begins.
###
### `seed_independent` also seeds each instrument that has
### its own mean or seasonality -- here MUR's mu_mur_<i> and
### c<k>_mur_<i> -- from that instrument's own rows. Those
### parameters otherwise start at zero, which for MUR means
### starting with no seasonal cycle at all.
###
### Sites with fewer than `min_obs` observations are left at
### the default start rather than seeded from a fit that the
### record cannot support.
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
### Fitting.
###
### `method` selects the optimiser for the final fit:
###
###   "kem"   the EM algorithm. Slow but robust, and the only
###           method that can estimate a day effect.
###   "BFGS"  quasi-Newton maximisation of the likelihood.
###           Estimates B, but like TMB cannot fit the
###           off-diagonal covariance a day effect adds.
###   "TMB"   direct maximisation through automatic
###           differentiation. Much the fastest, but see the
###           init block below: it does not estimate B at all.
###
### MARSS also accepts "BFGS_TMB" and "nlminb_TMB", which are
### subject to the same constraint on R.
###
### A method that cannot fit a configured day effect is
### refused rather than quietly fitting a reduced model:
### dropping a term the config asks for is a modelling
### decision. Set `day_effect = FALSE`, or pass
### --no-day-effects, to make it explicit.
###
### `controls` is keyed by method because MARSS rejects a
### control belonging to a different one: "minit", "abstol",
### "safe" and "conv.test.slope.tol" are EM only, while
### "eval.max", "iter.max" and the tolerances below are
### specific to the nlminb call TMB uses. See ?MARSS.
###
### Under EM the fit is run as a series of short chunks
### rather than one long call. Each chunk is saved, so a run
### that is interrupted -- or that is still improving after
### days of iterations -- leaves a usable fit behind, and the
### log likelihood trace across chunks shows whether it is
### still climbing. `chunks` caps the total; a run stops
### early as soon as MARSS reports convergence, which under
### TMB is normally at the first chunk.
###
### `warmup_controls` is the throwaway fit used only to
### recover the parameter vector ordering that the warm start
### writes into, so it runs for as few iterations as MARSS
### allows. It is always run under EM, whichever method the
### real fit uses, since the ordering is a property of the
### model rather than of the optimiser.
############################################################
method <- "BFGS"

controls <- list(
  
  kem = list(
    trace               = 1,
    maxit               = 20,
    minit               = 5,
    abstol              = 0.001,   # log likelihood convergence tolerance
    conv.test.slope.tol = 0.5,
    safe                = TRUE     # slower, more robust
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
### Initialisation stage.
###
### TMB reaches a good fit in one pass of a few minutes,
### where EM is still climbing after days. It cannot fit the
### whole model, so it is used here to find starting values
### rather than answers: the run fits a reduced model with a
### fast method first, then hands those estimates to the
### method above, which fits the model this file actually
### describes.
###
### Two restrictions define the reduced model.
###
### TMB does not estimate B. It returns whatever value B
### started at, and reports it as though it had been
### estimated -- on a test series with a true rho of 0.6 it
### returns exactly 1.0000. So B is held at `B_values` and
### removed from the parameter set, which at least makes what
### was and was not estimated legible from the fit. That
### makes B_values a modelling choice, not a tuning knob:
### 0.9 is strongly persistent day to day but stationary, so
### the states of this stage stay well behaved. A named list,
### e.g. list(rho_chi = 0.9, rho_eta = 0.95), sets the site
### and factor blocks separately.
###
### Day effects are dropped, for the same reason the direct
### optimisers refuse them: they put a shared covariance in
### the off-diagonal of R. Dropping them is safe here in a
### way it would not be for the final fit, because B and the
### day effect covariances are then simply not among the
### parameters transferred -- they keep their own starting
### values in the second stage.
###
### The two stages fit different models, so estimates are
### carried across by parameter name rather than by position.
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
### Reconstruction of temperatures from the fit.
###
### The fit leaves its states as scaled anomalies, with the
### seasonal cycle and the site means held separately as
### parameters. The reconstruction step puts them back
### together and returns the result to degrees C, in long
### format, one row per site and date. It runs at the end of
### every fitting run, and can be re-run on a finished run
### directory with R/run_reconstruct.R.
###
### The site view is always written, to states.csv/.rds. It
### is the fused reconstruction: the site's own latent
### temperature, informed by every instrument that observes
### it, built from the site state, the shared factors, the
### shared seasonal terms c<k>_<i> and the site means mu_<i>.
###
### `instruments` names variables to additionally reconstruct
### as they are seen by that instrument alone, each written
### to states_<tag>_only.csv/.rds. MUR is the interesting one
### and the default: it has `site_state = FALSE`, so its view
### is built from the shared factors and its own regional
### seasonal cycle with the site level anomaly left out. The
### difference between states.csv and states_mur_only.csv is
### therefore what the higher resolution instruments add to
### the interpolated product at a site. Set to NULL for the
### site view alone.
###
### The standard errors written alongside the values are
### state uncertainty at the parameter estimates: how well
### the Kalman smoother pins down the states, given the
### fitted parameters. They are approximate in two respects.
### They exclude parameter uncertainty entirely -- that is
### what `parameter_uncertainty` records, and FALSE is the
### only value implemented; a config asking for TRUE is
### rejected before the fit starts rather than after it. And
### they are propagated as if the states were independent,
### which ignores the smoother's covariance between a site
### state and the shared factors.
############################################################
enabled               <- TRUE
reconstruct_from      <- "mur_sst"
formats               <- c("csv", "rds")
parameter_uncertainty <- FALSE

reconstruction_params <- list(
  enabled               = enabled,
  instruments           = reconstruct_from,
  formats               = formats,
  parameter_uncertainty = parameter_uncertainty
)


############################################################
### Where a run is written.
###
### A run directory holds every chunk of the fit, the
### observation matrix, scaling constants, dates and design
### it was fit to, and a copy of this file, so a fit can
### always be traced back to the assumptions that produced it
### and handed to the reconstruction step.
###
### `output_root` is taken relative to the project root
### unless it is an absolute path. `run_name` names the
### subdirectory within that root; NULL derives it from this
### file's name and the date.
###
### The command line flags --outdir and --run-name override
### these values.
############################################################
output_root <- "examples/models/east_storm_bay"
run_name    <- NULL

output_params <- list(output_root = output_root, run_name = run_name)


model_config <- list(
  data           = data_params,
  scaling        = scaling_params,
  structure      = structure_params,
  inits          = inits_params,
  fitting        = fitting_params,
  reconstruction = reconstruction_params,
  output         = output_params
)
