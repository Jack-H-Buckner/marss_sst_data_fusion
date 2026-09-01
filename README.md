# Sea surface tempeaure data fusion with state space models 

This package uses state space models to fuse multiple sources of ocean temperature data into on longterm time series. The primary goal of the model is to leverage high resolution thermal infared sensors to detect high temperature events in hear shore waters where standard ocean temperaute data products can be baised. 

## Preprocessing

The first step takes a raw long format extract from `coastal_sst_data` and applies the
quality control the model depends on. It produces two things: a cleaned dataset, and a
complete record of what was removed and why, so every screening decision can be
inspected rather than taken on trust.

### Running it

```
Rscript R/run_preprocessing.R --input data/raw/tasi_salmon_20_yrs.csv
```

The script finds the project root from its own location, so it can be run from any
working directory. For interactive work the same pipeline is available directly:

```r
source("R/00_preprocessing.R")
source("R/00_diagnostics.R")

dat <- read.csv("data/raw/tasi_salmon_20_yrs.csv")
res <- run_pipeline(dat, "parameters/preprocessing.R")

res$data             # cleaned long format data
res$removed          # every removed row, tagged with the stage that removed it
res$removal_summary  # rows in / out / removed per stage
res$diagnostics$tables$variable_summary

write_pipeline_outputs(res, "outputs/my_run", res$config)   # optional
```

On the full 20-year extract this reads 9.9M rows, formats 4.24M of them, and removes
about 9,000 — roughly 2.5 minutes and 5 GB of memory.

### Options

| Flag | Meaning |
|---|---|
| `--input PATH` | Input CSV. Default `data/raw/tasi_salmon_20_yrs.csv` |
| `--config PATH` | Config file. Default `parameters/preprocessing.R` |
| `--outdir PATH` | Output directory, used verbatim. Overrides everything below |
| `--run-name NAME` | Run subdirectory name. Overrides the config's `run_name` |
| `--group-id NAME` | Override the config `group_id`. Pass `none` for pooled statistics only |
| `--no-diagnostics` | Skip the diagnostic tables and figures |
| `--no-variance` | Skip the match-up variance summary |
| `--overwrite` | Allow writing into a non-empty output directory |
| `--quiet` | Suppress progress messages |
| `--help` | Full option list |

A run refuses to write into a directory that already has files in it unless
`--overwrite` is given, so one run cannot be scattered through the artifacts of another.

### What gets written

```
outputs/<run-name>/
  clean_data.rds            cleaned long format data -- the input to the MARSS step
  removed_all.csv           every removed observation, tagged with `stage` and the
                            diagnostic that condemned it (residual, wind speed,
                            solar time, reference value, reason)
  tables/
    variable_summary.csv           distribution of each variable, kept and removed
    variable_summary_by_group.csv  the same, broken down by `group_id`
    removal_by_stage.csv           which stage removed what, per variable and group
    stage_counts.csv               rows in / out / removed, in pipeline order
    stage_NN_<stage>_<var>_*.csv   each stage's own summary, and the fitted
                                   coefficients for the match-up filter
    match_up_variance.csv          target vs reference statistics around the 1:1
                                   line -- the source of the observation error
                                   estimates for the MARSS R matrix
  diagnostics/
    dist_<var>.png                 full-sample distribution, removed values in red
    by_group_<var>.png             the same across groups; faceted histograms for a
                                   few groups, ordered boxplots for many
    by_group_removed_<var>.png     percent removed per group
    removal_cascade.png            how much each stage took out of each variable
    stage_NN_<stage>_<var>_*.png   the plot each screening step produced
    match_up_variance_*.png        match-ups against the 1:1 line
  diagnostics_report.pdf    every figure above, in one document
  config_used.R             a copy of the configuration this run used
  run_log.txt               the console log and sessionInfo()
```

Counts in the tables carry two denominators, because the extract is a padded grid: most
variables have a row at every key but a value at only a few. `n_input` counts rows,
`n_present` counts rows that actually carry a value, and `pct_removed` uses the latter —
the same denominator the figures show.

### Configuring the run

`parameters/preprocessing.R` is the single source of truth. It defines
`preprocess_config`, whose blocks are:

| Block | Controls |
|---|---|
| `group_id` | the column diagnostics break down by (`"point_id"`, `"aoi"`, or `NULL`) |
| `value_var` | which extracted variables to keep, and the name each is given |
| `convert_units` | Kelvin to Celsius conversion |
| `physics_filter` | plausible temperature bounds |
| `wind_filter` | rejection of skin temperatures under calm daytime conditions |
| `outliers` | seasonal-residual outlier passes |
| `matchup_filter` | screening against a co-located reference |
| `diagnostics` | figure sizes, histogram bins, when to switch to boxplots |
| `match_up_variance` | which targets and references to compare, and how |
| `output` | `output_root` and `run_name` — where the run is written |

The configuration is validated against the data before any work is done, and every
problem is reported at once: a variable name that does not exist, bounds the wrong way
round, a `group_id` that is not a column. A misspelled block name is warned about
immediately, since it would otherwise leave a stage silently running on nothing.

Where a run lands is resolved in this order:

1. `--outdir` — used verbatim
2. `--run-name` — a subdirectory of the config's `output_root`
3. the config's `output$run_name`
4. the input file name plus the date

A relative `output_root` is taken from the project root, not the working directory, so
the same config writes to the same place wherever it is launched from. Each run copies
its configuration to `config_used.R`, so a result can always be traced back to the
assumptions that produced it.

### Using the results

```r
clean <- readRDS("outputs/<run-name>/clean_data.rds")
```

`R/01_format_data_for_marss.R` turns that into the MARSS observation matrix. Note that
its `site.var` and `time.var` defaults (`"site_id"`, `"date"`) do not match the column
names the pipeline produces, so pass them explicitly:

```r
source("R/01_format_data_for_marss.R")
obs <- convert_long_format_marss(clean, time_step = 1,
                                 site.var = "point_id", time.var = "time")
```

## Model fitting

The second step fits the MARSS model to the observation matrix the preprocessing step
wrote. Every choice that used to be hard coded in a script — which sites, which dates,
which instruments share a mean or a seasonal cycle, how the states are structured — is a
value in `parameters/example_model.R`, so a fit is fully described by its config file.

### Running it

```
Rscript R/run_model.R --config parameters/example_model.R
```

Like the preprocessing script it finds the project root from its own location, so it can
be run from any working directory. For interactive work the pieces are available
directly:

```r
source("R/marss_matrix_functions.R")
source("R/02_model_specification.R")

cfg <- read_model_config("parameters/example_model.R")
dat <- readRDS(cfg$data$input)
validate_model_config(cfg, dat)

md    <- build_model_data(dat, cfg)     # subset, scale, index rows by instrument
model <- build_marss_model(md, cfg)     # Z, A, R, B, U, Q, D, x0, V0, d
```

### The model

The observation equation is `y_t = Z x_t + a + D d_t + v_t`, with `v_t ~ N(0, R)`, and
the state equation is `x_t = B x_{t-1} + w_t`, with `w_t ~ N(0, Q)`. `d_t` holds the
seasonal harmonics built by the preprocessing step.

Two state structures are available, selected by `structure$state_structure`:

| Value | States | Q on the factors |
|---|---|---|
| `site_plus_factors` | one AR(1) per site, plus `m_factors` shared dynamic factors | fixed at 1, which identifies the factor scale |
| `factors_only` | the shared factors alone | estimated, so the loadings are no longer separately identified |

`parameters/example_model.R` and `parameters/model_seasonal_diffs.R` are the same configuration
under the two structures.

Each instrument's place in the observation model is one entry in
`structure$instruments`:

| Field | Values |
|---|---|
| `intercept` | `site` → `mu_<i>`; `site+instrument` → `mu_<tag>+mu_<i>`; `independent` → `mu_<tag>_<i>`; `zero` → fixed at 0 |
| `seasonality` | `shared` → `c<k>_<i>`; `independent` → `c<k>_<tag>_<i>` |
| `error` | name of the observation variance on the diagonal of `R` |
| `day_effect` | `TRUE` adds `cov_<tag>`, one covariance shared between every pair of that instrument's rows |
| `site_state` | `TRUE` if the instrument observes the site's own latent anomaly; `FALSE` puts zeros in the site block of `Z` |

MUR is the reference the observations are scaled by, so its intercept is **fixed at zero**:
after per-site scaling its series has mean 0 at every site by construction, and an
estimated intercept would be a free parameter with nothing left to explain. Every other
intercept then reads as a difference from MUR at that site.

It keeps an independent seasonality and gets **no site state**, because it is an
interpolated product: its value at a site is a smoothed regional field rather than a
measurement of that site. Its seasonal cycle is regional, so forcing it to share the site
terms would bias them; and it cannot see a site-level anomaly, so loading it on the site
state would make it evidence about something it does not observe.

The high resolution TIR instruments are given a day effect because a cloud or wind event
affects a whole scene at once: on any given day their errors shift together across sites
rather than independently.

A site state that no instrument loads on has no observations and is not identified, so
`build_observation_matrix()` refuses that combination by name — it is easy to arrive at by
narrowing `data$sites` until a site's only remaining instrument is one excluded from the
site block.

### Fitting method

`fitting$method` selects the optimiser for the final fit. No method fits everything:

| | `kem` | `BFGS` | `TMB` |
|---|---|---|---|
| estimates `B`? | yes | yes | **no** |
| off-diagonal `R` (day effects)? | yes | **no** | **no** |
| speed on this model | does not converge in 75 chunks | — | converges in one ~3 min call |

`BFGS_TMB` and `nlminb_TMB` are also accepted and behave like `TMB` here.

**TMB does not estimate `B`.** It returns whatever value `B` started at and reports it as
though it had been estimated — on a test series with a true `rho` of 0.6, TMB returns
exactly `1.0000` while `kem` gives 0.544 and `BFGS` 0.539. This is why the initialisation
stage below holds `B` at a stated value and removes it from the parameter set, rather than
leaving it free and trusting the optimiser to move it.

**A day effect can only be fit by `kem`.** It puts a shared covariance in the off-diagonal
of `R`, and the direct optimisers require each block of a variance-covariance matrix to be
fixed, diagonal, or wholly unconstrained. `validate_model_config()` refuses the
combination up front, naming the instruments involved, rather than letting MARSS report it
as an opaque model specification error minutes into a run. To fit with TMB or BFGS as the
*final* method, either set `day_effect = FALSE` on those instruments or pass
`--no-day-effects`.

`fitting$controls` is keyed by method, because MARSS rejects a control belonging to a
different one — `minit`, `abstol`, `safe` and `conv.test.slope.tol` are EM only, `REPORT`
and `reltol` belong to BFGS, and `eval.max`, `iter.max` and the nlminb tolerances to TMB.

### Two-stage fitting

Because TMB is far the fastest but cannot fit the whole model, `fitting$init` uses it to
find *starting values* rather than answers. The run then has two stages:

1. **Initialisation** — a reduced model, `B` held at `init$B_values` and the day effects
   dropped, fit with `init$method` (TMB by default). Saved as `fit_init.rds`.
2. **The full fit** — the model the config actually describes, with `fitting$method`,
   started from stage 1's estimates.

The two stages fit different models, so estimates are carried across **by parameter name**,
not by position. Whatever stage 1 could not estimate simply has no name to match and keeps
its own starting value:

```
stage 1 (TMB, B fixed, no day effects):  68 params
stage 2 (kem, full model):               72 params
transferred by name:                     68
left at defaults:  R.cov_lst, R.cov_eco, B.rho_chi, B.rho_eta
```

Those are exactly the parameters stage 1 is structurally unable to estimate, which is what
makes dropping them safe here and not for the final fit.

### Where B starts

`init$B_values` governs **both** stages: stage 1 holds `B` fixed there, and the final fit,
where `B` is free again, starts from it — whichever method that fit uses.

That second part matters. `B` is the one parameter nothing else sets: stage 1 holds it at a
fixed value, so it is not a parameter there and has no name to transfer, and the warm start
seeds only the seasonal terms. Left alone, the final fit would begin `B` wherever the
throwaway ordering fit happened to leave it after two EM iterations from MARSS's own
default of 1 — a random walk, and a starting point that shifts whenever `warmup_controls`
is touched. Under TMB, which does not estimate `B`, that arbitrary value is not just the
starting point but the reported answer.

It is a single number applied to every free entry of `B`, or a named list for per-block
control:

```r
B_values <- 0.975                                    # every free entry
B_values <- list(rho_chi = 0.975, rho_eta = 0.99)    # site and factor blocks
```

0.975 is strongly persistent yet stationary: near-shore temperature anomalies decay on the
scale of weeks rather than days, and a value near 1 reflects that without letting the
states wander as a random walk would. Because TMB holds `B` at whatever it is given, this
is a modelling choice rather than a tuning knob.

Pass `--no-init` to skip the initialisation stage entirely — `B` still starts at
`init$B_values` — or `--init-method` / `--init-B` to override it.

Note that `--no-init` (skip the initialisation stage) and `--no-warm-start` (skip the
least-squares seeding of the seasonal terms) are different things; the warm start feeds
whichever fit runs first.

### Restarting a run from a saved fit

A run directory keeps every stage of its fit, so neither stage has to be paid for twice.

```
# skip stage 1, hand its estimates to a different final optimiser
Rscript R/run_model.R --run-name v3_kem --method kem --init-from models/v2_bfgs

# continue a run that stopped at the chunk limit while still climbing
Rscript R/run_model.R --run-name v2_kem_more --resume-from models/v2_kem

# hand a converged BFGS fit to EM so the day effects can be estimated
Rscript R/run_model.R --run-name v3_kem_day --method kem \
  --resume-from models/v2_bfgs/fit_chunk_01.rds
```

`--init-from` replaces the initialisation stage with one already on disk. This is safe
across final methods because **stage 1 does not depend on one**: it always holds `B`
fixed, always drops the day effects, and always uses `init$method`, whatever `method`
says. A saved stage 1 is reusable as long as `data`, `scaling`, `structure` and
`fitting$init` are unchanged.

`--resume-from` starts the final fit from a saved one. Either flag accepts a run
directory — resolving to `fit_init.rds` or to the highest-numbered `fit_chunk_NN.rds` —
or a path to the fit itself. Only one of the two may be given.

Both transfer estimates **by parameter name**, the same mechanism the two-stage fit uses,
so the model being fit does not have to be the one that produced them. The case worth
understanding is the third example: a fit produced under BFGS had to have the day effects
dropped, and continuing it under EM with them on is the natural next step. Those two
covariances have no counterpart in the saved fit, so they are left at their defaults and
reported as such:

```
resume | chunk NA | iters    78 | logLik 19291.724
Stage 2: kem.
Transferred 65 of 67 parameters from the resumed fit.
  left at defaults (not estimable in that stage): R.cov_lst, R.cov_eco
fit   | chunk  1 | iters     3 | logLik 19302.032
```

The loaded fit is copied into the new run directory (as `fit_init.rds` or
`fit_resumed.rds`) and recorded in `loglik.csv` as an `init` or `resume` row, so the new
directory stays self-contained and the trace shows the log likelihood the run started
from. A run refuses to start if the saved fit was fit to different observations — the
parameter names would still match, so the transfer would otherwise be silently wrong.

### Options

| Flag | Meaning |
|---|---|
| `--config PATH` | Config file. Default `parameters/example_model.R` |
| `--input PATH` | `marss_inputs.rds`, or the preprocessing run directory holding it |
| `--outdir PATH` | Output directory, used verbatim. Overrides everything below |
| `--run-name NAME` | Run subdirectory name. Overrides the config's `run_name` |
| `--sites LIST` | Comma separated site keys. Pass `all` for every site |
| `--start-date DATE`, `--end-date DATE` | Trim the time axis. Pass `none` to remove a bound |
| `--structure NAME` | `site_plus_factors` or `factors_only` |
| `--m-factors N` | Number of shared dynamic factors |
| `--method NAME` | `kem`, `TMB`, `BFGS`, `BFGS_TMB` or `nlminb_TMB` |
| `--no-day-effects` | Drop the per-instrument day effects. Required by every method other than `kem` |
| `--chunks N` | Maximum number of restarted fitting chunks |
| `--maxit N` | Iterations per chunk, for the selected method |
| `--no-scaling` | Fit in degrees C rather than scaling the observations |
| `--init-method NAME` | Method for the initialisation stage. Default `TMB` |
| `--init-B VALUE` | Value `B` is held at during initialisation, and where the final fit starts it. Default `0.975` |
| `--no-init` | Skip the initialisation stage and fit directly |
| `--init-from PATH` | Start from a saved initialisation fit: a run directory, or a `fit_init.rds` |
| `--resume-from PATH` | Start the final fit from a saved one: a run directory, or a `fit_chunk_NN.rds` |
| `--no-warm-start` | Skip the least squares seeding of the seasonal terms |
| `--no-figures` | Skip the coverage figure |
| `--overwrite` | Allow writing into a non-empty output directory |
| `--quiet` | Suppress progress messages |
| `--help` | Full option list |

Command line values are folded into the config before anything reads it, so
`config_used.R` plus `run_log.txt` is a complete record of the run.

### What gets written

```
models/<run-name>/
  fit_init.rds         the initialisation stage fit, if one ran
  fit_chunk_NN.rds     the fitted model after each chunk
  loglik.csv           stage, method, log likelihood and iteration count
  observations.rds     the observation matrix and its row keys
  scales.rds           the constants that return predictions to degrees C
  dates.rds            the date range the columns span
  design.rds           n_site, m_factors, state_structure, n_harmonics,
                       the instrument specification, method, day_effects
  states.csv/.rds      the fused reconstruction, in degrees C, long format
  states_mur_only.*    the same as MUR alone sees it, one file per instrument
                       named in `reconstruction$instruments`
  coverage.png         which cells of the observation matrix are observed
  config_used.R        a copy of the configuration this run used
  run_log.txt          the console log and sessionInfo()
```

`design.rds` carries the *resolved* instrument specification rather than a pointer at
`config_used.R`, because the config file is copied before the command line overrides are
applied: a run launched with `--structure` or `--m-factors` is described by `design.rds`
and not by the file sitting next to it.

Under EM the fit is run as a series of short restarted chunks rather than one long call.
Each is saved as it completes, so a run that is interrupted — or that is still climbing
after days of iterations — leaves a usable fit behind, and `loglik.csv` shows whether it
is still improving. The run stops as soon as MARSS reports convergence, or at
`fitting$chunks`; under TMB that is normally the first chunk.

### Scaling and the warm start

Two config blocks are worth understanding before changing them.

`scaling` is **per site**: every instrument at a site is centred and scaled by the mean and
standard deviation of one reference instrument — `scaling$variable` — at that same site.

That choice is what makes the rest of the model readable. Within a site every series is on
a common scale, so the estimated biases are departures from the reference; and the
reference itself becomes mean 0, variance 1 at every site, which is why its intercept is
fixed at zero rather than estimated. Every other parameter then describes a difference
from the reference, at the site and instrument level.

`scaling$variable` must therefore cover **every** selected site densely enough to give a
mean and a standard deviation — MUR, a gap-free interpolated product, is the only
instrument here that does, and the run refuses to start if the reference falls short at
any site. The constants go into `scales.rds` as vectors ordered by `scales$sites`; with
scaling disabled they are `mu = 0, sigma = 1` throughout, so downstream code does not have
to know which was used.

`inits` seeds the seasonal parameters from ordinary least squares, which saves the EM
algorithm from discovering the annual cycle one iteration at a time. `shared_from` names
the rows used for the shared site terms `mu_<i>` and `c<k>_<i>` — any instrument may be
named, and a gap free product covering every site is usually a better starting point
than an exact record covering two of them. `seed_independent` additionally seeds each
instrument that has its own mean or seasonality from its own rows; those parameters
otherwise start at zero, which for MUR means starting with no seasonal cycle at all.

## Reconstruction

The fit leaves its states as scaled anomalies, with the seasonal cycle and the site means
held separately as parameters. The third step puts them back together and returns the
result to degrees C. It runs automatically at the end of every fit, controlled by the
`reconstruction` block of the config, and writes long format tables next to the fit:

| Column | |
|---|---|
| `date` | the date the column spans |
| `site` | the site key |
| `variable` | `sst` for the fused reconstruction, the instrument's variable name otherwise |
| `value` | degrees C |
| `se` | standard error, degrees C |

### The two views

The **site view**, always written to `states.csv`/`states.rds`, is the fused
reconstruction: the site's own latent temperature, informed by every instrument that
observes it. It is built from the site state, the shared factors, the shared seasonal
terms `c<k>_<i>` and the site means `mu_<i>`.

An **instrument view**, written to `states_<tag>_only.*` for each variable named in
`reconstruction$instruments`, is what the model says that instrument alone sees at each
site, using that instrument's own intercept and seasonal cycle as
`structure$instruments` defines them, and — when it has `site_state = FALSE` — only the
shared factors rather than the site's own anomaly.

MUR is the default and the interesting one. Because it has no site state, its view is the
regional signal with the site-level anomaly left out, so the difference between
`states.csv` and `states_mur_only.csv` is what the higher resolution instruments add to
the interpolated product at a site. Set `instruments` to `NULL` for the site view alone.

An instrument whose parameterisation is the shared one — in situ, here — reconstructs to
exactly the site view, which is a useful check that a view is being built from the
parameters you think it is.

### What the standard errors are, and are not

They are state uncertainty at the parameter estimates: how well the Kalman smoother pins
down the states, given the fitted parameters. They are approximate in two respects.

They exclude **parameter uncertainty** entirely — that is what
`reconstruction$parameter_uncertainty` records, and `FALSE` is the only value implemented.
A config asking for `TRUE` is rejected by `validate_model_config()` before the fit starts
rather than hours later.

And they are propagated as `sqrt(Z^2 %*% states.se^2)`, which treats the states as
independent and so ignores the smoother's covariance between a site state and the shared
factors. The full covariance is available as `fit$kf$VtT` if this is tightened up.

### Re-running it on a finished fit

A fit can take hours, so the step also stands alone and can be re-run against a run
directory without refitting:

```
Rscript R/run_reconstruct.R --run models/v2_bfgs
Rscript R/run_reconstruct.R --run models/v2_bfgs --instruments mur_sst,modis_sst
Rscript R/run_reconstruct.R --run models/v2_bfgs --fit fit_chunk_03.rds --formats csv
```

| Flag | Meaning |
|---|---|
| `--run PATH` | The run directory to reconstruct. Required |
| `--fit NAME` | Which fit in it to use. Default: the highest numbered `fit_chunk_NN.rds`, or `fit_init.rds` if the run has no chunk |
| `--instruments LIST` | Comma separated variables to write an instrument view for. Pass `none` for the site view alone |
| `--formats LIST` | `csv`, `rds`, or both |
| `--outdir PATH` | Where to write. Default: the run directory itself |
| `--config PATH` | Take the `reconstruction` block from this config instead of the run's `design.rds` |

A run directory written before this step existed carries a shorter `design.rds`; the
instrument specification is then taken from the `config_used.R` beside it, with a message
saying so. For interactive work:

```r
source("R/marss_matrix_functions.R")
source("R/02_model_specification.R")
source("R/03_reconstruct_states.R")

run    <- read_run_dir("models/v2_bfgs")
states <- reconstruct_states(run$fit, run$scales, run$observations, run$design,
                             instruments = "mur_sst")

states$site      # the fused reconstruction, long format
states$mur_sst   # MUR's own view of the same sites and dates
```

`R/reconstricting_temperatures.R` and `R/reconstructing_temps_mur_fe.R` are the
exploratory predecessors of this step and are superseded by it. Both assume there are two
harmonic rows where these fits have four, so neither runs correctly against a current fit.
