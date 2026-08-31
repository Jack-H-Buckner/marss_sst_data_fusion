# Formalizing the MARSS fitting step into a configurable CLI stage

Status: **implemented and verified**, 2026-08-30. Both structures reproduce the
exploratory scripts' matrices exactly, and both configs run end to end.

Generated with Claude Code.

---

## Context

`EDA/base_marss_model.R` and `EDA/seasonal_diffs_marss_model.R` are the prototypes of the
model-fitting stage. They are ~94% byte-identical: only `Z`, `B` and `Q` differ (the base
model carries per-site random-walk states plus 2 DFA factors; the seasonal-diffs model drops
the site states and keeps the factors only, with a free `Q`).

The repo already has a finished template for a pipeline stage: `R/run_preprocessing.R`
(base-R CLI, `.script_path()` root resolution, run-directory resolution, tee'd logging,
`config_used.R` copy) driven by `parameters/preprocessing.R` (a plain-R `preprocess_config`
list). `R/02_model_specification.R` was a one-line stub — the placeholder for this work.

**Goal:** one entry point `Rscript R/run_model.R --config parameters/model.R`, constants in
a plain-R config, model construction in roxygen'd builder functions, both experiments
reproducible as two config files.

## Design decisions (agreed)

- Config is a **plain R list file** matching `parameters/preprocessing.R` — no new dependency.
- **One** entry script; a `state_structure` option (`"site_plus_factors"` | `"factors_only"`)
  selects the variant. Two config files capture the two experiments.
- **Fix all four** known bugs (below).
- Scaling is **configurable**: when on, a **single** `(mu, sigma)` pair is applied to the whole
  matrix, computed from a user-specified **instrument + site** (generalizing the old `y[1,]`,
  which was implicitly the first in-situ site).

---

## FINDING 1: `SF0017` is not a site in the data

**Resolved by preserving behaviour** — both configs now list the five sites that were
actually fit, with a comment pointing at `SF017`. Flip it if the sixth site was intended.

Both EDA scripts declare:

```r
sites <- c("CB001","CB008","SF001","SF0017","SF025","SF050")
```

`SF0017` is **not** a site key in `outputs/tasi_salmon_20_yrs_2026-08-29/marss_inputs.rds`.
The data contain `SF001` … `SF058`, including `SF017`, but no `SF0017`. Because the original
selects rows with `%in%`, the missing key matched nothing and was silently dropped — **every
fit produced by these scripts used five sites, not six.**

The new `validate_model_config()` catches this and refuses to run, which is how the finding
surfaced.

---

## FINDING 2: the seasonal half of the warm start never did anything

The exploratory scripts seed starting values with:

```r
init_vec[paste0("A.c1_", i)] <- pars$c1[i]
init_vec[paste0("A.c2_", i)] <- pars$c2[i]
init_vec[paste0("A.mu_",  i)] <- pars$mu[i]
```

In a fit of this model MARSS names the parameters `A.mu_1` … but `D.c1_1`, `D.c2_1`, …
— the means live in `A`, the harmonic coefficients in `D`. Confirmed directly on a fit:

```
A.* : A.mu_1, A.mu_2, A.mu_lst, A.mu_3, ...
D.* : D.c1_1, D.c1_2, ..., D.c1_mur_1, ...
"A.c1_1" in the parameter vector? FALSE     "D.c1_1"? TRUE
```

Assigning to a name that is not in the vector adds a new, ignored element rather than
erroring, so **only the site means were ever seeded**; every seasonal coefficient started
at zero regardless. `apply_inits()` now looks for each name under both the `A.` and `D.`
prefixes, which is why the run log reports `seeded 50 of 50`.

---

## FINDING 3: MUR's own parameters were never seeded

`seasonal_inits(..., mur_rows)` fits MUR's series and writes into `mu_<i>` / `c<k>_<i>`,
the *shared site* terms. Using MUR as the proxy is reasonable — it is gap free at every
site, whereas in situ exists at only two of five — but MUR's own `mu_mur_<i>` and
`c<k>_mur_<i>` were left at zero, i.e. starting with no seasonal cycle at all for the one
instrument given an independent one.

The `inits` block now separates the two: `shared_from` names the source for the shared
terms (MUR, deliberately), and `seed_independent` also seeds each instrument that has its
own terms from its own rows.

---

## FINDING 4: the reconstruction step cannot read these fits

`get_states_mur_factors()` in `R/reconstricting_temperatures.R` does
`model$par$A[paste0("c1_", i), ]` and `matrix(model$model$fixed$d, nrow = 2)`. Neither
holds: `fit$par$A` has no dimnames at all (name indexing errors outright, tested), the
harmonic coefficients are in `D` not `A`, and there are four harmonic rows, not two.

This is **pre-existing** — fits from the exploratory scripts have the same layout — and was
left alone as out of scope. `design.rds` now carries `n_harmonics` and `state_structure`
so the reconstruction can be made general when it is updated.

---

## What was written

| File | State |
|---|---|
| `parameters/model.R` | new — the base config |
| `parameters/model_seasonal_diffs.R` | new — the `factors_only` variant |
| `R/02_model_specification.R` | replaces the one-line `library(MARSS)` stub |
| `R/run_model.R` | new — the CLI entry point |
| `README.md` | new *Model fitting* section |

**Nothing in this project is under version control** — the git repo at
`/Users/johnbuckner/github` has zero tracked files. The two exploratory scripts were
therefore *not* deleted, contrary to the original plan: deleting them would be permanent.
They are superseded by the CLI and can be removed once there is a commit to fall back on.

### `parameters/model.R`

Blocks assembled into `model_config`: `data` (input, sites, start_date, end_date),
`scaling` (enabled, variable, site), `structure` (state_structure, m_factors, instruments,
init_var_x0), `inits` (enabled, shared_from, seed_independent, min_obs), `fitting`
(controls, warmup_controls, chunks), `output` (output_root, run_name).

The key generalization is the `instruments` list, which replaces five copy-pasted
per-instrument blocks in each of `a`, `D` and `R` with one loop, while reproducing the
existing parameter names exactly:

```r
instruments <- list(
  insitu_sst         = list(tag="insitu", intercept="site",            seasonality="shared",      error="sigma_2_insitu", cov=NULL),
  lst_sst_clean      = list(tag="lst",    intercept="site+instrument", seasonality="shared",      error="sigma_2_lst",    cov="cov_lst"),
  eco_sst_v002_clean = list(tag="eco",    intercept="site+instrument", seasonality="shared",      error="sigma_2_eco",    cov="cov_eco"),
  modis_sst          = list(tag="modis",  intercept="site+instrument", seasonality="shared",      error="sigma_2_modis",  cov=NULL),
  mur_sst            = list(tag="mur",    intercept="independent",     seasonality="independent", error="sigma_2_mur",    cov=NULL))
```

`intercept`: `"site"` → `mu_<i>`; `"site+instrument"` → `mu_<tag>+mu_<i>`;
`"independent"` → `mu_<tag>_<i>`. `seasonality`: `"shared"` → `c<k>_<i>`;
`"independent"` → `c<k>_<tag>_<i>`. `day_effect = TRUE` → `cov_<tag>` in every
off-diagonal cell between that instrument's rows (reproducing `cov_lst` and `cov_eco`).

### Fitting method

`fitting$method` is one of `kem`, `BFGS`, `TMB`, `BFGS_TMB`, `nlminb_TMB`, and
`fitting$controls` is a list keyed by method — MARSS rejects a control belonging to a
different one. Verified directly: passing the EM control list to `method = "TMB"` gives

```
elements minit is not allowed in arg control for method TMB
elements abstol is not allowed in arg control for method TMB
elements conv.test.slope.tol is not allowed in arg control for method TMB
elements safe is not allowed in arg control for method TMB
```

**Only `kem` can fit a day effect.** Verified on the real model:

```
WITH day effects    -> ERROR: The variance-covariance matrix R is not properly
                       constrained. when optim() or nlminb() are used, blocks in
                       var-cov matrices must either be fixed, diagonal (shared
                       values allowed) or unconstrained (no shared values) at t=1
WITHOUT day effects -> ok, logLik 24152.64, convergence 0
```

`validate_model_config()` catches this before the run starts and names the offending
instruments. The warm-up fit always runs under `kem` regardless of `method`, since the
parameter ordering it recovers is a property of the model, not the optimiser.

### FINDING 5: TMB does not estimate B

It returns whatever value `B` started at, and reports it as an estimated parameter without
saying it never moved. On a two-state test series with a true `rho` of 0.6:

```
kem    logLik  -865.665 | B.rho1=0.5439, B.rho2=0.5206
BFGS   logLik  -837.347 | B.rho1=0.5389, B.rho2=0.4837
TMB    logLik  -921.281 | B.rho1=1.0000, B.rho2=1.0000   <- MARSS's default init
```

This is what motivates the two-stage design: TMB is by far the fastest optimiser here, but
it can only be trusted for the parameters it actually moves. Fixing `B` numerically for the
init stage — rather than leaving it free — makes the fit honest: `B` then disappears from
the parameter set entirely instead of appearing as an estimate that is really an input.

### Two-stage fitting

`fitting$init` runs a reduced model first (`B` fixed at `init$B_values`, day effects
dropped) with a fast method, and hands its estimates to the final fit **by parameter name**,
since the two stages fit different models. Measured on the real model:

```
stage 1 (TMB, B fixed at 0.9, no day effects):  68 params, logLik 24029.52, conv 0, 515 iters
stage 2 (kem, full model):                      72 params
transferred by name:                            68   <- all of them, no orphans
left at defaults:  R.cov_lst, R.cov_eco, B.rho_chi, B.rho_eta
```

Exactly the four parameters stage 1 structurally cannot estimate are the four it does not
supply. That is what makes dropping the day effects safe for the init stage and not for the
final fit — hence the asymmetry: the init stage drops them automatically, while a final
method that cannot fit them is refused.

`transfer_params()` implements the handoff and is deliberately separate from
`apply_inits()`: the latter takes unprefixed names from the least-squares warm start and
has to guess between the `A.` and `D.` forms (Finding 2), whereas here both sides carry
full MARSS names and must match exactly.

---

## Model corrections (from the author, 2026-08-30)

Three changes to what the model says, made after the CLI was working. The first came from
`EDA/base_marss_model.R`; the other two follow from it.

### MUR no longer loads on the site states

Per-instrument field `site_state`. `FALSE` puts zeros in the site block of `Z`, so the
instrument sees the seasonal terms and the shared factors but not the site's own anomaly.
MUR is an interpolated product — its value at a site is a smoothed regional field, not a
measurement of that site — so loading it on the site state made it evidence about something
it cannot see.

Implemented by zeroing rows in place rather than assembling `chi_obs` from row subsets, so
the result does not depend on the instruments happening to fall in a particular order down
the matrix (the EDA version `rbind`s the MUR block last, which is only correct while MUR is
the final row block).

`build_observation_matrix()` now also refuses a site state that no instrument loads on: it
has no observations and is not identified, and it is easy to reach by narrowing
`data$sites` until a site's only remaining instrument is one excluded from the site block.

### Scaling is per site, not one global pair

Each site's observations are centred and scaled by `scaling$variable` at *that* site.
This replaces the single global `(mu, sigma)`; `scales.rds` now holds vectors ordered by
`scales$sites`. The reference must cover every selected site, which in practice means MUR.

Verified: `y_fit` is `all.equal()` to the scaled `y` the EDA script produces.

### MUR's intercepts are fixed at zero

New intercept kind `"zero"`. After per-site scaling by MUR, MUR's series has mean 0 and sd
1 at every site by construction — confirmed on the real data, row means ~1e-16 and row sds
exactly 1.0000 — so an estimated `mu_mur_<i>` would be a free parameter with nothing to
explain. Fixing it makes every other parameter read as a difference from MUR at the site
and instrument level.

This required `a` to become a **list** rather than a character vector, so fixed and
estimated entries can coexist; a character matrix would turn the fixed `0` into a parameter
named `"0"`. Verified that MARSS treats an all-character list matrix identically to a
character matrix (same logLik to 6 dp, same parameter set), so this costs nothing where
nothing is fixed.

Free `A` parameters: 21 → 16. Every matrix except `A` remains `identical()` to the EDA
script's construction.

### `R/02_model_specification.R`

| Function | Responsibility |
|---|---|
| `read_model_config(path)` | source + return `model_config`; mirrors `read_config()` |
| `validate_model_config(cfg, marss_inputs)` | reports **all** problems at once, as `validate_config()` does |
| `build_model_data(marss_inputs, cfg)` | row/column subsetting, `.dense_rank()`, per-instrument row and site-index blocks, scaling |
| `build_observation_matrix(md, cfg)` | `Z`; reuses `dfa_loadings()` from `R/marss_matrix_functions.R` |
| `build_observation_intercepts(md, cfg)` | `a` |
| `build_covariate_effects(md, cfg)` | `D`, sized from `nrow(md$d)` not a hardcoded 4 |
| `build_observation_covariance(md, cfg)` | `R` diagonal + within-instrument covariance blocks |
| `build_state_model(md, cfg)` | `B`, `Q`, `U`, `x0`, `V0` |
| `build_marss_model(md, cfg)` | assembles the model list |
| `seasonal_inits(md, cfg, variable)` | lm warm start for one instrument; loops harmonic rows instead of `d[1..4, ]` |
| `build_inits_seeds(md, cfg)` | shared site terms from `inits$shared_from`, plus each independent-term instrument from its own rows |
| `apply_inits(fit, seeds)` | `MARSSvectorizeparam()` round-trip |
| `fit_marss_chunks(md, model_list, cfg, outdir)` | chunk-restart loop, per-chunk `saveRDS`, `loglik.csv` |
| `write_model_outputs(md, cfg, outdir)` | `observations.rds`, `scales.rds`, `dates.rds`, `design.rds` |

Two robustness details worth keeping if this is ever rewritten: `apply_inits()` looks for
each parameter under **both** `A.<name>` and `D.<name>` (see Finding 2 — this is what makes
the seasonal warm start actually take effect), and `seasonal_inits()` skips rank-deficient
fits rather than writing `NA` into the parameter vector.

### `R/run_model.R`

Structural clone of `R/run_preprocessing.R`. Flags:
`--config --input --outdir --run-name --sites --start-date --end-date --structure
--m-factors --chunks --maxit --no-scaling --no-inits --no-figures --overwrite --quiet --help`,
with the same precedence rule (`--outdir` > `--run-name` > config `run_name` > derived).
`--input` accepts either a `marss_inputs.rds` path or a preprocessing run directory.
CLI overrides are folded into the config **before** anything downstream reads it, so
`config_used.R` plus the log is a complete record of the run.

---

## Bug fixes carried by the refactor

1. **Undefined objects.** `obs_mat_y`, `scales`, `dates` and (in seasonal_diffs) `y_scaled`
   are never defined — leftovers from `R/site_models_only.R:87-92`. Every `saveRDS` at the end
   of both scripts errors. `build_model_data()` defines all four. With scaling off, `scales` is
   written as `list(mu = 0, sigma = 1)` so the reconstruction path is unconditional.
2. **Missing output directory.** `models/marss_base/` does not exist and nothing creates it.
   `run_model.R` does `dir.create(outdir, recursive = TRUE)` after a non-empty-dir guard.
3. **`design.rds` written to the wrong run.** Both scripts write it into
   `models/marss_test_sites/`, clobbering an earlier run. It now goes to the run's own outdir.
4. **Warm start seeded from the wrong rows / into names that did not exist.** See Findings
   2 and 3 — the fix is `build_inits_seeds()` plus the dual `A.`/`D.` lookup. `c3`/`c4` are
   now seeded too; the original computed and discarded them.

Incidental: dead duplicate `sites` assignment (line 13); unreachable `if (is.null(fit))` in the
chunk loop; unused `ggplot2`/`dplyr`/`reshape2` imports (the local `dense_rank` exists only to
shadow `dplyr::dense_rank`, so it becomes `.dense_rank()`); `harminics_constants` typo.

---

## Verification performed

1. **Structural equivalence.** The original EDA construction, pasted verbatim, versus
   `build_marss_model()`: `identical()` on `Z`, `A`, `R`, `B`, `Q`, `D`, `U`, `x0`, `V0`, `d`
   — **all ten identical, for both structures.** The refactor is behaviour preserving.
   The configured scaling reference also reproduces the old implicit `y[1, ]` exactly
   (mu = 14.873, sigma = 2.163).
2. **Smoke runs.** Both configs, `--chunks 1 --maxit 3`. Each wrote `fit_chunk_01.rds`,
   the four record objects, `coverage.png`, `config_used.R`, `loglik.csv` and `run_log.txt`,
   and seeded 50 of 50 starting values. logLik −4241.49 (`site_plus_factors`, 7 states) and
   16033.76 (`factors_only`, 2 states); ~3 min per run on the 21 × 4959 matrix.
   After the method option was added, the EM run reproduced −4241.487 exactly.
2b. **TMB.** `--method TMB --no-day-effects` converged in one chunk, 497 iterations,
   logLik 24160.456, ~4.5 min. `--method TMB` without dropping the day effects is refused
   by name; `--method EM` is refused with the list of valid methods.
3. **CLI.** `--help` exits 0; bad `--config`, unknown flag and `--m-factors 0` each exit 1
   with a readable `ERROR:`; a non-empty output directory is refused without `--overwrite`;
   running from `/tmp` resolves every path identically.
4. **Config validation.** Unknown `state_structure`, `m_factors` above the site count, an
   absent `scaling$variable`, an `inits$shared_from` not in the model, and a variable with no
   `instruments` entry are each caught by name — and four bad fields at once are reported
   together in one error.
5. **Downstream contract.** Tested and **failing** — see Finding 4. Pre-existing.

## Remaining work

1. **Update `get_states_mur_factors()`** for the actual parameter layout (Finding 4): read
   the harmonics from `D` rather than `A`, index by position or restore dimnames, and take
   the harmonic count from `design$n_harmonics` instead of hardcoding 2.
2. **Decide whether `SF017` was meant** (Finding 1).
3. **Reconcile with the primary MARSS orchestrator refactor** — in particular whether
   `R/02_model_specification.R` is the right home for the builders.
4. **Get this project under version control**, then retire
   `EDA/base_marss_model.R` and `EDA/seasonal_diffs_marss_model.R`.
5. **Full runs**, one per config at the default `--chunks 75`, checking `loglik.csv`
   increases monotonically and the run stops early on convergence.
