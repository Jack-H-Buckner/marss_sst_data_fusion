# Speeding up the MARSS fits

## Context

The MARSS fits in this project take hours. A 75-chunk run of
`EDA/seasonal_diffs_marss_model.R` or `EDA/base_marss_model.R` is a ~4–5 hour
job, and both runs started this morning have yet to write a single chunk file.
The question that prompted this plan was whether multiple cores can help.

**They can't — not for one fit.** Two independent reasons:

1. The EM algorithm in `MARSSkem` is strictly sequential: iteration *k+1*
   consumes iteration *k*'s parameters. The `for (chunk in 1:75)` loop is the
   same dependency at a coarser grain — `inits = coef(fit)` chains each chunk
   to the last. There is no parallel axis inside a fit.
2. The matrices are small. `n = 21` observation rows, `m = 7` states
   (`base`) or `m = 2` (`seasonal_diffs`). Threaded BLAS needs matrices an
   order of magnitude larger before thread dispatch pays for itself.

The good news is that the fit is slow for a fixable reason that has nothing to
do with cores, and the fix is larger than anything 6 cores could have bought.

### Evidence

Measured, not assumed:

- Chunk files in `models/marss_test_mur_fe/` are written 3–4 minutes apart.
  At `maxit = 20` that is **~10 s per EM iteration**; 1500 iterations ≈ 4.5 h.
- Each `fit_chunk_*.rds` is **35 MB**; 75 of them is 2.6 GB per run.
- `sessionInfo()` reports reference BLAS (`libRblas.0.dylib`), R 4.5.1
  x86_64, 6 physical / 12 logical cores.
- MARSS 3.11.9, KFAS 1.6.0, TMB 1.9.21 installed. **marssTMB is not.**

### Root cause: `safe = TRUE`

Every call site sets it. From the installed `MARSS:::MARSSkem` source, that
flag calls `rerun.kf()` — a full Kalman filter/smoother pass — after **each of
eight parameter-block updates**:

| line | block | guarded by `!fixed[[…]]`? |
|------|-------|---------------------------|
| 244 | `R`  | yes |
| 361 | `Q`  | yes |
| 614 | `x0` | yes |
| 847 | `U`  | yes |
| 651 | `V0` | **no — reruns even when fixed** |
| 704 | `A`  | **no** |
| 922 | `B`  | **no** |
| 993 | `Z`  | **no** |

`R`, `Q`, `x0` and `U` are all free in these models, so all eight fire. Each
EM iteration therefore costs **~9 Kalman passes over T ≈ 4959 time steps
instead of 1.**

The comment next to the flag — `# slower, more robust — useful while
debugging` — is accurate. It was a debugging setting that became the default.

### Secondary cause: `trace = 1`

Also set at every call site. `trace` accepts any whole number `>= -1`
(`is.marssMLE:139-147`), but MARSS only ever branches on **two** thresholds —
`trace != -1` and `trace > 0` — so for `method = "kem"` there are exactly three
effective levels, and `trace = 1` is identical to `trace = 2` or `trace = 99`:

| level | behaviour |
|-------|-----------|
| `-1` | skips input checking (`MARSS():71`), `MARSSkemcheck()` model validation (`MARSSkem:6`), and the post-fit checks at `MARSS():98,105,147,168`. Fastest; the Rd notes these checks "are time expensive so this can speed up model fitting." |
| `0` | **the MARSS default.** Basic error checking, brief error messages. |
| `> 0` | everything below. |

What `trace > 0` costs, specifically:

- `MARSSkem:124` does `iter.record$par <- rbind(iter.record$par, coef(...))`
  every iteration — reallocating and copying a growing matrix each time. This
  is why each saved fit is 35 MB.
- `MARSS():256-284` attaches `$kf`, `$Ey`, `$Innov`, `$Sigma`, `$J`, `$Kt`,
  `$J0`. Because `fun.kf` is `MARSSkfas`, populating those calls
  **`MARSSkfss()` up to twice** (lines 258 and 280) — and `MARSSkfss` is the
  *slow pure-R* filter, not the fast KFAS one. So this is two slow extra
  Kalman passes per chunk, not one.
- Full per-iteration error context (`iter=N ...`) instead of a one-line
  summary.

**Nothing downstream consumes any of it.** A grep for `$kf`, `$Ey`, `$Innov`,
`$Sigma`, `$Kt`, `$J0` and `iter.record` across all `.R` files returns zero
hits. The scripts use only `fit$logLik`, `fit$convergence`, `fit$numIter`,
`fit$pass`, and the reconstruction scripts read `model$model`, `model$par`,
`model$states` — all present at any `trace` level. Dropping to `trace = 0`
loses nothing that is actually used.

Keep `trace = 0` rather than `-1`: `-1` also disables `MARSSkemcheck()`, and
the error reporting it suppresses is worth more than the validation time on a
run this long. Use `-1` only if profiling shows the checks are material.

`fun.kf` is already the fast path — it defaults to `"MARSSkfas"` via
`match.arg`, so there is nothing to win there.

---

## Step 1 — Fix the control list (expected ~5–9×, zero new dependencies)

The single highest-value change. In each of these files, replace the
`controls` list:

- [EDA/seasonal_diffs_marss_model.R:217-224](EDA/seasonal_diffs_marss_model.R#L217-L224)
- [EDA/base_marss_model.R:228-234](EDA/base_marss_model.R#L228-L234)
- [parameters/model_seasonal_diffs.R:126-135](parameters/model_seasonal_diffs.R#L126-L135)
- [parameters/model.R:207-214](parameters/model.R#L207-L214)

```r
controls <- list(
  trace               = 0,      # was 1: no per-iteration rbind, no extra kf pass
  maxit               = 100,    # was 20: amortize per-call setup over more iterations
  minit               = 5,
  abstol              = 0.001,
  conv.test.slope.tol = 0.5,
  safe                = FALSE   # was TRUE: 1 Kalman pass per iteration, not 9
)

warmup_controls <- list(minit = 1, maxit = 2, safe = FALSE)
```

`abstol` and `conv.test.slope.tol` already match the MARSS defaults; leave
them.

**Risk and mitigation.** `safe = FALSE` is the MARSS default, but it exists
because an unlucky parameter update can make the Kalman filter fail or the
log-likelihood drop. That is a real possibility here — `R` is a dense 21×21
with `cov_lst`/`cov_eco` blocks. Mitigate by wrapping the chunk call rather
than reverting globally:

```r
fit_chunk <- function(y, mod, inits, controls) {
  out <- tryCatch(
    MARSS(y, mod, inits = inits, control = controls),
    error = function(e) NULL
  )
  if (is.null(out) || !is.finite(out$logLik)) {
    message("  chunk failed with safe=FALSE; retrying with safe=TRUE")
    out <- MARSS(y, mod, inits = inits,
                 control = modifyList(controls, list(safe = TRUE)))
  }
  out
}
```

This keeps the fast path as the default and pays the 9× cost only on the rare
chunk that needs it.

### Also in the chunk loop

At [EDA/seasonal_diffs_marss_model.R:248-258](EDA/seasonal_diffs_marss_model.R#L248-L258)
(and the equivalent block in every other script):

- Raising `maxit` to 100 means ~15 chunks instead of 75, so 60 fewer
  `MARSS()` setup/teardown cycles. Checkpoints get coarser; that is the
  trade, and it is a good one at 10 s/iteration.
- Move `plot(plot_loglik)` out of the loop — it appends to `Rplots.pdf` on
  every chunk.
- With `trace = 0` the saved fits shrink substantially on their own. If
  they are still large, save `coef(fit, type = "list")` plus `fit$logLik`
  for the resume path and keep the full object only on the final chunk.

---

## Step 2 — Evaluate `marssTMB` (the real algorithmic win, verify before adopting)

MARSS 3.11.9's `method` formal already accepts `"TMB"`, `"BFGS_TMB"` and
`"nlminb_TMB"`, and TMB 1.9.21 is installed — only the `marssTMB` glue package
is missing. It replaces ~1500 slow EM iterations with a compiled objective and
autodiff gradients driven by a quasi-Newton optimizer, which typically needs
one to two orders of magnitude fewer iterations.

This is a **spike, not a commitment.** The model here is not vanilla: a dense
`R` with shared `cov_lst`/`cov_eco` covariance blocks, a `D`/`d` covariate
term, free `B`, `U`, `A` and `x0`. `marssTMB` does not support every
constraint form base MARSS does, so the first question is whether it accepts
this model at all.

```r
install.packages("marssTMB", repos = "https://atsa-es.r-universe.dev")
```

Spike script (write to the scratchpad, not the repo). Reuse the existing model
construction by sourcing the setup portion of
[EDA/seasonal_diffs_marss_model.R](EDA/seasonal_diffs_marss_model.R) up to
line 213, which builds `y`, `marss_model_list`, `d` and `pars`:

1. **Does it accept the model?**
   `MARSS(y, marss_model_list, method = "TMB", control = list(maxit = 5))` —
   if this errors on the `R` structure or on `D`/`d`, stop; the answer is no
   and Step 1 is the whole win.
2. **Does it agree with EM?** Fit both from the same `inits` and compare
   `logLik`. They should converge to the same optimum. A materially different
   log-likelihood means one of them is not fitting the model you think.
3. **Is it actually faster?** `system.time()` both to a fixed `abstol`, not
   to a fixed iteration count — TMB iterations are not comparable to EM
   iterations.
4. **Count the free parameters** with
   `length(coef(MARSS(y, marss_model_list, fit = FALSE), type = "vector"))`.
   This also tells you whether plain `method = "BFGS"` is worth a look — see
   below.

Adopt only if all of 1–3 pass. Keep `method = "kem"` reachable as a fallback.

### What TMB does and does not change about scaling

Measured on `MARSSkfas` directly (T = 200, n = 30, diagonal `Q` and `B` as in
this spec), one Kalman pass:

| m | npar | sec |
|---|------|-----|
| 2 | 7 | 0.029 |
| 4 | 13 | 0.044 |
| 8 | 25 | 0.139 |
| 12 | 37 | 0.349 |
| 16 | 49 | 0.776 |

A second sweep over a wider range (T = 400, n = 40) pins the asymptote:

| m | 2 | 4 | 8 | 16 | 24 | 32 | 48 |
|---|---|---|---|----|----|----|----|
| sec | 0.059 | 0.087 | 0.252 | 1.550 | 5.182 | 12.418 | 41.233 |

Overall exponent 2.85 for m ≥ 8, and the *local* slopes converge on exactly
cubic: 2.62 (8→16), 2.98 (16→24), 3.04 (24→32), 2.96 (32→48). So the cost is
**O(T·m³)**, not quadratic — the O(m³) covariance products dominate once m is
large enough that the O(m²) and O(nm) terms stop masking them (which is why
the narrower m ≤ 16 sweep reads ≈2.5).
Scaling in **n** is only ≈1.16 (near linear), because KFAS uses univariate
sequential filtering and skips the O(n³) innovation-covariance inversion —
relevant if the site set ever grows from 21 rows toward the full 267.
Parameter count is linear in m (exponent 0.97) because `Q` and `B` are
diagonal here; a full `Q` would make it quadratic.

Note the m³ is an implementation artifact, not a necessity: `B` and `Q` are
diagonal here, so `P ← B P B' + Q` is mathematically O(m²), but KFAS treats
`T` as dense and calls general matrix multiply regardless. Verified —
declaring `B` diagonal vs unconstrained changes filter time by nothing
(ratio 0.99–1.32 across m = 8/16/24). Exploiting it would mean hand-writing
the recursion; not worth it at m = 7, but it means added states cost more
than the math requires.

The dense `R` is free *in the filter*. KFAS requires a diagonal observation
covariance and handles a dense one by LDL-decomposing it (verified: state
dimension unchanged, `H` diagonal afterwards) rather than by augmenting the
state vector, which would have pushed m from 7 to 7 + 21. Measured at
n = 20, m = 5: `R = "unconstrained"` (225 free params) and
`R = "diagonal and unequal"` (35 free params) both filter in 0.045 s. The
extra parameters still cost in the EM M-step and in any gradient, just not in
the Kalman recursion.

**TMB does not change these exponents** — it evaluates the same likelihood, and
whether marssTMB uses a C++ Kalman filter or a Laplace approximation on the
joint density, both are O(T·m³) worst case. What it changes is (a) the
constant factor, compiled vs interpreted; (b) gradient cost relative to
*parameter* count — reverse-mode AD is O(1) function evals rather than the
p+1 of finite differencing, which is a p-axis win, not an m-axis one; and (c)
iteration count, where EM's linear convergence rate degrades as latent
dimension grows relative to the data while quasi-Newton's superlinear rate
does not. Total EM time therefore scales worse in m than per-iteration cost
implies, and TMB removes that penalty — as a better constant, not a better
exponent.

At the current m = 7 / m = 2 this is mostly academic; the measured m = 2 → 8
gap is ~4.8×, which is the real basis for the "runs faster" note at
`parameters/model_seasonal_diffs.R:122-124`. The 9× from `safe = TRUE` is
larger and is m-independent.

### Tracking progress during a TMB fit

The chunk loop currently doubles as the progress indicator — each chunk prints
`logLik` and appends to `plot_loglik`. TMB converges in far fewer iterations,
so chunking becomes largely unnecessary, and progress tracking has to come
from the optimizer instead.

`MARSS:::.onLoad` sets the per-method control defaults, and the TMB entries are
verbatim `stats::nlminb` arguments:

```r
alldefaults[["TMB"]]$control <- list(
  maxit = 5000, tmb.silent = TRUE, eval.max = 5000, iter.max = 5000,
  trace = 0, abs.tol = NULL, rel.tol = NULL, x.tol = NULL,
  xf.tol = NULL, step.min = NULL, step.max = NULL, sing.tol = NULL)
alldefaults[["nlminb_TMB"]] <- alldefaults[["TMB"]]
alldefaults[["BFGS_TMB"]]   <- alldefaults[["BFGS"]]   # + tmb.silent = TRUE
```

So `method = "TMB"` / `"nlminb_TMB"` drive `nlminb`, and `"BFGS_TMB"` drives
`optim(method = "BFGS")`. Two levers:

| control | effect |
|---------|--------|
| `trace = n` | `nlminb` prints objective + parameters every *n*th iteration (`Rd`: "printed every trace'th iteration"). For `BFGS_TMB` this is `optim`'s trace, paired with `REPORT = n`. |
| `tmb.silent = FALSE` | passed to `TMB::MakeADFun()`; prints every objective/gradient evaluation — finer than per-iteration, and very verbose. |

**`trace` means something different here than under `kem`.** For EM it is
MARSS's internal error-checking and recording level (see above); for the
TMB/BFGS methods it is forwarded to the optimizer as a print frequency. The
`trace = 0` recommendation in Step 1 applies only to `method = "kem"`.

Use `trace = 10`, not `trace = 1`: `nlminb` prints the whole parameter vector
on every traced line, and this model has enough free parameters to make
per-iteration output unreadable.

To recover a plottable series equivalent to `plot_loglik`, capture and parse
stdout — the format is `iter:  objective:  par1 par2 ...`:

```r
log <- capture.output(
  fit <- MARSS(y, marss_model_list, method = "TMB",
               control = list(trace = 10, iter.max = 5000))
)
prog <- read.table(text = grep("^ *[0-9]+:", log, value = TRUE),
                   sep = ":", strip.white = TRUE)[, 1:2]
names(prog) <- c("iter", "negLogLik")
plot(prog$iter, -prog$negLogLik, type = "l", ylab = "logLik")
```

`nlminb` **minimizes**, so the traced objective is the negative
log-likelihood — flip the sign before comparing against `fit$logLik` or the
EM `plot_loglik` series.

Caveat for the spike: base MARSS defines the `MARSSfit` generic with only
`.kem`, `.BFGS` and `.default` methods — `MARSSfit.TMB` is registered by
`marssTMB` itself. Since that package is not installed I could not verify the
control forwarding end to end; the defaults above are strong evidence but
confirm it in spike step 1.

### On plain `method = "BFGS"` (no new dependency)

Available today, but likely **slower** here, not faster: base MARSS BFGS uses
numerically-differentiated gradients, so one gradient costs one likelihood
evaluation per free parameter. With a dense `R` and per-site `A`/`D` terms the
parameter count is high enough that this is a bad trade. The conventional
MARSS advice — run EM to get close, then BFGS from that solution to polish —
is still reasonable as a *final* step once EM has stopped making progress. Do
not use it as the primary optimizer without the parameter count from spike
step 4.

---

## Optional — Swap the BLAS to Apple Accelerate

This is the only sense in which "use multiple cores" applies to a single fit,
which is why it is documented here even though I expect little from it.

**What it is.** R ships with a reference BLAS/LAPACK — a plain, correct,
single-threaded implementation of the matrix primitives (`dgemm`, `dpotrf`,
…). Optimized replacements use blocking, SIMD, and multiple threads for the
same operations. On macOS, Apple's Accelerate framework is such a
replacement, and CRAN's R builds ship a prebuilt shim for it. Swapping is a
symlink: R loads whichever library `libRblas.dylib` points at, so no package
or code changes are involved.

```sh
cd "$(Rscript -e 'cat(R.home("lib"))')"
ls -l libRblas*                                    # see what is available
ln -sf libRblas.vecLib.dylib libRblas.dylib        # switch to Accelerate
# revert with:  ln -sf libRblas.0.dylib libRblas.dylib
```

Confirm with `sessionInfo()$BLAS`.

**Why I expect a modest gain at best.** The Kalman recursion is the hot loop,
and it runs 21×21 and 7×7 operations sequentially across ~4959 time steps.
At that size, per-call overhead dominates and threading is counterproductive.
Any benefit would come from better single-threaded kernels in the EM M-step,
where the vectorized updates for the dense `R` and the `D`/`d` covariate term
build genuinely larger intermediates. Call it 1.2–2×, mostly not from
parallelism.

**Caveats.** It modifies your R installation, not this project. Results can
differ in the last bits of floating-point precision, so a fit resumed across a
BLAS swap is not bit-reproducible. Benchmark before and after; revert if it
does not help.

---

## Verification

Benchmark on a short run rather than a full 75-chunk job. Do this from the
scratchpad so nothing in `models/` is disturbed:

1. **Baseline.** Source `EDA/seasonal_diffs_marss_model.R` through line 213,
   then time one 20-iteration chunk with the current `controls` —
   `safe = TRUE, trace = 1`. Record `system.time()`, `fit$logLik`,
   `fit$numIter`, and `object.size(fit)`.
2. **After Step 1.** Same, with `safe = FALSE, trace = 0`. Expect roughly
   5–9× less wall time for the same iteration count, and a much smaller fit
   object. **Confirm `logLik` after 20 iterations is essentially unchanged** —
   `safe` affects how the E-step is refreshed, not what is being optimized, so
   a materially different likelihood means something went wrong rather than
   faster.
3. **Convergence, not just speed.** Run ~5 chunks and check `plot_loglik` is
   still monotone increasing. A drop is the failure mode `safe = TRUE` guards
   against, and the signal to keep the `tryCatch` fallback.
4. **TMB spike.** Steps 1–4 of the spike above, comparing `logLik` at
   convergence against the EM fit from the same `inits`.
5. **Full run.** Only once 1–4 look right, relaunch
   `EDA/seasonal_diffs_marss_model.R` end to end and confirm chunk files
   appear in `models/marss_season_diffs/` at the improved cadence.

## Files to change

- [EDA/seasonal_diffs_marss_model.R](EDA/seasonal_diffs_marss_model.R) — `controls` (217-224), warmup call (228), chunk loop (248-258)
- [EDA/base_marss_model.R](EDA/base_marss_model.R) — same pattern at 228-234, 238, 257-267
- [parameters/model_seasonal_diffs.R](parameters/model_seasonal_diffs.R) — `controls`/`warmup_controls`/`chunks` (126-136)
- [parameters/model.R](parameters/model.R) — same block at 207-214, and the `# slower, more robust` comment at 213

The same `controls` block is duplicated in seven more scripts under
[R/](R/) (`mur_only_model.R`, `site_models_only.R`,
`site_only_macquaray_bay.R`, `sites_only_mur_fe.R`,
`model_exploratory_analysis.R`, `tir_only_model.R`,
`site_only_tir_slope.R`). Apply the same edit there only if you intend to
rerun them; they are not on the critical path for the current two runs.

## Out of scope, noted

- **No parallel driver.** Per your call, this plan targets single-fit latency.
  If you later want the ~9 script variants running concurrently, that is a
  separate change: `future`/`future.apply` are already installed, and the
  natural home is the empty
  [R/02_model_specification.R](R/02_model_specification.R), which currently
  contains only `library(MARSS)`.
- **`model_config` is dead code.** Nothing reads the `model_config` lists in
  `parameters/model*.R` — the EDA scripts hardcode the same values. So the
  `controls` edit must be made in *both* places to take effect, and the
  duplication is a standing hazard.
- **`SF0017` matches zero rows** in `outputs/.../tables/marss_row_metadata.csv`
  (likely a typo for `SF017`), so `n_site` is 5, not the 6 sites listed at
  line 14 of both scripts. Unrelated to speed, but it means the model you are
  waiting hours on may not be the one you intended.
