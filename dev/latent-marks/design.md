# Latent marks: heritable unobserved lineage states in emphasis

Branch `latent-marks`. The design of an extension that lets a lineage carry a state the tree does
not show, inherited at speciation and changing along branches, which multiplies its rates beside
the endogenous covariates. Three processes are its targets: a latent rate (ClaDS), geography
(GeoSSE) and protracted speciation (PBD).

## The model

For lineage $s$ with mark $m_s(t) \in \{1, \dots, S\}$, under the exponential link,

$$\lambda_s(t) = e^{\eta^\lambda_s(t)}\, c^\lambda_{m_s(t)}, \qquad
  \mu_s(t) = e^{\eta^\mu_s(t)}\, c^\mu_{m_s(t)},$$

where $\eta$ is the linear predictor emphasis has now (intercept, $N$, $D$, $\ED$, $K$) and
$c^\lambda, c^\mu$ are positive multipliers per mark, with one reference mark fixed at 1. A mark
changes along a branch by a rate matrix $Q$ (anagenesis) and at a speciation by a kernel
$C(m \to m', m'')$, the joint law of the two daughters' marks given the parent's (cladogenesis).
The two crown lineages arise from a root mark drawn from $\pi$ by $C$. Given the complete tree the
marks are a Markov process on its branches, independent of everything except through the rates.

The three targets are parameterisations of $(c, Q, C, \pi)$:

| Process | Marks | $c$ | $Q$ | $C$ | Parameters |
|---|---|---|---|---|---|
| Latent rate (ClaDS2) | a grid of $S$ log-rate offsets $u_1 < \dots < u_S$ | $c^\lambda_m = c^\mu_m = e^{u_m}$ (constant turnover) | 0 | each daughter $u' = u + \mathcal N(\log\alpha, \sigma^2)$ independently, discretised on the grid | $\sigma$, $\alpha$ |
| Geography (GeoSSE) | A, B, AB | within-A, within-B and between-region speciation; extinction by region | range expansion A$\to$AB, B$\to$AB; extirpation AB$\to$A, AB$\to$B | AB splits into (A, AB), (B, AB) or (A, B) in proportion to those rates; A and B split into two copies | $s_A, s_B, s_{AB}, x_A, x_B, d_A, d_B$ (one rate fixed by $\eta$) |
| Protracted speciation (PBD) | good, incipient | initiation from good and from incipient; extinction of each | incipient $\to$ good at the completion rate | a split gives the parent's mark and an incipient daughter | $b_1/b_2$, $\mu_1/\mu_2$, $\lambda_c$ |

The general form (any $S$, any $Q$, any $C$) is the multitype model; the three are presets on it.

## Why it fits the inference emphasis has

The E-step draws hidden lineages $z$ and weighs each complete tree by $f(y,z\mid\theta)/g(z)$; the
M-step maximises the weighted $\log f$. With marks, $f(y,z\mid\theta)$ becomes the mark-marginal
density

$$f(y,z\mid\theta) = \sum_{\text{mark histories}} f(y,z,m\mid\theta),$$

and on a complete tree that sum is exact and cheap by pruning (Felsenstein's algorithm on the
lineages of `forest_t`): no mark is sampled, so the Monte Carlo error and the effective sample of
the E-step are not made worse by the marks themselves. The proposal $g$ is unchanged in form; it
reads mark-averaged rates, which affects only its efficiency.

**The pruning.** A lineage in emphasis continues through its own splits and starts a daughter at
each, so the recursion runs on lineages in decreasing order of birth (every daughter before its
parent). For lineage $k$, from its end backwards, a vector $v_k(m)$:

- at its end: $\mu$'s event factor $e^{\eta^\mu} c^\mu_m$ if it dies, 1 if it is alive at the
  present (sampling $\rho$ does not depend on the mark, unless tip marks are observed, below);
- across a stretch between two of its splits: with $Q = 0$, multiplied by
  $\exp(-c^\lambda_m I^\lambda - c^\mu_m I^\mu)$, where $I^\lambda, I^\mu$ are the stretch's
  integrals of $e^{\eta}$, the same `exp_integral` terms `loglik_ed` sums now, collected per
  lineage; with $Q \neq 0$, by the product over the stretch's segments of
  $\exp\{(Q - \mathrm{diag}(c^\lambda e^{\eta^\lambda} + c^\mu e^{\eta^\mu}))\,\Delta t\}$;
- at its split at $s$ with daughter $d$:
  $v_k(m) \leftarrow e^{\eta^\lambda(s)} c^\lambda_m \sum_{m',m''} \bar C(m\to m',m'')\, v_k(m')\, v_d(m'')$,
  with $\bar C$ the kernel symmetrised in the two daughters, since which daughter continues the
  parent is a convention of the lineage representation, not of the process.

At the crown, $\log f = \log \sum_m \pi_m \sum_{m',m''} \bar C(m\to m',m'')\, v_a(m')\, v_b(m'')$,
with every vector carried in log-scale. With one mark ($S = 1$) this is `loglik_ed` exactly.

**Cost.** With $Q = 0$ (the latent rate) the multipliers factor out of every stretch, so the
integrals are computed once per tree as now and the marks add $O(nS^2)$ per evaluation for $n$
lineages: with $S = 21$ and 300 lineages, about $10^5$ operations beside the $10^5$ of
`loglik_ed`. With $Q \neq 0$ each segment of each lineage needs a matrix exponential of size $S$;
for GeoSSE ($S = 3$) and PBD ($S = 2$) that is cheap.

## What changes in the package

1. **Likelihood** (`inst/include/model.hpp`): `loglik_marks`, beside `loglik_ed`, on the same
   table with the lineage index of every row added; dispatched when a model has marks.
2. **Simulator** (`general_tree.hpp`, `simulate_tree()`): marks drawn forward with the tree; the
   mark history is returned with the tree and is what the tests check against.
3. **Proposal** (`augment_tree.cpp`): rates averaged over the mark distribution at the current
   parameters; no change to the proposal's density code.
4. **Parameters** (`param_t`, the R side): the mark parameters appended after slot 11, so a
   12-element vector still means "no marks"; `model_bin` gains no flag, marks are a separate
   argument (`marks = latent_rate()`, `geography()`, `protracted()`, `multitype(S)`).
5. **Survival conditioning** (`cond = "surface"`): by direct simulation, which needs only (2).
6. **Observed tip marks** (geography): optional; an observed tip's end vector is the indicator of
   its mark (times $\rho$), which turns latent geography into trait-dependent diversification with
   endogenous covariates.

## Scope of the first version

- Exponential link only. Under the linear link $\max(0, \eta + \cdot)$ a multiplier does not
  factor out; the linear link refuses marks.
- $Q \neq 0$ together with the raw $\ED$ (whose $\eta$ has a slope within a segment) needs a
  time-ordered exponential; the first version refuses that combination. Centred $\ED$, $N$, $D$
  and $K$ have constant $\eta$ on a segment and combine with every mark model.
- Protracted speciation changes the observation map: an incipient lineage alive at the present is
  not a species of its own, so which tips the reconstructed tree shows depends on the marks. It is
  the last of the three and needs its own design of the observation step before any code.

## Order of work and checks

| Phase | What | Check before the next |
|---|---|---|
| 1 | `loglik_marks` for $Q = 0$, latent-rate preset, on complete trees | $S = 1$ equals `loglik_ed` to rounding; on a small complete tree, equals brute-force enumeration of all mark histories; the R entry point for tests |
| 2 | Simulator with marks; ML on complete trees | recovery of $\sigma$ on simulated complete trees; with $\sigma = 0$ the trees match `simulate_tree()` in distribution |
| 3 | Full MCEM on reconstructed trees, latent rate | effective sample at the generating values against the same trees without marks; recovery at 100 tips |
| 4 | $Q \neq 0$: multitype and GeoSSE, optional tip marks | brute force on small trees; GeoSSE with no endogenous covariate against `diversitree`'s GeoSSE likelihood fitted on the same reconstructed trees |
| 5 | Protracted speciation: the observation step | its own design note first |

Model choice between `nk` and a latent-rate model fitted on the same trees is what the extension
is for; its estimability is a study of its own, after phase 3.

## Expected difficulties

- Latent multitype models are weakly identified from reconstructed trees (Louca and Pennell 2020
  apply to them as to the endogenous models); the estimability study, not the code, decides what
  can be read.
- The grid for the latent rate: its span must cover the spread the random walk reaches after the
  largest $K$ in the tree ($\sigma\sqrt{K_{\max}}$), and its spacing must be fine against $\sigma$.
  The first version fixes $S$ and scales the grid with $\sigma$; the check against brute force
  and a finer grid measures the discretisation error.
- The proposal ignores the marks; where marks spread the rates widely the effective sample of the
  E-step may fall, as it already does for `nedc` and `ned` at turnover 0.40 (E46's gates).
