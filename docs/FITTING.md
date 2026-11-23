# Fitting

A fit finds the parameter values that bring the model closest to one or several test logs at once,
and then says how well those logs determine each value.

This file covers the fit file format, the hold-out and profile checks, how the search works, and
how to read the report. A fit minimises the score defined in [STUDIES.md](STUDIES.md), so read
the scoring section there first. `myfit.toml` below stands for a fit file of your own.

```
julia -t 8 --project bin/n2o-fill-sim.jl fit myfit.toml [--write-scenes]
julia -t 8 --project bin/n2o-fill-sim.jl fit myfit.toml --holdout coldflow      # leave one log out, then predict it
julia -t 8 --project bin/n2o-fill-sim.jl fit myfit.toml --profile VENT_ORIFICE.Cd [--profile-n 9]
julia --project bin/n2o-fill-sim.jl fit --init scenes/hf4.toml      # a draft fit file from the scene's [fit] tags
```

Two things about it matter more than the search itself.

Fit several logs together. A value fitted to one log quietly absorbs whatever that log's other
unknowns were. Fitting several logs with one set of shared values is what shows whether the model
can actually satisfy them all, and where it cannot.

Read the report's section "What these data do and do not determine". A fit will always return a number. Whether the data
determine that number is a separate question, and the report answers it: it names every pair of
parameters the data cannot tell apart, and every combination of parameters that could change
without changing any reading.

## Fit file

```toml
name = "wall-model"
explore = 80                    # optional: runs spread over the bounds first; the search starts from the best
seed = 1
max_evals = 600                 # evaluations for the search (default 60 per parameter)

[[dataset]]                     # one per log. The keys are those of a study (docs/STUDIES.md)
name = "hf4"
scene = "../scenes/hf4.toml"
replay = "../test/data/hf4_fill_10hz.csv"       # or manifest = "..."
t0 = 9020
until = 990
weight = 1.0                    # optional: how much this dataset counts against the others
set = { "TANK.wall_split" = 1 } # optional: fixed changes to this dataset's scene
[dataset.score]
channels = [{ log = "PT202_psig", model = "TANK.P" }]
phases = { stall = [140, 236], vent1 = [236, 310] }
[dataset.targets]               # optional: features scored beside the residuals (docs/STUDIES.md)
vent1_min = { fn = "min", log = "PT202_psig", model = "TANK.P", phase = "vent1", tol = 5 }

[[dataset]]
name = "coldflow"
scene = "../scenes/hf3cf.toml"
# ...

[[parameter]]                   # a parameter applies to every dataset (a shared value) ...
name = "VENT_ORIFICE.Cd"
from = 0.5
to = 1.2
start = 0.7                     # optional: default is the scene's present value

[[parameter]]
name = "settings.interphase_E"
from = 100
to = 2000
log = true                      # search in even ratios: for a quantity known only to an order of magnitude

[[parameter]]                   # ... unless it names the datasets it belongs to. The same name listed
name = "BOTTLE.fill_kg"         # again for another dataset is fitted separately there.
from = 10
to = 29
datasets = ["hf4"]
```

What is minimised is the same normalised score a study reports: every scored number (a channel's
RMS difference over a phase, or its mean difference over a time bin with `bin_s`), each divided
by how far that reading can be trusted (`sigma`), and every target over its tolerance, over all
datasets together (`weight` scales a whole dataset). See [STUDIES.md](STUDIES.md).

## Setup

A fit always returns numbers. Whether they mean anything is decided before it runs, by these
five choices.

1. **Fix what was measured.** A bottle temperature read from its pressure, a wall mass from a
   drawing, the weather on the day: put them in the scene or in `set` and do not fit them. Every
   value left free is one the fit can use to hide an error somewhere else.
2. **Count the scored numbers against the parameters.** A fit with about as many parameters as
   scored numbers can match anything and determines nothing. Two channels over three phases are
   six numbers. `bin_s` under `[dataset.score]` turns each long phase into several, and is the
   first thing to reach for when a fit has too little to go on.
3. **Choose which parameters are shared.** Physics and hardware that did not change between
   tests (heat-transfer coefficients, a vent's discharge coefficient, a wall's mass) should be
   shared. Things that differ from day to day and were not measured (how full the bottle was)
   should be per dataset. Getting that split wrong in either direction is the commonest way to
   fit nonsense.
4. **State what the model cannot do.** If the model is known to be wrong in part of a log, say
   so with a `discrepancy` on that phase ([STUDIES.md](STUDIES.md), "Stated model error").
   Without it the fit spends its parameters on an error they cannot remove, and returns biased
   values with a small spread, which is the worst combination.
5. **Keep a log back.** Fit with one log held out and see whether it is predicted better than
   before (next section). A fit that improves its own logs and worsens the one it did not see
   has fitted those logs, not the hardware.

## Hold-out

```
julia -t 8 --project bin/n2o-fill-sim.jl fit myfit.toml --holdout coldflow
```

This leaves the dataset called `coldflow` out of the fit entirely, fits the rest, and then runs
`coldflow` once with the fitted values. A fit can always be made to match the logs it was given.
Whether its values are properties of the hardware, or only of those logs, shows on a log it was
not given. Do it once for each dataset in turn: the scores of the held-out logs, side by side,
are the honest measure of the fit.

A parameter shared with the fitted datasets is fitted on them and used for the held-out one. A
parameter belonging only to the held-out dataset (`datasets = ["coldflow"]`) is not fitted: it
stays at the `start` given for it, or at the scene's value. So give such a parameter a measured
start value, or the prediction is being judged with a guess in it.

The report and the calibration record give each fitted dataset's own score, and for the held-out
one its score with the fitted values beside its score with every parameter at its start value.
If the second is the better of the two, the fit made the model worse at predicting.

## Profiles

```
julia -t 8 --project bin/n2o-fill-sim.jl fit myfit.toml --profile VENT_ORIFICE.Cd --profile-n 9
```

After the fit, the named parameter (by its label in the report, so `BOTTLE.fill_kg@hf4` for one
that belongs to a dataset) is held at each of 9 values evenly across its bounds, and at its fitted
value, and at each one all the other parameters are fitted again, starting from the best point.
The best score reached at each value is the parameter's profile, written to
`profile_<label>.csv` with the other parameters' refitted values.

- A profile that rises steeply either side of the fitted value: the data pin the parameter.
- A flat profile: whatever value this parameter is given, the others can make up for it. The
  data do not determine it. The rule used is stated in the report: the score rises by less than
  5 % of its lowest value across the whole range. Fix such a parameter at a measured or published
  value and take it out of the fit.
- A profile flat on one side only gives a bound, not a value.

The other columns of the file show which parameters did the making up: one that swings across
its range as the profiled one moves is its partner.

This costs a full fit per value, so use it on the parameters the report is doubtful about. It is
worth the cost because it assumes nothing: the spread and the correlations in the report come
from straight-line sensitivities at one point, which is a poor description of a long curved
valley. If the model fails at one of the values, no verdict is given.

## Search

Optionally a Latin hypercube over the bounds (`explore`) to find a good place to start. Then, by
default, the Levenberg-Marquardt method (`method = "lm"`): at each iteration it takes the
sensitivity of every residual to every parameter (one extra evaluation per parameter, differences
over 2 % of each range, well above solver noise), solves for a step, and tries four step lengths
at once. All of those evaluations are independent, so they run side by side: start Julia with as
many threads as there are cores (`-t 8`). Values are kept inside their bounds, and a parameter the
data push against a bound is held there.

`method = "nelder-mead"` is the simplex search instead: it uses no sensitivities at all, which
makes it the fallback if a model is too rough for them, but it needs ten or more times as many
evaluations and can use spare cores only to look one step ahead.

Either finds a nearby minimum, not necessarily the best one. `evals.csv` records every evaluation;
if the exploratory runs show several separate good regions, the fit is telling you the problem has
more than one answer.

### Failed runs

A parameter set for which a dataset cannot be simulated to the end (the
solver fails, or it takes more than `work_limit` part-flow evaluations) is not scored as a bad
fit: the search simply does not step there, and it is listed in `evals.csv` with the reason. The
limit counts work done, not seconds, so the result does not depend on how busy the machine is.
If the best points sit next to failed ones, the model's ability to run is deciding the answer,
not the data: look at why they fail.

## Report

`<fit file>_results/report.txt`, also printed, has these sections. With `--holdout` and `--profile`
it also has the held-out dataset's scores and the profile.

**Parameters.** The fitted values with an approximate spread, and a flag on any that ended on a
bound (which means the bound, not the data, chose the value).

**How close.** The RMS per dataset, channel and phase, in the log's own units; the offset taken
out of any channel whose zero is free; and each target's value in the model and in the log.

**What these data do and do not determine.** This comes from the sensitivity of every residual
to every parameter, taken at a tighter solver tolerance with steps of 1 % of each parameter's
range so that the differences are signal and not solver noise. It makes statements of these
kinds:

| Statement | Meaning |
|---|---|
| "A and B are not separately determined (correlation -0.97)" | The logs constrain some combination of the two, not each. A valve's bore and its discharge coefficient are the textbook case. Fix one from a measurement, or accept that only the combination is known. |
| "The data do not determine the combination ..." | Changing the parameters together in those proportions changes no reading at all. |
| "X has no measurable effect" | These logs say nothing about X. |

The collinearity index of the fitted set (Brun, Reichert and Künsch, 2001) is always printed.
Each parameter's sensitivities are scaled to unit length, so only the pattern counts (which
readings it moves, and in what proportion); the index is one over the square root of the
smallest eigenvalue of the scaled S'S. It reads as: a change in one parameter can be made up by
the others to within 1/index of its effect. 1 means the parameters act in entirely separate
ways. Above about 15 to 20 the set cannot be estimated together, and the report says so. A
correlation catches two parameters that cancel; the index also catches three or more that only
cancel as a group, where every pairwise correlation can look harmless. The cure is the same: fix
some, and fit fewer. (`collinearity(J)` on chosen columns asks the question of a subset.)

The spread treats each scored number as one independent observation, which is a rough stand-in
for the real number of independent readings in a time series (and too generous if bins are short:
see "Time bins" in [STUDIES.md](STUDIES.md)). Use it to compare parameters with each other, not as a
confidence interval.

## Output files

Everything goes in `<fit file>_results/`.

| File | Contents |
|---|---|
| `evals.csv` | Every evaluation of the search: parameter values, score, each dataset's score, failures. |
| `report.txt` | The report above. |
| `<date>-<name>.toml` | The calibration record (below). |
| `profile_<label>.csv` | With `--profile`: the profile, with the other parameters' refitted values. |
| `<dataset>.toml` | With `--write-scenes`: a copy of each fitted dataset's scene with the fitted values in place and a header saying which fit produced them. The scenes you gave are never modified. |

The calibration record holds the values, their bounds and spread, the datasets, channels,
phases, scored numbers and trust they were fitted against, the RMS reached, the commit of the
fit file's own folder (`commit`) and of the model that ran it (`model_commit`), the collinearity
index and the identifiability statements. With `--holdout` it has the held-out dataset's scores
under `[holdout]`, and with `--profile` the profile and its verdict under `[profile]`. A fitted
number without this record is a number nobody can check.
