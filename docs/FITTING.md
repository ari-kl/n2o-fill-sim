# Fitting

A fit finds the parameter values that bring the model closest to one or several test logs at once,
and then says how well those logs determine each value.

```
julia -t 8 --project bin/n2o-fill-sim.jl fit myfit.toml [--write-scenes]
julia --project bin/n2o-fill-sim.jl fit --init scenes/hf4.toml      # a draft fit file from the scene's [fit] tags
```

Two things about it matter more than the search itself.

**Fit several logs together.** A value fitted to one log quietly absorbs whatever that log's other
unknowns were. Fitting several logs with one set of shared values is what shows whether the model
can actually satisfy them all, and where it cannot.

**Read the last section of the report.** A fit will always return a number. Whether the data
determine that number is a separate question, and the report answers it: it names every pair of
parameters the data cannot tell apart, and every combination of parameters that could change
without changing any reading.

## The file

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
datasets together (`weight` scales a
whole dataset). See docs/STUDIES.md.

**Count the scored numbers against the parameters.** A fit with about as many parameters as
scored numbers can match anything and determines nothing. Two channels over three phases are six
numbers. `bin_s` under `[dataset.score]` turns each long phase into several, and is the first
thing to reach for when a fit has too little to go on.

Shared and per-dataset parameters are the point of the format. Physics and hardware that did not
change between tests (heat-transfer coefficients, a vent's discharge coefficient, a wall's mass)
should be shared. Things that differ from day to day and were not measured (how full the bottle
was, the air temperature) should be per dataset. Getting that split wrong in either direction is
the commonest way to fit nonsense.

**State what the model cannot do.** If the model is known to be wrong in part of a log, say so
with a `discrepancy` on that phase (docs/STUDIES.md, "Stated model error"). Without it the fit
spends its parameters on an error they cannot remove, and returns biased values with a small
spread, which is the worst combination. With it, the values that the rest of the log determines
come out as they should.

## How it searches

Optionally a Latin hypercube over the bounds (`explore`) to find a good place to start. Then, by
default, the **Levenberg-Marquardt** method (`method = "lm"`): at each iteration it takes the
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

**Runs that fail.** A parameter set for which a dataset cannot be simulated to the end (the
solver fails, or it takes more than `work_limit` part-flow evaluations) is not scored as a bad
fit: the search simply does not step there, and it is listed in `evals.csv` with the reason. The
limit counts work done, not seconds, so the result does not depend on how busy the machine is.
If the best points sit next to failed ones, the model's ability to run is deciding the answer,
not the data: look at why they fail.

## The report

`<fit file>_results/report.txt`, also printed:

- **Parameters**: the fitted values with an approximate spread, and a flag on any that ended on a
  bound (which means the bound, not the data, chose the value).
- **How close**: the RMS per dataset, channel and phase, in the log's own units; the offset taken
  out of any channel whose zero is free; and each target's value in the model and in the log.
- **What these data do and do not determine**: from the sensitivity of every residual to every
  parameter, taken at a tighter solver tolerance with steps of 1 % of each parameter's range so
  that the differences are signal and not solver noise.
  - "A and B are not separately determined (correlation -0.97)": the logs constrain some
    combination of the two, not each. A valve's bore and its discharge coefficient are the
    textbook case. Fix one from a measurement, or accept that only the combination is known.
  - "The data do not determine the combination ...": changing the parameters together in those
    proportions changes no reading at all.
  - "X has no measurable effect": these logs say nothing about X.

The spread treats each scored number as one independent observation, which is a rough stand-in
for the real number of independent readings in a time series (and too generous if bins are short:
see "Time bins" in docs/STUDIES.md). Use it to compare parameters with each other, not as a
confidence interval.

## What it writes

- `evals.csv`: every evaluation (parameter values, score, each dataset's score, failures).
- `report.txt`.
- `<date>-<name>.toml`: the **calibration record**. The values, their bounds and spread, the
  datasets, channels, phases, scored numbers and trust they were fitted against, the RMS reached, the commit of the fit file's own
  folder (`commit`) and of the model that ran it (`model_commit`), and the identifiability
  statements. A fitted number without this is a number nobody can check.
- with `--write-scenes`, `<dataset>.toml`: a copy of each dataset's scene with the fitted values
  in place and a header saying which fit produced them. The scenes you gave are never modified.
