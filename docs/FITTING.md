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

What is minimised is the same normalised score a study reports: each channel's RMS difference
divided by how far that reading can be trusted (`sigma`), every channel and phase counted once,
over all datasets together (`weight` scales a whole dataset). See docs/STUDIES.md.

Shared and per-dataset parameters are the point of the format. Physics and hardware that did not
change between tests (heat-transfer coefficients, a vent's discharge coefficient, a wall's mass)
should be shared. Things that differ from day to day and were not measured (how full the bottle
was, the air temperature) should be per dataset. Getting that split wrong in either direction is
the commonest way to fit nonsense.

## How it searches

Optionally a Latin hypercube over the bounds (`explore`), then the Nelder-Mead simplex method from
the best point, run twice (the second pass restarts from where the first stopped, with a smaller
simplex). It uses no derivatives, which suits a model whose output carries solver noise. Values
are kept inside their bounds.

Start Julia with as many threads as there are cores (`-t 8`). Each evaluation runs its datasets
side by side; the exploratory runs and the sensitivity runs go several at a time; and when there
are more threads than datasets the search looks ahead, evaluating together the four points its
next step could need. That takes the same path to the same answer as doing them one at a time.

Nelder-Mead finds a nearby minimum, not necessarily the best one, and slows down beyond a dozen or
so parameters. `evals.csv` records every evaluation; if the exploratory runs show several separate
good regions, the fit is telling you the problem has more than one answer.

## The report

`<fit file>_results/report.txt`, also printed:

- **Parameters**: the fitted values with an approximate spread, and a flag on any that ended on a
  bound (which means the bound, not the data, chose the value).
- **How close**: the RMS per dataset, channel and phase, in the log's own units.
- **What these data do and do not determine**: from the sensitivity of every residual to every
  parameter, taken at a tighter solver tolerance with steps of 1 % of each parameter's range so
  that the differences are signal and not solver noise.
  - "A and B are not separately determined (correlation -0.97)": the logs constrain some
    combination of the two, not each. A valve's bore and its discharge coefficient are the
    textbook case. Fix one from a measurement, or accept that only the combination is known.
  - "The data do not determine the combination ...": changing the parameters together in those
    proportions changes no reading at all.
  - "X has no measurable effect": these logs say nothing about X.

The spread treats each channel-and-phase as one independent observation, which is a rough stand-in
for the real number of independent readings in a time series. Use it to compare parameters with
each other, not as a confidence interval.

## What it writes

- `evals.csv`: every evaluation (parameter values, score, each dataset's score, failures).
- `report.txt`.
- `<date>-<name>.toml`: the **calibration record**. The values, their bounds and spread, the
  datasets, channels, phases and trust they were fitted against, the RMS reached, the commit, and
  the identifiability statements. A fitted number without this is a number nobody can check.
- with `--write-scenes`, `<dataset>.toml`: a copy of each dataset's scene with the fitted values
  in place and a header saying which fit produced them. The scenes you gave are never modified.
