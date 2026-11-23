# Studies

A study runs one scene many times with parameters varied, scores each run against a test log if
there is one, and reduces each run to a row of numbers. It is how to ask "which vent coefficient
fits this log best", "how sensitive is the loaded mass to the bottle temperature", or, with no log
at all, "how long does the tank stay above 800 psig for each of these fill plans".

```
julia -t 8 --project bin/n2o-fill-sim.jl study studies/hf4_vent.toml
```

Runs are independent and use as many threads as Julia was started with (`-t`). Results are the
same to the last digit however many threads there are. Each finished run is appended at once to
`<study file>_results/results.csv`; `--resume` skips the runs already there, and `--keep-traces`
also writes each run's scored and measured signals to `traces/<n>.csv`.

## The file

```toml
scene  = "../scenes/hf4.toml"                 # relative to the study file
replay = "../test/data/hf4_fill_10hz.csv"     # optional: a log to replay and score against,
# manifest = "log.toml"                       #   or a manifest (docs/LOGS.md)
t0 = 9020                                     # log time that is simulation time zero
until = 990                                   # seconds to run (default: the end of the log)
sample_dt = 0.5                               # shortest interval between recorded samples
reltol = 1e-5                                 # solver tolerance; 1e-4 is fine for a coarse first pass
work_limit = 5000000                          # give up on a run after this many part-flow evaluations
                                              #   (the HF4 replay takes 0.3 million); the same on any machine
wall_limit = 600                              # and, as a last resort, after this many seconds
mode = "grid"                                 # "grid": every combination
                                              # "oat": the scene as it is, then one parameter at a time
                                              # "lhs": n runs spread over the ranges (needs n; seed optional)

[set]                                         # optional: changes applied to every run
"settings.T_amb_C" = 28

[vary]
"VENT_ORIFICE.Cd" = { from = 0.6, to = 1.0, n = 5 }      # evenly spaced; log = true for even ratios
"settings.interphase_E" = [300, 700]                      # or a list

[[case]]                                      # optional variants, each crossed with every parameter set
name = "vent_stuck"
actions = [{ t = 250.0, part = "AV204", action = "stuck_open" }]

[score]                                       # needs a log
channels = [{ log = "PT202_psig", model = "TANK.P" },
            { log = "TH_tank_C", model = "TANK.T", sigma = 1.0, weight = 1.0 },   # sigma, weight: optional
            { log = "LC_kg", model = "TANK.m", offset = "free" }]                 # a reading whose zero is unknown
                                                                                  #   ("per_phase": and jumps between phases)
phases = { stall = [140, 236], vent1 = [236, 310] }      # seconds after t0, or
                                                          # { from = "AV201 open", to = "AV201 close" }
bin_s = 30                                    # optional: one scored number per 30 s of each phase
                                              #   (a phase table can carry its own bin_s)
# A phase written as a table can also state how far the model itself is known to be off there:
# phases = { fill_start = { from = 140, to = 160, discrepancy = 40 }, stall = [160, 236] }

[targets]                                     # optional: features of the log scored against the same
                                              #   features of the model, each with a tolerance
vent1_min = { fn = "min", log = "PT202_psig", model = "TANK.P", from = 236, to = 310, tol = 5 }
recovery = { fn = "slope", log = "PT202_psig", model = "TANK.P", from = "AV204 close", to = 400, tol = 0.05 }
hold_end = { fn = "final", log = "PT202_psig", model = "TANK.P", phase = "stall", tol = 5 }

[metrics]                                     # optional numbers taken from each run
loaded_kg = { fn = "final", signal = "TANK.m" }
peak_psig = { fn = "max", signal = "TANK.P", unit = "psig" }
seconds_above_800_psig = { fn = "time_in_band", signal = "TANK.P", lo_psig = 800, hi_psig = 2000 }
```

Parameter names are the same `PART.param` and `settings.key` that `--set` takes, and are checked
before anything runs, as are metric functions and log columns.

**Scoring.** The model is interpolated onto the log's sample times; times with no log reading are
skipped. Each channel and phase gives an RMS difference in the log's own units (psi, °C, kg), in a
column named `rms.<log column>.<phase>`. Those are the numbers to read.

`score` is the one number to sort by, and it has no units. It is made of **scored numbers**, each
a difference divided by how far that reading can be trusted, its `sigma`. Without `bin_s` there is
one scored number per channel and phase, weight · rms / sigma, and the score is their root mean
square:

    score = sqrt( mean over channels and phases of (weight · rms / sigma)² )

So a score of 1 means the model is, on average, as far from the log as the readings themselves
can be trusted, and 3 means three times that. Two things follow. A pressure and a temperature can
be scored together, because 5 psi and 1 K each count as one. And every phase counts once, however
long it is, so a ten-second vent weighs as much as a twenty-minute hold.

Every scored number is also a column, `res.<log column>.<phase>` (and `res.target.<name>`), so a row can be checked number
by number (for example, to throw out every run in which any one of them is over 3).

**Time bins.** One number for a twenty-minute hold says how far off the model was on average, and
nothing about whether it drifted the right way. With `bin_s`, each phase is cut into equal bins as
near to that many seconds as fit it, and each bin is one scored number:

    weight · (mean of model − log over the bin) / sigma

It keeps its sign, so a fitter can see that the model is high early and low late. The columns are
`res.<log column>.<phase>.<bin number>`; the per-phase `rms` columns stay as they are. The score is
still the root mean square of the scored numbers, so a phase now counts once per bin: a long phase
weighs more than a short one. A bin with no log reading is left out.

**A reading whose zero is unknown.** A load cell that was not tared, or was tared with plumbing
hanging on it, reads the right changes around the wrong zero. `offset = "free"` on a channel says
so. For each run, the mean of model − log over all of that channel's scored samples is taken out
first, so only changes are scored: how much went in during the fill, how much left in a vent. The
amount taken out is the column `offset.<log column>` (what to add to the log to lay it on the
model). Read it: if it differs widely between runs that fit equally well, the changes alone are
not telling those runs apart. (This is not the `offset` of a manifest channel, which is a known
number added when the log is read. See docs/LOGS.md.)

`offset = "per_phase"` is for a reading whose zero also jumps: a load cell that steps when a
valve slams or a line stiffens. The mean is then taken out of each phase by itself (columns
`offset.<log column>.<phase>`), so put the phase boundaries at the events where it steps. The
price is stated plainly: **only changes inside a phase are scored**. What the reading did from
one phase to the next, which is where a step hides but also where real mass moved, tells the
score nothing. A phase scored whole then says only how well the shape inside it matches; with
`bin_s` each bin says how far it sits from the phase's mean, and a phase of one bin says nothing
at all. If the offsets found for neighbouring phases differ by more than the reading could have
stepped, the model moved a different mass between them than the test did.

**Stated model error.** Sometimes the model is known to be wrong in one stretch of a test: it
lacks a mechanism that matters for the first seconds of a fill, say. Scored as it stands, that
stretch dominates the score, and a fitter will bend every parameter it has to make the error
smaller, spoiling values that the rest of the log pins down well. Leaving the stretch out hides
the problem instead of recording it.

`discrepancy` on a phase says, in the log's units, how far off the model is expected to be there
whatever its parameters. It is added to `sigma` in quadrature for that phase's scored numbers:

    trusted to = sqrt( sigma² + discrepancy² )

So `discrepancy = 40` on a pressure with a sigma of 5 psi means a 40 psi miss in that phase counts
about as one, not as eight. Give one number (`discrepancy = 40`), or one per log column when
several channels are scored (`discrepancy = { PT202_psig = 40, TH_tank_C = 3 }`; columns not
named get none). The phase must be written as a table, with `from` and `to` as times or events.
The `rms` columns are not changed by it: they still show the miss in full. A stated error is a
written admission, with a size, of what the model does not do. It belongs in the file, where the
next reader sees it.

Make bins **30 s or longer** for a slow process such as a hold or a fill. A fitter treats every
scored number as a separate piece of information. Neighbouring samples of a slowly changing
pressure are not: if the model is 10 psi high now it will be about 10 psi high two seconds from
now, whatever the parameters. Short bins therefore count the same information many times, and
the fit comes out more confident than the data justify. A bin should be about as long as it takes
for the model's error to change of its own accord.

`sigma` comes from the channel in `[score]` if given, else from the manifest's channel, else from
the kind of reading: **5 psi for a pressure, 1 K for a temperature, 0.25 kg for a mass**. Set it
from what you know of the sensor (two transducers that disagree by 5 psi cannot be trusted to
better than that). `weight` (default 1) is a further factor for when one channel or dataset should
simply matter more.

With a manifest and no `[score] channels`, every reading the manifest maps to a model signal is
scored, over the manifest's phases.

**Targets.** A point-by-point difference asks the model to be right everywhere at once. Often
what matters, and what a parameter actually controls, is a feature: how low the pressure went in
a vent, how fast it recovered afterwards, what it was at the end of a hold. A target names one
such feature with a metric function (below), takes it from the log column and from the model
signal over the same window, and adds one scored number:

    (feature of the model − feature of the log) / tol

`tol` is how close counts as agreement, in the feature's own units: psi for a minimum, psi per
second for a slope, seconds for a crossing time. The window is `from` and `to` (seconds or
events, as for a phase) or `phase = "<name>"`; without one it is the whole run. Both sides are in
the **log column's units** (psi gauge, °C, kg), so `lo`, `hi` and `value` are given plainly in
those and there is no `unit`. The model's signal is read at the log's own sample times, so both
sides are worked out from the same instants; times with no log reading are skipped on both. The
columns are `target.<name>.model`, `target.<name>.log` and the scored number `res.target.<name>`.
A feature the log does not have (a level it never reaches) is an error when the file is read. One
the model does not have in a given run counts as 100 tolerances out.

A slope or a change does not depend on a reading's zero, so a target is also a way to use a load
cell that cannot be trusted in absolute terms.

**Metrics.** `final`, `max`, `min`, `mean_over`, `slope` (per second, by least squares),
`time_in_band` (with `lo`, `hi`) and `first_crossing` (with `value`). With `from` and `to` in
seconds, only that window of the run is looked at (so `final` with `to = 300` is the value at
300 s); without them, the whole run. Signals are the names in an exported run (`TANK.P`, `TANK.m`,
`AV203.mdot`, `vented.m`, ...) and are in SI units. An argument can be given as a gauge pressure
or in Celsius by its name (`lo_psig`, `value_C`), and `unit = "psig"` or `"C"` converts a result
that is a pressure or a temperature (not a slope: that stays in SI per second).

**Every row** also carries how far the run got (`t_end`), why it failed if it did (`failed`), how
many junction solves did not settle, the mass and energy balance errors, the mass vented and the
wall time. A run that cannot be built or does not finish is still a row.

## Reading the result

`results.csv` opens in any spreadsheet. Two things worth doing with it:

- Sort by `score`, then look at the per-phase columns: a parameter that improves one phase and
  spoils another is telling you about the model, not about its value.
- Plot `score` against two varied parameters. A long valley means the log cannot tell those two
  apart: any point along it fits equally well.

The traces from `--keep-traces` load into the window through "Load test log..." to lay runs over
each other and over the log.
