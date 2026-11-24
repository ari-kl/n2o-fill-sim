# Studies

A study runs one scene many times with parameters varied, scores each run against a test log if
there is one, and reduces each run to a row of numbers. It is how to ask "which vent coefficient
fits this log best", "how sensitive is the loaded mass to the bottle temperature", or, with no log
at all, "how long does the tank stay above 800 psig for each of these fill plans".

This file covers the study file format, how runs are scored, control policies (a controller
inside the runs), screening for the parameters that matter, and prediction bands (the spread of a prediction over every parameter set a log
could not rule out). The scoring section also applies to fits ([FITTING.md](FITTING.md)), which minimise
the same score.

```
julia -t 8 --project bin/n2o-fill-sim.jl study studies/hf4_vent.toml
```

Runs are independent and use as many threads as Julia was started with (`-t`). Results are the
same to the last digit however many threads there are. Each finished run is appended at once to
`<study file>_results/results.csv`; `--resume` skips the runs already there, and `--keep-traces`
also writes each run's scored and measured signals to `traces/<n>.csv`.

## Study file

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
                                              # "morris": screening, n trajectories (see below)
                                              # "sets": the rows of a file (see "Prediction bands")
# sets = "kept.csv"                           #   that file, for mode = "sets"

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

[bands]                                       # optional: the spread of signals over the runs, to bands.csv
signals = ["TANK.P", "TANK.m"]                #   (see "Prediction bands")
```

Parameter names are the same `PART.param` and `settings.key` that `--set` takes, and are checked
before anything runs, as are metric functions and log columns.

## Scoring

The model is interpolated onto the log's sample times; times with no log reading are
skipped. Each channel and phase gives an RMS difference in the log's own units (psi, °C, kg), in a
column named `rms.<log column>.<phase>`. Those are the numbers to read.

`score` is the one number to sort by, and it has no units. It is made of scored numbers, each
a difference divided by how far that reading can be trusted, its `sigma`. Without `bin_s` there is
one scored number per channel and phase, weight · rms / sigma, and the score is their root mean
square:

    score = sqrt( mean over channels and phases of (weight · rms / sigma)² )

So a score of 1 means the model is, on average, as far from the log as the readings themselves
can be trusted, and 3 means three times that. Two things follow. A pressure and a temperature can
be scored together, because 5 psi and 1 K each count as one. And every phase counts once, however
long it is, so a ten-second vent weighs as much as a twenty-minute hold.

Every scored number is also a column, `res.<log column>.<phase>` (and `res.target.<name>`), so a
row can be checked number by number (for example, to throw out every run in which any one of them
is over 3).

With a manifest and no `[score] channels`, every reading the manifest maps to a model signal is
scored, over the manifest's phases.

### Sigma and weight

`sigma` comes from the channel in `[score]` if given, else from the manifest's channel, else from
the kind of reading: 5 psi for a pressure, 1 K for a temperature, 0.25 kg for a mass. Set it
from what you know of the sensor (two transducers that disagree by 5 psi cannot be trusted to
better than that). `weight` (default 1) is a further factor for when one channel or dataset should
simply matter more.

### Time bins

One number for a twenty-minute hold says how far off the model was on average, and
nothing about whether it drifted the right way. With `bin_s`, each phase is cut into equal bins as
near to that many seconds as fit it, and each bin is one scored number:

    weight · (mean of model − log over the bin) / sigma

It keeps its sign, so a fitter can see that the model is high early and low late. The columns are
`res.<log column>.<phase>.<bin number>`; the per-phase `rms` columns stay as they are. The score is
still the root mean square of the scored numbers, so a phase now counts once per bin: a long phase
weighs more than a short one. A bin with no log reading is left out.

Make bins 30 s or longer for a slow process such as a hold or a fill. A fitter treats every
scored number as a separate piece of information. Neighbouring samples of a slowly changing
pressure are not: if the model is 10 psi high now it will be about 10 psi high two seconds from
now, whatever the parameters. Short bins therefore count the same information many times, and
the fit comes out more confident than the data justify. A bin should be about as long as it takes
for the model's error to change of its own accord.

### Unknown offsets

A load cell that was not tared, or was tared with plumbing
hanging on it, reads the right changes around the wrong zero. `offset = "free"` on a channel says
so. For each run, the mean of model − log over all of that channel's scored samples is taken out
first, so only changes are scored: how much went in during the fill, how much left in a vent. The
amount taken out is the column `offset.<log column>` (what to add to the log to lay it on the
model). Read it: if it differs widely between runs that fit equally well, the changes alone are
not telling those runs apart. (This is not the `offset` of a manifest channel, which is a known
number added when the log is read. See [LOGS.md](LOGS.md).)

`offset = "per_phase"` is for a reading whose zero also jumps: a load cell that steps when a
valve slams or a line stiffens. The mean is then taken out of each phase by itself (columns
`offset.<log column>.<phase>`), so put the phase boundaries at the events where it steps. The
price is that only changes inside a phase are scored. What the reading did from
one phase to the next, which is where a step hides but also where real mass moved, tells the
score nothing. A phase scored whole then says only how well the shape inside it matches; with
`bin_s` each bin says how far it sits from the phase's mean, and a phase of one bin says nothing
at all. If the offsets found for neighbouring phases differ by more than the reading could have
stepped, the model moved a different mass between them than the test did.

### Sensor readings

A channel or a target can name an instrument's reading as its model signal: `model =
"PT202.P_read"`, `"PT202.T_read"`, or a load cell's `"LC.m_read"`. The reading carries the
sensor parameters set on that Gauge or LoadCell in the scene (offset, gain, lag, range, the hose
step; [MODEL.md](MODEL.md) section 8), which can be set per run with `set`, swept, or fitted like
any other parameter. A gauge with no sensor parameter set records no reading; in `[score]` and
`[targets]` its `P_read` and `T_read` then mean its true `P` and `T`, so a fit may start a
sensor parameter from its default. Metrics, `[bands]` and `--columns` take recorded signals
only.

Which tool for an instrument that reads off:

- A fixed offset or scale that is known (from a calibration, or from a stretch of the log where
  the truth is known) belongs in the manifest (`gain`, `offset`: the log is corrected as it is
  read) or on the gauge (`offset_psi`, `offset_kg`: the model reports what the instrument would).
  Use one, not both. The manifest is the place when the number belongs to the log file; the
  gauge when it should show in the window and the exported run, or when it is to be fitted.
- An offset that is not known belongs in the score: `offset = "free"`, or `"per_phase"` if it
  also jumps. Nothing is fitted and the channel then only scores changes. Fitting `offset_psi`
  instead spends a parameter to learn the same number and keeps the absolute level in the
  score, which is only worth it if something else in the log pins that level down.
- A lag, a range limit or a hose step is not an offset and no score setting stands in for it:
  put it on the instrument and score against the reading.

A lagged reading depends a little on `sample_dt`; see the note on the lag in MODEL.md section 8.

### Stated model error

Sometimes the model is known to be wrong in one stretch of a test: it
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

## Targets

A point-by-point difference asks the model to be right everywhere at once. Often
what matters, and what a parameter actually controls, is a feature: how low the pressure went in
a vent, how fast it recovered afterwards, what it was at the end of a hold. A target names one
such feature with a metric function (below), takes it from the log column and from the model
signal over the same window, and adds one scored number:

    (feature of the model − feature of the log) / tol

`tol` is how close counts as agreement, in the feature's own units: psi for a minimum, psi per
second for a slope, seconds for a crossing time. The window is `from` and `to` (seconds or
events, as for a phase) or `phase = "<name>"`; without one it is the whole run. Both sides are in
the log column's units (psi gauge, °C, kg), so `lo`, `hi` and `value` are given plainly in
those and there is no `unit`. The model's signal is read at the log's own sample times, so both
sides are worked out from the same instants; times with no log reading are skipped on both. The
columns are `target.<name>.model`, `target.<name>.log` and the scored number `res.target.<name>`.
A feature the log does not have (a level it never reaches) is an error when the file is read. One
the model does not have in a given run counts as 100 tolerances out.

A slope or a change does not depend on a reading's zero, so a target is also a way to use a load
cell that cannot be trusted in absolute terms.

## Metrics

A metric is a number taken from each run. The functions are `final`, `max`, `min`, `mean_over`,
`slope` (per second, by least squares), `time_in_band` (with `lo`, `hi`) and `first_crossing`
(with `value`). With `from` and `to` in
seconds, only that window of the run is looked at (so `final` with `to = 300` is the value at
300 s); without them, the whole run. Signals are the names in an exported run (`TANK.P`, `TANK.m`,
`AV203.mdot`, `vented.m`, ...) and are in SI units. An argument can be given as a gauge pressure
or in Celsius by its name (`lo_psig`, `value_C`), and `unit = "psig"` or `"C"` converts a result
that is a pressure or a temperature (a slope stays in SI per second: `unit` with `slope` is an
error).

Every row also carries how far the run got (`t_end`), why it failed if it did (`failed`), how
many junction solves did not settle, the mass and energy balance errors, the mass vented and the
wall time. A run that cannot be built or does not finish is still a row.

## Control policies

A study's runs normally follow a schedule fixed in advance: the scene's own, a replayed log's and
each case's actions. A `[controller]` table adds a controller that reads the instruments at a
fixed tick and commands valves and heaters as the run goes on, the way a vent controller or an
operator would. It is how to ask "how many times does the vent cycle to hold this band", or "how
much is loaded by the time the tank has settled", where the answer depends on what was commanded
and the commands depend on what the run does.

```toml
scene = "../scenes/hf4.toml"
until = 900

[[case]]
name = "fill"
actions = [{ t = 2, part = "AV201", action = "open" }, { t = 2, part = "AV203", action = "open" }]

[controller]
dt = 0.2                                      # the tick, seconds (default 0.2)
latency_s = 0.0                               # a reading is what the instrument showed this long before the tick
actuation_delay_s = 0.0                       # a command reaches its part this long after the tick that sent it
stand_ins = []                                # true signals a policy may read as if they were instruments
noise = { "PT202.P_read" = { sigma_psi = 2 } } # optional: noise added to a reading (sigma in SI, sigma_psi or sigma_K)

[[controller.policy]]
name = "vent"
kind = "hysteresis"
reading = "PT202.P_read"
actuator = "AV204"
open_above_psig = 800
close_below_psig = 760
min_open_s = 1                                # it stays open at least this long once opened
min_shut_s = 10                               #   and shut at least this long once shut

[vary]
"controller.vent.close_below_psig" = [740, 760, 780]
```

Without a `[controller]` a study runs exactly as it always has: none of what follows is used.

### Ticks

The controller ticks at t = 0, dt, 2 dt, and so on, to the end of the run. At each tick the run is
brought to the tick, its readings are taken, and the policies act in the order they are listed.
The solver is not stopped at the ticks: it runs ahead as far as the physics allows and the state
at the tick is read off its interpolant, so a tick at which nothing is commanded costs a sample,
not a restart. A command goes into the run's schedule like any other action and is applied
`actuation_delay_s` after its tick. A policy sends a command only when what it wants changes, so
a valve held open is not told to open again at every tick.

A sample is recorded at every tick, so a run with a controller is recorded at least that densely,
whatever `sample_dt` says.

### Readings and stand-ins

A policy sees only the readings it declares, taken as an instrument would report them: a Gauge's
`P_read` or `T_read`, a LoadCell's `m_read` (with the sensor parameters set on that part; MODEL.md
section 8). A gauge with no sensor parameter set reads the truth, so its `P_read` is its `P`, as
in scoring. Each reading is the recorded signal at `latency_s` before the tick, plus the noise
given for it, if any.

Sometimes a policy has to act on something no instrument on the stand measures, such as the true
mass in the tank. Name that signal in `stand_ins` and a policy may read it. The study then says so:
every row has a `stand_ins` column listing them, so a result that leaned on a quantity nobody
could have measured is marked as such. A true signal not in `stand_ins` is an error when the study
is read, as is a stand-in nothing reads.

The noise on a reading is drawn from the same small generator a seeded study uses, seeded from
the run's parameters other than the controller's and its case. Two runs that differ only in
controller settings therefore see the same noise, so a difference between them is the
controller's doing and not the luck of the draw, and the same run gives the same noise on any
number of threads and with `--resume`.

### Units

Settings that are levels (`open_above`, `close_below`, `above`, `lo`, `hi`) are given in SI, as an absolute
pressure in Pa or a temperature in K, or by their name with `_psig` or `_C` for a gauge pressure or
a temperature in Celsius, as metric arguments are. Times are in seconds. A noise sigma is a
difference, so it takes `sigma` (SI), `sigma_psi` or `sigma_K`, never a gauge zero. The scene has
no display units of its own (those are a window setting), so a study file says its units by name.

### Policies

| `kind` | What it does | Settings |
|---|---|---|
| `hysteresis` | opens a valve when a reading goes above one level and shuts it when the reading falls below a lower one, each time only once the valve has been open (or shut) for its minimum | `reading`, `actuator`, `open_above`, `close_below`, `min_open_s`, `min_shut_s` |
| `pulse_wait` | when a reading is above a level, opens a valve for `pulse_s`, shuts it, and waits `wait_s` before it may open again | `reading`, `actuator`, `above`, `pulse_s`, `wait_s` |
| `band_watch` | marks the moment a reading has stayed inside a band for `hold_s` (see "Marks") | `reading`, `lo`, `hi`, `hold_s`, `mark` |
| `heater_hold` | switches a vessel's heater on (`on = 1`, the default) or off. While on, the heater runs as its own `heater_ctl` and `heater_set` say | `heater`, `on` |
| `julia` | a policy written in Julia (below) | `file`, `reads`, `actuators`, `settings` |

Every policy has a `name` (used in `[vary]` and in its signals) and may start disabled
(`enabled = false`). A hysteresis or pulse-and-wait policy starts from where its valve is when the run starts: shut, or
open for a normally-open valve, which a run starts with de-energised. A
heater's set-point is a parameter of the vessel (`TANK.heater_set`), so it is varied or set per
run like any other.

### Varying controller settings

`[vary]` takes controller settings by name: `controller.<policy>.<setting>`, for example
`"controller.vent.open_above_psig" = [790, 800, 810]`, and `controller.dt`,
`controller.latency_s` and `controller.actuation_delay_s`. Every run's controller is built and
checked when the study is read.

### Signals

Each policy adds signals to the run's history, held from one tick to the next: `ctl.<name>.enabled`
(1 or 0), `ctl.<name>.cmd` (its last command: 1 for open or on, 0 for close or off) and, for a
policy with a reading, `ctl.<name>.obs`, the reading it acted on, noise and all. They can be
metrics, `[bands]` signals and `--keep-traces` columns like any other.

### Sequences

Most procedures are a series of steps: get ready, fill, hold until ready, fire. `[[controller.step]]`
lists them in order. The sequence starts with the first step at t = 0, and each step ends when its
`duration_s` has passed, when its `until` condition is met (and then `operator_delay_s` later,
for a person who has to see it and act), or at its `timeout_s`, whichever is first. The next step
starts at the same tick.

```toml
[[controller.step]]
name = "fill"
on_enter = [{ part = "AV201", action = "open" }]     # commands sent on entering the step
on_exit = [{ part = "AV201", action = "close" }]     #   and on leaving it
until = { reading = "LC1.m_read", ge = 9.0 }
operator_delay_s = 2
timeout_s = 600                                     # it ends here whatever; `fill.timed_out` says so
enable = ["vent"]                                   # policies switched on and off on entering it
disable = []

[[controller.step]]
name = "hold"
duration_s = 600
enable = ["ready"]
```

A step's commands are logged and delayed like a policy's. A policy switched on by a step that
starts at a tick acts at that tick, with that step's settings. Once the last step has ended the
policies keep its settings.

A setting can differ from step to step: `open_above_psig = { fill = 820, hold = 800 }`, with
`default = ...` for the steps not named. Every step must have a value, and every step named must
exist. Varied, it is `controller.vent.open_above_psig.hold`.

A run can end at a step's edge rather than at a time: `until = { step = "hold", edge = "end",
max_s = 3600 }` at the top of the study file ends each run when its hold ends, and gives up on it
at 3600 s (a run stopped by `max_s` is a failed row).

### Conditions

A step's `until` is a condition, written as a table:

| Condition | True when |
|---|---|
| `{ reading = "PT202.P_read", ge_psig = 800 }` | the reading compares so (`gt`, `ge`, `lt`, `le`; SI, or `_psig`, `_C`) |
| `{ rate = "LC1.m_read", window_s = 30, lt = 0.002 }` | its rate of change, from the earliest tick inside the window to now, compares so (SI per second, or `_psi` for psi per second) |
| `{ elapsed_ge = 60 }` | the step has run this many seconds |
| `{ held = { ... }, for_s = 30 }` | the condition inside has been true at every tick for this long |
| `{ mark = "ready", plus_s = 0 }` | the mark has happened (and this long ago) |
| `{ any = [ ... ] }`, `{ all = [ ... ] }` | any or all of the conditions in the list |

The readings are held to the same rule as a policy's: an instrument's reading, or a stand-in. A
rate and a hold keep what they have seen since their step began. A condition is looked at only
while its step runs.

### Marks

A `band_watch` policy marks the moment a reading has stayed inside a band for a while:

```toml
[[controller.policy]]
name = "ready"
kind = "band_watch"
enabled = false                                      # switched on by the hold step
reading = "PT202.P_read"
lo_psig = 780
hi_psig = 820
hold_s = 30                                          # in the band (ends included) at every tick for 30 s
mark = "ready"                                       # the default is the policy's name
```

It counts only while it is enabled, starts again from zero whenever the reading leaves the band,
and marks once. A mark can end a step (`{ mark = "ready" }`), place a metric's window, align the
bands, and key a case's actions.

### Result columns

A study with a controller adds columns to `results.csv`: `stand_ins`, `escape_hatch`, for each
step `t_<step>.start`, `t_<step>.end` and `<step>.timed_out` (1 if its timeout ended it), and for
each mark `t_<mark>`, all in seconds from the start of the run. A step that never started or a mark
that never happened leaves the cell empty, and that alone does not fail the run.

A metric's window can be given by events instead of seconds: `from = { mark = "ready" }`,
`to = { mark = "ready", plus_s = 300 }`, `from = { step = "fill", edge = "end" }`. A metric whose
window needs an event that did not happen, or ends after the run did, is empty: the value a
signal held at the end of the run is not presented as its value at a later time.

```toml
[metrics]
m_fill_end = { fn = "final", signal = "TANK.m", to = { step = "fill", edge = "end" } }
m_ready_300 = { fn = "final", signal = "TANK.m", to = { mark = "ready", plus_s = 300 } }
```

`[bands] align = "ready"` puts each run's grid on the time since that mark: `t_s` in `bands.csv` is
then seconds after it, so runs that got ready at different times are compared at the same point
in their procedure. A run that never made the mark is left out of the bands, and the number left
out is printed.

The controller's own signals include `ctl.step` (the number of the step running, 0 before the
first and after the last) and `ctl.mark.<name>` (0, then 1 from the mark on).

### Cases with a controller

A case can be a whole scenario: `set` changes parameters for its runs only, as `[set]` does for
the study, and `group` is a label carried into a `group` column for sorting and filtering (cases
are still run, and banded, one by one).

```toml
[[case]]
name = "hot_day"
group = "summer"
set = { "settings.T_amb_C" = 32, "BOTTLE.T_C" = 32, "controller.heat.on" = 0 }
actions = [
    { t = 300, part = "AV204", action = "stuck_closed" },                    # at a time, as always
    { mark = "ready", plus_s = 60, part = "AV204", action = "stuck_closed" }, # a minute after the mark
    { step = "hold", edge = "start", policy = "vent", action = "disable" },  # on a policy
    { t = 400, reading = "PT202.P_read", fault = "freeze" },                 # on a reading
]
```

A case's `set` is applied before the run's own varied parameters, and a parameter cannot be both
set by a case and varied. Its controller settings (`controller.<policy>.<setting>`) count as the
controller's: they do not change the noise. Its other parameters do, like a varied one.

An action with a time and a part goes into the schedule as before. The controller carries out the
rest at the first tick at or after it falls due: one keyed to a mark or a step's edge (with
`plus_s` after it), one on a policy (`enable`, `disable`), and a fault on a reading:

- `freeze`: the reading stays at what it was at that tick, whatever the truth does, and the
  policies go on acting on it, as they would on a reading that has stopped updating.
- `offset`: `value` (SI), `value_psi` or `value_K` is added to it.
- `drop`: the reading is missing. A missing reading is NaN, and every comparison with it is false,
  so a policy on it sends nothing and a condition on it is not met: absence is never taken for a
  value.
- `clear`: the reading is restored.

A policy or reading named in a case must exist, and a keyed action needs a `[controller]`.

### Metrics for control

These come beside the metrics above, in `[metrics]`, and take the same `from` and `to` (seconds,
or events as above).

Regulation, on a `signal` with a band `lo`, `hi` (SI, or `_psig`, `_C`):

| `fn` | |
|---|---|
| `time_in_band_frac` | the share of the window spent in the band |
| `max_above`, `max_below` | how far the signal went above `hi` or below `lo` at most (0 if never); in Pa or K, or psi with `unit = "psi"` |
| `time_to_band` | seconds from the start of the window to first being in the band |
| `band_exits` | how many times it left the band (counted on the recorded samples) |
| `peak` | the same as `max` |
| `cold_exposure` | kelvin-seconds below `below` (or `below_C`): how long, and how far, a temperature sat under it |

Margins and mass, on a part:

| `fn` | |
|---|---|
| `vented_kg` | the mass through `part`, flowing a to b, over the window: one vent path |
| `relief_kg`, `relief_lifts` | the same through a relief valve `part`, and how many times it began to flow |
| `prv_margin_psi` | the smallest margin between a relief valve's set pressure and `signal`, in psi (negative: it lifted) |
| `liquid_frac_max` | the highest liquid fraction of `vessel` |
| `m_at` | `signal` at `at` (seconds, or an event); empty if the run did not get there |

Switching, from the controller's own log of what it commanded (`actuator`, a part): what was sent,
not where the valve was.

| `fn` | |
|---|---|
| `cmds`, `cycles` | commands sent to it; those that opened it (or switched it on) |
| `cmds_per_hour`, `min_interval_s` | their rate over the window; the shortest time between two |
| `max_cmds_in` | the most commands in any `window_s` seconds |
| `recharge_violations` | commands less than `recharge_s` after the one before (a pneumatic actuator's supply that has not recovered) |
| `n2_used` | `per_cmd` for each command, in whatever unit `per_cmd` is |
| `open_s`, `duty` | time commanded open, and its share of the window |
| `longest_open_s`, `longest_shut_s` | the longest stretch commanded open, and shut |

Procedure:

| `fn` | |
|---|---|
| `step_s` | how long `step` ran |
| `t_mark` | when `mark` happened, in seconds after `from` (default the start of the run) |
| `stalled` | 1 if `signal` was below `below` when `step` ended (it stopped short), else 0 |
| `m_at_fire` | `signal` where the run ended on its step's edge (`until = { step, edge }`) |

These are numbers to compare candidates by. They are not features of a log, so targets do not
take the ones on parts or on the controller.

### A policy written in Julia

When none of the policies above will do, `kind = "julia"` names a file whose last expression is a
function of `(obs, state)`:

```julia
(obs, state) -> begin
    x = obs.readings["PT202.P_read"]
    want = x > obs.settings["above"]
    want == get(state, "open", false) && return []
    state["open"] = want
    [("AV204", want ? "open" : "close")]
end
```

`obs` has `t`, `step`, `readings` (only those listed in the policy's `reads`) and `settings` (the
policy's `settings` table); `state` is a Dict of its own, kept from tick to tick. It returns the
commands to send, each `(part, action)` or `(part, action, value)`, on parts listed in its
`actuators`. The same checks apply as for any policy: a reading must be an instrument's or a
stand-in. The file is code and runs as such, so a study with one runs only with
`--allow-julia-policy` on the command line, and its rows have `escape_hatch` = 1.

A controller outside the simulator, talking to it over a line protocol, is not supported.

### The HF4 example

`studies/hf4_control.toml` runs the HF4 stand as a procedure: fill to 9 kg (on the true mass, a
stand-in, since the stand has no load cell), with a vent that opens above 800 psig and shuts below
770 or 785 psig, then a 600 s hold that starts a 30 s band watch (760 to 810 psig). A second case
freezes the tank gauge's reading at the mark, so the vent goes on acting on the value it had then.

```
julia -t 2 --project bin/n2o-fill-sim.jl study studies/hf4_control.toml
```

| case | shut below | fill_s | vent cycles in hold | in band after the mark | PRV margin, psi | m_fire, kg |
|---|---:|---:|---:|---:|---:|---:|
| nominal | 770 | 120.0 | 0 | 1.00 | 201.8 | 9.01 |
| nominal | 785 | 127.4 | 1 | 1.00 | 200.8 | 8.96 |
| frozen_gauge | 770 | 120.4 | 0 | 1.00 | 199.4 | 9.01 |
| frozen_gauge | 785 | 127.2 | 0 | 1.00 | 199.4 | 9.00 |

The vent opens during the fill and keeps the back pressure down, so the fill is quick; the narrower
band (shut below 785) costs a few seconds and opens the vent once more in the hold. With the gauge
frozen the vent never opens in the hold, and the pressure creeps up by a couple of psi more than it
would have: a 600 s hold in this weather is too short for the frozen reading to matter, and a
longer or sunnier hold is where it would. Nothing in the frozen reading shows either way.
(Each case draws its own noise, so the two fills differ slightly before the fault.)

## Screening

Before fitting ten parameters it is worth knowing which of them the scored numbers respond to at
all. `mode = "morris"` answers that cheaply, by the method of elementary effects (Morris, 1991).

```toml
mode = "morris"
n = 20                                        # trajectories (default 20); seed optional
[vary]
"VENT_ORIFICE.Cd" = { from = 0.5, to = 1.0 }  # the range to screen over; log = true for even ratios
"settings.interphase_E" = { from = 100, to = 2000, log = true }
```

Each parameter's range is cut into four levels (0, 1/3, 2/3 and all of the way across). A
trajectory starts at a random point on that grid and then moves one parameter at a time, in a
random order, by two thirds of its range. With k parameters that is k + 1 runs, and each step is
one elementary effect: the change in an output, divided by the step. So `n = 20` costs
20 × (k + 1) runs and gives every parameter 20 effects, each taken with the other parameters
somewhere different. A seed gives the same runs on any machine. A morris study takes one `[[case]]`
at most.

Besides `results.csv` (one row per run, ids beginning `m<trajectory>.<step>;`), the study writes
`screening.csv`: for every varied parameter and every output (the score, each scored number
`res.*`, each target's model value and each metric),

- `mu_star`: the mean of the absolute effects. How much the output moves when the parameter
  crosses its whole range. Because it is per whole range, parameters can be compared with each
  other, and a wider range gives a larger number: the answer is about the ranges you gave.
- `sigma`: the spread of the effects. Small beside `mu_star` means the parameter acts the same
  way wherever the others are. Large means its effect depends on the others, or bends.
- `n`: how many effects went into it. A step with a run at either end that failed, or stopped
  partway, is left out.

The ranking by `mu_star` on the score is printed. Read the table by output as well: a parameter
with a small effect on the score can still be the only one that moves a particular scored number,
and that number is then what pins it. A parameter whose `mu_star` is small on every scored number
cannot be found from this log, whatever a fit reports: fix it at a measured or published value.
Screening says what matters, not what the values are, and not which parameters can be told apart
from each other (two can both matter and still only act together: the "What these data do and do not determine" section of the fit
report, described in [FITTING.md](FITTING.md), is for that).

## Prediction bands

A fit returns one value for each parameter. Often the log does not determine that value: the
parameter can sit anywhere in a wide range, with the others shifted to suit, and the model still
matches. One run with the fitted values then shows a single line where the honest answer is a
spread. A prediction band is that spread: the model run once for every parameter set the log could
not rule out, and each signal's percentiles over those runs at every time.

It takes two studies and one command between them.

1. Sample. Run a study with `mode = "lhs"` over the ranges the parameters could plausibly lie
   in, scored against the log, with any model error you know of stated as a `discrepancy`.

2. Keep the runs the log could not rule out.

```
julia --project bin/n2o-fill-sim.jl study keep studies/hf4_sample_results/results.csv [--cutoff 3] [--nth 2]
```

Every `res.*` column is a difference between model and log divided by how far it is trusted,
sqrt(sigma² + stated discrepancy²). A run is ruled out when its `--nth` largest |res| is over
`--cutoff`: by default, when two or more of its scored numbers are more than 3 of those out. One is let through because, among many scored numbers, one will stray without
the parameters being wrong; with dozens of scored numbers (time bins), raise `--nth` in proportion.
A run that failed or stopped early is ruled out. The command prints how many runs were kept and,
for each parameter, the range left among the kept runs beside the range sampled, and writes
`kept.csv` beside the results: one row per kept run, one column per parameter.

3. Run the kept sets on the scenario to predict.

```toml
scene = "../scenes/hf4.toml"                  # need not be the scene or the log the sets came from
mode = "sets"
sets = "hf4_sample_results/kept.csv"          # relative to the study file

[bands]
signals = ["TANK.P", "TANK.m"]                # the names in an exported run
percentiles = [5, 50, 95]                     # the default
grid_s = 1.0                                  # the default: one row per second
```

`mode = "sets"` takes its parameter sets from the file instead of making them: a header of
`PART.param` and `settings.key` names, then one set per row, each crossed with every `[[case]]` as
usual. An empty cell leaves that parameter as the scene has it. The scene and the schedule are
this study's own, which is the point: sets kept on old logs, run on a fill that has not happened.
A column that is not a parameter of this scene is an error that names it. Run ids begin `s<row>;`.

`[bands]` works with any mode. The study then also writes

- `bands.csv`: a row per grid time, with `t_s`, `n` and, for every signal and percentile, a column
  such as `TANK.P.p5`, in SI units like an exported run. Each run's recorded history is
  interpolated onto the grid. `n` is how many runs reached that time: a run that failed or stopped
  early counts up to where it stopped and not after, so a falling `n` is worth a look at
  `results.csv`.
- `metric_bands.csv`: a row per metric, with `case`, `metric`, `n` and a column per percentile
  (`p5`, `p50`, `p95`), over the runs that finished.

With several cases there is a `bands_<case>.csv` for each, because cases are different scenarios
and their runs are not pooled. A study with band signals cannot be resumed (`--resume`): the runs
already done left no history to take a band from.

4. Load the bands in the window. In run mode, "Load bands..." (or `--bands bands.csv` on the command line) draws each
band on the plot of its kind: shaded between its lowest and highest percentile, with the median as
a thin line, behind the live lines and in the colour of its signal's line.

### Interpretation

A band is the spread of the model's output over the parameter sets that the stated tolerances
(each reading's `sigma` and each phase's `discrepancy`) could not rule out, among those that were
sampled. Read it as: given this model, and given how far I said the log can be trusted, the data
do not tell these runs apart, and this is how far they disagree about the signal.

It is not a confidence interval. Three things it does not contain:

- It does not include model error outside the stated allowances. Every run in the band is the same model. What the
  model leaves out moves all of them together, and the band does not widen to show it. A real test
  can land outside a narrow band for that reason alone.
- It does not include anything that was not varied. A parameter held fixed contributes no width, however poorly it is known.
- It does not give probabilities. The percentiles are over the kept sets, which were spread evenly over the
  sampled ranges and then cut by a pass/fail rule. Wider ranges, a different cutoff or a different
  `sigma` give a different band. "5-95 %" says where nine tenths of the kept runs lie, not that
  the truth is inside with nine chances in ten.

A band is most useful when it is wide: it then shows, in the units of the plot, what the logs have
not pinned down, and which signal a new measurement would have to be taken on to pin it.

### Example

`studies/hf4_sample.toml` spreads 400 sets over four numbers of the HF4 scene (the vent's discharge
coefficient, the boiling enhancement, the inner heat-transfer multiplier and the bottle's fill),
each scored against the logged tank pressure in four phases. `studies/hf4_bands.toml` runs the
kept sets on the same fill.

```
julia -t 8 --project bin/n2o-fill-sim.jl study studies/hf4_sample.toml
julia --project bin/n2o-fill-sim.jl study keep studies/hf4_sample_results/results.csv
julia -t 8 --project bin/n2o-fill-sim.jl study studies/hf4_bands.toml
julia --project bin/n2o-fill-sim.jl scenes/hf4.toml --replay test/data/hf4_fill_10hz.csv --t0 9020 --run --bands studies/hf4_bands_results/bands.csv
```

The second command prints (the first takes about 20 s on eight threads):

```
400 runs, 235 kept (ruled out: failed, or |res| over 3.0 in 2 or more of the 4 scored numbers)
  parameter                         kept                      sampled
  BOTTLE.fill_kg                    15.05 to 28.74            15.0 to 28.99
  VENT_ORIFICE.Cd                   0.5016 to 0.9976          0.5005 to 0.9998
  settings.h_in_mult                0.1032 to 0.998           0.102 to 0.998
  settings.interphase_E             100.9 to 1974.0           100.5 to 1994.0
```

Two runs in five are ruled out, and yet every parameter keeps its whole range: it is
combinations that the log excludes, and no one of these four numbers can be read off it. The bands
show what that leaves open. Tank pressure in psig and tank mass in kg, 5th to 95th percentile
(`hf4_sample.toml` has `[bands]` too, so the spread over everything sampled can be set beside it):

| t [s] | logged psig | all 400 sampled: psig | 235 kept: psig | all 400 sampled: kg | 235 kept: kg |
|------:|------------:|----------------------:|---------------:|--------------------:|-------------:|
| 200 | 800.9 | 786.9 to 818.7 (31.8) | 788.6 to 813.4 (24.7) | 7.83 to 8.85 (1.01) | 7.91 to 8.85 (0.94) |
| 500 | 788.3 | 763.4 to 807.1 (43.6) | 767.5 to 797.8 (30.4) | 8.37 to 10.33 (1.96) | 8.57 to 10.33 (1.75) |
| 900 | 849.2 | 823.8 to 864.4 (40.6) | 826.3 to 857.2 (30.9) | 8.37 to 10.33 (1.96) | 8.57 to 10.33 (1.75) |

Scoring against the pressure narrows the pressure band by a fifth to a third, to 25 to 31 psi (a
run passes with an RMS miss of up to 15 psi in all but one phase: 3 times a `sigma` of 5 psi). It
barely narrows the mass band: the loaded mass is 8.57 to 10.33 kg with a median of 9.96
(`metric_bands.csv`), against 8.37 to 10.33 kg before any run was thrown out. This pressure log, scored this way, says almost nothing about how much
nitrous is in the tank. That is what a single fitted run cannot show, and the reason to draw the
band.

## Results

`results.csv` opens in any spreadsheet. Two things worth doing with it:

- Sort by `score`, then look at the per-phase columns: a parameter that improves one phase and
  spoils another is telling you about the model, not about its value.
- Plot `score` against two varied parameters. A long valley means the log cannot tell those two
  apart: any point along it fits equally well.

The traces from `--keep-traces` load into the window through "Load test log..." to lay runs over
each other and over the log.
