# Test logs and manifests

The simulator reads a test log for two things: to replay its valve and heater sequence, and to
draw its readings next to the model's. There are two ways to tell it what the columns are.

**By guessing** (`--replay log.csv`, or "Load test log..." with a CSV). The first column is time.
Columns are recognised by name: `[PSI]` or `_psig` is a gauge pressure, `[kg]` or `_kg` a mass,
`[C]` or `_C` a temperature, and 0/1 columns are matched to valves and heaters in the scene
(`valve_AV-201` and `solenoid_AV201` both match a part called `AV201`). `valve_` columns are taken
as 1 = open; `solenoid_` and `relay_` columns are never replayed until you say what a 1 means.
This is enough for a single-board log with tidy names.

**By a manifest** (`--manifest log.toml`, or "Load test log..." with the manifest). A manifest
states everything, so nothing is guessed and one test's conventions do not end up in code.

```toml
file = "raw_log.csv"            # relative to the manifest
time = "device_timestamp"       # the time column (default: the first)
t0 = 21950.0                    # log time that is simulation time zero (default: the first reading)
window = [21940, 23170]         # optional: read only this span of log time
sample_dt = 1.0                 # the result is put on an even grid of this spacing (default 0.1 s)
debounce_s = 0.5                # a valve or heater change counts only if it then lasts this long

# Only for logs written by several boards, each stamping rows with its own clock:
source_column = "source"        # the column naming the board that wrote the row
reference_source = "BOARD-A"    # every other board's clock is shifted onto this one
clock_offsets = "auto"          # or given: { "BOARD-B" = 100.18, "BOARD-C" = 1869.08 }

[[channel]]
column = "PT202 [PSI]"
role = "pressure"               # pressure | mass | temperature | valve | heater
maps_to = "TANK.P"              # the recorded signal to compare it with (used when scoring a run)
name = "PT202_psig"             # optional: what to call it (default: the column name)
gain = 1.0                      # reading = gain * value + offset, which must come out as
offset = 0.0                    #   psi gauge, kg or °C
sigma = 5.0                     # optional: how far this reading can be trusted, in those units
                                #   (used when scoring; default 5 psi, 1 K, 0.25 kg by role)

[[channel]]
column = "solenoid_AV201"
role = "valve"
maps_to = "AV201"               # the part it drives
one_means = "closed"            # "open" or "closed": required, never assumed

[[channel]]
column = "TH2 [C]"
role = "temperature"
group = "TH_tank"               # channels in a group are averaged into one reading
maps_to = "TANK.T"

[phases]                        # optional named spans, used when scoring
fill = { from = "AV201 open", to = "AV201 close" }
hold = { from = "AV204 close", to = "AV201 open#2" }
soak = [170, 1210]              # or plain seconds after t0
```

Rules worth knowing:

- **Clocks.** With `clock_offsets = "auto"`, each board's offset is estimated from the order the
  rows arrived in: for every row of that board, the last time seen on a reference-board row minus
  the row's own time, and then the median of those. It is good to about one reference sampling
  period. If the offset over the first third of the file differs from the last third by more than
  0.5 s, a warning says that board's clock drifts.
- **The grid.** Each grid point takes the median of the readings within half a `sample_dt` either
  side (for valves and heaters, the last state in that interval). **An interval with no reading is
  left empty.** A gap is never filled in, so a silent sensor cannot pass for a steady one.
- **Bounce.** With `debounce_s`, a change of a 0/1 column is kept only if the new value then holds
  for that long; a kept change is therefore reported up to that much late.
- **Phase events.** `"AV201 open"` is the first time the replayed schedule opens AV201, and
  `"AV201 open#2"` the second. A part's state at the start of the log is not an event. A `to`
  without `#n` is the first such event after `from`.
- **Mistakes are errors.** An unknown key, a missing column, a valve without `one_means` or a
  phase event that never happens stops the load with a message; nothing is skipped quietly.

`n2o-fill-sim manifest suggest log.csv scene.toml` prints a draft from the name-based guesses, with
`"?"` wherever it will not guess. It is a starting point to read and correct.

`n2o-fill-sim manifest extract log.toml out.csv` writes the log as the manifest describes it: one
clock, the even grid, readings under their manifest names, and valve columns as 1 = open whatever
the original meant. That file is small, needs no manifest to replay, and is how a raw log of
hundreds of megabytes becomes something to keep next to a scene.

## Saturation checks

In a tank holding liquid, pressure and temperature are tied together. The log dialog's "also draw
what saturation implies" overlays, for each logged pressure, the saturation temperature it implies
(on the temperature plot), and for each logged temperature the saturation pressure (on the pressure
plot). A wall thermometer that follows the first is following the liquid. A pressure well above
the second points at a warm surface layer or a sensor fault.
