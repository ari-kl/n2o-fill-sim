# n2o-fill-sim

> [!WARNING]
> **Experimental and unvalidated.** This is a research and learning tool, not an engineering
> one. Its results are model predictions that have been compared with two fill logs from one
> test stand and with no measured temperature at all. Do not use it to set fill quantities,
> pressures, relief settings, procedures or operating limits, or as the basis for any safety
> decision. Nitrous oxide is hazardous; work with it under your organisation's safety rules
> and qualified supervision.

A simulator for filling a hybrid rocket's nitrous oxide run tank. You draw the plumbing (bottle,
hoses, valves, orifices, regulator, tank, vent), press Simulate, and click valves while it runs.
It shows pressure, temperature, mass and flow everywhere, and how the tank behaves when you vent,
heat, wait in the sun, or leave a valve stuck. Think of it as LTspice for a nitrous fill.

![A recorded fill replayed: the plumbing coloured by pressure while the tank vents, with the logged pressures and load cell overlaid on the plots](docs/screenshot.png)

It is an analysis tool. It does not talk to any control system or hardware, and its numbers are
model predictions: see [How far to trust it](#how-far-to-trust-it).

The same model runs without the window from the command line, many times over with parameters
varied (a study), or inside a fitter that adjusts parameters to match test logs. Four logs come
with it, as 10 Hz or 1 Hz extracts in `test/data`. "HF3" and "HF4" are two static hot-fire tests
of a student hybrid motor (June and July 2026), and `scenes/hf4.toml` is that test stand's fill
plumbing. The other two are from the same hardware: an HF3 cold flow with nitrous
(`scenes/hf3cf.toml`), and a fill and sealed hold from a launch campaign, "LC26 day 4", which has
thermistors on the tank (`scenes/lc26_day4.toml`).

## Installing and running

Julia 1.13 is needed. It runs natively on Windows and Linux, and on WSL through WSLg.

```
git clone https://github.com/ari-kl/n2o-fill-sim
cd n2o-fill-sim
julia --project -e "using Pkg; Pkg.instantiate()"
```

On Windows, install Julia first (`winget install julia -s msstore`, then `juliaup add 1.13.4`).
After the lines above, double-click `n2o-fill-sim.cmd`, or from a terminal run
`n2o-fill-sim.cmd scenes\hf4.toml`.

On Linux or WSL, install [juliaup](https://github.com/JuliaLang/juliaup) first. Then:

```
julia --project bin/n2o-fill-sim.jl scenes/hf4.toml             # open the window on a scene
julia --project bin/n2o-fill-sim.jl scenes/hf4.toml --replay test/data/hf4_fill_10hz.csv --t0 9020 --run
julia --project bin/n2o-fill-sim.jl --help
```

The window appears about 3 s after launch. The first launch after installing or changing the code
compiles first, which can take a minute.

Nothing is fetched at run time. To use it somewhere without a network, first run
`Pkg.instantiate()` for the project and for `test/`, start the program once, and set
`JULIA_PKG_OFFLINE=true`. `Manifest.toml` is committed so every machine gets the same versions.

## The window

In the editor, pick a part from the palette and click to place it, then click one pin and another
to wire them. Select a part to edit its parameters; hovering a parameter, or a part on the
schematic, says what it means. The Settings tab has ambient conditions, an optional target
pressure window and the model's tuning coefficients. The Schedule tab is a table of scripted
actions (time, part, action). Problems that would stop a run are listed in red at the bottom of
the Build tab.

| Editing | |
|---|---|
| Select | click; shift-click to add or remove; drag a box over empty space; `Ctrl+A` for everything |
| Move, rotate, delete | drag the selection; `R`; `Delete` |
| Copy, paste, duplicate | `Ctrl+C`; `Ctrl+V` at the pointer; `Ctrl+D` in place |
| Wire corners | click a wire to add a corner and drag it; right-click a corner to remove it, a wire to cut it |
| Undo, redo | `Ctrl+Z`, `Ctrl+Y` |
| View | right-drag pans, the wheel zooms, `F` fits the schematic to the window |
| Files | `Ctrl+O`, `Ctrl+S` |

Simulate starts the run paused; press Run or `Space`. Keys `1` to `5` set the speed (×1, ×5, ×20,
×100, max). Left-click a valve to open or close it. Right-click any part to stick a valve open or
closed, switch a heater or change a set-point. Hover any part for its pressures, flow and
temperature. The plumbing cannot be changed while running; Stop / edit returns to the editor.

A shut valve is drawn filled, with a bar across it; an open one as an outline only; one part-way
as an outline with part of the fill. The colours (red, green, amber) say the same again. Check
valves, relief valves and regulators are filled while nothing passes and outlines while something
does.

A normally-open valve is open until something commands it shut, as on a powered-down panel. The
panel says which valves are open for that reason.

The Tanks tab shows each tank's pressure, liquid share and temperatures, and what the present
load would do if the tank were closed and warmed to 30 °C and 35 °C. The Signals tab plots any
recorded signal. The View tab sets units and colours the plumbing by pressure, temperature or
vapour share. Each wire is one colour, because a wire has no pressure drop in this model; the
step in colour across a part is the drop it causes.

If the simulation cannot keep up with the speed asked for, it falls behind rather than
stuttering, and the panel shows the speed achieved. On WSL the window is software-rendered at
10-20 ms a frame; on Windows it draws on the GPU.

"Export CSV" writes two files next to the scene:

- `NAME_run.csv`: every signal, in SI units. Per vessel: pressure, liquid/ullage/wall temperature,
  mass, level, liquid share by volume and by mass, vapour pressure of the liquid, heater on/off.
  Per flow part: mass flow. Per valve: `position` (0 shut to 1 open), `open_cmd` (1 = commanded
  open) and `stuck`. Regulator and relief set-points. Ambient temperature, sun and wind.
- `NAME_run_events.csv`: every valve, heater and set-point action with its time, whether scripted
  or clicked. Loaded as a schedule it repeats the run. "Keep my clicks as the schedule" does the
  same without leaving the window.

## Headless runs

`--headless` runs the schedule or a replay with no window and prints a summary. Changes are given
on the command line, so a scene file does not have to be edited to try something:

```
julia --project bin/n2o-fill-sim.jl scenes/regulated_fill.toml --headless --until 1800 --csv out.csv

julia --project bin/n2o-fill-sim.jl scenes/hf4.toml --headless \
    --replay test/data/hf4_fill_10hz.csv --t0 9020 \
    --set VENT_ORIFICE.Cd=0.85 --set settings.interphase_E=300 \
    --csv out.csv --columns TANK.P,TANK.m
```

`--set` changes a part's parameter or a setting by name and refuses a name that does not exist.
`--schedule file.csv` adds actions from a file with columns `t_s,part,action,value`, the format
of the exported `_events.csv`. `--sample-dt 0.5` sets the shortest interval between recorded samples (default 0.1 s). A replay runs to the end of its log unless `--until` says
otherwise. Warnings (for example a heater that nothing switches on) are printed before the run.
`--help` lists every option. The HF4 fill, 16 minutes of real time, simulates in about a third of
a second.

## Test logs

"Load test log..." or `--replay log.csv` reads a telemetry CSV. Pressure, mass and temperature
columns can be overlaid on the plots, together with what saturation implies (the saturation
temperature of a logged pressure, the saturation pressure of a logged temperature). Valve and
heater columns are matched to parts by name and can be replayed as the schedule.

Check what a 1 means before replaying. Older logs (`valve_` columns, HF3 and HF4) use 1 = open.
Some later logs write the same valves as `solenoid_` columns with 1 for a reported state of
CLOSED. `solenoid_` and `relay_` columns are therefore never replayed until you choose.

A log that needs more than guessing (several boards on different clocks, valve columns with
inverted meaning, bouncing states, offsets) is described once in a manifest and loaded with
`--manifest log.toml` or through the same dialog. [docs/LOGS.md](docs/LOGS.md) has the rules and
the manifest format.

## Studies and fitting

A study runs one scene many times with parameters varied, scores each run against a log, and
writes one row per run. See [docs/STUDIES.md](docs/STUDIES.md).

```
julia -t 8 --project bin/n2o-fill-sim.jl study studies/hf4_vent.toml
```

A fit adjusts parameters to match one or several logs at once and reports which of them the logs
cannot tell apart. It takes a fit file that you write, or draft with `fit --init`. See
[docs/FITTING.md](docs/FITTING.md).

```
julia --project bin/n2o-fill-sim.jl fit --init scenes/hf4.toml
julia -t 8 --project bin/n2o-fill-sim.jl fit myfit.toml
```

## Scene files

A scene is a TOML file you can also edit by hand. The unit of every number is in its name.

```toml
[settings]
T_amb_C = 28.0            # also solar_W_m2, wind_m_s, P_atm_kPa, target_lo_psig, target_hi_psig,
                          # meop_psig, interphase_E, cond_tau_s, h_in_mult, h_ext_mult

[[part]]
kind = "Valve"            # Bottle, Tank, Pipe, Valve, Orifice, Regulator, CheckValve, ReliefValve,
name = "AV203"            #   Atmosphere, Junction, Gauge
x = 21                    # position on the grid, and quarter turns clockwise
y = 4
rot = 0
d_mm = 9.4                # parameters of that kind; omitted ones take their defaults
Cd = 0.5

[[wire]]
from = "AV203.b"          # part.pin; pins are a/b on two-ended parts, out on a Bottle,
to = "TANK.bottom"        #   top/vent/bottom on a Tank, a on Atmosphere, Junction and Gauge
via = [[24, 7]]           # optional: corners the drawn wire passes through (drawing only)

[[schedule]]              # optional scripted actions
t = 63.1
part = "AV203"
action = "open"           # open, close, stuck_open, stuck_closed, unstick (valves);
                          # on, off (heater, named by its vessel); set + value (regulator, relief)
```

Wires are ideal, so everything wired together is one point. A Bottle or Tank cannot be wired
straight to another vessel or to Atmosphere; put a valve, orifice or pipe between them.
Regulator, CheckValve and ReliefValve pass flow from `a` to `b` only. A Tank's `bottom` pin draws
liquid, `top` draws ullage, and `vent` draws ullage until the liquid reaches `vent_height_frac`,
then liquid. A Bottle draws liquid if `draw_height_frac` is 0 (siphon tube or tipped), vapour if 1.

| Part | What it models |
|---|---|
| Bottle, Tank | Liquid zone, ullage and wall at separate temperatures, exchanging heat and mass; cool as they deliver; can go liquid-full |
| Pipe | Friction for a length and bore, plus minor losses `K` |
| Valve | On/off, bore and `Cd` or a `Cv`; normally-open or normally-closed; stroke time; can be stuck |
| Orifice | Fixed restriction: orifice plate, quick-disconnect, a fitting |
| Regulator | Closes as its outlet reaches the set pressure; saturated liquid flashes across it |
| CheckValve, ReliefValve | One-way; cracking pressure / set pressure |
| Atmosphere | Where vented nitrous goes; totals what was lost |
| Junction, Gauge | A tee; a pressure and temperature reading |

Each part kind's parameters, with their meanings, are in one table at the top of `src/scene.jl`;
the property panel shows the same text.

## How far to trust it

[docs/MODEL.md](docs/MODEL.md) has the equations, their sources, what is left out, the
calibration, and background numbers for planning a fill. In short:

- Properties are the reference equation of state and match CoolProp to round-off.
- The plumbing in `scenes/hf4.toml` was fitted to the HF4 hot-fire log (five numbers adjusted).
  After fitting, tank pressure follows the log within about 5 psi RMS at the stall, 16 psi through
  the vents and 3 psi through the heater period, and the final load is 9.5 kg against 9.5-9.8 kg
  measured.
- Run unchanged against the HF3 log (a different day and procedure), it gets the stall pressure
  within 6 psi and the final load within the load cell's uncertainty, but reaches the stall about
  20 s too early.
- Two further logs, again with nothing in the model re-fitted: an HF3 cold flow (tank pressure
  within 2.5 psi RMS through a 100 s sealed period; 4.6 kg loaded against 4.3 kg on the load cell)
  and LC26 day 4 (a 1000 s sealed hold followed within 5 psi RMS once the outside heat input was
  fitted; the vented fill only within about 33 psi). The LC26 thermistors show the tank wall within
  about 1 K of the liquid, where the model's single wall temperature sits several kelvin above it:
  see MODEL.md, "Further checks".
- Known disagreements: the middle of the fast fill is too quick in both tests; at HF4's first stall
  the model holds 8.3 kg where the load cell reads 6.7 kg (the load cell may have been hanging up,
  see MODEL.md section 8, but that is not established).
- Several inputs are guesses nobody has measured (bottle size and contents, hose lengths, wall
  masses, heater power, vent discharge coefficient). They are marked `[guess]` or `[fit]` in the
  scene file. MODEL.md section 9 lists what to measure.
- The first ~20 s of a fill and the contents of the lines are not modelled (see MODEL.md,
  "What is left out").
- It has not been checked against a regulated fill or a fill to a lower pressure, and against a
  measured temperature and a long hold only once (LC26 day 4). Treat those results as showing
  which way things go and roughly how fast.

## Tests

```
julia --project=test test/runtests.jl     # model, scene, log and replay tests; no display needed
julia --project test/gui_input.jl         # drives the real window with injected mouse and keys
julia --project test/bench.jl [--gui]     # speed numbers and work counts
julia --project test/stress.jl            # worst cases: a supercritical bottle, random plumbing
```

`test/` needs its own `Pkg.instantiate()` once (it adds CoolProp, to check the property code).
Keep the mouse away from the window while `gui_input.jl` runs, because real input reaches it too.

## Layout

```
n2o-fill-sim.cmd      Windows launcher
bin/n2o-fill-sim.jl   command line entry
src/n2o.jl            N2O properties
src/flow.jl           mass flux through a restriction
src/vessel.jl         bottle / tank model and heat transfer
src/scene.jl          part catalogue, scene file, edit operations
src/network.jl        scene -> equations, junction solve, time stepping, recording
src/logs.jl           reading test logs, column meaning, manifests, overlay and replay
src/study.jl          studies: parameter sets, scoring, metrics
src/fit.jl            fitting to several logs, identifiability report, calibration records
src/gui*.jl           the window: state and loop, schematic, panels, plots, dialogs
scenes/               hf4.toml (fitted to the HF4 log), hf3cf.toml and lc26_day4.toml (further
                      logs), regulated_fill.toml (regulated hot-day fill)
studies/              an example study file
assets/               the window icon
test/                 tests, benchmark, stress cases, and extracts of the fill logs
docs/MODEL.md         equations, sources, omissions, calibration, sensor notes, background numbers
docs/LOGS.md          how test logs are read; the manifest format
docs/STUDIES.md       running a scene many times with parameters varied, scored against a log
docs/FITTING.md       fitting parameters to several logs at once; what the logs do and do not determine
```

## Licence

MIT; see `LICENSE`. The equation-of-state coefficients in `src/n2o.jl` were copied from
[CoolProp](https://github.com/CoolProp/CoolProp) (MIT licence), as were the surface tension
coefficients; the viscosity and thermal conductivity correlations are those of ESDU 91022 (1991).
