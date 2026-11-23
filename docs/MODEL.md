# n2o-fill-sim model

What the simulator computes, where each piece comes from, what it leaves out, and how well it
reproduced the test logs it has been compared with. It is for anyone who wants to trust, check or
change a number the simulator gives: read it before relying on one.

Sections 1 to 5 describe the model, section 6 what it leaves out, sections 7 and 8 the comparison
with test logs, and section 9 gives property numbers for planning a fill. Sources are in
section 10.

The code uses SI units throughout, and pressures are absolute. Gauge values in the window and in
scene files are relative to the `P_atm_kPa` setting.

## 1. Structure

A scene compiles to a graph.

| Element | Role |
|---|---|
| Vessels (Bottle, Tank) | Hold nitrous and carry all the state. |
| Flow parts (Pipe, Valve, Orifice, Regulator, CheckValve, ReliefValve) | Edges. Each gives a mass flow from the pressure at each end and the enthalpy of whatever is upstream. They hold no fluid. |
| Junctions (where wires meet without a vessel) | Hold no fluid either. Their pressures are whatever makes inflow equal outflow. |
| Atmosphere | A fixed pressure. Mass and enthalpy sent to it are totalled. |

### Junction pressures

Junction pressures are solved every time the equations are evaluated, the way a circuit simulator
solves node voltages. Streams meeting at a junction mix (enthalpy is carried with the flow), so
flashing across a regulator or along a hose appears without special handling.

Dead ends are removed before the solve. A junction with only one open path cannot pass flow, so
it takes the pressure of what it hangs off and the part leading to it carries exactly zero. With
every valve shut there is nothing left to solve. Where the part leading to the dead end shuts on
the dead end's own pressure, it stops earlier: behind a shut valve a regulator's outlet rests at
min(supply, set-point) (lock-up) and never falls with the supply (it is one-way); a check valve
stops one cracking pressure short; a relief valve reseats at its set pressure.

The solve is Newton's method on the junction pressures, to a mass imbalance of 1e-9 of the
largest flow. It also stops once no pressure moves by more than 0.01 Pa, because the flow laws
carry round-off of their own (largest with supercritical fluid upstream) and near it that
tolerance cannot be met; an imbalance under 1e-6 of the largest flow is then accepted. If Newton
does not settle (a check valve cracking, a restriction starting to choke), each junction is
balanced in turn with its neighbours held, repeatedly. A solve that still fails is counted, and
the window shows the count in amber. A junction that is in balance and that no flow responds to
(every part on it at rest) is held where it is while the others are solved; left in, its empty
row makes the whole system singular.

### Time stepping

The vessel states are integrated with a stiff ODE solver (Rodas5P, adaptive step). Valve strokes,
schedule actions and switches between vessel modes happen between steps.

The solver chooses its own step size (milliseconds as a valve opens, seconds when nothing is
changing) and is allowed to run ahead of the time being displayed; what is shown and recorded is
read off its interpolant. When a valve is clicked the solver is wound back to the displayed time
first, so the click takes effect exactly then. The alternative, forcing a step to end on every
screen frame, costs 60 steps per simulated second and makes the window lag.

The solver needs the Jacobian of the right-hand side at every step. Each vessel's own equations
are differenced, one column per state. The junction solve is not solved again for each column.
A change in a vessel's state reaches the rest of the system only through that vessel's pressure
and the enthalpy at its ports, and to first order the junctions answer those linearly. So once per
Jacobian the network is linearised about the solution just found: each part's flow is differenced
with respect to the pressure at each end and the enthalpy arriving from upstream, and the mass
balance and mixture enthalpy of every junction are written out to first order and factorised.
Each column is then one evaluation of the perturbed vessel, a back-substitution for the change in
every flow and junction enthalpy, and the vessels' derivatives with those flows. A wall
temperature changes neither the vessel's pressure nor what leaves through its ports, and a vessel
with every connection shut has no flow to change: those columns need only that vessel's
derivative. If the junction solve did not converge there is nothing sound to linearise about, and
every column is differenced through the whole right-hand side instead. On the HF4 replay a
Jacobian costs 42 to 54 part-flow evaluations this way, against 250 to 520 for the plain
differences, and the two agree to about 3e-4 of the largest entry in each row (the plain
differences are themselves no closer than that to a central difference).

### Balances and recording

A running total of mass and energy (vessel contents + wall heat + vented − heat taken in) is kept
and shown; it should stay at zero. In the tests it holds to about 1e-6 (mass) and 1e-5 (energy).

A run records, at up to ten samples a second, every vessel and junction state, every flow, each
valve's position and command, heater states, set-points and the ambient conditions, plus a
timestamped log of every action applied. The README lists the exported columns.

## 2. Nitrous oxide properties (`src/n2o.jl`)

- Equation of state: Lemmon & Span (2006), the short Helmholtz-energy form used by REFPROP and
  CoolProp, implemented directly. Coefficients were copied from CoolProp's `NitrousOxide.json`.
- The saturation line is solved from the equation of state at load (equal pressure and Gibbs
  energy) and tabulated up to 309.40 K; above that the fluid is treated as single-phase. The
  table has 2048 rows and is read as a cubic with fourth-order slopes, which reproduces the
  equation of state to about 1e-11. That accuracy is not for its own sake: where an expansion
  enters the dome the table takes over from the equation of state, and any mismatch is a step in
  the flux that the junction solve can cycle around (a 512-row table, good to 1e-7 inside and
  1e-4 at its ends, does exactly that).
- States are found from (ρ, u), (P, h) or (P, s) by Newton's method; inside the dome they are
  equilibrium liquid/vapour mixtures.
- Checked in the tests against CoolProp (saturation line to 1e-5 or better, single-phase states to
  1e-8) and against the tables in section 9.
- Viscosity and thermal conductivity: ESDU 91022 (1991) saturated-liquid and saturated-vapour
  correlations. ESDU gives viscosity for −90 to 30 °C and conductivity for −90 to 10 °C and says
  not to extrapolate, so outside those ranges the nearest in-range value is used. ESDU flags most
  of its own transport values as not backed by measurement. They only enter through heat-transfer
  coefficients and pipe Reynolds numbers.
- No solid phase: below the triple point (−90.8 °C) temperature is clamped.

## 3. Flow through a restriction (`src/flow.jl`)

Every flow part is an effective area times one mass-flux function G(P₁, h₁, P₂).

- Liquid or two-phase upstream: the Dyer non-homogeneous non-equilibrium model as written in
  Zimmerman et al. (2013), eq. 7:
  G = (κ·G_SPI + G_HEM)/(1 + κ), κ = √((P₁ − P₂)/(P_sat(T₁) − P₂)),
  with G_SPI = √(2ρ₁(P₁ − P₂)) and G_HEM = ρ₂√(2(h₁ − h₂)) evaluated along the isentrope and capped
  at its maximum (choking). Subcooled liquid discharging above its vapour pressure is pure G_SPI.
  The maximum is searched for along the saturation table when the expansion is inside the dome.
  Coming from single phase, the point where the expansion enters the dome is evaluated directly
  (the speed of sound drops there and the flux often peaks exactly at it), and the single-phase
  stretch above it is searched only if the flow is already supersonic at that point.
  With vapour already present upstream the result is faded towards G_HEM by quality.
- Vapour or supercritical upstream: G_HEM alone, i.e. real-gas isentropic nozzle flow, choked when
  the back pressure is low enough.
- Within 2 K of the critical temperature the G_SPI share is faded out linearly, so that the two
  cases above meet. This is a numerical choice, not a measured one: liquid and vapour become one
  fluid there, and without the fade the flux stepped by about 15 % where a supercritical supply
  crossed into the two-phase region. Nothing below 34.2 °C is affected (HF4's bottle was 32.2 °C).
- Supercritical and near-critical conditions (a bottle left in the sun past 36.4 °C) run, but are
  unvalidated and slow: every flux evaluation there needs iterative property solves instead of
  table look-ups (measured: 51 µs against 3.6 µs for saturated liquid at 25 °C), and the window
  can still hitch for up to a second when such a bottle is opened into an empty tank.
- Effective area: `Cd·π·d²/4`; a valve may give `Cv` instead (1.698e-5 m² per unit Cv); a valve's
  area scales with how far open it is.
- Pipe: area × (1 + K + f·L/D)^-½, with the Churchill friction factor. Below Re = 4000 the value
  at 4000 is used (see section 6).
- Regulator: area × clamp((set − P_outlet)/droop, 0, 1), flow a → b only.
- Relief valve: opens linearly from its set pressure to fully open 10 % above. Check valve: opens
  over a band above its cracking pressure.

`Cd` on every part is a calibration knob. Typical sharp-edged values are 0.6-0.8.

## 4. Bottles and tanks (`src/vessel.jl`)

### Stratified mode

The normal mode. States: liquid mass and temperature, ullage N2O mass and
temperature, wall temperature, air mass. Following Zilliac & Karabeyoglu (2005) as set out in
Zimmerman et al. (2013):

- Liquid has saturated-liquid density and energy at its own temperature. Ullage is real gas, which
  may be superheated or slightly supersaturated. Tank pressure is the ullage N2O pressure plus the
  air partial pressure (ideal gas).
- The surface is saturated at the N2O partial pressure. Heat reaches it from the liquid and leaves
  it to the ullage by natural convection (Nu = 0.15·Ra^⅓). Evaporation is the net heat arriving
  divided by the latent heat. When the liquid is hotter than the surface (it wants to boil) the
  liquid-side coefficient is multiplied by `interphase_E`; the review reports fitted values from
  240 to 1300 and no way to predict it. Default 700.
- Ullage denser than saturated vapour condenses with time constant `cond_tau_s` and the
  condensate joins the liquid.
- Wall: one temperature. Natural convection to liquid and to ullage (Nu = 0.021·Ra^0.4) over the
  wetted and dry areas, scaled by `h_in_mult`. Outside: natural convection (0.59·Ra^¼) combined
  with wind over a cylinder (Hilpert), scaled by `h_ext_mult`; absorbed sunlight on the projected
  area; infrared exchange with ambient; heater power.
- Optional wall model, off by default. Each part of it is a parameter of the vessel, and with
  all three left off the equations and results are those of the single wall temperature above:
  - `wall_split = 1`: two wall temperatures instead of one. The tube (a share `wall_tube_frac` of
    the wall mass) exchanges heat with the liquid over its wetted part and with the ullage over
    its dry part, and takes the sun and the heater; the ends take the bottom to the liquid and the
    top to the ullage; `wall_link_W_K` conducts between them. A thin tube follows the contents in
    seconds while heavy bulkheads lag for minutes, which one temperature cannot do.
    `wall_ends_W_K`, if set, puts a conduction path in series with the ends' contact with the
    contents: a thick bulkhead gives up its heat only as fast as the metal conducts it to the
    wetted face, however hard the liquid boils there.
  - `h_fill_mult`: adds the stirring of the liquid by the incoming jet to the wall-to-liquid
    coefficient (Dittus-Boelter on the tank diameter with the jet velocity decayed as 6·d/D,
    combined with natural convection as the cube root of the sum of cubes).
  - `h_boil_mult`: a wetted wall hotter than saturation boils the liquid against it. The flux is
    the pool-boiling data of Gorenflo and Kotthoff (2005) for carbon dioxide as Zimmerman (2015,
    eqs 5.76-5.77) fitted it for tank models, carried over to nitrous oxide by reduced pressure:
    ln h20 = 0.3092 (ln Pr)³ + 1.649 (ln Pr)² + 3.641 ln Pr + 5.272 (kW/m²K at 20 kW/m²) with
    h ∝ q^n, n = 1 − 0.3·Pr^0.3. It rises about as ΔT⁴, so a wetted wall cannot stay more than a
    kelvin or two above saturation. Capped at Zuber's critical heat flux; drying-out beyond it
    is not modelled. The data stop at Pr = 0.625 (about 14 °C for nitrous), so warm tanks are an
    extrapolation, and Zimmerman expected the relation to need adjusting: hence the multiplier.
  These exist because thermistors on a tank wall were seen within about 1 K of the liquid through
  a long hold, and falling 28 K within seconds of a fill starting, where the single wall
  temperature with HF4's fitted coupling sits several kelvin above the liquid. They are not yet
  fitted: with the other fitted values left alone, no setting tried improves every log at once
  (stronger coupling helps the thermistor log and hurts HF4's vents and heater period).
- Ports draw liquid or ullage by height: fully liquid once the level is 1 cm above the port.
- Inflow through a port below the liquid level joins the liquid with all its enthalpy (flash
  vapour bubbling up through the liquid heats it; if that leaves it superheated it boils at the
  surface). Inflow through a port in the ullage is flashed at tank pressure and split.
- Geometry is an upright cylinder of the given volume and diameter.
- Air in the ullage is a second, ideal gas with its own mass, heat capacity and partial pressure.
  It leaves with vented ullage gas. Section 6 gives the limits of this.

### Homogeneous mode

One mixed volume at (m, U) with the full equation of state. A vessel switches
to it when ullage falls below 0.5 % of the volume, or the liquid passes 35.3 °C, and switches back
when a real ullage reappears. Mass and energy carry over exactly. This is what makes an overfilled
tank warming up show a steep pressure rise and lift its relief valve instead of failing.

## 5. Warm-soak check

The run panel shows, for the mass now in each tank, the pressure and liquid fraction if the tank
were closed and brought to 30 °C and 35 °C. That is a property calculation on (mass, volume,
temperature), independent of the dynamic model, and it reproduces the liquid-full thresholds in
section 9.

## 6. What is left out, and numerical compromises

- **The lines hold no fluid.** Giving every junction a volume and a fluid state would show a
  hose charging and liquid trapped between closed valves, but it could not be solved reliably:
  a few millilitres of nearly incompressible liquid sitting on the edge of boiling make the
  equations so stiff and kinked that every solver tried took microsecond steps or stalled.
  Consequences: no hose-charging transient; the slug of liquid in a pre-charged hose does not
  arrive at the start of a fill; no heat picked up in the hose; an isolated section of line just
  holds its last pressure. In HF4 the hose was charged before the tank valve opened and the model
  is about 100 psi low for the first 10-20 s; in HF3 the hose was not pre-charged and the first
  20 s agree to 27 psi, which fits that explanation without proving it. To study trapped liquid,
  put a small Tank where the trapped volume is (there is a test doing exactly that).
- **Flow near zero pressure difference** is smoothed within about 300 Pa so it passes through zero
  with a finite slope; the true √ΔP law has infinite slope there.
- **Laminar pipe flow** is not modelled (friction factor held below Re = 4000). It only matters
  for trickles of a gram per second or less.
- No hydrostatic head or elevation differences (about 1 psi per metre of liquid).
- No nitrogen pressurisation or dissolved nitrogen. Pure N2O plus inert air in a vessel's ullage.
- **Air in a vessel is a second gas in the ullage only.** It has its own heat capacity and
  partial pressure, and leaves with the ullage gas in proportion to how much of the gas it is.
  Once in the lines the stream is counted as N2O (the lines carry one fluid, sized by the flow
  law as if it were N2O vapour at the same pressure and temperature, which overstates the mass
  flow of nearly pure air by about a fifth), and anything that flows back in from the atmosphere
  arrives as N2O. A tank's worth of air is tens of grams, so this matters only until the first
  N2O has displaced it.
- One wall temperature per vessel (two with `wall_split`, section 4); no conduction along the
  wall beyond that link, or into stands and plumbing.
- A single liquid temperature: no warm surface layer. Vent rebound comes only from the lag between
  liquid and ullage.
- No boiling swell. In a real tank, bubbles formed in the liquid while venting raise the level for
  as long as it boils, and can put liquid at the dip tube before the settled level gets there. Here
  the liquid level is always the settled level, so "liquid %" can fall during a vent (the liquid
  boils off and contracts as it cools, about 1 % in volume per kelvin near 23 °C) and liquid
  reaches the vent later than it would in practice.
- Pipe friction is treated as an equivalent orifice, which is rough for a long line that flashes
  along its length.
- The relief valve has no pop and reseat hysteresis.

## 7. Calibration against HF4

HF4 was a static hot fire on 2026-07-11. `scenes/hf4.toml` is the fill panel and tank as drawn
in the QRET Chimera Final Design Report (2026), Fig. 35. The logged valve and heater states
(`test/data/hf4_fill_10hz.csv`, a 10 Hz slice of the hot-fire log) drive the model from device
time 9020 s to 10010 s. Values fitted, in the order they were pinned down:

| Parameter | Value | Pinned by |
|---|---|---|
| Bottle temperature | 32.2 °C | PT201 before the fill (934-952 psig) |
| Bottle contents / volume | 20 kg / 40 L | how far bottle pressure falls as it delivers (stall pressure) |
| Restriction bottle → panel | 2.5 mm at Cd 0.7 | fill rate; the procedure lists a 1/4" hose with a CGA 326 nipple here |
| `h_in_mult` | 0.3 | vent drop and rebound, and the heater slope |
| Heater | 300 W | heater period, 775 → 843 psig in 469 s |

Left at defaults or guesses: `interphase_E` (HF4 does not constrain it: 100 to 2000 changes the
result by a few psi), vent `Cd` 0.7, tank wall 10 kg, ambient 28 °C, sun 600 W/m².

Tank pressure, model against PT202 (RMS over each period):

| Period (device s) | What happens | RMS error |
|---|---|---|
| 9084-9104 | first 20 s of fill | 109 psi (probably the pre-charged hose; see section 6) |
| 9104-9160 | rest of the fast fill | 28 psi |
| 9160-9256 | stall | 5 psi |
| 9256-9330 | vent 1 and rebound | 16 psi (model drops to 723, log to 697) |
| 9330-9372 | vent 2 and close | 8 psi |
| 9387-9856 | heater | 3 psi |
| PT201, 9100-10000 | bottle side | 9 psi |

Mass: 9.5 kg at the end of fill; the load cell, tared with the hose charged, reads 9.5 kg at fill
close and 9.8 kg before firing. **Open discrepancy:** at the stall the model has 8.3 kg where the
load cell reads 6.7 kg. Section 8 gives a reason to doubt the load cell there, but it is a
hypothesis; if the load cell is right, the model overstates how much goes in before the first vent.

This is a fit to one fill, with five adjusted numbers. It shows the model has the right structure
(bottle cooling sets the stall; the vent sets a floor near the liquid's vapour pressure; heat
through the wall sets the slow drift). It is not a validation.

### Independent check: HF3

HF3 (2026-06-07) was a different procedure on the same hardware: AV-203 open first, then one unvented fill
through AV-201 for 288 s, bottle at 849 psig. The HF4 scene was run against it with nothing
re-fitted. Only the bottle's starting pressure (measured) and the ambient temperature (22 °C, a
guess; the log has no temperature) were changed. `test/data/hf3_fill_10hz.csv`.

| Period (device s) | What happens | RMS error, tank pressure |
|---|---|---|
| 4528-4548 | first 20 s | 27 psi |
| 4548-4620 | fast fill | 53 psi (model reaches the stall about 20 s early) |
| 4620-4814 | stall | 6 psi (model 722-725 psig, log 711-723) |

Final load 9.0 kg against a load cell reading of 9.3 kg that includes roughly 0.6 kg of hose load.
So the stall pressure and the load carried over to a different day and bottle state; the fill rate
in the middle of the fill is too high in both tests and is the model's clearest error. Moving the
guessed ambient temperature by ±6 °C moves the final load by about ±0.5 kg, so that guess matters.
After the fill, HF3's log shows the run heater switched on for 150 s while tank pressure fell
about 5 psi; the model with HF4's 300 W rises 20 psi. The log records that the heater was switched
on, not what it delivered.

### Speed

The 990 s HF4 replay takes about 0.3 s (559 solver steps, about 230,000 equation-of-state
evaluations). In the window at ×1 the solver uses a few hundredths of a millisecond of a typical
frame and 3.5 ms in the worst single step. `test/stress.jl` drives scenes nobody tuned (a bottle
heated until supercritical, and random plumbing from default parts); chains of check valves and
regulators are the remaining weak spot, with junction solves that fail to settle and frames of up
to about 0.6 s. `test/bench.jl` prints these numbers.

### Further checks: HF3 cold flow and LC26 day 4

Two more logs, replayable with `scenes/hf3cf.toml` and `scenes/lc26_day4.toml` (data in `test/data`).

- **HF3 cold flow (nitrous, unvented fill, 4.3 kg).** HF4 fit unchanged: tank pressure through the
  100 s sealed period within 5.1 psi RMS, bottle recovery within about 4 psi, load 4.6 kg against 4.3 kg on
  the stand load cell. Post-fill warming is 0.08 psi/s in the model against 0.13 logged (ambient is a
  guess). Open: the model empties the tank through the injector, while the load cell shows ~2.3 kg
  left at 0 psig.
- **HF4 vents prefer a larger vent.** Vent `Cd` 0.85 with `interphase_E` 300 takes vent 1 from 16 to
  8 psi RMS (minimum 696 against 696 logged) and vent 2 from 8 to 3, but the heater period from 3 to
  7 psi. Final load 9.7 kg. The scenes keep the published values until the wall model can carry
  both.
- **LC26 day 4 has tank temperatures.** The two heater thermistors sit within about 1 K of the
  saturation temperature of PT202 through the 1000 s sealed hold, and 1-2 K above it during the
  vented fill. With `h_in_mult` 0.3 the model's wall runs ~6 K warmer than its liquid after the fill,
  so the single wall node couples too weakly for the thin tube wall while HF4 still wants 0.3: a sign
  the wall wants two nodes (thin wetted tube, heavy bulkheads and fittings), not a different number.
- **Not identifiable from LC26:** the tower load cell is unusable, so the vent size, fill flow and
  bottle contents trade off; the scene uses vent `Cd` 1.0 and a 28 kg bottle, and follows the tank
  through the hold within 5 psi RMS but the vented fill only within about 33 psi. The hold
  needs 2-4x the model's outside heat input, but ambient, sun and `h_ext_mult` are degenerate.

### Open findings

- **Measured inputs.** The HF4 and HF3 cold flow scenes carry the test days' weather from a
  historical archive (air temperature good to about 2 °C, sun, wind reduced from 10 m to stand
  height by estimate), the tank wall mass from the vehicle mass budget (9.06 kg: tube 7.18 kg,
  bulkheads 0.88 and 1.00 kg) and bottle start temperatures from the bottle-side pressure before
  the fill. Tank start temperature on the stand was not measured and is set to ambient.
- **The fill line is probably tighter than the scenes say.** The load cell's fill rate against the
  logged line pressures puts the restriction between bottle and panel at about 2.2-2.4 mm bore
  (at `Cd` 0.7) on HF3, its cold flow and HF4; the scenes use 2.5 mm. With 2.3 mm the cold flow
  loads 4.0 kg and its sealed pressure comes out 36 psi low. Raising wall and outside heat
  transfer about tenfold and fourfold restores the three stand logs but makes LC26 day 4, held
  out, worse (9.2 against 4.8 in the scoring of `docs/STUDIES.md`). So the published 2.5 mm and
  wall coefficient 0.3 agree with the stand logs partly by compensating each other, and the
  scenes keep them until the tank's heat balance during a fill is understood.
- **What the logs do not determine.** With measured inputs fixed, sampling 4,000 parameter sets
  per log leaves the boiling factor `interphase_E`, the wall coefficient, the outside heat
  multiplier and the vent `Cd` free across nearly their whole plausible ranges on every log; the
  bottle wall coefficient has no visible effect at all. They are published values, not
  calibrated ones. What moves the result most is what was not recorded: how full the bottle was.
- **LC26 day 4 is not reproduced by any setting tried**: the bottle-side pressure in the first
  20 s of the fill, tank pressure and wall temperature through the fill, and tank pressure late in
  the 1000 s hold each miss by more than their tolerance whatever the parameters.

## 8. Reading the two logs

Two things about the HF4 and HF3 logs matter when comparing the model with them.

- **The two pressure transducers disagree by about ±5 psi.** With the fill valve open and almost
  no flow, PT201 (fill panel) reads 3-5 psi *above* PT202 (tank) in HF4 and 4-7 psi *below* it in
  HF3. Flow was into the tank both times, so the HF3 sign is physically impossible: the
  difference is sensor mismatch, and it changes between tests. A "stall 3-6 psi below the fill
  line" read off the HF4 log is therefore inside that mismatch and says nothing about the
  pressure difference driving the flow. Before the HF4 fill PT201 read 936 rising to 953 psig,
  which puts the bottle at about 33 °C.
- **The load cell** is usable for the end state, not for the fill history.
  - Charging the fill hose adds 0.8 kg with nothing in the tank.
  - Every vent opening and closing sets it ringing at about 5 Hz with ±2 kg amplitude for several
    seconds. Through the ringing at the first vent opening its mean rises about 1 kg within 0.3 s
    and stays up after the vent closes. Nitrous cannot arrive that fast (the fill was passing
    under 0.1 kg/s) and vent thrust would have gone away at closure. The second vent opening and
    the first closing ring the same way with no shift in the mean. One reading of this is that
    the stand was hanging up and the first jolt freed it, in which case readings before that vent,
    including the 6.7 kg at the stall, are about a kilogram low. The log cannot confirm it. HF3's
    load cell, with no vent event, climbs smoothly and agrees with the model.
  - With the tank isolated during heating it climbs 0.36 kg with no flow at all.

## 9. Background numbers for planning a fill

Property values are for saturated N2O; gauge pressures use 14.4 psi for the atmosphere. The tests
check the tables against the property code.

- **A pressure target is a temperature target.** A self-pressurising tank sits at the vapour
  pressure of its liquid. Reaching a higher pressure means adding heat (HF4: 469 s on the heater
  to reach 850 psig); holding a lower one on a warm day means removing heat and keeping it out.

  | Tank (psig) | Tsat (°C) | ρ liquid (kg/m³) | ρ vapour (kg/m³) | h_fg (kJ/kg) | cp liquid (kJ/kg·K) | dP/dT (psi/K) | Boil-off to cool 1 K (% of liquid) |
  |---|---|---|---|---|---|---|---|
  | 550 | 8.8 | 859 | 111 | 208 | 2.52 | 13.7 | 1.21 |
  | 700 | 18.9 | 794 | 152 | 174 | 3.09 | 16.3 | 1.77 |
  | 850 | 27.4 | 719 | 208 | 135 | 4.54 | 19.0 | 3.37 |
  | 950 | 32.4 | 651 | 264 | 99 | 8.36 | 20.9 | 8.45 |

  The other way: 0 °C = 438 psig, 10 °C = 566, 20 °C = 718, 25 °C = 805, 30 °C = 900,
  35 °C = 1005. A colder load is denser: the 19.16 L HF4 tank at 90 % liquid by volume holds
  15.0 kg at 550 psig and 12.8 kg at 850 psig.
- **Flash on entry.** Saturated liquid from a warm bottle partly flashes on entering a
  lower-pressure tank. That vapour has to be vented and is lost.

  | Supply liquid | into 550 psig | into 700 psig | into 850 psig |
  |---|---|---|---|
  | 25 °C (805 psig) | 20.8 % | 10.3 % | no flow |
  | 30 °C (900 psig) | 29.2 % | 20.3 % | 7.1 % |
  | 33 °C (962 psig, about HF4's bottle) | 35.7 % | 28.1 % | 17.2 % |

- **Venting to cool the load** costs the boil-off column above for every kelvin removed.
- **Cooling the supply bottle** reduces flash loss but also lowers the pressure difference that
  drives the fill. Padding the bottle with nitrogen restores it, but nitrogen dissolves in
  nitrous and breaks "pressure = vapour pressure of the liquid"; this model has no nitrogen.
- **An unvented fill stalls at equalisation** (HF4: 0.089 kg/s, then 0.008 kg/s for about 107 s).
  In a vented fill the vent sets the tank pressure, which settles where vent outflow balances
  flash vapour plus displaced ullage; the vent has to carry roughly the flash fraction times the
  fill flow. Vent flow and duration are also the icing exposure of the vent valve.
- **A vent pulse cools the surface, then the pressure rebounds** (HF4 vent 1: 803 to 698 psi in
  13.8 s, back to about 775 within a minute). A single-temperature tank model cannot do that,
  which is why liquid, ullage and wall are separate here.
- **Going liquid-full is the main overfill hazard.** A cold, dense load expands as it warms. If
  the tank is closed and soaks to ambient, a tank filled too full goes liquid-full, and pressure
  then rises hydraulically, far faster than the vapour-pressure curve.

  | Filled at | Liquid-full on warming to 30 °C if liquid volume fraction above | to 35 °C if above |
  |---|---|---|
  | 550 psig | 77 % | 64 % |
  | 700 psig | 84 % | 68 % |
  | 850 psig | 94 % | 75 % |

  The allowed fill fraction should come from the worst-case warm soak (sun, vent failed closed),
  not from the operating pressure: 35 °C is 1005 psig whatever the fill pressure was. HF4 for
  reference: about 9.3 kg in 19.16 L at about 775 psig, 53 % liquid by volume.
- **Reading the sensors.** Pressure transducers are gauge, so a saturation temperature needs the
  local atmospheric pressure (a setting here). A thermocouple needs better than ±1 °C to be
  useful: 1 °C is 14-21 psi. Pressure above the vapour pressure of the bulk liquid means a warm
  surface layer or a sensor problem.

### Inputs worth measuring

Several inputs in the scene files are guesses (marked `[guess]` or `[fit]`). Measuring these would
replace them:

- Run tank and bottle dimensions, wall material, thickness and surface area.
- Vent orifice size; vent and fill valve Cv.
- Heater power (it turns the HF4 heater period into a calibration of the tank's thermal mass).
- The design day at the launch site: maximum ambient temperature and sun, and the expected wait
  between disconnect and fire.
- What bottle temperature control is available on site (shade, ice, a cooler).
- Temperatures: of the four logs only LC26 day 4 has any (thermistors on the tank). Thermocouples
  on the fill and vent lines would show flash and liquid at the vent directly.

## 10. Sources

- Lemmon, E. W. and Span, R., "Short fundamental equations of state for 20 industrial fluids",
  J. Chem. Eng. Data 51 (2006), doi:10.1021/je050186n. Coefficients via CoolProp,
  `dev/fluids/NitrousOxide.json`.
- ESDU 91022, "Thermophysical properties of nitrous oxide" (1991), eqs 4.9-4.12.
- Zimmerman, J. E., "Self-Pressurizing Propellant Tank Dynamics", PhD dissertation, Stanford
  University (2015): evaluation of the equilibrium, Zilliac-Karabeyoglu and Casalino-Pastrone
  models (E between 240 and 1300 with no pattern, mean 693), the wall-boiling relation used here,
  and the observation that boiling dominates heat and mass transfer in these tanks.
- Gorenflo, D. and Kotthoff, S., "Review on pool boiling heat transfer of carbon dioxide",
  Int. J. Refrigeration 28 (2005) 1169-1185.
- Zuber, N., "On the stability of boiling heat transfer", Trans. ASME 80 (1958) (critical heat flux).
- Mulero, A., Cachadiña, I. and Parra, M. I., "Recommended correlations for the surface tension
  of common fluids", J. Phys. Chem. Ref. Data 41 (2012); coefficients via CoolProp.
- Zimmerman, Waxman, Cantwell and Zilliac, "Review and Evaluation of Models for Self-Pressurizing
  Propellant Tank Dynamics", AIAA 2013-4045 (tank model, Dyer flow model, convection constants).
- Zilliac and Karabeyoglu, "Modeling of Propellant Tank Pressurization", AIAA 2005-3549.
- Churchill, S. W., "Friction-factor equation spans all fluid-flow regimes", Chem. Eng. 84 (1977).
- Bergman et al., "Fundamentals of Heat and Mass Transfer" (natural convection on a vertical
  surface; Hilpert correlation for a cylinder in cross flow). The constants for these two and for
  air properties were written from the textbook values, not checked against a copy.
- QRET, Chimera Final Design Report (2026), and the team's static fire procedure, for the plumbing.
