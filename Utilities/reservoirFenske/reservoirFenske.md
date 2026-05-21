# reservoirFenske — Function Object Documentation

Compiled function object in
`Utilities/reservoirFenske/` (library `libreservoirFenske.so`).
Computes HETP, Fenske minimum stages, and global species conservation for total-reflux
rectification simulations using the `inletOutletCondenser` / `inletOutletBoiler`
boundary conditions.

Replaces the previous `coded`-type inline block in `system/functions`. Being a
pre-compiled library it avoids dynamic compilation at solver startup and can be reused
across all cases without copy-pasting code.

---

## Purpose

At every execute interval the function object:

1. Reads the reservoir states (ND, xD, NB, xB, yB) directly from the active BC objects.
2. Integrates molar density fields over the column volume (NColumn, N1Column, N2Column).
3. Computes the system-wide mole totals (column + both reservoirs).
4. Tracks drift of NSystem, N1System, N2System relative to t = 0 as a conservation check.
5. Computes the volume-weighted geometric mean of A12 over the column.
6. Applies the Fenske equation to obtain Nmin and HETP.

All results are appended to `postProcessing/<name>/0/values.dat`, where `<name>` is
the key used for this entry in `system/functions` (not necessarily the type name).

---

## Build

```bash
cd ~/OpenFOAM/elhaw-13/run/meineSolvers/incompressibleVoFTC
wmake Utilities/reservoirFenske
```

Or via the top-level build script:

```bash
./Allwmake
```

Output: `$FOAM_USER_LIBBIN/libreservoirFenske.so`

Depends on `librectificationCompositionBCs.so` (must be built first — the `Allwmake`
script ensures the correct order).

---

## Configuration in system/functions

```cpp
reservoirFenske             // key = name used in the output path
{
    type            reservoirFenske;
    libs            ("libreservoirFenske.so");

    executeControl  timeStep;
    executeInterval 1;
    writeControl    none;           // write() is a no-op; output happens in execute()

    speciesField    cyclohexane.X;  // mole-fraction field (X1)
    A12Field        A12;            // relative volatility field

    columnTotalMolesField    nMolesTotal;
    columnSpecies1MolesField nMoles1;
    columnSpecies2MolesField nMoles2;

    // One representative patch per apparatus is sufficient.
    // All patches in a group share the same reservoir state.
    topReservoirPatches    ( bottom_Top_on_top    );
    bottomReservoirPatches ( bottom_Top_on_bottom );

    packingLength   0.05;           // [m] total packing height — must match mesh
}
```

### Multi-patch group coupling (e.g. 13 patches per apparatus)

Only ONE patch from each apparatus group needs to be listed. The function reads
the scalar reservoir state (ND, xD, NB, xB, yB) from that single BC object. Since
`advanceReservoir` runs identically on all patches in a group and produces the same
result via `reduce`, any representative patch returns the correct system-wide values.

Listing all 13 patches would not change the result but is unnecessary.

---

## Data sources

### Reservoir states — read from BC objects

`readCondenserReservoir` and `readBoilerReservoir` iterate over the configured
`topReservoirPatches` / `bottomReservoirPatches` lists and cast the first matching
patch BC to `inletOutletCondenserFvPatchScalarField` or
`inletOutletBoilerFvPatchScalarField`. They then return the scalar members:

| Symbol | BC member       | Meaning                                        |
|--------|-----------------|------------------------------------------------|
| ND     | reservoirMoles_ | Total moles in condenser reservoir [mol]       |
| xD     | reservoirX_     | Liquid mole fraction of species 1 in condenser |
| NB     | reservoirMoles_ | Total moles in boiler reservoir [mol]          |
| xB     | reservoirX_     | Liquid mole fraction of species 1 in boiler    |
| yB     | reservoirY_     | Vapour mole fraction in boiler (VLE with xB)   |

These are scalar values, not face-averaged quantities. No mass balance error is
introduced by reading them — the function object is purely diagnostic.

### Column molar inventory — volume integrals

`volumeIntegral` sums `f[celli] * V[celli]` over all local cells and reduces across
processors:

| Quantity  | Field             | Meaning                                   |
|-----------|-------------------|-------------------------------------------|
| NColumn   | nMolesTotal       | Total moles in column volume [mol]        |
| N1Column  | nMoles1           | Moles of species 1 in column volume [mol] |
| N2Column  | nMoles2           | Moles of species 2 in column volume [mol] |

Only `internalField()` is used — no boundary contribution.

### Geometric mean A12 — volume-weighted

`volumeGeoMean` computes:

```
A12_geo = exp( sum_i(V_i * ln(A12_i)) / sum_i(V_i) )
```

where the sum runs over all cells with `A12 > SMALL`. This is the correct mean for a
multiplicative quantity and avoids the arithmetic mean overestimating A12 in columns
with composition gradients.

---

## System-wide mole balance

```
NSystem  = NColumn + ND + NB
N1System = N1Column + ND*xD + NB*xB
N2System = N2Column + ND*(1-xD) + NB*(1-xB)
```

The reference values are captured in the class member variables `NSystem0_`,
`N1System0_`, `N2System0_` at the first time step where `NSystem > VSMALL`.
The `initialized_` flag guards the capture. These members are initialised to zero in
the constructor and reset whenever `foamRun` starts, so `dN` is always relative to the
beginning of the current run (not the original t = 0 when restarting a simulation).

Drift quantities written to file:

| Column | Formula              | Meaning                        |
|--------|----------------------|--------------------------------|
| dN     | NSystem  − NSystem0  | Absolute total mole drift [mol] |
| dN1    | N1System − N1System0 | Absolute species-1 drift [mol] |
| dN2    | N2System − N2System0 | Absolute species-2 drift [mol] |
| relN   | dN  / NSystem0       | Relative total drift [-]       |
| relN1  | dN1 / N1System0      | Relative species-1 drift [-]   |
| relN2  | dN2 / N2System0      | Relative species-2 drift [-]   |

A nonzero drift means moles are created or destroyed numerically. Acceptable causes are
VoF compressibility corrections and the explicit-Euler temporal splitting of the
reservoir ODE. Values below ~0.5 % are typical for physical diffusivities.

---

## HETP and Fenske calculation

### Fenske minimum stages

```
Nmin = ln( xD/(1-xD) * (1-xB)/xB ) / ln(A12_geo)
```

Valid only when: `xD > 0`, `xB > 0`, `xD < 1`, `xB < 1`, `xD > xB`, `A12_geo > 1`.
If any condition fails, `nan` is written.

**Convention:** The boiler BC acts as a partial reboiler (equilibrium stage): it
computes `yB = A12·xB / (1 + (A12−1)·xB)` from the reservoir liquid `xB`. The
standard Fenske equation with these `xD` and `xB` values therefore counts the boiler
as one of the Nmin stages:

```
Nmin = N_packing + 1
```

where `N_packing` is the number of equilibrium stages provided by the packed column
itself. When `Nmin = 1` the column contributes zero stages (only the boiler's VLE
counts), which is equivalent to `xD ≈ yB`.

### Packing stages and HETP

```
N_packing = Nmin − 1        (column stages only, boiler excluded)
HETP      = packingLength / N_packing
```

`packingLength` [m] is the total active packing height of the simulated domain.
It must equal the column height in the mesh (verify with blockMeshDict).

`HETP` is written as `nan` when `Nmin ≤ 1` (i.e. `N_packing ≤ 0`), because the
column provides no measurable separation and HETP is physically undefined in that case.

For the 2D validation case: `packingLength = 0.05 m` (blockMeshDict: `y = 50 mm`).

---

## Output file format

`postProcessing/<name>/0/values.dat` — space-separated, one row appended per
`executeInterval`. The subfolder `<name>` is the key used in `system/functions`.

```
# Time  ND xD  NB xB yB  NColumn N1Column N2Column
        NSystem N1System N2System
        dN dN1 dN2  relN relN1 relN2
        A12ColumnGeo  NminColumnGeo  HETPColumnGeo
```

| Column        | Formula                            | Meaning                                    |
|---------------|------------------------------------|--------------------------------------------|
| NminColumnGeo | Fenske eq.                         | Total stages incl. boiler [-]              |
| HETPColumnGeo | packingLength / (NminColumnGeo − 1)| Height equivalent to one column stage [m]  |

`nan` is written for `A12ColumnGeo`, `NminColumnGeo`, `HETPColumnGeo` when inputs
are outside the valid range (e.g. during the initial transient when xD ≤ xB, or when
Nmin ≤ 1 for HETP).

---

## Parallel correctness

| Operation                        | Method                                          | Status |
|----------------------------------|-------------------------------------------------|--------|
| volumeIntegral (NColumn etc.)    | Foam::reduce after local cell loop              | ✓      |
| volumeGeoMean (A12)              | Foam::reduce after local cell loop              | ✓      |
| readCondenserReservoir (scalars) | All procs hold same state via advanceReservoir reduce | ✓ |
| File write                       | Pstream::master() guard                         | ✓      |
| Baseline initialisation          | NSystem from reduce → same on all procs → initialized_ set consistently | ✓ |

---

## Known limitations

### Baseline resets on foamRun restart

`initialized_`, `NSystem0_`, `N1System0_`, `N2System0_` are class member variables.
They are zero-initialised in the constructor and re-initialised from the first time
step of each `foamRun` call. When continuing a simulation (`startFrom latestTime`),
`dN` measures drift from that restart point, not from the original t = 0. The Python
postprocessing script `HETP.py` handles this transparently by reading the raw values
from the file without assuming a common baseline.

### HETP undefined when Nmin ≤ 1

When the packed column provides no separation (only the boiler's VLE counts), Nmin = 1
and `HETPColumnGeo` is written as `nan`. This is the physically correct behaviour:
HETP = packingLength / 0 is undefined. A `nan` HETP indicates that the column mass
transfer is insufficient to produce measurable enrichment above the boiler stage.

### A12 geometric mean excludes cells with A12 ≤ SMALL

Cells where A12 ≈ 0 (unphysical; relative volatility must be positive) are excluded.
Under normal operation all cells have A12 > 1 and no cells are excluded.

### Fenske equation assumes constant A12 through the column

The geometric mean A12 is the standard approximation for a column with a spatially
varying relative volatility. The resulting Nmin (and therefore HETP) is a
column-averaged estimate.

---

## Dependency on other utilities

| Dependency                          | Required for                                       |
|-------------------------------------|----------------------------------------------------|
| librectificationCompositionBCs.so   | Casting BC to inletOutletCondenser/Boiler, reading reservoir state |
| nMolesTotal, nMoles1, nMoles2 fields | Computed by incompressibleVoFTC solver            |
| A12 field                           | Computed by compositionPredictor                   |
| cyclohexane.X (or configured X field) | Solved by compositionPredictor                   |
