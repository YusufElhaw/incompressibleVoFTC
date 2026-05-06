# adjustContinuity

OpenFOAM Foundation v13 utility to correct mapped `U` boundary values so that the global boundary flux balance satisfies continuity.

This variant corrects both local inflow and local outflow on the selected patches. It supports any number of patch names. A single selected patch may contain both inflow faces and outflow faces.

## Build

```bash
cd adjustContinuity
wmake
```

If you previously built an older version, clean first:

```bash
wclean
wmake
```

## Install dictionary

Copy the example into your case:

```bash
cp system/adjustContinuityDict.example /path/to/case/system/adjustContinuityDict
```

Edit the patch list, for example:

```cpp
patches (top bottom inletLiquid outletGas);
time 0;
split 0.5;
```

## Run

From the case directory:

```bash
adjustContinuity
```

For another case:

```bash
adjustContinuity -case /path/to/case
```

## What is changed

For every selected boundary face, the utility computes:

```text
phi = U & Sf
```

- `phi > 0`: local outflow
- `phi < 0`: local inflow

The global imbalance is:

```text
netFlux = sum over all boundary faces(phi)
```

The utility scales selected outflow and selected inflow with two factors:

```text
outScale = 1 - split * netFlux / selectedOut
inScale  = 1 + (1 - split) * netFlux / selectedIn
```

With `split 0.5`, half of the correction is applied to selected outflow and half to selected inflow. Only the normal component of `U` is changed; tangential velocity components are preserved.

The corrected velocity is written back to the configured time directory, e.g. `0/U`, unless `dryRun true` is set.

## OpenFOAM v13 note

This version uses `mesh.boundaryMesh().findIndex(patchName)` instead of the older/non-existent `findPatchID(...)` call.

## Important limitation

This fixes the total volumetric flux balance of `U`. It does not independently enforce liquid and gas phase balances (`alphaPhi`). For a VoF case this is usually a restart/initialisation correction before the solver recomputes fluxes.
