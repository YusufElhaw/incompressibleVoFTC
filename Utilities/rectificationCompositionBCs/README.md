# rectificationCompositionBCs

Dynamic scalar boundary conditions for total-reflux rectification simulations with `incompressibleVoFTC`.

The library provides two OpenFOAM-v13 `fvPatchScalarField` boundary conditions for mole-fraction fields such as:

```foam
cyclohexane.X
nHeptane.X
```

Available types:

```foam
inletOutletCondensator
inletOutletBoiler
```

## Phase detection

The phase order in `constant/phaseProperties` is **not** used to decide which phase is liquid.

The liquid phase is detected from the single existing/registered field:

```text
alpha.<phase>
```

In the intended `incompressibleVoFTC` setup there is only one `alpha.<phase>` field, and this field is always the liquid volume fraction.

Example:

```text
0/alpha.liquid
```

means:

```text
liquid phase = liquid
gas phase    = the other phase in constant/phaseProperties
```

If `constant/phaseProperties` contains:

```foam
phases (gas liquid);
```

and the case contains:

```text
0/alpha.liquid
```

then the BC still detects:

```text
liquid phase = liquid
gas phase    = gas
```

So the order in `phaseProperties` does not matter for choosing the liquid phase. The `phaseProperties` file is only used to find the other phase name.

The user does **not** specify `phi`, `samplePhi` in the boundary file. These names are determined internally.

## Build

```bash
cd rectificationCompositionBCs
./Allwmake
```

Add the library to `system/controlDict`:

```foam
libs
(
    "librectificationCompositionBCs.so"
);
```

## General local switching

Both boundary conditions inherit from `mixedFvPatchScalarField` and behave locally like an `inletOutlet` condition.

For each boundary face:

```text
selected phi < 0 -> inflow  -> dynamic fixedValue
selected phi > 0 -> outflow -> zeroGradient
```

OpenFOAM's usual sign convention is used: positive boundary flux leaves the CFD domain, negative boundary flux enters the CFD domain.

## Condensator

The `inletOutletCondensator` condition represents a total-reflux condenser.

Physical meaning:

```text
outgoing vapour/gas composition -> incoming liquid composition
```

Internally it uses:

```text
phi       = alphaPhi.<liquidPhase>
samplePhi = alphaPhi.<gasPhase>
```

Example with:

```text
0/alpha.liquid
```

and:

```foam
phases (gas liquid);
```

becomes internally:

```text
phi       = alphaPhi.liquid
samplePhi = alphaPhi.gas
```

The local switching is therefore:

```text
alphaPhi.<liquidPhase> < 0 -> liquid inflow  -> dynamic fixedValue
alphaPhi.<liquidPhase> > 0 -> liquid outflow -> zeroGradient
```

The measured condenser value is calculated from the outgoing gas flux:

```text
X_in = sum(max(alphaPhi.<gasPhase>, 0)*X_patchInternal)
     / sum(max(alphaPhi.<gasPhase>, 0))
```

The summation is performed over all selected patches and over all MPI ranks. The measured value uses `patchInternalField()`, not the imposed boundary value.

Minimal example:

```foam
top
{
    type        inletOutletCondensator;
    value       uniform 0;
}
```

Recommended multi-patch example:

```foam
top1
{
    type        inletOutletCondensator;
    group       topCondensator;
    value       uniform 0;
}

top2
{
    type        inletOutletCondensator;
    group       topCondensator;
    value       uniform 0;
}
```

All patches with the same `group` and the same boundary-condition type are coupled and receive the same dynamic inlet value.

## Boiler

The `inletOutletBoiler` condition represents a total-reflux boiler/reboiler.

Physical meaning in Variant A:

```text
outgoing liquid composition -> incoming gas composition
```

This is a total-reflux reboiler mass-balance closure without an explicit
reboiler reservoir state. The BC does **not** apply VLE directly to the
incoming/outgoing bottom liquid. This avoids the artificial algebraic
self-enrichment loop `x_liquid -> y_eq -> condenser -> liquid -> ...`.


Internally it uses:

```text
phi       = alphaPhi.<gasPhase>
samplePhi = alphaPhi.<liquidPhase>
```

Example with:

```text
0/alpha.liquid
```

and:

```foam
phases (gas liquid);
```

becomes internally:

```text
phi       = alphaPhi.gas
samplePhi = alphaPhi.liquid
```

The local switching is therefore:

```text
alphaPhi.<gasPhase> < 0 -> gas inflow  -> dynamic fixedValue
alphaPhi.<gasPhase> > 0 -> gas outflow -> zeroGradient
```

The boiler first measures the outgoing liquid mole fraction:

```text
x = sum(max(alphaPhi.<liquidPhase>, 0)*X_patchInternal)
  / sum(max(alphaPhi.<liquidPhase>, 0))
```

The incoming gas mole fraction is then set from the total-reflux
mass balance:

```text
y_in,gas = x
```


Minimal example:

```foam
bottom
{
    type        inletOutletBoiler;
    value       uniform 0;
}
```

Recommended multi-patch example:

```foam
bottom1
{
    type        inletOutletBoiler;
    group       bottomBoiler;
    value       uniform 0;
}

bottom2
{
    type        inletOutletBoiler;
    group       bottomBoiler;
    value       uniform 0;
}
```



## Multi-patch coupling

The recommended syntax is `group`.

All patches with the same BC type and the same group name are coupled together.

Example:

```foam
group topCondensator;
```

or:

```foam
group bottomBoiler;
```

## Explicit patch lists

Instead of `group`, an explicit patch list can be used:

```foam
top1
{
    type        inletOutletCondensator;
    patches     (top1 top2 top3);
    value       uniform 0;
}
```

The current patch is added automatically if it is missing.

Legacy entries are still accepted:

```foam
otherCondensatorPatches (top2 top3);
otherBoilerPatches      (bottom2 bottom3);
```

## Optional entries

```foam
group       topCondensator;
patches     (top1 top2 top3);
relax       1;
inletValue  0;
minValue    0;
maxValue    1;
minFlux     1e-300;
```

`relax` is numerical damping for the dynamic inlet value. `relax 1` means no damping.

`inletValue` is a restart/fallback value only. In normal operation, the value is calculated dynamically.

`value` should be present for OpenFOAM restart/read consistency.

## Notes

- Do not specify `phi`.
- Do not specify `samplePhi`.
- The liquid phase is detected from the existing/registered `alpha.<phase>` field.
- The gas phase is the other phase listed in `constant/phaseProperties`.
- The condensator samples outgoing gas flux and imposes incoming liquid composition.
- The boiler samples outgoing liquid flux and imposes the same value as incoming gas composition: `y_in,gas = x_out,liquid`.
- Positive flux means outflow from the CFD domain.
- Negative flux means inflow into the CFD domain.
- The measured composition uses `patchInternalField()`, not the imposed boundary value.


## Fenske post-processing with Variant A

For Fenske evaluation, use:

```text
x_D = outgoing top vapour composition
x_B = equivalent reboiler reservoir liquid composition
```

The boiler BC gives the entering gas composition from the bottom total-reflux
mass balance:

```text
y_B,gas = x_out,liquid
```

The equivalent reservoir liquid composition is then inferred from the inverse
relative-volatility relation:

```text
x_B,reservoir = y_B,gas/(A12 - (A12 - 1)*y_B,gas)
```

The file `examples/fenskePatchData.massBalance.controlDict` contains a coded
functionObject that writes these values over time and also computes

```text
Nmin = ln((x_D/(1-x_D))*((1-x_B)/x_B))/ln(A12_columnGeo)
```
