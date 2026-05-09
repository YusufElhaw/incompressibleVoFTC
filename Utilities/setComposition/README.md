# setComposition

OpenFOAM utility for initializing two-species, two-phase composition fields for `incompressibleVoFTC`.

It reads:

- phase names from `constant/phaseProperties`
- phase1 from the existing `0/alpha.<phase>` field
- species names from `constant/physicalProperties.<phase1>/species`
- `defaultSpecie` from each `constant/physicalProperties.<phase>`
- `molWeight` from `<species>/specie/molWeight`
- phase density `rho` from each `constant/physicalProperties.<phase>` when needed
- for `perfectGas`: pressure from `0/p_rgh`, `0/p`, or `.orig` fallback and temperature from `0/T` or `0/T.orig`

It writes:

- `<species1>.X`
- `<species2>.X`
- `<species1>.<phase1>`
- `<species2>.<phase1>`
- `<species1>.<phase2>`
- `<species2>.<phase2>`

## Mixed-cell mole-fraction weighting

For interface cells, the total mole fraction field `species.X` is weighted with phase molar concentration:

```text
X_i = (alpha1*C1*x_i1 + alpha2*C2*x_i2)/(alpha1*C1 + alpha2*C2)
```

For incompressible phases:

```text
C = rho/Mmix
Mmix = x1*W1 + x2*W2
```

For phases whose `thermoType/equationOfState` is `perfectGas`:

```text
C = p/(R*T)
```

The same alpha*C logic is used for the written patch `value` of the `*.X` fields.
For `perfectGas`, patch values of `p` and `T` are used there as well.
The `inletValue` itself remains the prescribed inlet mole fraction from `setCompositionDict`.

The default gas constant is `8314.46261815324` J/(kmol K), consistent with OpenFOAM `molWeight` values in kg/kmol.

If `p_rgh` is relative/gauge pressure, add for example:

```foam
pressureOffset 101325;
```

or phase-specific:

```foam
pressureOffset
{
    gas 101325;
}
```

## Build

```bash
wmake
```

## Run

From the case directory:

```bash
setComposition
```
