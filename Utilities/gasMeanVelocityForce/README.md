# gasMeanVelocityForce for OpenFOAM Foundation v13

This user library adds an `fvConstraint` named `gasMeanVelocityForce`.
It is based on the OpenFOAM Foundation v13 `meanVelocityForce` model, but
uses the gas volume fraction `alphaG = 1 - alphaLiquid` as a weight.

## Purpose

It is intended for periodic/NCC-cyclic VoF cases where you want to impose a
prescribed **gas loading / superficial gas velocity** without introducing an
inlet/outlet pair.

The controlled quantity is

\[
j_G = \frac{1}{V_{zone}} \sum_{cells \in zone} \alpha_G (\hat e \cdot U) V
\]

where `e` is the direction of `Ubar`.

For a periodic packing segment with constant cross-section this corresponds to
superficial gas velocity in m/s, i.e. the same unit as `m^3/(m^2 s)`.

## Dictionary example

Put this into `system/fvConstraints`:

```foam
gasDrive
{
    type            gasMeanVelocityForce;

    cellZone        all;

    U               U;
    alpha           alpha.liquid;

    // target superficial gas velocity, e.g. upward in z
    Ubar            (0 0.8 0);

    relaxation      0.1;

    // Optional: ignore cells with very small gas fraction
    useCutoff       false;
    alphaCutoff     1e-3;
}
```

## Build

From the library directory:

```bash
wmake libso
```

This writes the shared library to:

```bash
$FOAM_USER_LIBBIN/libgasMeanVelocityForce.so
```

## Load the library

In `system/controlDict`:

```foam
libs
(
    "libgasMeanVelocityForce.so"
);
```

## Notes

1. `alpha` must be the **liquid** VoF field. The model computes `alphaG = 1 - alpha`.
2. The force is weighted with `alphaG`, so it acts mainly in gas-dominated cells.
3. The regulator targets superficial gas velocity, not gas-only interstitial velocity.
4. Start with a small `relaxation` value, for example `0.05 ... 0.2`.
5. This was prepared against the OpenFOAM Foundation v13 `meanVelocityForce`
   structure and may need small adjustments if your local tree differs.



## Note for incompressibleVoF

This library expects the **liquid** VoF field name in the fvConstraints entry, for example `alpha alpha.liquid;`.
It does **not** expect `alpha1` or `alpha2` from the solver internals.
Internally the forcing always uses `alphaGas = max(0, min(1, 1 - alphaLiquid))`, so the momentum forcing acts only through the gas-weighted term.
