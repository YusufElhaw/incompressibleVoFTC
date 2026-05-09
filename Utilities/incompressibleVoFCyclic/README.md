# incompressibleVoFCyclic

`incompressibleVoFCyclic` is an OpenFOAM Foundation v13 solver module derived from `incompressibleVoF`.
It keeps the VOF/two-phase transport behaviour of `incompressibleVoF`, but overrides the momentum predictor and pressure corrector so that the solved pressure is the static pressure `p`, not the buoyant pressure `p_rgh`.

## Mathematical change

Original `p_rgh` force form in the OpenFOAM `VoFSolver` momentum predictor:

```cpp
surfaceTensionForce()
- buoyancy.ghf*fvc::snGrad(rho)
- fvc::snGrad(p_rgh)
```

Using

```text
p = p_rgh + rho*gh
```

with constant gravity gives

```text
-grad(p_rgh) - gh*grad(rho) = -grad(p) + rho*g
```

Therefore the direct-`p` solver uses

```cpp
surfaceTensionForce()
+ fvc::interpolate(rho)*(buoyancy.g & mesh.Sf())/mesh.magSf()
- fvc::snGrad(p)
```

The same replacement is applied in the pressure corrector through `phig`.

## Important inheritance note

OpenFOAM v13 `VoFSolver` constructs and reads `p_rgh` through the `buoyancy` object. This module deliberately does **not** change `VoFSolver` or any parent class. Consequently, a `p_rgh` file is still required by the inherited constructor, but the overridden equations solve and correct with `p`. The field `p_rgh` is updated as

```cpp
p_rgh = p - rho*gh;
```

for compatibility/output only.

## Build

Place the folder somewhere under your user OpenFOAM directory, for example:

```bash
mkdir -p $WM_PROJECT_USER_DIR/modules
cp -r incompressibleVoFCyclic $WM_PROJECT_USER_DIR/modules/
cd $WM_PROJECT_USER_DIR/modules/incompressibleVoFCyclic
./Allwmake
```

or simply:

```bash
wmake 
```

## Case setup

In `system/controlDict` load the library:

```foam
libs
(
    "libincompressibleVoFCyclic.so"
);
```

Run with:

```bash
foamRun -solver incompressibleVoFCyclic
```

The case must contain both:

- `0/p` — the actual pressure solved by this module
- `0/p_rgh` — still needed by the inherited `VoFSolver`/`buoyancy` construction path

For cyclic cases, put the intended cyclic pressure boundary condition on `p`. `p_rgh` can be kept consistent as a compatibility field.

## Files changed relative to `incompressibleVoF`

- `momentumPredictor.C`: direct `p` form with explicit `rho*g`
- `pressureCorrector.C`: pressure equation solves `p`, not `p_rgh`
- `incompressibleVoFCyclic.C`: reads `p`, installs a `p`-based pressure reference, and keeps `p_rgh` synchronized for output
