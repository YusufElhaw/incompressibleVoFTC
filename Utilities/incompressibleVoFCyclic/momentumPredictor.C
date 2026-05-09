/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Derived solver module
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derived from OpenFOAM and released under the GNU General
    Public License, version 3 or later.
\*---------------------------------------------------------------------------*/

#include "incompressibleVoFCyclic.H"

#include "fvmDiv.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleVoFCyclic::momentumPredictor()
{
    volVectorField& U = U_;

    tUEqn =
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + divDevTau(U)
     ==
        fvModels().source(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        // Direct-p form corresponding to
        //   -grad(p_rgh) - gh*grad(rho)  ==  -grad(p) + rho*g
        // for p = p_rgh + rho*gh.
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    surfaceTensionForce()
                  + fvc::interpolate(rho)*(buoyancy.g & mesh.Sf())/mesh.magSf()
                  - fvc::snGrad(p)
                )*mesh.magSf()
            )
        );

        fvConstraints().constrain(U);
    }
}


// ************************************************************************* //
