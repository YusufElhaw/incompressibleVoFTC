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

#include "incompressibleVoFTC.H"

#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "findRefCell.H"
#include "fvcMeshPhi.H"
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "surfaceInterpolate.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleVoFTC::pressureCorrector()
{
    volVectorField& U = U_;
    surfaceScalarField& phi = phi_;

    fvVectorMatrix& UEqn = tUEqn.ref();
    setrAU(UEqn);

    const surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));

    while (pimple.correct())
    {
        volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p));

        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::flux(HbyA)
          + fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf)
        );

        MRF.makeRelative(phiHbyA);

        if (p.needReference())
        {
            fvc::makeRelative(phiHbyA, U);
            adjustPhi(phiHbyA, U, p);
            fvc::makeAbsolute(phiHbyA, U);
        }

        // Direct-p body/interface force contribution:
        //   p = p_rgh + rho*gh
        //   -grad(p_rgh) - gh*grad(rho) = -grad(p) + rho*g
        // Hence phig contains surface tension plus rho*g, not -gh*grad(rho).
        surfaceScalarField phig
        (
            (
                surfaceTensionForce()
              + fvc::interpolate(rho)*(buoyancy.g & mesh.Sf())/mesh.magSf()
            )*rAUf*mesh.magSf()
        );

        phiHbyA += phig;

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAUf, MRF);

        // Cache any sources using the solved pressure field p
        fvScalarMatrix pEqnSource
        (
            fvModels().sourceProxy(alpha1, p)
          + fvModels().sourceProxy(alpha2, p)
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvc::div(phiHbyA) - fvm::laplacian(rAUf, p)
             ==
                pEqnSource
            );

            pEqn.setReference
            (
                pressureReference().refCell(),
                getRefCellValue(p, pressureReference().refCell())
            );

            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();

                p.relax();

                U = HbyA
                  + rAU()*fvc::reconstruct((phig + pEqn.flux())/rAUf);
                U.correctBoundaryConditions();
                fvConstraints().constrain(U);
            }
        }

        continuityErrors();

        // Correct Uf if the mesh is moving
        fvc::correctUf(Uf, U, phi, MRF);

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (p.needReference())
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pressureReference().refValue()
              - getRefCellValue(p, pressureReference().refCell())
            );
        }

        // Compatibility/output field only.  The equations above use p directly.
        p_rgh_ = p - rho*buoyancy.gh;
    }

    clearrAU();
    tUEqn.clear();
}


// ************************************************************************* //
