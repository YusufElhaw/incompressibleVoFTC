/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "incompressibleVoFTC.H"
#include "localEulerDdtScheme.H"
#include "fvCorrectPhi.H"
#include "geometricZeroField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDdt.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleVoFTC, 0);
    addToRunTimeSelectionTable(solver, incompressibleVoFTC, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleVoFTC::incompressibleVoFTC(fvMesh& mesh)
:
    twoPhaseVoFSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new incompressibleTwoPhaseVoFMixtureTC(mesh))
    ),

    mixture
    (
        refCast<incompressibleTwoPhaseVoFMixtureTC>(twoPhaseVoFSolver::mixture)
    ),

    p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*buoyancy.gh
    ),
   // H(
   //     IOobject
   //     (
   //       "H",
   //       runTime.name(),
   //       mesh,
   //       IOobject::NO_READ,
   //       IOobject::NO_WRITE
   //     ),
   //     mesh,
   //     dimensionedScalar("H", dimEnergy/dimMass, 0 ),
   //     mixture.T().boundaryField().types()
   //   
   //   ),
    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict()
    ),
   alphaRhoPhi1
    (
        IOobject::groupName("alphaRhoPhi", alpha1.group()),
        fvc::interpolate(mixture.thermo1().rho())*alphaPhi1
    ),

    alphaRhoPhi2
    (
        IOobject::groupName("alphaRhoPhi", alpha2.group()),
        fvc::interpolate(mixture.thermo2().rho())*alphaPhi2
    ),
    K("K", 0.5*magSqr(U)),

    momentumTransport
    (
        U,
        phi,
        alphaPhi1,
        alphaPhi2,
        mixture
    ),
    
    thermophysicalTransport(momentumTransport),
    Etherm   // Test energy conservation
    (
        IOobject
        (
            "Etherm",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Etherm", dimEnergy/dimVolume, 0)
    ),
    nMoles1 // Test mass conservation
    (
        IOobject
        (
            "nMoles1",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "nMoles1",dimMoles/dimVolume,0
        )
    ),

    nMoles2 // Test mass conservation
    (
        IOobject
        (
            "nMoles2",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "nMoles2",dimMoles/dimVolume,0
        )
    ),

    nMolesTotal // Test mass conservation
    (
        IOobject
        (
            "nMolesTotal",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "nMolesTotal",
            dimMoles/dimVolume,
            0
        )
    )
{
    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }

    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p_rgh,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleVoFTC::~incompressibleVoFTC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleVoFTC::prePredictor()
{
    twoPhaseVoFSolver::prePredictor();
    //     original with constant density
        //   const dimensionedScalar& rho1 = mixture.rho1();
        //   const dimensionedScalar& rho2 = mixture.rho2();
  
    const volScalarField& rho1 = mixture.thermo1().rho();
    const volScalarField& rho2 = mixture.thermo2().rho();
    alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;
    alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;

    // Calculate the mass-flux
    rhoPhi = alphaRhoPhi1 + alphaRhoPhi2;

    contErr1 =
    (
        fvc::ddt(alpha1, rho1)()()
      + fvc::div(alphaRhoPhi1)()()
    );

    contErr2 =
    (
        fvc::ddt(alpha2, rho2)()()
      + fvc::div(alphaRhoPhi2)()()
    );
}


void Foam::solvers::incompressibleVoFTC::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::incompressibleVoFTC::thermophysicalTransportPredictor()
{
    thermophysicalTransport.predict();
}


void Foam::solvers::incompressibleVoFTC::pressureCorrector()
{
    incompressiblePressureCorrector(p);
}



void Foam::solvers::incompressibleVoFTC::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::incompressibleVoFTC::thermophysicalTransportCorrector()
{
    thermophysicalTransport.correct();
}


// ************************************************************************* //
