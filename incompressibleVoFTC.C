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

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict()
    ),
///////////////////
// vorhanden in compressible Modell // bei Bedarf hinzufügen
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
////////////////////
    momentumTransport
    (
        U,
        phi,
        alphaPhi1,
        alphaPhi2,
        mixture
    ),
    
    thermophysicalTransport(momentumTransport)
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

    alphaRhoPhi1 = fvc::interpolate((mixture.thermo1().rho()))*alphaPhi1;
    alphaRhoPhi2 = fvc::interpolate((mixture.thermo2().rho()))*alphaPhi2; 

    // Calculate the mass-flux
    rhoPhi =
        alphaPhi1*fvc::interpolate(mixture.thermo1().rho())
      + alphaPhi2*fvc::interpolate(mixture.thermo2().rho());
    
    
    contErr1 =
    (
       ( fvc::ddt(alpha1, mixture.thermo1().rho())()() + fvc::div(alphaRhoPhi1)()()
      - (fvModels().source(alpha1, mixture.thermo1().rho())&mixture.thermo1().rho())()
    ));

    contErr2 =
    (
         (fvc::ddt(alpha2, mixture.thermo2().rho())()() + fvc::div(alphaRhoPhi2)()()
      - (fvModels().source(alpha2, mixture.thermo2().rho())&mixture.thermo2().rho())()
    ));
   
    compositionPredictor();   
}


void Foam::solvers::incompressibleVoFTC::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::incompressibleVoFTC::thermophysicalTransportPredictor()
{
    thermophysicalTransport.predict();// hinzugefügt
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
    thermophysicalTransport.correct(); // hinzugefügt
}


// ************************************************************************* //
