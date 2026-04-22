/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
#include "fvcMeshPhi.H"
#include "fvcDdt.H"             
#include "fvmDiv.H"             
#include "fvmSup.H"             
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleVoFTC::thermophysicalPredictor()
{
  compositionPredictor();
  //T_0  298.15; // reference temperature of 298.15 K   $FOAM_ETC/controlDic
  Info <<nl<< "thermophysicalPredictor läuft"<<endl;
  
  const volScalarField& rho1(mixture.thermo1().rho()); 
  const volScalarField& rho2(mixture.thermo2().rho());
  //const volScalarField& rho(mixture.rho());

  const volScalarField& Cpv1(mixture.thermo1().Cpv());
  const volScalarField& Cpv2(mixture.thermo2().Cpv());

  // Energy of the phases (in physicalProperties Choose Enthalpy or InternalEnergy)
  const volScalarField& e1(mixture.thermo1().he());   //he Enthalpy or InternalEnergy
  const volScalarField& e2(mixture.thermo2().he());   //he Enthalpy or InternalEnergy

  const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
  const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));

  volScalarField& T = mixture.T();
   /*
   const volScalarField rhoCpv("rhoCpv", alpha1*rho1*Cpv1+alpha2*rho2*Cpv2);
   
     
   const surfaceScalarField rhoCpvPhi
   (
       "rhoCpvPhi",
       // alphaRhoCpvPhi1                  
       fvc::interpolate(Cpv1)*fvc::interpolate((mixture.thermo1().rho()))*alphaPhi1
       // + alphaRhoCpvPhi2
     + fvc::interpolate(Cpv2)*fvc::interpolate((mixture.thermo2().rho()))*alphaPhi2
   );
    
   fvScalarMatrix TEqn 
   (
     fvm::ddt(rhoCpv, T)
   + fvm::div(rhoCpvPhi, T)
   - fvm::Sp(fvc::ddt(rhoCpv) + fvc::div(rhoCpvPhi), T)                       //added for the stability by changes in rhoCp
   - fvm::laplacian(thermophysicalTransport.kappaEff(), T)
   ==
     (e1Source&e1)
   + (e2Source&e2)
   );
   
   TEqn.relax();
   fvConstraints().constrain(TEqn);
   TEqn.solve();
   fvConstraints().constrain(T);
   mixture.correctThermo();
   mixture.correct();*/
    const dimensionedScalar T0("T0", dimTemperature, 298.15);


    const volScalarField limAlpha1( min(max(alpha1, scalar(0)), scalar(1)) );
    const wordList TpatchTypes = T.boundaryField().types();
    surfaceScalarField alphaEffRho
    (
        IOobject
        (
            "alphaEffRho",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)
      *fvc::interpolate(thermophysicalTransport.alphaEff())
    );
//Info<<"Cpv1 beträgt"<<Cpv1<<endl;

    H = ( (T - T0)*(limAlpha1*rho1*Cpv1 + (1.0 - limAlpha1)*rho2*Cpv2) )/rho;

    H.correctBoundaryConditions();

    const scalar RelaxFac = 1.0;
    fvScalarMatrix EEqn
      (
          fvm::ddt(rho, H)
        + fvm::div(rhoPhi, H)
        ==
         fvc::laplacian(thermophysicalTransport.kappaEff(), T)
        + RelaxFac*
        (   
          fvm::laplacian(alphaEffRho, H) 
        - fvc::laplacian(alphaEffRho, H) 
        ) 
              
      );
        EEqn.relax();
        fvConstraints().constrain(EEqn);
        EEqn.solve();
        //Now reevaluate T for the updated enthalpy fields
        T = T0 + rho*H/(limAlpha1*rho1*Cpv1 + (1.0 - limAlpha1)*rho2*Cpv2);
        mixture.correctThermo();
        mixture.correct();
}

//void Foam::solvers::incompressibleVoFTC::energyPredictor()

// ************************************************************************* //
