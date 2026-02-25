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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleVoFTC::thermophysicalPredictor()
{
  compositionPredictor();

  Info <<nl<< "thermophysicalPredictor läuft"<<endl;
  const volScalarField& rho1(mixture.thermo1().rho()); 
  const volScalarField& rho2(mixture.thermo2().rho());
  
  const volScalarField& Cpv1(mixture.thermo1().Cpv());
  const volScalarField& Cpv2(mixture.thermo2().Cpv());

  // Energy of the phases (in physicalProperties Choose Enthalpy or InternalEnergy)
  const volScalarField& e1(mixture.thermo1().he());   //he Enthalpy or InternalEnergy
  const volScalarField& e2(mixture.thermo2().he());   //he Enthalpy or InternalEnergy

  const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
  const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));

   volScalarField& T = mixture.T();
  ////////////////////////////////
    
    const volScalarField rhoCpv("rhoCpv", alpha1*rho1*Cpv1+alpha2*rho2*Cpv2);

    const surfaceScalarField alphaRhoCpvPhi1
    (
      "alphaRhoCpvPhi1", fvc::interpolate(Cpv1)*alphaRhoPhi1
    );

    const surfaceScalarField alphaRhoCpvPhi2
    (
      "alphaRhoCpvPhi2", fvc::interpolate(Cpv2)*alphaRhoPhi2
    );

    const surfaceScalarField rhoCpvPhi
    (   
        "rhoCpvPhi",
        alphaRhoCpvPhi1 + alphaRhoCpvPhi2
    );
     
    fvScalarMatrix TEqn
    (
      fvm::ddt(rhoCpv, T)
    + fvm::div(rhoCpvPhi, T)
    - fvm::Sp(fvc::ddt(rhoCpv) + fvc::div(rhoCpvPhi), T)
    - fvm::laplacian(thermophysicalTransport.kappaEff(), T)
    ==
      (e1Source&e1)
    + (e2Source&e2)
    );


 /////////////////////////////

/*
////////////////////////////////
  // getrennte Energie
 
  tmp<volScalarField> tAlphaRhoCpv1(alpha1*rho1*mixture.thermo1().Cpv() );
  const volScalarField& alphaRhoCpv1=tAlphaRhoCpv1();
  tmp<volScalarField> tAlphaRhoCpv2(alpha2*rho2*mixture.thermo2().Cpv() );
  const volScalarField& alphaRhoCpv2=tAlphaRhoCpv2();
  
  const surfaceScalarField alphaRhoCpvPhi1 = fvc::interpolate(alphaRhoCpv1)*phi;
  const surfaceScalarField alphaRhoCpvPhi2 = fvc::interpolate(alphaRhoCpv2)*phi;
  
  volScalarField& T = mixture.T();

  fvScalarMatrix TEqn
 (
  
  fvm::ddt(alphaRhoCpv1, T) +fvm::ddt(alphaRhoCpv2, T) +
  fvm::div(alphaRhoCpvPhi1, T) +  fvm::div(alphaRhoCpvPhi2, T) 

  ==
  //-
  fvm::laplacian(thermophysicalTransport.kappaEff(), T) 
  
 // ==
 //   (e1Source&e1)
 // + (e2Source&e2)
 );  
 
 /////////////////////////////
*/

// 
// 
// 
//
//  const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
//  const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));
//  const dimensionedScalar ZeroA ("ZeroA", rho1.dimensions()/dimTime, 0);
//
//  volScalarField& T = mixture.T();
//
//
////////////////////////
// const surfaceScalarField alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;
// const surfaceScalarField alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;
//  /*
// fvScalarMatrix TEqn
//  (
//     correction
//     (
//         mixture.thermo1().Cv()()
//         *(
//             fvm::ddt(alpha1, rho1, T) + fvm::div(alphaRhoPhi1, T)
//           - (
//                e1Source.hasDiag()
//               ? fvm::Sp(e1Source.A(), T)
//               : fvm::Sp(ZeroA, T)
//             )
//         )
//       + mixture.thermo2().Cv()()
//         *(
//             fvm::ddt(alpha2, rho2, T) + fvm::div(alphaRhoPhi2, T)
//           - (
//                 e2Source.hasDiag()
//                 ? fvm::Sp(e2Source.A(), T)
//                 : fvm::Sp(ZeroA, T)
//             )
//         )
//     )
//
//    + fvc::ddt(alpha1, rho1, e1) + fvc::div(alphaRhoPhi1, e1)
//    
//    + fvc::ddt(alpha2, rho2, e2) + fvc::div(alphaRhoPhi2, e2)
//    
//    - fvm::laplacian(thermophysicalTransport.kappaEff(), T)
//
//    + (
//          mixture.totalInternalEnergy()
//        ?
//          fvc::div(fvc::absolute(phi, U), p)()()
//        + (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))()()
//        - (U()&(fvModels().source(rho, U)&U)())
//        :
//          p*fvc::div(fvc::absolute(phi, U))()()
//      )
//   ==
//      (e1Source&e1)
//    + (e2Source&e2)
//  );
//*/
 
 

    TEqn.relax();

    fvConstraints().constrain(TEqn);

    TEqn.solve();

    fvConstraints().constrain(T);

    mixture.correctThermo();
    mixture.correct();

}

//void Foam::solvers::incompressibleVoFTC::energyPredictor()

// ************************************************************************* //
