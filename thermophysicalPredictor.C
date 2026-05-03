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
  Info <<nl<< "thermophysicalPredictor läuft"<<endl;

  const volScalarField& rho1(mixture.thermo1().rho()); 
  const volScalarField& rho2(mixture.thermo2().rho());
  //const volScalarField& rho(mixture.rho());

  const volScalarField& Cpv1(mixture.thermo1().Cpv());
  const volScalarField& Cpv2(mixture.thermo2().Cpv());
   //Info<< "Cpv1 beträgt: " << Cpv1.internalField()<<nl<<nl<< " Cpv2 beträgt: " << Cpv2.internalField() << endl;
  // Energy of the phases (in physicalProperties Choose Enthalpy or InternalEnergy)
  const volScalarField& e1(mixture.thermo1().he());   //he Enthalpy or InternalEnergy
  const volScalarField& e2(mixture.thermo2().he());   //he Enthalpy or InternalEnergy

  const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
  const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));

  volScalarField& T = mixture.T();
  K = 0.5*magSqr(U);// Update the kinetic energy field with the current velocity field


 // const volScalarField rhoCpv1("rhoCpv1", rho1*Cpv1); 
 // const volScalarField rhoCpv2("rhoCpv2", rho2*Cpv2);
 // const surfaceScalarField alphaRhoCpvPhi1("alphaRhoCpv1", fvc::interpolate(rhoCpv1)*alphaPhi1);
 // const surfaceScalarField alphaRhoCpvPhi2("alphaRhoCpv2", fvc::interpolate(rhoCpv2)*alphaPhi2);

  fvScalarMatrix TEqn
(
    correction
    (
        Cpv1*
        (
            fvm::ddt(alpha1, rho1, T)
          + fvm::div(alphaRhoPhi1, T)
          - (
                e1Source.hasDiag()
              ? fvm::Sp(contErr1(), T) + fvm::Sp(e1Source.A(), T)
              : fvm::Sp(contErr1(), T)
            )
        )

      + Cpv2*
        (
            fvm::ddt(alpha2, rho2, T)
          + fvm::div(alphaRhoPhi2, T)
          - (
                e2Source.hasDiag()
              ? fvm::Sp(contErr2(), T) + fvm::Sp(e2Source.A(), T)
              : fvm::Sp(contErr2(), T)
            )
        )
    )

  + fvc::ddt(alpha1, rho1, e1)
  + fvc::div(alphaRhoPhi1, e1)
  - contErr1()*e1

  + fvc::ddt(alpha2, rho2, e2)
  + fvc::div(alphaRhoPhi2, e2)
  - contErr2()*e2

  - fvm::laplacian(thermophysicalTransport.kappaEff(), T)

  + (
        mixture.totalInternalEnergy()
      ?
        fvc::div(fvc::absolute(phi, U), p)()()
      + (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))()()
      - (U()&(fvModels().source(rho, U)&U)())
      - (contErr1() + contErr2())*K
      :
        p*fvc::div(fvc::absolute(phi, U))()()
    )
 ==
      (e1Source & e1)
    + (e2Source & e2)
);

  TEqn.relax();

  fvConstraints().constrain(TEqn);

  TEqn.solve();

  fvConstraints().constrain(T);

  mixture.correctThermo();
  mixture.correct();

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*  // one Field Methode // works but subtracts some energy for stabilisation
  
  volScalarField rhoCp
  (
      IOobject
      (
          "rhoCp",
          runTime.name(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
      alpha1*rho1*Cpv1 + alpha2*rho2*Cpv2
  );

  surfaceScalarField CpvRhoPhi
  (
      IOobject
      (
          "CpvRhoPhi",
          runTime.name(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
      fvc::interpolate(rhoCp/rho)*rhoPhi
  );

  fvScalarMatrix TEqn
  (
      fvm::ddt(rhoCp, T)        // =T*ddt(rhoCp) + rhoCp*ddt(T)
    + fvm::div(CpvRhoPhi, T)    // =T*fvc::div(CpvRhoPhi) + CpvRhoPhi&fvc::grad(T)
    - fvm::laplacian(thermophysicalTransport.kappaEff(), T)
    - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(CpvRhoPhi), T) // subtract T*(ddt(rhoCp) + div(CpvRhoPhi)) for stabilisation
  );


  TEqn.relax();
  fvConstraints().constrain(TEqn);
  TEqn.solve();
  fvConstraints().constrain(T);

  mixture.correctThermo();
  mixture.correct();
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Nach Mehdi Nabil
//    const dimensionedScalar T0("T0", dimTemperature, 298.15);
//
//
//  const volScalarField limAlpha1( min(max(alpha1, scalar(0)), scalar(1)) );
//  const wordList TpatchTypes = T.boundaryField().types();
//  surfaceScalarField alphaEffRho
//  (
//      IOobject
//      (
//          "alphaEffRho",
//          mesh.time().name(),
//          mesh,
//          IOobject::NO_READ,
//          IOobject::NO_WRITE
//      ),
//      fvc::interpolate(rho1)*
//      (
//        fvc::interpolate(thermophysicalTransport.kappaEff()/(Cpv*rho))/*+fvc::interpolate(turbulence_.turbulence_->nut())*/  
//      )
//    + fvc::interpolate(rho2)*
//      (
//        fvc::interpolate(thermophysicalTransport.kappaEff()/(Cpv*rho))/*+fvc::interpolate(turbulence_.turbulence_->nut())*/ 
//      )
//  );
//  //    
//
//  volScalarField H
//    (
//      IOobject
//      (
//        "H",
//        runTime.name(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//      ),
//      mesh,
//      dimensionedScalar("H", dimEnergy/dimMass,0),
//      T.boundaryField().types()
//    );
//  
//  H = ( (T - T0)*(limAlpha1*rho1*Cpv1 + (1.0 - limAlpha1)*rho2*Cpv2) )/rho;
//
//  H.correctBoundaryConditions();
//
//  const scalar RelaxFac = 1.0;
//  fvScalarMatrix EEqn
//    (
//        fvm::ddt(rho, H)
//      + fvm::div(rhoPhi, H)
//      ==
//        fvc::laplacian(thermophysicalTransport.kappaEff(), T)
//      + RelaxFac*
//      (   
//        fvm::laplacian(alphaEffRho, H) 
//      - fvc::laplacian(alphaEffRho, H) 
//      ) 
//            
//    );
//      EEqn.relax();
//      fvConstraints().constrain(EEqn);
//      EEqn.solve();
//      //Now reevaluate T for the updated enthalpy fields
//      T = T0 + rho*H/(limAlpha1*rho1*Cpv1 + (1.0 - limAlpha1)*rho2*Cpv2);
//      mixture.correctThermo();
//      mixture.correct();
//
//
//

  // Test energy conservation
  // Update total thermal energy field 
    Etherm = alpha1*mixture.thermo1().rho()*mixture.thermo1().he()
           + alpha2*mixture.thermo2().rho()*mixture.thermo2().he();
    Etherm.correctBoundaryConditions();
}

//void Foam::solvers::incompressibleVoFTC::energyPredictor()

// ************************************************************************* //
