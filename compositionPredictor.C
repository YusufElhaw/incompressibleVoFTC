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

#include "fvMatrix.H"  
#include "incompressibleVoFTC.H"
#include "fvcMeshPhi.H"
#include "fvcDdt.H"                    //#include "fvmDdt.H" in multiphaseEuler
#include "fvmDiv.H"                    //#include "fvcDiv.H" in multiphaseEuler
#include "fvmSup.H"                    //#include "fvcSup.H" in multiphaseEuler
#include "fvmLaplacian.H"
#include "zeroGradientFvPatchFields.H"
#include "VLEConstant.H"                

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleVoFTC::compositionPredictor()
{
 Info<< "compositionPredictor() called, Meldung from compositionPredictor.C " <<  nl;

 // --- Phase densities from thermo (rhoConst -> still a volScalarField) ---
   const volScalarField& rho1(mixture.thermo1().rho());
   const volScalarField& rho2(mixture.thermo2().rho());
 
 // --- Species mass-fraction fields (4 fields in 0/ and time/) ---
   volScalarField& Y1p1 = mixture.s1Phase1();  // species1.phase1
   volScalarField& Y2p1 = mixture.s2Phase1();  // species2.phase1
   volScalarField& Y1p2 = mixture.s1Phase2();  // species1.phase2
   volScalarField& Y2p2 = mixture.s2Phase2();  // species2.phase2
  // --- Species mol-weight ---
  const dimensionedScalar W1 = mixture.thermo1().Wi(0)/1e3;   // mol mass in kg
  const dimensionedScalar W2 = mixture.thermo1().Wi(1)/1e3;   // mol mass in kg
  
  // bounded alphas (ALPHA is 0<=alpha<=1)
    // alphas smaller than aTol are set to 0   
    const scalar aTol = 1e-4; 
    volScalarField a("a", min(max(alpha1, scalar(0)), scalar(1))); // Alpha1 only between 0 and 1
    volScalarField ALPHA1("ALPHA1",(scalar(1) - pos(aTol-a) - pos(a-(scalar(1)-aTol)))*a + pos(a-(scalar(1)-aTol)));
    volScalarField ALPHA2("ALPHA2", scalar(1) - ALPHA1);

  // Reading/calculating the concentration fields of the species  
    const word C1name = mixture.species1Name() + ".concentration";
    const word C2name = mixture.species2Name() + ".concentration";
    
    // taking the boundarys from the Mass fraction of Species1 phase 1 
    const wordList CpatchTypes = Y1p1.boundaryField().types();
    
    // creating the concentration fields
    if (!mesh.foundObject<volScalarField>(C1name))
    {
      mesh.objectRegistry::store
      (
        new volScalarField
        (
          IOobject
          (
            C1name,
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionedScalar("C1", rho1.dimensions()/W1.dimensions(), 0.0)
          ,CpatchTypes
        )
      );
    }

    if (!mesh.foundObject<volScalarField>(C2name))
    {
      mesh.objectRegistry::store
      (
        new volScalarField
        (
          IOobject
          (
            C2name,
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionedScalar("C2", rho2.dimensions()/W2.dimensions(), 0.0)
          ,CpatchTypes
        )
      );
    }
  
    volScalarField& C1 = const_cast<volScalarField&>(mesh.lookupObject<volScalarField>(C1name));
    volScalarField& C2 = const_cast<volScalarField&>(mesh.lookupObject<volScalarField>(C2name));
    
      // saving the BC for fixed Value with the value calculated from the mass fractions 
      volScalarField C1p1("C1p1", ALPHA1*(rho1/W1)*Y1p1);
      volScalarField C1p2("C1p2", ALPHA2*(rho2/W1)*Y1p2);
      volScalarField C2p1("C2p1", ALPHA1*(rho1/W2)*Y2p1);
      volScalarField C2p2("C2p2", ALPHA2*(rho2/W2)*Y2p2);

      volScalarField C1calc ("C1calc", C1p1 + C1p2);
      volScalarField C2calc ("C2calc", C2p1 + C2p2);

      // internal immer setzen
      C1.internalFieldRef() = C1calc.internalField();
      C2.internalFieldRef() = C2calc.internalField();
      // NUR fixedValue-Patches setzen (Wert aus Ccalc übernehmen)
      forAll(C1.boundaryField(), patchi)
      {
        if (C1.boundaryField()[patchi].type() == "fixedValue")
        {
          C1.boundaryFieldRef()[patchi] == C1calc.boundaryField()[patchi];
        }
      }

      forAll(C2.boundaryField(), patchi)
      {
        if (C2.boundaryField()[patchi].type() == "fixedValue")
        {
          C2.boundaryFieldRef()[patchi] == C2calc.boundaryField()[patchi];
        }
      }
      C1.correctBoundaryConditions();
      C2.correctBoundaryConditions();

    // intialisation of C Fields only once 
      // read C Fields, if found in the Time-Folder to be started with
      // if C Fields not found, calculate them based on mass fractions
      static bool CinitDone = false;  // intialisation is not done

      if (!CinitDone)
      {      
        static word startTimeName;  
        const scalar tStart = runTime.startTime().value();   // Start time from controlDict
        startTimeName = Time::timeName(tStart, 6);           // Precision: 6

        IOobject C1h(C1name, startTimeName, mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE);
        IOobject C2h(C2name, startTimeName, mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE);

        if (C1h.headerOk())
        {
          volScalarField C1read
          (
            IOobject(C1name, startTimeName, mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
            mesh
          );
          C1 = C1read;   // copy values into existing field
          Info<<C1name<<" read from the "<< startTimeName << " folder."<<endl;
        }
        else
        {
          C1 = C1p1 + C1p2;
          Info<<C1name<<" calculated from the mass fractions of "<<startTimeName <<" folder."<<endl;
        }

        if (C2h.headerOk())
        {
          volScalarField C2read
          (
              IOobject(C2name, startTimeName, mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
              mesh
          );
          C2 = C2read;
          Info<<C2name<<" read from the "<< startTimeName << " folder."<<endl;
        }
        else
        {
          C2 = C2p1 + C2p2;
          Info<<C2name<<" calculated from the mass fractions of "<<startTimeName <<" folder."<<endl;
        }
        CinitDone = true;   // intialisation is done
      }

  // VLE constant 
    const volScalarField& T = mixture.T();

    static autoPtr<Foam::VLEConstant> vle1Ptr;
    static autoPtr<Foam::VLEConstant> vle2Ptr;

    if (!vle1Ptr.valid())
    {
        vle1Ptr.reset
        (
            new Foam::VLEConstant
            (
              momentumTransport,
              mixture,
              mixture.species1Name(),
              "VLEProperties"
            )
        );
    }

    if (!vle2Ptr.valid())
    {
        vle2Ptr.reset
        (
            new Foam::VLEConstant
            (
              momentumTransport,
              mixture,
              mixture.species2Name(),
              "VLEProperties"
            )
        );
    }
    // K1 = K_Raoult * (C_gas/C_liquid) 
      //if the liquid phase is absent, C_liq=0 => dont devide by 0 
      // Dummy value of division by 1e3 causes no Problem because K1 & K2 only work on mixed cells
      const dimensionedScalar noDivBy0("noDivBy0", C1p1.dimensions(), 1 );
     
      volScalarField K1("K1", vle1Ptr->K(p, T) /* *  (C1p2+C2p2)/max((C1p1+C2p1), noDivBy0)*/);
      volScalarField K2("K2", vle2Ptr->K(p, T) /* *  (C1p2+C2p2)/max((C1p1+C2p1), noDivBy0)*/);

      Info<<"Relative volatility min / max = "<<min(K1/K2).value()<<" / "<<max(K1/K2).value()<<endl;

  // Diffusivity D_j (Start: alpha-gewichtet -erstamliger Wert) // to be changed
 
    volScalarField tDEff =  thermophysicalTransport.DEff();
 
  // Solver for C-Equations
    // Eq.10: Phi_j = -D * [ C(1-K)/(alpha1 + K*alpha2) ] grad(alpha1) // doi:10.1016/j.ces.2010.01.012
      volScalarField fraction1(C1*(scalar(1) - K1)/(ALPHA1 + K1*ALPHA2));
      volScalarField fraction2(C2*(scalar(1) - K2)/(ALPHA1 + K2*ALPHA2));
 
      volVectorField Phi1("Phi1", -tDEff*fraction1*fvc::grad(ALPHA1));
      volVectorField Phi2("Phi2", -tDEff*fraction2*fvc::grad(ALPHA1));  
    //  Info<<"Der K1Raoult beträgt :"<<vle1Ptr->K(p, T)()<<endl;
    //  Info<<"Der K1 beträgt :"<<K1()<<endl;
    //  Info<<"Der tDEff beträgt :"<<fvm::laplacian(tDEff, C2)<<endl;
    //  Info<<"Der Phi1 beträgt :"<<fvc::div(Phi1)<<endl;
 
    // Eq.9 :
    // ddt(C) + div(phi,C)  =  laplacian(D,C)+ div(Phi) // doi:10.1016/j.ces.2010.01.012
      fvScalarMatrix C1Eqn
      (
          fvm::ddt(C1)
        + fvm::div(phi, C1)
        ==
          fvm::laplacian(tDEff, C1)
        + fvc::div(Phi1) 
      );

      fvScalarMatrix C2Eqn
      (
          fvm::ddt(C2)
        + fvm::div(phi, C2)
        ==
          fvm::laplacian(tDEff, C2)
        + fvc::div(Phi2)
      );

      C1Eqn.relax(); fvConstraints().constrain(C1Eqn); C1Eqn.solve(/*mesh.solution().solverDict(C1.name())*/);
      C2Eqn.relax(); fvConstraints().constrain(C2Eqn); C2Eqn.solve(/*mesh.solution().solverDict(C2.name())*/);
      if (runTime.writeTime()) // writing the C-Fields for the writeTime
      {
          C1.write();
          C2.write();
      }
    
  // back substitution in Y with CiWi/sum(CjWj)
        
    // seperate between cells with pure phases and cells with mixed phases  
      volScalarField alpha1_pure("alpha1_pure", pos(ALPHA1 - (scalar(1) - aTol)));  // only cells with pure phases, mixed cells 0
      volScalarField alpha2_pure("alpha2_pure", pos(ALPHA2 - (scalar(1) - aTol)));  // only cells with pure phases, mixed cells 0 
      volScalarField alpha1_mix("alpha1_mix"  , pos(ALPHA1 - alpha1_pure)*(ALPHA1 - alpha1_pure)); // alpha1 for mixed cells, 0 for pure cells
      volScalarField alpha2_mix("alpha2_mix"  , pos(ALPHA2 - alpha2_pure)*(ALPHA2 - alpha2_pure)); // alpha2 for mixed cells, 0 for pure cells
      volScalarField Is_mix_Cell("Is_mix_Cell", pos(ALPHA1 - aTol)*pos(ALPHA2 - aTol)); // gives 1 for mix cells,
    
    // cells of pure phases  
      volScalarField Y1p1_pure("Y1p1_pure", alpha1_pure*C1*W1/(C1*W1+C2*W2));   // gives Y1p1 where alpha1 =1, everywhere else (0 as Dummy)
      volScalarField Y2p1_pure("Y2p1_pure", 1-Y1p1_pure);                       // gives Y2p1 where alpha1 =1, everywhere else (1 as Dummy)
      volScalarField Y2p2_pure("Y2p2_pure", alpha2_pure*C2*W2/(C1*W1+C2*W2));   // gives Y2p2 where alpha2 =1, everywhere else (0 as Dummy)
      volScalarField Y1p2_pure("Y1p2_pure", 1-Y2p2_pure);                       // gives Y1p2 where alpha2 =1, everywhere else (1 as Dummy)
    
    // cells of mixed phases (2-film Model -> Vapour liquid equilibrium on these cells)
    // concentration of the species for each phase
      volScalarField C1p1_mix("C1p1_mix",   Is_mix_Cell*(C1/max(alpha1_mix+K1*alpha2_mix, SMALL))); // Protection devision by 0 for the pure cells 
      volScalarField C1p2_mix("C1p2_mix",   K1*C1p1_mix);
      volScalarField C2p1_mix("C2p1_mix",   Is_mix_Cell*(C2/max(alpha1_mix+K2*alpha2_mix, SMALL))); // Protection devision by 0 for the pure cells
      volScalarField C2p2_mix("C2p2_mix",   K2*C2p1_mix);
    // Mass fraction of the mixed phases
      volScalarField Y1p1_mix("Y1p1_mix",   alpha1_mix*C1p1_mix*W1/(C1*W1+C2*W2));
      volScalarField Y1p2_mix("Y1p2_mix",   alpha2_mix*C1p2_mix*W1/(C1*W1+C2*W2));
      volScalarField Y2p1_mix("Y2p1_mix",   alpha1_mix*C2p1_mix*W2/(C1*W1+C2*W2));
      volScalarField Y2p2_mix("Y2p2_mix",   alpha2_mix*C2p2_mix*W2/(C1*W1+C2*W2));

//
    // Piecewise assembly
      // Y1p1 and Y2p2 have all the Field their real Values
      // Y2p1 and Y1p2 have 1 in cells, where their phases does not exist, as a Dummy-value
      // Y2p1 and Y1p2 have in the rest of the cells their real values
      Y1p1 = Y1p1_pure + Y1p1_mix ;                                 
      Y2p1 = Y2p1_pure + Y2p1_mix - Is_mix_Cell;      // subtraction of the Dummy 1 in the Y2p1_pure in the mix cells 
      Y1p2 = Y1p2_pure + Y1p2_mix - Is_mix_Cell;      // subtraction of the Dummy 1 in the Y1p2_pure in the mix cells
      Y2p2 = Y2p2_pure + Y2p2_mix;                                  
  
     // Info<<"the massfractions of Y1p1 are : "<<Y1p1<<nl<<endl;
     // Info<<"the massfractions of Y2p1 are : "<<Y2p1<<nl<<endl;
     // Info<<"the massfractions of Y1p2 are : "<<Y1p2<<nl<<endl;
     // Info<<"the massfractions of Y2p2 are : "<<Y2p2<<nl<<endl;
      // correct boundary conditions
    Y1p1.correctBoundaryConditions();
    Y2p1.correctBoundaryConditions();
    Y1p2.correctBoundaryConditions();
    Y2p2.correctBoundaryConditions();

   
  // Push mass-fractions into thermo objects + update derived properties
    mixture.correctComposition();
   // mixture.correctThermo();
   // mixture.correct();
}



///////////////////////////////////////////////

//void Foam::solvers::incompressibleVoFTC::energyPredictor()

// ************************************************************************* //
