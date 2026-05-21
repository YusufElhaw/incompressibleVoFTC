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
#include "fvcDdt.H"                  
#include "fvmDiv.H"                  
#include "fvmSup.H"                  
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "zeroGradientFvPatchFields.H"
#include "VLEConstant.H"          

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleVoFTC::compositionPredictor()
{
  // --- Phase densities from thermo (rhoConst -> still a volScalarField) ---
    const volScalarField& rho1(mixture.thermo1().rho());
    const volScalarField& rho2(mixture.thermo2().rho());
    //const volScalarField& rho(mixture.rho());
   
  // --- Species mass-fraction fields (4 fields in 0/ and time/) ---
    volScalarField& Y1L = mixture.s1Phase1();  // species1.Liquid
    volScalarField& Y2L = mixture.s2Phase1();  // species2.Liquid
    volScalarField& Y1G = mixture.s1Phase2();  // species1.Gas
    volScalarField& Y2G = mixture.s2Phase2();  // species2.Gas
  // --- Species mol-weight ---
    const dimensionedScalar W1 = mixture.thermo1().Wi(0)/1e3;   // mol mass in kg
    const dimensionedScalar W2 = mixture.thermo1().Wi(1)/1e3;   // mol mass in kg

  // bounded alphas (ALPHA is 0<=alpha<=1)
    // alphas smaller than aTol are set to 0   
    const scalar aTol = 1e-5; 
    volScalarField a("a", min(max(alpha1, scalar(0)), scalar(1))); // Alpha1 only between 0 and 1
    volScalarField ALPHA1("ALPHA1",(1 - pos(aTol-a) - pos(a-(1-aTol)))*a + pos(a-(1-aTol)));
    volScalarField ALPHA2("ALPHA2", 1 - ALPHA1);
    volScalarField Is_mix_Cell("Is_mix_Cell", pos(ALPHA1 - aTol)*pos(ALPHA2 - aTol)); // gives 1 for mix cells,

  // Concentrations
    volScalarField C1L("C1L", pos(alpha1-aTol) * (rho1/W1)*Y1L);  //multiplied with alpha in the concentration equation
    volScalarField C1G("C1G", pos(alpha2-aTol) * (rho2/W1)*Y1G);  //multiplied with alpha in the concentration equation
    volScalarField C2L("C2L", pos(alpha1-aTol) * (rho1/W2)*Y2L);  //multiplied with alpha in the concentration equation
    volScalarField C2G("C2G", pos(alpha2-aTol) * (rho2/W2)*Y2G);  //multiplied with alpha in the concentration equation
    
    volScalarField C1 ("C1", ALPHA1*C1L + ALPHA2*C1G); // Concentration of Species 1  
    volScalarField C2 ("C2", ALPHA1*C2L + ALPHA2*C2G); // Concentration of Species 2

    volScalarField CL ("CL", C1L + C2L); // Concentration of the liquid phase
    volScalarField CG ("CG", C1G + C2G); // Concentration of the gas phase

  
  // Reading/calculating the mol fractions fields of the species      
    const word X1name = mixture.species1Name() + ".X";
    const word X2name = mixture.species2Name() + ".X";

    // Use the actual start-time folder, not the current time step
    const word startTimeName = Time::timeName(runTime.startTime().value(), 6);

    // Default fallback BC only if no X-field exists in the start-time folder.
    // Do NOT copy inletOutlet only by type name, because inletOutlet needs dictionary entries.
    const wordList zeroGradientPatchTypes
    (
        mesh.boundary().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    auto readOrCreateX =
    [
        this,
        &startTimeName,
        &zeroGradientPatchTypes
    ]
    (
        const word& Xname,
        const volScalarField& Xcalc,
        const word& Xlabel
    )
    {
        if (mesh.foundObject<volScalarField>(Xname))
        {
            return;
        }

        IOobject Xheader
        (
            Xname,
            startTimeName,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );

        if (Xheader.headerOk())
        {
            mesh.objectRegistry::store
            (
                new volScalarField
                (
                    IOobject
                    (
                        Xname,
                        startTimeName,
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );

            Info<< Xname << " read from the " << startTimeName
                << " folder." << endl;
        }
        else
        {
            mesh.objectRegistry::store
            (
                new volScalarField
                (
                    IOobject
                    (
                        Xname,
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(Xlabel, dimless, 0.0),
                    zeroGradientPatchTypes
                )
            );

            volScalarField& X =
                const_cast<volScalarField&>
                (
                    mesh.lookupObject<volScalarField>(Xname)
                );

            X.internalFieldRef() = Xcalc.internalField();
            X.correctBoundaryConditions();

            Info<< Xname << " calculated from the mass fractions of "
                << startTimeName << " folder." << endl;
        }
    };
    // mol fractions
    volScalarField C("C", C1 + C2); // total concentration Field
    volScalarField X1Calc ("X1Calc", C1/C);     // mole fraction of Species 1  
    volScalarField X2Calc ("X2Calc", 1-X1Calc); // mole fraction of Species 2  

    readOrCreateX(X1name, X1Calc, "X1");
    readOrCreateX(X2name, X2Calc, "X2");

    volScalarField& X1 = const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>(X1name)
        );

    volScalarField& X2 = const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>(X2name)
        );
  // VLE constant 
    const volScalarField& T = mixture.T();
    const volScalarField& pVLE = mixture.p();

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

    volScalarField K1eq("K1eq", vle1Ptr->K(pVLE, T));
    volScalarField K2eq("K2eq", vle2Ptr->K(pVLE, T));

    // relative volatility A_ij = K1eq/K2eq
      tmp<volScalarField> tA12(K1eq/K2eq);

      if (mesh.foundObject<volScalarField>("A12"))
      {
          volScalarField& A12 =
              const_cast<volScalarField&>
              (
                  mesh.lookupObject<volScalarField>("A12")
              );

          A12 = tA12();
          A12.correctBoundaryConditions();
      }
      else
      {
          mesh.objectRegistry::store
          (
              new volScalarField
              (
                  IOobject
                  (
                      "A12",
                      runTime.name(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::AUTO_WRITE
                  ),
                  tA12()
              )
          );
      }

      volScalarField& A12 = const_cast<volScalarField&>
      (
        mesh.lookupObject<volScalarField>("A12")
      );
      
      Info<<"Relative volatility min / max = "<< min(A12).value()
          << " / " << max(A12).value() << endl;
    volScalarField X1L_lag("X1L_lag", pos(ALPHA1)*C1L/max(CL, dimensionedScalar("SMALL", CL.dimensions(), SMALL))); // lagged value of the mole fraction of species 1 in the liquid phase, used for the first guess of Kgcst
    //volScalarField X1G_lag("X1G_lag", pos(ALPHA2)*C1G/max(CG, dimensionedScalar("SMALL", CG.dimensions(), SMALL))); // lagged value of the mole fraction of species 1 in the gas phase, used for the first guess of Kgcst

    volScalarField Kgcst("Kgcst",A12/(1+(A12-1)*X1L_lag));

  // ---------------------------------------------------------------------
  // Concentration Equation GCST 
  // ---------------------------------------------------------------------  
    // Number of Picard Iterations
      const label nCompositionPicard
      (
        max(scalar(1), pimple.dict().lookupOrDefault<label>("compositionPicardIterations", 2))
      );

    volScalarField fractionL((ALPHA1*CL + ALPHA2*CG)/(ALPHA1*CL + Kgcst*ALPHA2*CG));
    volScalarField fractionG(Kgcst*(ALPHA1*CL + ALPHA2*CG)/(ALPHA1*CL + Kgcst*ALPHA2*CG));

    volScalarField D1Eff("D1Eff", ALPHA1*CL*thermophysicalTransport.D1Eff());
    volScalarField D2Eff("D2Eff", ALPHA2*CG*thermophysicalTransport.D2Eff());

    volScalarField D1fL("D1fL", D1Eff*fractionL);
    volScalarField D2fG("D2fG", D2Eff*fractionG);

    surfaceScalarField alphaCLPhi1
    (
        "alphaCLPhi1", alphaPhi1*fvc::interpolate(CL)
    );

    surfaceScalarField alphaCGPhi2
    (
        "alphaCGPhi2", alphaPhi2*fvc::interpolate(CG)
    );
    // correction of the numerical speed instabilties, which violates continuity
      volScalarField contErrCL1 
      (
          "contErrCL1", fvc::ddt(alpha1, CL) + fvc::div(alphaCLPhi1)
      );

      volScalarField contErrCG2 
      (
          "contErrCG2", fvc::ddt(alpha2, CG) + fvc::div(alphaCGPhi2)
      );

    for (label picardIter = 0; picardIter < nCompositionPicard; ++picardIter)
      {
          Info<< "compositionPredictor Picard iteration "
              << picardIter + 1 << "/" << nCompositionPicard << endl;
        
       // Kgcst=A12/(1+(A12-1)*X1);
        Kgcst = A12/(1 + (A12 - 1)*X1L_lag);
        fractionL =  (ALPHA1*CL + ALPHA2*CG)
                    /(ALPHA1*CL + Kgcst*ALPHA2*CG);
        fractionG =  Kgcst*(ALPHA1*CL + ALPHA2*CG)
                    /(ALPHA1*CL + Kgcst*ALPHA2*CG);

        D1fL = D1Eff*fractionL;
        D2fG = D2Eff*fractionG;


        volVectorField J1CorrCoeff
        (
            "J1CorrCoeff", D1Eff*fvc::grad(fractionL)
        );
        
        volVectorField J2CorrCoeff
        (
            "J2CorrCoeff", D2Eff*fvc::grad(fractionG)
        );
        
        surfaceScalarField J1CorrPhi
        (
            "J1CorrPhi", fvc::interpolate(J1CorrCoeff) & mesh.Sf()
        );
        
        surfaceScalarField J2CorrPhi
        (
            "J2CorrPhi", fvc::interpolate(J2CorrCoeff) & mesh.Sf()
        );


        fvScalarMatrix X1Eqn
        (
            correction
            (
              fvm::ddt(alpha1, CL, X1) + fvm::div(alphaCLPhi1, X1)
            - fvm::Sp(contErrCL1, X1)
            
            + fvm::ddt(alpha2, CG, X1) + fvm::div(alphaCGPhi2, X1)
            - fvm::Sp(contErrCG2, X1)
            )
            +  fvc::ddt(alpha1, CL, X1) + fvc::div(alphaCLPhi1, X1) 
            - contErrCL1*X1
            + fvc::ddt(alpha2, CG, X1) + fvc::div(alphaCGPhi2, X1) 
            - contErrCG2*X1
          ==  
            correction
            (
                  fvm::laplacian(D1fL, X1)
                + fvm::laplacian(D2fG, X1)
                + fvm::div(J1CorrPhi, X1)
                + fvm::div(J2CorrPhi, X1)
            )
          + fvc::laplacian(D1fL, X1)
          + fvc::laplacian(D2fG, X1)
          + fvc::div(J1CorrPhi, X1)
          + fvc::div(J2CorrPhi, X1)
        );

        X1Eqn.relax();
        fvConstraints().constrain(X1Eqn);
        X1Eqn.solve();
        fvConstraints().constrain(X1);


        X1 = min(max(X1, scalar(0)), scalar(1));
        X1.correctBoundaryConditions();
        volScalarField X1L_new("X1L_new",X1*fractionL);

        X1L_lag =
          Is_mix_Cell*X1L_new
        + (1 - Is_mix_Cell)*X1L_lag;

        X1L_lag.correctBoundaryConditions();
        
      }

      X2 = 1 - X1;
      X2.correctBoundaryConditions();
      
    if (runTime.writeTime())
    {
        X1.write();
        X2.write();
    }

    Kgcst = A12/(1 + (A12 - 1)*X1L_lag);

    fractionL =
        (ALPHA1*CL + ALPHA2*CG)
      /(ALPHA1*CL + Kgcst*ALPHA2*CG);

    fractionG =
        Kgcst*(ALPHA1*CL + ALPHA2*CG)
      /(ALPHA1*CL + Kgcst*ALPHA2*CG);
  // ---------------------------------------------------------------------
  // back-substitution to  mass fractions
  // ---------------------------------------------------------------------
   
    volScalarField ALPHA1_pure("ALPHA1_pure", pos(ALPHA1 - (1 - aTol)));  // only cells with pure phases, mixed cells 0
    volScalarField ALPHA2_pure("ALPHA2_pure", pos(ALPHA2 - (1 - aTol)));  // only cells with pure phases, mixed cells 0 
    volScalarField ALPHA1_mix("ALPHA1_mix"  , pos(ALPHA1 - ALPHA1_pure)*(ALPHA1 - ALPHA1_pure)); // alpha1 for mixed cells 
    volScalarField ALPHA2_mix("ALPHA2_mix"  , pos(ALPHA2 - ALPHA2_pure)*(ALPHA2 - ALPHA2_pure)); // alpha2 for mixed cells 
    
    const dimensionedScalar denomMinW ("denomMinW",W1.dimensions(), SMALL);
    
    volScalarField Y1L_pure("Y1L_pure", ALPHA1_pure*X1*W1/max((X1*W1+X2*W2),denomMinW));   // gives Y1L where alpha1 =1, everywhere else (0 as Dummy)
    volScalarField Y2L_pure("Y2L_pure", 1-Y1L_pure);                       // gives Y2L where alpha1 =1, everywhere else (1 as Dummy)
    volScalarField Y2G_pure("Y2G_pure", ALPHA2_pure*X2*W2/max((X1*W1+X2*W2),denomMinW));   // gives Y2G where alpha2 =1, everywhere else (0 as Dummy)
    volScalarField Y1G_pure("Y1G_pure", 1-Y2G_pure);                       // gives Y1G where alpha2 =1, everywhere else (1 as Dummy)

    const dimensionedScalar denomMinC ("denomMinC",CL.dimensions(), SMALL);
    

    // cells of mixed phases (2-film Model -> Vapour liquid equilibrium on these cells)
    // concentration of the species for each phase
      volScalarField X1L=X1*(ALPHA1_mix*CL+ALPHA2_mix*CG)/max((ALPHA1_mix*CL+Kgcst*ALPHA2_mix*CG), denomMinC);
      volScalarField X1G=X1L*Kgcst;
      volScalarField X2L=Is_mix_Cell*(1-X1L);
      volScalarField X2G=Is_mix_Cell*(1-X1G);

    // Mass fraction of the mixed phases
      volScalarField Y1L_mix("Y1L_mix",   Is_mix_Cell*X1L*W1/max((X1L*W1+X2L*W2),denomMinW));
      volScalarField Y1G_mix("Y1G_mix",   Is_mix_Cell*X1G*W1/max((X1G*W1+X2G*W2),denomMinW));
      volScalarField Y2L_mix("Y2L_mix",   Is_mix_Cell*(1-Y1L_mix));
      volScalarField Y2G_mix("Y2G_mix",   Is_mix_Cell*(1-Y1G_mix));


    // Piecewise assembly
      // Y1L and Y2G have all the Field their real Values
      // Y2L and Y1G have 1 in cells, where their phases does not exist, as a Dummy-value
      // Y2L and Y1G have in the rest of the cells their real values
      Y1L = Y1L_pure + Y1L_mix ;                   // has a dummy of 1 in the pure Gas cells              
      Y2L = Y2L_pure + Y2L_mix - Is_mix_Cell;      // subtraction of the Dummy 1 in the Y2L_pure in the mix cells 
      Y1G = Y1G_pure + Y1G_mix - Is_mix_Cell;      // subtraction of the Dummy 1 in the Y1G_pure in the mix cells
      Y2G = Y2G_pure + Y2G_mix;                    // has a dummy of 1 in the pure Liquid cells   


          
      Y1L.correctBoundaryConditions();
      Y2L.correctBoundaryConditions();
      Y1G.correctBoundaryConditions();
      Y2G.correctBoundaryConditions();

    // Test mole conservation for post processing
      nMoles1 =
          alpha1*mixture.thermo1().rho()*Y1L/W1
        + alpha2*mixture.thermo2().rho()*Y1G/W1;

      nMoles2 =
          alpha1*mixture.thermo1().rho()*Y2L/W2
        + alpha2*mixture.thermo2().rho()*Y2G/W2;

      nMolesTotal = nMoles1 + nMoles2;

      nMoles1.correctBoundaryConditions();
      nMoles2.correctBoundaryConditions();
      nMolesTotal.correctBoundaryConditions();
  
  // Push mass-fractions into thermo objects + update derived properties
    mixture.correctComposition();
    mixture.correctThermo(p);
    mixture.correct();
    p_rgh_ = p - rho*buoyancy.gh;
}

// ************************************************************************* //
