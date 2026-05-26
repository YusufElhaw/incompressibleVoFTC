/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    GCST-style composition predictor for incompressibleVoFTC.

    Core structure:

        LHS:
            old stable phase-separated GCST advection

            ddt(alpha1, CL, X)
          + div(alphaCLPhi1, X)
          - Sp(contErrCL1, X)

          + ddt(alpha2, CG, X)
          + div(alphaCGPhi2, X)
          - Sp(contErrCG2, X)

        RHS:
            new algebraic VLE diffusion

            div(D1Base grad(X1L_alg))
          + div(D2Base grad(X1G_alg))

        with algebraic split in mixed cells:

            q = M*X
            q = aL*X1L + aG*X1G
            X1G = A12*X1L/(1 + (A12 - 1)*X1L)


\*---------------------------------------------------------------------------*/

#include "fvMatrix.H"
#include "incompressibleVoFTC.H"

#include "fvcMeshPhi.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "zeroGradientFvPatchFields.H"

#include "VLEConstant.H"
#include "inletOutletCondenserFvPatchScalarField.H"
#include "inletOutletBoilerFvPatchScalarField.H"
#include "nonConformalProcessorCyclicFvPatch.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleVoFTC::compositionPredictor()
{
    // ---------------------------------------------------------------------
    // Phase densities and species mass-fraction fields
    // ---------------------------------------------------------------------

        const volScalarField& rho1 = mixture.thermo1().rho();
        const volScalarField& rho2 = mixture.thermo2().rho();

        volScalarField& Y1L = mixture.s1Phase1();
        volScalarField& Y2L = mixture.s2Phase1();
        volScalarField& Y1G = mixture.s1Phase2();
        volScalarField& Y2G = mixture.s2Phase2();

        const dimensionedScalar W1 = mixture.thermo1().Wi(0)/1e3;
        const dimensionedScalar W2 = mixture.thermo1().Wi(1)/1e3;

    // ---------------------------------------------------------------------
    // Alpha masks
    // ---------------------------------------------------------------------

        const scalar aTol = 1e-9;

        volScalarField a("a",min(max(alpha1, scalar(0)), scalar(1)));
        volScalarField ALPHA1("ALPHA1", (scalar(1) - pos(aTol - a) - pos(a - (scalar(1) - aTol)))*a+ pos(a - (scalar(1) - aTol)));
        volScalarField ALPHA2("ALPHA2", scalar(1) - ALPHA1);

        volScalarField Is_mix_Cell("Is_mix_Cell",pos(ALPHA1 - aTol)*pos(ALPHA2 - aTol));

    // ---------------------------------------------------------------------
    // Phase molar concentrations
    // ---------------------------------------------------------------------

        volScalarField C1L("C1L", pos(alpha1 - aTol)*(rho1/W1)*Y1L);
        volScalarField C1G("C1G", pos(alpha2 - aTol)*(rho2/W1)*Y1G);
        volScalarField C2L("C2L", pos(alpha1 - aTol)*(rho1/W2)*Y2L);
        volScalarField C2G("C2G", pos(alpha2 - aTol)*(rho2/W2)*Y2G);

        //volScalarField CL("CL", C1L + C2L);
        //volScalarField CG("CG", C1G + C2G);
        if (!mesh.foundObject<volScalarField>("CL"))
        {
            mesh.objectRegistry::store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "CL",
                        runTime.name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    C1L + C2L
                )
            );
        }
        else
        {
            volScalarField& CLstored =
                const_cast<volScalarField&>
                (
                    mesh.lookupObject<volScalarField>("CL")
                );

            CLstored = C1L + C2L;
            CLstored.correctBoundaryConditions();
        }

        if (!mesh.foundObject<volScalarField>("CG"))
        {
            mesh.objectRegistry::store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "CG",
                        runTime.name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    C1G + C2G
                )
            );
        }
        else
        {
            volScalarField& CGstored =
                const_cast<volScalarField&>
                (
                    mesh.lookupObject<volScalarField>("CG")
                );

            CGstored = C1G + C2G;
            CGstored.correctBoundaryConditions();
        }

        const volScalarField& CL = mesh.lookupObject<volScalarField>("CL");
        const volScalarField& CG = mesh.lookupObject<volScalarField>("CG");
        // Molar inventory coefficients
        volScalarField aL("aL", ALPHA1*CL);
        volScalarField aG("aG", ALPHA2*CG);

        volScalarField M("M", aL + aG);

        const dimensionedScalar Mmin("Mmin", M.dimensions(),SMALL);

        // Bulk species inventory for initial X calculation
        volScalarField C1("C1",ALPHA1*C1L + ALPHA2*C1G);
        volScalarField X1Calc("X1Calc",C1/max(M, Mmin));
        volScalarField X2Calc("X2Calc",scalar(1) - X1Calc);

    // ---------------------------------------------------------------------
    // Read/create mole-fraction fields
    // ---------------------------------------------------------------------

        const word X1name = mixture.species1Name() + ".X";
        const word X2name = mixture.species2Name() + ".X";

        const word startTimeName = Time::timeName(runTime.startTime().value(), 6);

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

        readOrCreateX(X1name, X1Calc, "X1");
        readOrCreateX(X2name, X2Calc, "X2");

        volScalarField& X1 =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(X1name)
            );

        volScalarField& X2 =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(X2name)
            );

        X1.correctBoundaryConditions();

    // ---------------------------------------------------------------------
    // VLE and relative volatility A12
    // ---------------------------------------------------------------------

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

        if (!mesh.foundObject<volScalarField>("A12"))
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
                    K1eq/K2eq
                )
            );
        }

        volScalarField& A12 =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>("A12")
            );

        A12 = K1eq/K2eq;
        A12.correctBoundaryConditions();

        Info<< "Relative volatility min / max = "
            << min(A12).value()
            << " / "
            << max(A12).value()
            << endl;

    // ---------------------------------------------------------------------
    // Algebraic VLE potential for diffusion only
    // ---------------------------------------------------------------------
    //
    // X1 is still the old transported bulk/effective field.
    //
    // For diffusion, construct the liquid-side equilibrium coordinate xEq:
    //
    // Pure liquid:
    //     X1 = x
    //     xEq = X1
    //
    // Pure gas:
    //     X1 = y
    //     y = A*x/(1 + (A - 1)*x)
    //     xEq = inverseVLE(y)
    //
    // Mixed/interface cell:
    //     q = M*X1 = aL*x + aG*y
    //     y = A*x/(1 + (A - 1)*x)
    //
    // Diffusion flux:
    //     J = -[aL*D_L + aG*D_G*(dy/dx)] grad(xEq)
    //
    // Therefore steady pure diffusion gives uniform xEq,
    // and the VLE jump is independent of D.
    // ---------------------------------------------------------------------

        volScalarField xEq_alg("xEq_alg", X1);
        volScalarField yEq_alg("yEq_alg", X1);

        volScalarField dxEqdX("dxEqdX", X1*scalar(0) + scalar(1));
        volScalarField dYdXeq("dYdXeq", X1*scalar(0) + scalar(1));

        // Single-point VLE split: given bulk mole fraction X_val and phase
        // volume fractions aa (liquid), bb (gas), total moles mm, relative
        // volatility A, compute equilibrium liquid x, gas y, dy/dx, and
        // chain-rule factor dx/dX used for the diffusion Laplacian.
        auto vlePoint =
        [](
            scalar X_val, scalar aa, scalar bb, scalar mm, scalar A,
            scalar& xEq_out, scalar& yEq_out,
            scalar& dydx_out, scalar& dxDX_out
        )
        {
            if (aa > SMALL && bb > SMALL && mm > SMALL)
            {
                // Mixed/interface cell:
                // q = aa*x + bb*y,  y = A*x/(1+(A-1)*x)
                const scalar q  = mm*X_val;
                const scalar c2 = aa*(A - scalar(1));
                const scalar c1 = aa + bb*A - q*(A - scalar(1));

                scalar x = X_val;
                if (mag(c2) > VSMALL)
                {
                    const scalar disc = c1*c1 + scalar(4)*c2*q;
                    x = (-c1 + sqrt(max(disc, scalar(0))))/(scalar(2)*c2);
                }
                else
                {
                    x = q/(aa + bb*A);
                }

                const scalar den  = scalar(1) + (A - scalar(1))*x;
                const scalar y    = A*x/den;
                const scalar s    = A/sqr(den);   // dy/dx
                const scalar dqdx = aa + bb*s;

                xEq_out  = x;
                yEq_out  = y;
                dydx_out = s;
                dxDX_out = mag(dqdx) > VSMALL ? mm/dqdx : scalar(1);
            }
            else if (aa > SMALL && mm > SMALL)
            {
                // Pure liquid: X1 is liquid composition x
                const scalar den = scalar(1) + (A - scalar(1))*X_val;
                const scalar s   = A/sqr(den);

                xEq_out  = X_val;
                yEq_out  = A*X_val/den;
                dydx_out = s;
                dxDX_out = scalar(1);
            }
            else if (bb > SMALL && mm > SMALL)
            {
                // Pure gas: X1 is gas composition y; invert VLE to get x
                const scalar denomInv = A - (A - scalar(1))*X_val;
                const scalar x        = X_val/max(denomInv, VSMALL);
                const scalar den      = scalar(1) + (A - scalar(1))*x;
                const scalar s        = A/sqr(den);   // dy/dx at x

                xEq_out  = x;
                yEq_out  = X_val;
                dydx_out = s;
                // x = inverseVLE(y), so dx/dX = dx/dy = 1/(dy/dx)
                dxDX_out = scalar(1)/max(s, VSMALL);
            }
            else
            {
                xEq_out  = X_val;
                yEq_out  = X_val;
                dydx_out = scalar(1);
                dxDX_out = scalar(1);
            }
        };

        auto reconstructVLEPotential =
        [
            &X1,
            &xEq_alg,
            &yEq_alg,
            &dxEqdX,
            &dYdXeq,
            &aL,
            &aG,
            &M,
            &A12,
            &vlePoint
        ]()
        {
            scalarField& xEq  = xEq_alg.primitiveFieldRef();
            scalarField& yEq  = yEq_alg.primitiveFieldRef();
            scalarField& dxDX = dxEqdX.primitiveFieldRef();
            scalarField& dydx = dYdXeq.primitiveFieldRef();

            const scalarField& X   = X1.primitiveField();
            const scalarField& aLI = aL.primitiveField();
            const scalarField& aGI = aG.primitiveField();
            const scalarField& MI  = M.primitiveField();
            const scalarField& AI  = A12.primitiveField();

            forAll(X, celli)
            {
                vlePoint
                (
                    X[celli], aLI[celli], aGI[celli], MI[celli], AI[celli],
                    xEq[celli], yEq[celli], dydx[celli], dxDX[celli]
                );
            }

            forAll(X1.boundaryField(), patchi)
            {
                scalarField& xEqp  = xEq_alg.boundaryFieldRef()[patchi];
                scalarField& yEqp  = yEq_alg.boundaryFieldRef()[patchi];
                scalarField& dxDXp = dxEqdX.boundaryFieldRef()[patchi];
                scalarField& dydxp = dYdXeq.boundaryFieldRef()[patchi];

                const fvPatchScalarField& Xp  = X1.boundaryField()[patchi];
                const fvPatchScalarField& aLp = aL.boundaryField()[patchi];
                const fvPatchScalarField& aGp = aG.boundaryField()[patchi];
                const fvPatchScalarField& Mp  = M.boundaryField()[patchi];
                const fvPatchScalarField& Ap  = A12.boundaryField()[patchi];

                forAll(Xp, facei)
                {
                    vlePoint
                    (
                        Xp[facei], aLp[facei], aGp[facei], Mp[facei], Ap[facei],
                        xEqp[facei], yEqp[facei], dydxp[facei], dxDXp[facei]
                    );
                }
            }

            xEq_alg.correctBoundaryConditions();
            yEq_alg.correctBoundaryConditions();
            dxEqdX.correctBoundaryConditions();
            dYdXeq.correctBoundaryConditions();
        };

        reconstructVLEPotential();

    // ---------------------------------------------------------------------
    // Phase molar fluxes
    // ---------------------------------------------------------------------

        surfaceScalarField alphaCLPhi1("alphaCLPhi1", alphaPhi1*fvc::interpolate(CL));
        surfaceScalarField alphaCGPhi2("alphaCGPhi2", alphaPhi2*fvc::interpolate(CG));
        // Make boundary molar fluxes consistent with the final corrected U.
        // This is essential for conservative rectification inletOutlet BCs.
        // For nonConformalProcessorCyclic (NCPC) patches: apply upwind-correct
        // molar density. alphaPhi*fvc::interpolate(CL/CG) uses the neighbour
        // boundary value (MPI-received), which is WRONG for OUTFLOW faces where
        // the exiting fluid carries the OWN cell's molar concentration.
        // Use own-cell CL/CG for outflow (phi>0), neighbour CL/CG for inflow
        // (phi<0 — incoming fluid comes from the cyclic neighbour side).
        // Do NOT use U & Sf here (wrong face-normal convention for NCC geometry).
        forAll(mesh.boundary(), patchi)
        {
            if (isA<nonConformalProcessorCyclicFvPatch>(mesh.boundary()[patchi]))
            {
                const fvsPatchScalarField& alphaPhi1p =
                    alphaPhi1.boundaryField()[patchi];
                const fvsPatchScalarField& alphaPhi2p =
                    alphaPhi2.boundaryField()[patchi];

                const scalarField CLp_own(CL.boundaryField()[patchi].patchInternalField());
                const scalarField CGp_own(CG.boundaryField()[patchi].patchInternalField());
                const fvPatchScalarField& CLp_nbr = CL.boundaryField()[patchi];
                const fvPatchScalarField& CGp_nbr = CG.boundaryField()[patchi];

                fvsPatchScalarField& alphaCLPhi1p =
                    alphaCLPhi1.boundaryFieldRef()[patchi];
                fvsPatchScalarField& alphaCGPhi2p =
                    alphaCGPhi2.boundaryFieldRef()[patchi];

                // NCPC patches come in two flavours:
                //
                //   Owner-side (boiler, proc 0/1):  standard convention
                //     phi1 > 0 = liquid OUTFLOW to boiler reservoir
                //     phi2 < 0 = gas   INFLOW  from boiler reservoir
                //
                //   Neighbour-side (condenser, proc 2/3):  INVERTED convention
                //     phi1 > 0 = liquid INFLOW  from condenser reservoir
                //     phi2 < 0 = gas   OUTFLOW  to condenser reservoir
                //
                // For the condenser side we must NEGATE alphaCL/CGPhi so that the
                // FVM assembler (which always uses phi>0=outflow) sees the correct
                // direction, and applyInletValuesOnNCPCPatches() can supply xD at
                // the now-negative (inflow) liquid faces.
                //
                // Global mole conservation requires:
                //   liquid CL: use boiler (NBR at condenser) so both sides of the
                // All NCPC patches: standard upwind (own for outflow, nbr for inflow).
                // alphaCLPhi1 keeps the raw phi sign so measureCompositionPredictorFluxes
                // in the BC classes can read it with the expected convention.
                forAll(alphaCLPhi1p, facei)
                {
                    const scalar phi1 = alphaPhi1p[facei];
                    const scalar phi2 = alphaPhi2p[facei];
                    alphaCLPhi1p[facei] =
                        max(phi1, scalar(0))*CLp_own[facei]
                    + min(phi1, scalar(0))*CLp_nbr[facei];
                    alphaCGPhi2p[facei] =
                        max(phi2, scalar(0))*CGp_own[facei]
                    + min(phi2, scalar(0))*CGp_nbr[facei];
                }
                continue;
            }

            const fvsPatchVectorField& Sfp =
                mesh.Sf().boundaryField()[patchi];

            const fvPatchVectorField& Up =
                U.boundaryField()[patchi];

            const fvPatchScalarField& alpha1p =
                alpha1.boundaryField()[patchi];

            const fvPatchScalarField& alpha2p =
                alpha2.boundaryField()[patchi];

            const fvPatchScalarField& CLp =
                CL.boundaryField()[patchi];

            const fvPatchScalarField& CGp =
                CG.boundaryField()[patchi];

            fvsPatchScalarField& alphaCLPhi1p =
                alphaCLPhi1.boundaryFieldRef()[patchi];

            fvsPatchScalarField& alphaCGPhi2p =
                alphaCGPhi2.boundaryFieldRef()[patchi];

            forAll(alphaCLPhi1p, facei)
            {
                const scalar phiPatch = Up[facei] & Sfp[facei];

                alphaCLPhi1p[facei] =
                    alpha1p[facei]*phiPatch*CLp[facei];

                alphaCGPhi2p[facei] =
                    alpha2p[facei]*phiPatch*CGp[facei];
            }
        }

        X1.correctBoundaryConditions();

    // In parallel, the original NCC patches (inletOutletBoiler/Condenser)
    // have size=0 so their updateCoeffs() is a no-op.  The actual NCC faces
    // live on nonConformalProcessorCyclic (NCPC) patches whose field values
    // are set by the MPI exchange (patchInternalField from the neighbour).
    // For incoming faces this gives the wrong composition; override them now
    // with the reservoir composition BEFORE fvMatrix assembly reads them.
    auto applyNCPCInletValues = [&]()
    {
        if (!Pstream::parRun()) return;
        forAll(X1.boundaryField(), patchi)
        {
            const fvPatchScalarField& pf = X1.boundaryField()[patchi];
            if (isA<inletOutletCondenserFvPatchScalarField>(pf))
            {
                refCast<const inletOutletCondenserFvPatchScalarField>(pf)
                    .applyInletValuesOnNCPCPatches(X1);
            }
            else if (isA<inletOutletBoilerFvPatchScalarField>(pf))
            {
                refCast<const inletOutletBoilerFvPatchScalarField>(pf)
                    .applyInletValuesOnNCPCPatches(X1);
            }
        }
    };

    applyNCPCInletValues();

        X2 = scalar(1) - X1;
        X2.correctBoundaryConditions();
    // ---------------------------------------------------------------------
    // new concentration equation
    // ---------------------------------------------------------------------

        const label nCompositionPicard
        (
            max(label(1),pimple.dict().lookupOrDefault<label>("compositionPicardIterations",2))
        );

        const bool boundComposition(pimple.dict().lookupOrDefault<bool>("compositionBoundX",false));

        // Molar-continuity errors: depend only on CL/CG and the fluxes,
        // which are fixed for the duration of the Picard loop.
        const volScalarField contErrCL1("contErrCL1", fvc::ddt(alpha1, CL) + fvc::div(alphaCLPhi1));
        const volScalarField contErrCG2("contErrCG2", fvc::ddt(alpha2, CG) + fvc::div(alphaCGPhi2));

        for (label picardIter = 0; picardIter < nCompositionPicard; ++picardIter)
        {
            Info<< "compositionPredictor Picard iteration "
                << picardIter + 1 << "/" << nCompositionPicard << endl;

            reconstructVLEPotential();

            // Diffusion coefficient in the common VLE potential xEq:
            //
            // Species flux:
            //     J = -aL*D_L*grad(x)
            //         -aG*D_G*grad(y)
            //
            // with:
            //     y = f(x)
            //     grad(y) = (dy/dx)*grad(x)
            //
            // Therefore:
            //     J = -[aL*D_L + aG*D_G*(dy/dx)] grad(x)
            volScalarField DeqBase
            (
                "DeqBase",
                ALPHA1*CL*thermophysicalTransport.D1Eff()
              + ALPHA2*CG*thermophysicalTransport.D2Eff()*dYdXeq
            );

            // Since equation variable is X1, but diffusion potential is xEq(X1):
            //
            //     grad(xEq) = (dxEq/dX1)*grad(X1)
            //
            // implicit slope:
            volScalarField DeqSlope("DeqSlope", DeqBase*dxEqdX);

            fvScalarMatrix X1Eqn
            (
                    correction
                    (
                        fvm::ddt(alpha1, CL, X1)
                        + fvm::div(alphaCLPhi1, X1)
                        - fvm::Sp(contErrCL1, X1)

                        + fvm::ddt(alpha2, CG, X1)
                        + fvm::div(alphaCGPhi2, X1)
                        - fvm::Sp(contErrCG2, X1)
                    )

                + fvc::ddt(alpha1, CL, X1)
                + fvc::div(alphaCLPhi1, X1)
                - contErrCL1*X1

                + fvc::ddt(alpha2, CG, X1)
                + fvc::div(alphaCGPhi2, X1)
                - contErrCG2*X1

                ==
                    // -------------------------------------------------------------
                    // RHS: new algebraic VLE diffusion
                    //
                    // target:
                    //     div(D1Base grad(X1L_alg))
                    //   + div(D2Base grad(X1G_alg))
                    //
                    // Picard form:
                    //     explicit fvc part with algebraic split
                    //     implicit correction with dX1L/dX1 and dX1G/dX1
                    // -------------------------------------------------------------
                    correction
                    (
                        fvm::laplacian(DeqSlope, X1)
                    )

                + fvc::laplacian(DeqBase, xEq_alg)
            );

            X1Eqn.relax();

            fvConstraints().constrain(X1Eqn);

            X1Eqn.solve();

            fvConstraints().constrain(X1);

            if (boundComposition)
            {
                Info<< "compositionBoundX before bounding X1 min/max = "
                    << min(X1).value() << " / " << max(X1).value()
                    << endl;

                X1 = min(max(X1, scalar(0)), scalar(1));
            }

            X1.correctBoundaryConditions();
            applyNCPCInletValues();
        }

    // ---------------------------------------------------------------------
    // Post-solve reservoir advance
    // Advance condenser/boiler reservoirs now that X1 is solved (X1_new).
    // alphaCLPhi1 and alphaCGPhi2 are still in scope and registered in the
    // mesh, so advanceReservoir() can look them up via the existing
    // compositionFluxFieldsAvailable() / lookupSurfaceFlux() mechanism.
    // ---------------------------------------------------------------------

        forAll(X1.boundaryField(), patchi)
        {
            fvPatchScalarField& pf =
                X1.boundaryFieldRef()[patchi];

            if (isA<inletOutletCondenserFvPatchScalarField>(pf))
            {
                refCast<inletOutletCondenserFvPatchScalarField>(pf)
                    .advanceReservoir();
            }
            else if (isA<inletOutletBoilerFvPatchScalarField>(pf))
            {
                refCast<inletOutletBoilerFvPatchScalarField>(pf)
                    .advanceReservoir();
            }
        }

    // ---------------------------------------------------------------------
    // NCPC molar conservation correction
    // ---------------------------------------------------------------------
    // The nonConformalProcessorCyclic (NCPC) patches at the column interior
    // cyclic interface are not perfectly conservative: NCC interpolation can
    // create or destroy moles there. Measure the molar flux imbalance at all
    // interior NCPC patches (i.e. those whose referPatch is NOT a boiler/
    // condenser patch) and redistribute the imbalance to the reservoirs so
    // that N_sys = N_col + N_D + N_B remains constant.
    //
    // Sign convention: alphaCLPhi1/alphaCGPhi2 > 0 = leaving the column.
    // Global sum over interior NCPC faces < 0 → NCPC creates moles in column.
    // Correction deltaN = flux_sum * dt: negative → remove moles from reservoir.
    // Liquid-phase imbalance is attributed to the boiler;
    // gas-phase imbalance to the condenser.
    if (Pstream::parRun())
    {
        // Identify reservoir patch names to exclude their NCPC counterparts
        wordHashSet reservoirPatchNames;
        forAll(X1.boundaryField(), patchi)
        {
            if
            (
                isA<inletOutletBoilerFvPatchScalarField>
                    (X1.boundaryField()[patchi])
             || isA<inletOutletCondenserFvPatchScalarField>
                    (X1.boundaryField()[patchi])
            )
            {
                reservoirPatchNames.insert(mesh.boundary()[patchi].name());
            }
        }

        scalar fluxLiqN  = 0;   // liquid-phase total molar flux at interior NCPC
        scalar fluxLiqN1 = 0;   // species-1 liquid-phase
        scalar fluxGasN  = 0;   // gas-phase total molar flux
        scalar fluxGasN1 = 0;   // species-1 gas-phase

        forAll(mesh.boundary(), patchi)
        {
            if (!isA<nonConformalProcessorCyclicFvPatch>(mesh.boundary()[patchi]))
                continue;

            const auto& ncpc =
                refCast<const nonConformalProcessorCyclicFvPatch>
                    (mesh.boundary()[patchi]);

            const word refName =
                ncpc.nonConformalProcessorCyclicPatch().referPatch().name();

            if (reservoirPatchNames.found(refName))
                continue;

            const scalarField& phiL  = alphaCLPhi1.boundaryField()[patchi];
            const scalarField& phiG  = alphaCGPhi2.boundaryField()[patchi];
            const scalarField  X1int =
                X1.boundaryField()[patchi].patchInternalField();

            forAll(phiL, facei)
            {
                fluxLiqN  += phiL[facei];
                fluxLiqN1 += phiL[facei] * X1int[facei];
                fluxGasN  += phiG[facei];
                fluxGasN1 += phiG[facei] * X1int[facei];
            }
        }

        reduce(fluxLiqN,  sumOp<scalar>());
        reduce(fluxLiqN1, sumOp<scalar>());
        reduce(fluxGasN,  sumOp<scalar>());
        reduce(fluxGasN1, sumOp<scalar>());

        const scalar dt = runTime.deltaT().value();

        const scalar deltaBoilerN  = fluxLiqN  * dt;
        const scalar deltaBoilerN1 = fluxLiqN1 * dt;
        const scalar deltaCondN    = fluxGasN  * dt;
        const scalar deltaCondN1   = fluxGasN1 * dt;

        if (Pstream::master())
        {
            Info<< "NCPC conservation: dN_boiler=" << deltaBoilerN
                << " dN_cond=" << deltaCondN
                << " dN_total=" << deltaBoilerN + deltaCondN
                << nl;
        }

        forAll(X1.boundaryFieldRef(), patchi)
        {
            fvPatchScalarField& pf = X1.boundaryFieldRef()[patchi];

            if (isA<inletOutletBoilerFvPatchScalarField>(pf))
            {
                refCast<inletOutletBoilerFvPatchScalarField>(pf)
                    .correctReservoir(deltaBoilerN, deltaBoilerN1);
            }
            else if (isA<inletOutletCondenserFvPatchScalarField>(pf))
            {
                refCast<inletOutletCondenserFvPatchScalarField>(pf)
                    .correctReservoir(deltaCondN, deltaCondN1);
            }
        }
    }

    // ---------------------------------------------------------------------
    // Final algebraic split for mixed-cell back-substitution
    // ---------------------------------------------------------------------

        reconstructVLEPotential();

        volScalarField X1L("X1L", (scalar(1) - Is_mix_Cell)*X1 + Is_mix_Cell*xEq_alg);
        volScalarField X1G("X1G", (scalar(1) - Is_mix_Cell)*X1 + Is_mix_Cell*yEq_alg);

        X1L.correctBoundaryConditions();
        X1G.correctBoundaryConditions();

        volScalarField X2L("X2L", scalar(1) - X1L);
        volScalarField X2G("X2G", scalar(1) - X1G);

        X2L.correctBoundaryConditions();
        X2G.correctBoundaryConditions();

        X2 = scalar(1) - X1;
        X2.correctBoundaryConditions();


        if (runTime.writeTime())
        {
            X1.write();
            X2.write();
            /*
            // Useful diagnostics.  Enable while debugging.
            A12.write();
            X1L.write();
            X1G.write(); 
            volScalarField VLEResidual ("VLEResidual",Is_mix_Cell*(X1G- A12*X1L/(scalar(1) + (A12 - scalar(1))*X1L)));
            VLEResidual.write();
            */
        }

    // ---------------------------------------------------------------------
    // Back-substitution to phase mass fractions
    // ---------------------------------------------------------------------

        volScalarField ALPHA1_pure("ALPHA1_pure", pos(ALPHA1 - (scalar(1) - aTol)));    // only cells with pure phases, mixed cells 0
        volScalarField ALPHA2_pure("ALPHA2_pure", pos(ALPHA2 - (scalar(1) - aTol)));    // only cells with pure phases, mixed cells 0
                
        const dimensionedScalar denomMinW("denomMinW", W1.dimensions(),SMALL);
        // Mass fraction of the pure phases

            volScalarField Y1L_pure("Y1L_pure", ALPHA1_pure*X1*W1/max(X1*W1 + X2*W2, denomMinW));       // gives Y1L where alpha1 =1, everywhere else (0 as Dummy)
            volScalarField Y2L_pure("Y2L_pure", scalar(1) - Y1L_pure        );                          // gives Y2L where alpha1 =1, everywhere else (1 as Dummy)
            volScalarField Y2G_pure("Y2G_pure", ALPHA2_pure*X2*W2/max(X1*W1 + X2*W2, denomMinW));       // gives Y2G where alpha2 =1, everywhere else (0 as Dummy)
            volScalarField Y1G_pure("Y1G_pure", scalar(1) - Y2G_pure);                                  // gives Y1G where alpha2 =1, everywhere else (1 as Dummy)
        // Mass fraction of the mixed phases

            volScalarField Y1L_mix("Y1L_mix", Is_mix_Cell*X1L*W1/max(X1L*W1 + X2L*W2, denomMinW));
            volScalarField Y1G_mix("Y1G_mix", Is_mix_Cell*X1G*W1/max(X1G*W1 + X2G*W2, denomMinW));
            volScalarField Y2L_mix("Y2L_mix", Is_mix_Cell*(scalar(1) - Y1L_mix));
            volScalarField Y2G_mix("Y2G_mix", Is_mix_Cell*(scalar(1) - Y1G_mix));

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

        nMoles1 = alpha1*mixture.thermo1().rho()*Y1L/W1
                + alpha2*mixture.thermo2().rho()*Y1G/W1;

        nMoles2 = alpha1*mixture.thermo1().rho()*Y2L/W2
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