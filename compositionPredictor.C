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

        D does not appear in the VLE jump.  D appears only in D1Base/D2Base.

    Important:
        If D1 = D2 = 0, the RHS becomes zero and the equation is exactly the
        old stable phase-separated advection equation.  No artificial VLE
        relaxation is introduced in pure advection.

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

    const scalar aTol = 1e-5;

    volScalarField a
    (
        "a",
        min(max(alpha1, scalar(0)), scalar(1))
    );

    volScalarField ALPHA1
    (
        "ALPHA1",
        (scalar(1) - pos(aTol - a) - pos(a - (scalar(1) - aTol)))*a
      + pos(a - (scalar(1) - aTol))
    );

    volScalarField ALPHA2
    (
        "ALPHA2",
        scalar(1) - ALPHA1
    );

    volScalarField Is_mix_Cell
    (
        "Is_mix_Cell",
        pos(ALPHA1 - aTol)*pos(ALPHA2 - aTol)
    );

    // ---------------------------------------------------------------------
    // Phase molar concentrations
    // ---------------------------------------------------------------------

    volScalarField C1L
    (
        "C1L",
        pos(alpha1 - aTol)*(rho1/W1)*Y1L
    );

    volScalarField C1G
    (
        "C1G",
        pos(alpha2 - aTol)*(rho2/W1)*Y1G
    );

    volScalarField C2L
    (
        "C2L",
        pos(alpha1 - aTol)*(rho1/W2)*Y2L
    );

    volScalarField C2G
    (
        "C2G",
        pos(alpha2 - aTol)*(rho2/W2)*Y2G
    );

    volScalarField CL
    (
        "CL",
        C1L + C2L
    );

    volScalarField CG
    (
        "CG",
        C1G + C2G
    );

    // Molar inventory coefficients
    volScalarField aL
    (
        "aL",
        ALPHA1*CL
    );

    volScalarField aG
    (
        "aG",
        ALPHA2*CG
    );

    volScalarField M
    (
        "M",
        aL + aG
    );

    const dimensionedScalar Mmin
    (
        "Mmin",
        M.dimensions(),
        SMALL
    );

    // Bulk species inventory for initial X calculation
    volScalarField C1
    (
        "C1",
        ALPHA1*C1L + ALPHA2*C1G
    );

    volScalarField X1Calc
    (
        "X1Calc",
        C1/max(M, Mmin)
    );

    volScalarField X2Calc
    (
        "X2Calc",
        scalar(1) - X1Calc
    );

    // ---------------------------------------------------------------------
    // Read/create mole-fraction fields
    // ---------------------------------------------------------------------

    const word X1name = mixture.species1Name() + ".X";
    const word X2name = mixture.species2Name() + ".X";

    const word startTimeName =
        Time::timeName(runTime.startTime().value(), 6);

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

    volScalarField K1eq
    (
        "K1eq",
        vle1Ptr->K(pVLE, T)
    );

    volScalarField K2eq
    (
        "K2eq",
        vle2Ptr->K(pVLE, T)
    );

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

volScalarField xEq_alg
(
    "xEq_alg",
    X1
);

volScalarField yEq_alg
(
    "yEq_alg",
    X1
);

volScalarField dxEqdX
(
    "dxEqdX",
    X1*scalar(0) + scalar(1)
);

volScalarField dYdXeq
(
    "dYdXeq",
    X1*scalar(0) + scalar(1)
);

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
    &A12
]()
{
    scalarField& xEq = xEq_alg.primitiveFieldRef();
    scalarField& yEq = yEq_alg.primitiveFieldRef();
    scalarField& dxDX = dxEqdX.primitiveFieldRef();
    scalarField& dydx = dYdXeq.primitiveFieldRef();

    const scalarField& X = X1.primitiveField();
    const scalarField& aLI = aL.primitiveField();
    const scalarField& aGI = aG.primitiveField();
    const scalarField& MI = M.primitiveField();
    const scalarField& AI = A12.primitiveField();

    forAll(X, celli)
    {
        const scalar aa = aLI[celli];
        const scalar bb = aGI[celli];
        const scalar mm = MI[celli];
        const scalar A = AI[celli];

        if (aa > SMALL && bb > SMALL && mm > SMALL)
        {
            // Mixed/interface cell:
            // q = aa*x + bb*y
            // y = A*x/(1+(A-1)*x)

            const scalar q = mm*X[celli];

            const scalar c2 = aa*(A - scalar(1));
            const scalar c1 = aa + bb*A - q*(A - scalar(1));

            scalar x = X[celli];

            if (mag(c2) > VSMALL)
            {
                const scalar disc =
                    c1*c1 + scalar(4)*c2*q;

                x =
                    (-c1 + sqrt(max(disc, scalar(0))))
                   /(scalar(2)*c2);
            }
            else
            {
                x = q/(aa + bb*A);
            }

            const scalar den = scalar(1) + (A - scalar(1))*x;
            const scalar y = A*x/den;
            const scalar s = A/sqr(den);   // dy/dx

            const scalar dqdx = aa + bb*s;

            xEq[celli] = x;
            yEq[celli] = y;
            dydx[celli] = s;

            if (mag(dqdx) > VSMALL)
            {
                dxDX[celli] = mm/dqdx;
            }
            else
            {
                dxDX[celli] = scalar(1);
            }
        }
        else if (aa > SMALL && mm > SMALL)
        {
            // Pure liquid:
            // X1 is liquid composition x.
            const scalar x = X[celli];

            const scalar den = scalar(1) + (A - scalar(1))*x;
            const scalar s = A/sqr(den);

            xEq[celli] = x;
            yEq[celli] = A*x/den;
            dydx[celli] = s;
            dxDX[celli] = scalar(1);
        }
        else if (bb > SMALL && mm > SMALL)
        {
            // Pure gas:
            // X1 is gas composition y.
            // For diffusion potential, convert y -> x via inverse VLE.
            const scalar y = X[celli];

            const scalar denomInv =
                A - (A - scalar(1))*y;

            scalar x = y/max(denomInv, VSMALL);

            const scalar den = scalar(1) + (A - scalar(1))*x;
            const scalar s = A/sqr(den);   // dy/dx at x

            xEq[celli] = x;
            yEq[celli] = y;
            dydx[celli] = s;

            // x = inverseVLE(y), so dx/dX = dx/dy = 1/(dy/dx)
            dxDX[celli] = scalar(1)/max(s, VSMALL);
        }
        else
        {
            xEq[celli] = X[celli];
            yEq[celli] = X[celli];
            dydx[celli] = scalar(1);
            dxDX[celli] = scalar(1);
        }
    }

    forAll(X1.boundaryField(), patchi)
    {
        scalarField& xEqp = xEq_alg.boundaryFieldRef()[patchi];
        scalarField& yEqp = yEq_alg.boundaryFieldRef()[patchi];
        scalarField& dxDXp = dxEqdX.boundaryFieldRef()[patchi];
        scalarField& dydxp = dYdXeq.boundaryFieldRef()[patchi];

        const fvPatchScalarField& Xp = X1.boundaryField()[patchi];
        const fvPatchScalarField& aLp = aL.boundaryField()[patchi];
        const fvPatchScalarField& aGp = aG.boundaryField()[patchi];
        const fvPatchScalarField& Mp = M.boundaryField()[patchi];
        const fvPatchScalarField& Ap = A12.boundaryField()[patchi];

        forAll(Xp, facei)
        {
            const scalar aa = aLp[facei];
            const scalar bb = aGp[facei];
            const scalar mm = Mp[facei];
            const scalar A = Ap[facei];

            if (aa > SMALL && bb > SMALL && mm > SMALL)
            {
                const scalar q = mm*Xp[facei];

                const scalar c2 = aa*(A - scalar(1));
                const scalar c1 = aa + bb*A - q*(A - scalar(1));

                scalar x = Xp[facei];

                if (mag(c2) > VSMALL)
                {
                    const scalar disc =
                        c1*c1 + scalar(4)*c2*q;

                    x =
                        (-c1 + sqrt(max(disc, scalar(0))))
                       /(scalar(2)*c2);
                }
                else
                {
                    x = q/(aa + bb*A);
                }

                const scalar den = scalar(1) + (A - scalar(1))*x;
                const scalar y = A*x/den;
                const scalar s = A/sqr(den);

                const scalar dqdx = aa + bb*s;

                xEqp[facei] = x;
                yEqp[facei] = y;
                dydxp[facei] = s;

                if (mag(dqdx) > VSMALL)
                {
                    dxDXp[facei] = mm/dqdx;
                }
                else
                {
                    dxDXp[facei] = scalar(1);
                }
            }
            else if (aa > SMALL && mm > SMALL)
            {
                const scalar x = Xp[facei];

                const scalar den = scalar(1) + (A - scalar(1))*x;
                const scalar s = A/sqr(den);

                xEqp[facei] = x;
                yEqp[facei] = A*x/den;
                dydxp[facei] = s;
                dxDXp[facei] = scalar(1);
            }
            else if (bb > SMALL && mm > SMALL)
            {
                const scalar y = Xp[facei];

                const scalar denomInv =
                    A - (A - scalar(1))*y;

                const scalar x = y/max(denomInv, VSMALL);

                const scalar den = scalar(1) + (A - scalar(1))*x;
                const scalar s = A/sqr(den);

                xEqp[facei] = x;
                yEqp[facei] = y;
                dydxp[facei] = s;
                dxDXp[facei] = scalar(1)/max(s, VSMALL);
            }
            else
            {
                xEqp[facei] = Xp[facei];
                yEqp[facei] = Xp[facei];
                dydxp[facei] = scalar(1);
                dxDXp[facei] = scalar(1);
            }
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

    surfaceScalarField alphaCLPhi1
    (
        "alphaCLPhi1",
        alphaPhi1*fvc::interpolate(CL)
    );

    surfaceScalarField alphaCGPhi2
    (
        "alphaCGPhi2",
        alphaPhi2*fvc::interpolate(CG)
    );

    // ---------------------------------------------------------------------
    // GCST equation
    // ---------------------------------------------------------------------

    const label nCompositionPicard
    (
        max
        (
            label(1),
            pimple.dict().lookupOrDefault<label>
            (
                "compositionPicardIterations",
                2
            )
        )
    );

    const bool boundComposition
    (
        pimple.dict().lookupOrDefault<bool>
        (
            "compositionBoundX",
            false
        )
    );

    for (label picardIter = 0; picardIter < nCompositionPicard; ++picardIter)
    {
        Info<< "compositionPredictor Picard iteration "
            << picardIter + 1 << "/" << nCompositionPicard << endl;

        reconstructVLEPotential();

        // Old stable phase-wise continuity errors
        volScalarField contErrCL1
        (
            "contErrCL1",
            fvc::ddt(alpha1, CL) + fvc::div(alphaCLPhi1)
        );

        volScalarField contErrCG2
        (
            "contErrCG2",
            fvc::ddt(alpha2, CG) + fvc::div(alphaCGPhi2)
        );

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
        volScalarField DeqSlope
        (
            "DeqSlope",
            DeqBase*dxEqdX
        );

        fvScalarMatrix X1Eqn
        (
            // -------------------------------------------------------------
            // LHS: old stable GCST advection, unchanged
            // -------------------------------------------------------------
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
            // RHS: algebraic VLE diffusion
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

        X2 = scalar(1) - X1;
        X2.correctBoundaryConditions();

        Info<< "X1 min/max = "
            << min(X1).value()
            << " / "
            << max(X1).value()
            << endl;
    }

    // ---------------------------------------------------------------------
    // Final algebraic split for mixed-cell back-substitution
    // ---------------------------------------------------------------------

    reconstructVLEPotential();

    volScalarField X1L
    (
        "X1L",
        (scalar(1) - Is_mix_Cell)*X1 + Is_mix_Cell*xEq_alg
    );

    volScalarField X1G
    (
        "X1G",
        (scalar(1) - Is_mix_Cell)*X1 + Is_mix_Cell*yEq_alg
    );

    X1L.correctBoundaryConditions();
    X1G.correctBoundaryConditions();

    volScalarField X2L
    (
        "X2L",
        scalar(1) - X1L
    );

    volScalarField X2G
    (
        "X2G",
        scalar(1) - X1G
    );

    X2L.correctBoundaryConditions();
    X2G.correctBoundaryConditions();

    X2 = scalar(1) - X1;
    X2.correctBoundaryConditions();


    if (runTime.writeTime())
    {
        X1.write();
        X2.write();

        // Useful diagnostics.  Enable while debugging.
        A12.write();
        X1L.write();
        X1G.write(); 
        
        volScalarField VLEResidual
        (
            "VLEResidual",
            Is_mix_Cell*
            (
                X1G
              - A12*X1L/(scalar(1) + (A12 - scalar(1))*X1L)
            )
        );

        VLEResidual.write();
    }

    // ---------------------------------------------------------------------
    // Back-substitution to phase mass fractions
    // ---------------------------------------------------------------------

    volScalarField ALPHA1_pure
    (
        "ALPHA1_pure",
        pos(ALPHA1 - (scalar(1) - aTol))
    );

    volScalarField ALPHA2_pure
    (
        "ALPHA2_pure",
        pos(ALPHA2 - (scalar(1) - aTol))
    );

    volScalarField ALPHA1_mix
    (
        "ALPHA1_mix",
        pos(ALPHA1 - ALPHA1_pure)*(ALPHA1 - ALPHA1_pure)
    );

    volScalarField ALPHA2_mix
    (
        "ALPHA2_mix",
        pos(ALPHA2 - ALPHA2_pure)*(ALPHA2 - ALPHA2_pure)
    );

    const dimensionedScalar denomMinW
    (
        "denomMinW",
        W1.dimensions(),
        SMALL
    );

    volScalarField Y1L_pure
    (
        "Y1L_pure",
        ALPHA1_pure*X1*W1/max(X1*W1 + X2*W2, denomMinW)
    );

    volScalarField Y2L_pure
    (
        "Y2L_pure",
        scalar(1) - Y1L_pure
    );

    volScalarField Y2G_pure
    (
        "Y2G_pure",
        ALPHA2_pure*X2*W2/max(X1*W1 + X2*W2, denomMinW)
    );

    volScalarField Y1G_pure
    (
        "Y1G_pure",
        scalar(1) - Y2G_pure
    );

    volScalarField Y1L_mix
    (
        "Y1L_mix",
        Is_mix_Cell*X1L*W1/max(X1L*W1 + X2L*W2, denomMinW)
    );

    volScalarField Y1G_mix
    (
        "Y1G_mix",
        Is_mix_Cell*X1G*W1/max(X1G*W1 + X2G*W2, denomMinW)
    );

    volScalarField Y2L_mix
    (
        "Y2L_mix",
        Is_mix_Cell*(scalar(1) - Y1L_mix)
    );

    volScalarField Y2G_mix
    (
        "Y2G_mix",
        Is_mix_Cell*(scalar(1) - Y1G_mix)
    );

    // Dummy-value convention retained
    Y1L = Y1L_pure + Y1L_mix;
    Y2L = Y2L_pure + Y2L_mix - Is_mix_Cell;
    Y1G = Y1G_pure + Y1G_mix - Is_mix_Cell;
    Y2G = Y2G_pure + Y2G_mix;

    Y1L.correctBoundaryConditions();
    Y2L.correctBoundaryConditions();
    Y1G.correctBoundaryConditions();
    Y2G.correctBoundaryConditions();

    // ---------------------------------------------------------------------
    // Mole conservation diagnostics
    // ---------------------------------------------------------------------

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