/*---------------------------------------------------------------------------*\
  Dynamic reboiler finite-holdup inletOutlet composition boundary condition
\*---------------------------------------------------------------------------*/

#include "inletOutletBoilerFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "typeInfo.H"
#include <fstream>
#include "OSspecific.H"
#include "nonConformalProcessorCyclicFvPatch.H"
#include "HashSet.H"

namespace Foam
{

scalar inletOutletBoilerFvPatchScalarField::clamp
(
    const scalar x
) const
{
    return min(max(x, minValue_), maxValue_);
}


scalar inletOutletBoilerFvPatchScalarField::yVLE
(
    const scalar x,
    const scalar A
) const
{
    const scalar xb = clamp(x);
    const scalar Ab = max(A, SMALL);

    return clamp(Ab*xb/max(scalar(1) + (Ab - scalar(1))*xb, SMALL));
}


wordList inletOutletBoilerFvPatchScalarField::selectedPatchNames() const
{
    if (usePatchList_)
    {
        return patchNames_;
    }

    if (groupName_.empty())
    {
        return wordList(1, patch().name());
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (!mesh.foundObject<volScalarField>(this->internalField().name()))
    {
        return wordList(1, patch().name());
    }

    const volScalarField& x =
        mesh.lookupObject<volScalarField>(this->internalField().name());

    wordList names(x.boundaryField().size());
    label n = 0;

    forAll(x.boundaryField(), patchi)
    {
        const fvPatchScalarField& pf = x.boundaryField()[patchi];

        if (isA<inletOutletBoilerFvPatchScalarField>(pf))
        {
            const inletOutletBoilerFvPatchScalarField& bpf =
                refCast<const inletOutletBoilerFvPatchScalarField>(pf);

            if (bpf.groupName_ == groupName_)
            {
                names[n++] = mesh.boundary()[patchi].name();
            }
        }
    }

    if (n == 0)
    {
        return wordList(1, patch().name());
    }

    names.setSize(n);
    return names;
}


label inletOutletBoilerFvPatchScalarField::patchID
(
    const word& patchName
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const label patchi = mesh.boundaryMesh().findIndex(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Cannot find patch " << patchName << nl
            << "Available patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return patchi;
}


const surfaceScalarField&
inletOutletBoilerFvPatchScalarField::lookupSurfaceFlux
(
    const word& name
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (!mesh.foundObject<surfaceScalarField>(name))
    {
        FatalErrorInFunction
            << "Cannot find surfaceScalarField " << name << nl
            << "Available surfaceScalarFields are "
            << mesh.toc<surfaceScalarField>() << nl
            << "For incompressibleVoFTC the expected defaults are:" << nl
            << "    liquidFlux alphaCLPhi1;" << nl
            << "    gasFlux    alphaCGPhi2;"
            << exit(FatalError);
    }

    return mesh.lookupObject<surfaceScalarField>(name);
}


bool inletOutletBoilerFvPatchScalarField::compositionFluxFieldsAvailable() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    return
        mesh.foundObject<surfaceScalarField>(liquidFluxName_)
     && mesh.foundObject<surfaceScalarField>(gasFluxName_)
     && mesh.foundObject<volScalarField>(this->internalField().name())
     && mesh.foundObject<volScalarField>(A12Name_);
}

void inletOutletBoilerFvPatchScalarField::measureCompositionPredictorFluxes
(
    const wordList& patchNames,
    scalar& nLtoReservoir,
    scalar& nGtoReservoir,
    scalar& nLfromReservoir,
    scalar& nGfromReservoir,
    scalar& n1LtoReservoir,
    scalar& n1GtoReservoir,
    scalar& n1LfromReservoir,
    scalar& n1GfromReservoir
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (!compositionFluxFieldsAvailable())
    {
        FatalErrorInFunction
            << "Cannot measure compositionPredictor fluxes. Missing fields:" << nl
            << "  liquidFlux=" << liquidFluxName_
            << " found=" << mesh.foundObject<surfaceScalarField>(liquidFluxName_) << nl
            << "  gasFlux=" << gasFluxName_
            << " found=" << mesh.foundObject<surfaceScalarField>(gasFluxName_) << nl
            << "  X=" << this->internalField().name()
            << " found=" << mesh.foundObject<volScalarField>(this->internalField().name()) << nl
            << "  A12=" << A12Name_
            << " found=" << mesh.foundObject<volScalarField>(A12Name_) << nl
            << exit(FatalError);
    }

    const surfaceScalarField& liquidFlux =
        lookupSurfaceFlux(liquidFluxName_);

    const surfaceScalarField& gasFlux =
        lookupSurfaceFlux(gasFluxName_);

    const volScalarField& X =
        mesh.lookupObject<volScalarField>(this->internalField().name());

    const scalar Amean = meanA12(patchNames);
    const scalar xB = reservoirX_;
    const scalar yB = yVLE(xB, Amean);

    nLtoReservoir = 0;
    nGtoReservoir = 0;
    nLfromReservoir = 0;
    nGfromReservoir = 0;

    n1LtoReservoir = 0;
    n1GtoReservoir = 0;
    n1LfromReservoir = 0;
    n1GfromReservoir = 0;

    // In parallel, original NCC patches have size=0; include procBoundary patches.
    DynamicList<label> measureIDs(patchNames.size());
    {
        forAll(patchNames, i) measureIDs.append(patchID(patchNames[i]));
        if (Pstream::parRun())
        {
            const wordHashSet nameSet(patchNames);
            const fvBoundaryMesh& bm = patch().boundaryMesh();
            forAll(bm, pi)
            {
                if (!isA<nonConformalProcessorCyclicFvPatch>(bm[pi])) continue;
                const auto& ncpc =
                    refCast<const nonConformalProcessorCyclicFvPatch>(bm[pi]);
                if (nameSet.found(
                        ncpc.nonConformalProcessorCyclicPatch()
                            .referPatch().name()))
                    measureIDs.append(pi);
            }
        }
    }

    forAll(measureIDs, _i)
    {
        const label patchi = measureIDs[_i];

        const fvsPatchScalarField& Lp =
            liquidFlux.boundaryField()[patchi];

        const fvsPatchScalarField& Gp =
            gasFlux.boundaryField()[patchi];

        const scalarField xPatchInternal
        (
            X.boundaryField()[patchi].patchInternalField()
        );

        forAll(Lp, facei)
        {
            // compositionPredictor convention for alphaCLPhi1/alphaCGPhi2:
            // positive boundary flux leaves the CFD domain into the apparatus.
            const scalar Lto = max( Lp[facei], scalar(0));
            const scalar Gto = max( Gp[facei], scalar(0));

            const scalar Lfrom = max(-Lp[facei], scalar(0));
            const scalar Gfrom = max(-Gp[facei], scalar(0));

            // Outflow from column into apparatus: VLE split for mixed cells.
            // Solve X_overall*(Lto+Gto) = Lto*xL + Gto*yVLE(xL,A12) for xL.
            const scalar xOut = clamp(xPatchInternal[facei]);
            const scalar S = Lto + Gto;

            scalar xL, yG;
            if (S < SMALL || mag(Amean - 1) < 1e-10)
            {
                xL = xOut;
                yG = xOut;
            }
            else if (Lto < SMALL)
            {
                // Pure gas: yG = xOut, invert VLE for xL
                xL = clamp(xOut / max(Amean - xOut*(Amean - 1), SMALL));
                yG = clamp(xOut);
            }
            else if (Gto < SMALL)
            {
                xL = xOut;
                yG = yVLE(xOut, Amean);
            }
            else
            {
                const scalar a = Lto*(Amean - 1);
                const scalar b = Lto + Gto*Amean - xOut*S*(Amean - 1);
                const scalar c = -xOut*S;
                xL = clamp((-b + sqrt(max(b*b - 4*a*c, scalar(0))))/(2*a));
                yG = yVLE(xL, Amean);
            }

            nLtoReservoir += Lto;
            nGtoReservoir += Gto;

            nLfromReservoir += Lfrom;
            nGfromReservoir += Gfrom;

            n1LtoReservoir += Lto*xL;
            n1GtoReservoir += Gto*yG;

            // Boiler return into column: liquid gets xB, gas gets yB.
            n1LfromReservoir += Lfrom*xB;
            n1GfromReservoir += Gfrom*yB;
        }
    }

    reduce(nLtoReservoir, sumOp<scalar>());
    reduce(nGtoReservoir, sumOp<scalar>());
    reduce(nLfromReservoir, sumOp<scalar>());
    reduce(nGfromReservoir, sumOp<scalar>());

    reduce(n1LtoReservoir, sumOp<scalar>());
    reduce(n1GtoReservoir, sumOp<scalar>());
    reduce(n1LfromReservoir, sumOp<scalar>());
    reduce(n1GfromReservoir, sumOp<scalar>());
}

void inletOutletBoilerFvPatchScalarField::readPatchSelection
(
    const dictionary& dict
)
{
    if (dict.found("patches"))
    {
        dict.lookup("patches") >> patchNames_;
        usePatchList_ = true;
    }
    else
    {
        patchNames_.clear();
        usePatchList_ = false;
    }

    if (usePatchList_)
    {
        bool hasThisPatch = false;

        forAll(patchNames_, i)
        {
            if (patchNames_[i] == patch().name())
            {
                hasThisPatch = true;
                break;
            }
        }

        if (!hasThisPatch)
        {
            FatalErrorInFunction
                << "Patch list for " << typeName
                << " must include the current patch " << patch().name()
                << ". Given patches = " << patchNames_
                << exit(FatalError);
        }

        forAll(patchNames_, i)
        {
            patchID(patchNames_[i]);
        }
    }
}


scalar inletOutletBoilerFvPatchScalarField::meanA12
(
    const wordList& patchNames
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (!mesh.foundObject<volScalarField>(A12Name_))
    {
        FatalErrorInFunction
            << "Cannot find volScalarField " << A12Name_ << nl
            << "compositionPredictor should create and update A12 before "
            << typeName << " is evaluated." << nl
            << "Available volScalarFields are "
            << mesh.toc<volScalarField>()
            << exit(FatalError);
    }

    const volScalarField& A12 =
        mesh.lookupObject<volScalarField>(A12Name_);

    scalar sumA = 0;
    scalar sumA12A = 0;

    // In parallel, original NCC patches have size=0; actual faces are on
    // nonConformalProcessorCyclic patches. Include both in the measurement.
    DynamicList<label> measureIDs(patchNames.size());
    {
        forAll(patchNames, i) measureIDs.append(patchID(patchNames[i]));
        if (Pstream::parRun())
        {
            const wordHashSet nameSet(patchNames);
            const fvBoundaryMesh& bm = patch().boundaryMesh();
            forAll(bm, pi)
            {
                if (!isA<nonConformalProcessorCyclicFvPatch>(bm[pi])) continue;
                const auto& ncpc =
                    refCast<const nonConformalProcessorCyclicFvPatch>(bm[pi]);
                if (nameSet.found(
                        ncpc.nonConformalProcessorCyclicPatch()
                            .referPatch().name()))
                    measureIDs.append(pi);
            }
        }
    }

    forAll(measureIDs, _i)
    {
        const label patchi = measureIDs[_i];

        const scalarField& magSf = mesh.boundary()[patchi].magSf();
        const fvPatchScalarField& A12p = A12.boundaryField()[patchi];

        forAll(magSf, facei)
        {
            const scalar A = magSf[facei];

            sumA += A;
            sumA12A += A12p[facei]*A;
        }
    }

    reduce(sumA, sumOp<scalar>());
    reduce(sumA12A, sumOp<scalar>());

    if (sumA <= SMALL)
    {
        // NCC patch area is 0 until stitcher_->connect() runs in the first
        // PIMPLE iteration. Return neutral A12 = 1 (no VLE selectivity) as
        // startup default; correct value is used from the next iteration on.
        return scalar(1);
    }

    return max(sumA12A/sumA, SMALL);
}


void inletOutletBoilerFvPatchScalarField::advanceReservoir()
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const label timeIndex = mesh.time().timeIndex();

    if (timeIndex == lastReservoirUpdateIndex_)
    {
        return;
    }

    if (reservoirMoles_ <= minReservoirMoles_)
    {
        reservoirMoles_ = minReservoirMoles_;
    }

    // During field construction/read, OpenFOAM may call updateCoeffs()
    // before compositionPredictor has created alphaCLPhi1/alphaCGPhi2/A12.
    // Do not abort and do not mark this time step as updated. Later in the
    // same time step, compositionPredictor creates the fields and calls
    // correctBoundaryConditions(); then the reservoir is advanced.
    if (!compositionFluxFieldsAvailable())
    {
        if (log_)
        {
            Info<< typeName << " " << patch().name()
                << ": waiting for compositionPredictor flux fields before reservoir update:"
                << " liquidFlux=" << liquidFluxName_
                << " found=" << mesh.foundObject<surfaceScalarField>(liquidFluxName_)
                << ", gasFlux=" << gasFluxName_
                << " found=" << mesh.foundObject<surfaceScalarField>(gasFluxName_)
                << ", A12=" << A12Name_
                << " found=" << mesh.foundObject<volScalarField>(A12Name_)
                << ", X=" << this->internalField().name()
                << " found=" << mesh.foundObject<volScalarField>(this->internalField().name())
                << endl;
        }

        return;
    }

    const scalar dt = mesh.time().deltaTValue();

    if (dt <= SMALL)
    {
        lastReservoirUpdateIndex_ = timeIndex;
        return;
    }

    const wordList patchNames = selectedPatchNames();

    const scalar Amean = meanA12(patchNames);

    scalar nLtoReservoir = 0;
    scalar nGtoReservoir = 0;
    scalar nLfromReservoir = 0;
    scalar nGfromReservoir = 0;

    scalar n1LtoReservoir = 0;
    scalar n1GtoReservoir = 0;
    scalar n1LfromReservoir = 0;
    scalar n1GfromReservoir = 0;

    measureCompositionPredictorFluxes
    (
        patchNames,
        nLtoReservoir,
        nGtoReservoir,
        nLfromReservoir,
        nGfromReservoir,
        n1LtoReservoir,
        n1GtoReservoir,
        n1LfromReservoir,
        n1GfromReservoir
    );

    const scalar n1ToReservoir =
        n1LtoReservoir + n1GtoReservoir;

    const scalar n1FromReservoir =
        n1LfromReservoir + n1GfromReservoir;


    const scalar NOld = reservoirMoles_;
    const scalar xOld = reservoirX_;
    const scalar yOld = yVLE(xOld, Amean);

    const scalar nToReservoir =
        nLtoReservoir + nGtoReservoir;

    const scalar nFromReservoir =
        nLfromReservoir + nGfromReservoir;


    // Complete finite-holdup reboiler balances:
    //   dN/dt       = L_toB + G_toB - L_fromB - G_fromB
    //   d(N*x)/dt   = L_toB*x_L + G_toB*y_G
    //               - L_fromB*x_B - G_fromB*y_VLE(x_B)
    const scalar dNdt =
        nToReservoir - nFromReservoir;

    const scalar dN1dt =
        n1ToReservoir - n1FromReservoir;

    scalar NEuler =
        NOld + dt*dNdt;

    if (NEuler < minReservoirMoles_)
    {
        if (log_)
        {
            WarningInFunction
                << typeName << " " << patch().name()
                << ": boiler reservoir would be depleted. "
                << "Clamping N from " << NEuler
                << " to minReservoirMoles=" << minReservoirMoles_
                << nl;
        }

        NEuler = minReservoirMoles_;
    }

    scalar N1Euler =
        NOld*xOld + dt*dN1dt;

    N1Euler =
        min(max(N1Euler, minValue_*NEuler), maxValue_*NEuler);

    // Relax both conserved variables consistently.
    // relax = 1 gives the physical explicit Euler update.
    const scalar N1Old = NOld*xOld;

    const scalar NRelaxed =
        max(NOld + relax_*(NEuler - NOld), minReservoirMoles_);

    scalar N1Relaxed =
        N1Old + relax_*(N1Euler - N1Old);

    N1Relaxed =
        min(max(N1Relaxed, minValue_*NRelaxed), maxValue_*NRelaxed);

    reservoirMoles_ = NRelaxed;
    reservoirX_ = clamp(N1Relaxed/reservoirMoles_);
    reservoirY_ = yVLE(reservoirX_, Amean);
    writeApparatusBalance
    (
        dt,
        NOld,
        xOld,
        reservoirMoles_,
        reservoirX_,
        nLtoReservoir,
        nGtoReservoir,
        nLfromReservoir,
        nGfromReservoir,
        n1ToReservoir,
        n1FromReservoir
    );
    lastReservoirUpdateIndex_ = timeIndex;

    lastA12_ = Amean;
    lastLiquidToReservoir_ = nLtoReservoir;
    lastGasToReservoir_ = nGtoReservoir;
    lastLiquidFromReservoir_ = nLfromReservoir;
    lastGasFromReservoir_ = nGfromReservoir;
    lastSpeciesToReservoir_ = n1ToReservoir;
    lastSpeciesFromReservoir_ = n1FromReservoir;
    lastMolesResidual_ = dNdt;
    lastSpeciesResidual_ = dN1dt;

    lastLiquidXToReservoir_ =
        nLtoReservoir > minFlux_
      ? n1LtoReservoir/nLtoReservoir
      : xOld;

    lastGasXToReservoir_ =
        nGtoReservoir > minFlux_
      ? n1GtoReservoir/nGtoReservoir
      : yOld;

    if (log_)
    {
        Info<< typeName << " " << patch().name()
            << ": NB=" << reservoirMoles_
            << ", xB=" << reservoirX_
            << ", yB=" << reservoirY_
            << ", A12=" << Amean
            << ", LtoB=" << nLtoReservoir
            << ", GtoB=" << nGtoReservoir
            << ", LfromB=" << nLfromReservoir
            << ", GfromB=" << nGfromReservoir
            << ", xLtoB=" << lastLiquidXToReservoir_
            << ", xGtoB=" << lastGasXToReservoir_
            << ", speciesToB=" << n1ToReservoir
            << ", speciesFromB=" << n1FromReservoir
            << ", dNdt=" << dNdt
            << ", dN1dt=" << dN1dt
            << endl;
    }
}

void Foam::inletOutletBoilerFvPatchScalarField::writeApparatusBalance
(
    const scalar dt,
    const scalar NOld,
    const scalar xOld,
    const scalar NNew,
    const scalar xNew,
    const scalar nLtoReservoir,
    const scalar nGtoReservoir,
    const scalar nLfromReservoir,
    const scalar nGfromReservoir,
    const scalar n1ToReservoir,
    const scalar n1FromReservoir
) const
{
    if (!writeBalance_ || !Pstream::master())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const scalar N1Old = NOld*xOld;
    const scalar N1New = NNew*xNew;

    const scalar Nto = nLtoReservoir + nGtoReservoir;
    const scalar Nfrom = nLfromReservoir + nGfromReservoir;

    const scalar dNdt = Nto - Nfrom;
    const scalar dN1dt = n1ToReservoir - n1FromReservoir;

    const scalar dNSoll = dt*dNdt;
    const scalar dNIst = NNew - NOld;

    const scalar dN1Soll = dt*dN1dt;
    const scalar dN1Ist = N1New - N1Old;

    const scalar NSoll = NOld + dNSoll;
    const scalar N1Soll = N1Old + dN1Soll;

    const scalar N2Old = NOld - N1Old;
    const scalar N2New = NNew - N1New;

    const scalar N2to = Nto - n1ToReservoir;
    const scalar N2from = Nfrom - n1FromReservoir;

    const scalar dN2Soll = dt*(N2to - N2from);
    const scalar dN2Ist = N2New - N2Old;

    const scalar xTo =
        Nto > minFlux_
      ? n1ToReservoir/Nto
      : xOld;

    const scalar xFrom =
        Nfrom > minFlux_
      ? n1FromReservoir/Nfrom
      : xOld;

    const fileName dir =
        mesh.time().globalPath()
       /"postProcessing"
       /"apparatusBalance"
       /mesh.time().startTime().name();

    mkDir(dir);

    const fileName file = dir/(patch().name() + ".dat");

    const bool newFile = !isFile(file);

    std::ofstream os(file.c_str(), std::ios::app);

    if (newFile)
    {
        os
            << "# Time dt "
            << "N_old x_old N1_old "
            << "N_to N_from N1_to N1_from "
            << "dN_soll dN_ist dN_error "
            << "N_soll N_ist "
            << "dN1_soll dN1_ist dN1_error "
            << "N1_soll N1_ist "
            << "N2_soll N2_ist dN2_error "
            << "L_to G_to L_from G_from "
            << "x_to x_from"
            << std::endl;
    }

    os
        << mesh.time().value() << " "
        << dt << " "
        << NOld << " "
        << xOld << " "
        << N1Old << " "
        << Nto << " "
        << Nfrom << " "
        << n1ToReservoir << " "
        << n1FromReservoir << " "
        << dNSoll << " "
        << dNIst << " "
        << dNIst - dNSoll << " "
        << NSoll << " "
        << NNew << " "
        << dN1Soll << " "
        << dN1Ist << " "
        << dN1Ist - dN1Soll << " "
        << N1Soll << " "
        << N1New << " "
        << N2Old + dN2Soll << " "
        << N2New << " "
        << dN2Ist - dN2Soll << " "
        << nLtoReservoir << " "
        << nGtoReservoir << " "
        << nLfromReservoir << " "
        << nGfromReservoir << " "
        << xTo << " "
        << xFrom
        << std::endl;
}

tmp<scalarField>
inletOutletBoilerFvPatchScalarField::reservoirInletFraction() const
{
    tmp<scalarField> tfrac(new scalarField(patch().size(), scalar(0)));
    scalarField& frac = tfrac.ref();

    if (!compositionFluxFieldsAvailable())
    {
        return tfrac;
    }

    const surfaceScalarField& liquidFlux =
        lookupSurfaceFlux(liquidFluxName_);

    const surfaceScalarField& gasFlux =
        lookupSurfaceFlux(gasFluxName_);

    const label patchi = patch().index();

    const fvsPatchScalarField& Lp =
        liquidFlux.boundaryField()[patchi];

    const fvsPatchScalarField& Gp =
        gasFlux.boundaryField()[patchi];

    forAll(frac, facei)
    {
        const scalar Lfrom = max(-Lp[facei], scalar(0));
        const scalar Gfrom = max(-Gp[facei], scalar(0));

        frac[facei] = ((Lfrom + Gfrom) > minFlux_) ? scalar(1) : scalar(0);
    }

    return tfrac;
}

tmp<scalarField>
inletOutletBoilerFvPatchScalarField::reservoirInletValue() const
{
    tmp<scalarField> tval(new scalarField(patch().size(), reservoirY_));
    scalarField& val = tval.ref();

    if (!compositionFluxFieldsAvailable())
    {
        return tval;
    }

    const surfaceScalarField& liquidFlux =
        lookupSurfaceFlux(liquidFluxName_);

    const surfaceScalarField& gasFlux =
        lookupSurfaceFlux(gasFluxName_);

    const label patchi = patch().index();

    const fvsPatchScalarField& Lp =
        liquidFlux.boundaryField()[patchi];

    const fvsPatchScalarField& Gp =
        gasFlux.boundaryField()[patchi];

    const scalar xB = reservoirX_;
    const scalar yB = reservoirY_;

    forAll(val, facei)
    {
        const scalar Lfrom = max(-Lp[facei], scalar(0));
        const scalar Gfrom = max(-Gp[facei], scalar(0));
        const scalar totalFrom = Lfrom + Gfrom;

        if (totalFrom > minFlux_)
        {
            val[facei] = clamp((Lfrom*xB + Gfrom*yB)/totalFrom);
        }
        else
        {
            val[facei] = yB;
        }
    }

    return tval;
}


void inletOutletBoilerFvPatchScalarField::correctReservoir
(
    const scalar deltaN,
    const scalar deltaN1
)
{
    const scalar N1old = reservoirMoles_ * reservoirX_;
    reservoirMoles_ = max(reservoirMoles_ + deltaN, minReservoirMoles_);
    if (reservoirMoles_ > SMALL)
    {
        reservoirX_ = clamp((N1old + deltaN1) / reservoirMoles_);
        reservoirY_ = yVLE(reservoirX_, lastA12_);
    }
}


void inletOutletBoilerFvPatchScalarField::applyInletValuesOnNCPCPatches
(
    volScalarField& X
) const
{
    if (!Pstream::parRun() || !compositionFluxFieldsAvailable())
    {
        return;
    }

    const wordList patchNames = selectedPatchNames();
    const wordHashSet nameSet(patchNames);

    const surfaceScalarField& liquidFlux = lookupSurfaceFlux(liquidFluxName_);
    const surfaceScalarField& gasFlux    = lookupSurfaceFlux(gasFluxName_);

    const scalar xB = reservoirX_;
    const scalar yB = reservoirY_;

    const fvBoundaryMesh& bm = patch().boundaryMesh();

    label nPatches = 0;
    label nFacesIn = 0;
    label nFacesOut = 0;

    forAll(bm, pi)
    {
        if (!isA<nonConformalProcessorCyclicFvPatch>(bm[pi]))
        {
            continue;
        }

        const auto& ncpc =
            refCast<const nonConformalProcessorCyclicFvPatch>(bm[pi]);

        const word refName =
            ncpc.nonConformalProcessorCyclicPatch().referPatch().name();

        if (!nameSet.found(refName))
        {
            continue;
        }

        ++nPatches;

        const fvsPatchScalarField& Lp = liquidFlux.boundaryField()[pi];
        const fvsPatchScalarField& Gp = gasFlux.boundaryField()[pi];

        const scalarField xInternal(X.boundaryField()[pi].patchInternalField());
        fvPatchField<scalar>& Xp = X.boundaryFieldRef()[pi];

        // const scalar LpMin = min(Lp);
        // const scalar LpMax = max(Lp);
        // const scalar GpMin = min(Gp);
        // const scalar GpMax = max(Gp);

        // Pout<< "inletOutletBoiler::applyNCPC patch=" << bm[pi].name()
        //     << " refPatch=" << refName
        //     << " size=" << Xp.size()
        //     << " Lp=[" << LpMin << "," << LpMax << "]"
        //     << " Gp=[" << GpMin << "," << GpMax << "]"
        //     << " xB=" << xB << " yB=" << yB << nl;

        forAll(Xp, facei)
        {
            const scalar Lfrom = max(-Lp[facei], scalar(0));
            const scalar Gfrom = max(-Gp[facei], scalar(0));
            const scalar totalFrom = Lfrom + Gfrom;

            if (totalFrom > minFlux_)
            {
                Xp[facei] = clamp((Lfrom*xB + Gfrom*yB)/totalFrom);
                ++nFacesIn;
            }
            else
            {
                Xp[facei] = xInternal[facei];
                ++nFacesOut;
            }
        }
    }

    // if (nPatches > 0 || true)
    // {
    //     Pout<< "inletOutletBoiler::applyNCPC " << patch().name()
    //         << " nPatches=" << nPatches
    //         << " nFacesIn=" << nFacesIn
    //         << " nFacesOut=" << nFacesOut
    //         << " xB=" << xB << " yB=" << yB << nl;
    // }
}


inletOutletBoilerFvPatchScalarField::
inletOutletBoilerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    liquidFluxName_("alphaCLPhi1"),
    gasFluxName_("alphaCGPhi2"),
    A12Name_("A12"),
    groupName_(word()),
    patchNames_(),
    usePatchList_(false),
    reservoirMoles_(1),
    reservoirX_(0.5),
    reservoirY_(0.5),
    relax_(1),
    minValue_(0),
    maxValue_(1),
    minReservoirMoles_(SMALL),
    minFlux_(VSMALL),
    log_(true),
    lastReservoirUpdateIndex_(-1),
    lastA12_(1),
    lastLiquidToReservoir_(0),
    lastGasToReservoir_(0),
    lastLiquidFromReservoir_(0),
    lastGasFromReservoir_(0),
    lastSpeciesToReservoir_(0),
    lastSpeciesFromReservoir_(0),
    lastMolesResidual_(0),
    lastSpeciesResidual_(0),
    lastLiquidXToReservoir_(0),
    lastGasXToReservoir_(0),
    writeBalance_(true),
    useRawSpeciesFlux_(true),
    UName_("U"),
    alphaLiquidName_("alpha.liquid"),
    rhoLiquidName_("rho.liquid"),
    rhoGasName_("rho.gas"),
    species1LiquidName_("cyclohexane.liquid"),
    species2LiquidName_("nHeptane.liquid"),
    species1GasName_("cyclohexane.gas"),
    species2GasName_("nHeptane.gas"),
    W1_(84.16/1000.0),
    W2_(100.2/1000.0)
{
    refValue() = reservoirY_;
    refGrad() = 0;
    valueFraction() = 0;
    fvPatchScalarField::operator=(reservoirY_);
}


inletOutletBoilerFvPatchScalarField::
inletOutletBoilerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    liquidFluxName_(dict.lookupOrDefault<word>("liquidFlux", "alphaCLPhi1")),
    gasFluxName_(dict.lookupOrDefault<word>("gasFlux", "alphaCGPhi2")),
    A12Name_(dict.lookupOrDefault<word>("A12", "A12")),
    groupName_(dict.lookupOrDefault<word>("group", word())),
    patchNames_(),
    usePatchList_(false),
    reservoirMoles_(dict.lookupOrDefault<scalar>("reservoirMoles", 1)),
    reservoirX_(dict.lookupOrDefault<scalar>("reservoirX", 0.5)),
    reservoirY_(dict.lookupOrDefault<scalar>("reservoirY", reservoirX_)),
    relax_(dict.lookupOrDefault<scalar>("relax", 1)),
    minValue_(dict.lookupOrDefault<scalar>("minValue", 0)),
    maxValue_(dict.lookupOrDefault<scalar>("maxValue", 1)),
    minReservoirMoles_(dict.lookupOrDefault<scalar>("minReservoirMoles", SMALL)),
    minFlux_(dict.lookupOrDefault<scalar>("minFlux", VSMALL)),
    log_(dict.lookupOrDefault<Switch>("log", true)),
    lastReservoirUpdateIndex_(-1),
    lastA12_(dict.lookupOrDefault<scalar>("lastA12", 1)),
    lastLiquidToReservoir_
    (
        dict.lookupOrDefault<scalar>("lastLiquidToReservoir", 0)
    ),
    lastGasToReservoir_
    (
        dict.lookupOrDefault<scalar>("lastGasToReservoir", 0)
    ),
    lastLiquidFromReservoir_
    (
        dict.lookupOrDefault<scalar>("lastLiquidFromReservoir", 0)
    ),
    lastGasFromReservoir_
    (
        dict.lookupOrDefault<scalar>("lastGasFromReservoir", 0)
    ),
    lastSpeciesToReservoir_
    (
        dict.lookupOrDefault<scalar>("lastSpeciesToReservoir", 0)
    ),
    lastSpeciesFromReservoir_
    (
        dict.lookupOrDefault<scalar>("lastSpeciesFromReservoir", 0)
    ),
    lastMolesResidual_
    (
        dict.lookupOrDefault<scalar>("lastMolesResidual", 0)
    ),
    lastSpeciesResidual_
    (
        dict.lookupOrDefault<scalar>("lastSpeciesResidual", 0)
    ),
    lastLiquidXToReservoir_
    (
        dict.lookupOrDefault<scalar>("lastLiquidXToReservoir", reservoirX_)
    ),
    lastGasXToReservoir_
    (
        dict.lookupOrDefault<scalar>("lastGasXToReservoir", reservoirY_)
    ),
    writeBalance_(dict.lookupOrDefault<Switch>("writeBalance", true)),
    useRawSpeciesFlux_(dict.lookupOrDefault<Switch>("useRawSpeciesFlux", true)),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    alphaLiquidName_(dict.lookupOrDefault<word>("alphaLiquid", "alpha.liquid")),
    rhoLiquidName_(dict.lookupOrDefault<word>("rhoLiquid", "rho.liquid")),
    rhoGasName_(dict.lookupOrDefault<word>("rhoGas", "rho.gas")),
    species1LiquidName_(dict.lookupOrDefault<word>("species1Liquid", "cyclohexane.liquid")),
    species2LiquidName_(dict.lookupOrDefault<word>("species2Liquid", "nHeptane.liquid")),
    species1GasName_(dict.lookupOrDefault<word>("species1Gas", "cyclohexane.gas")),
    species2GasName_(dict.lookupOrDefault<word>("species2Gas", "nHeptane.gas")),
    W1_(dict.lookupOrDefault<scalar>("W1", 84.16)/1000.0),
    W2_(dict.lookupOrDefault<scalar>("W2", 100.2)/1000.0)
{
    readPatchSelection(dict);

    reservoirX_ = clamp(reservoirX_);
    reservoirY_ = clamp(reservoirY_);
    relax_ = min(max(relax_, scalar(0)), scalar(1));
    minReservoirMoles_ = max(minReservoirMoles_, VSMALL);
    reservoirMoles_ = max(reservoirMoles_, minReservoirMoles_);

    refValue() = reservoirY_;
    refGrad() = 0;
    valueFraction() = 0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", iF.dimensions(), dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(reservoirY_);
    }
}


inletOutletBoilerFvPatchScalarField::
inletOutletBoilerFvPatchScalarField
(
    const inletOutletBoilerFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    liquidFluxName_(ptf.liquidFluxName_),
    gasFluxName_(ptf.gasFluxName_),
    A12Name_(ptf.A12Name_),
    groupName_(ptf.groupName_),
    patchNames_(ptf.patchNames_),
    usePatchList_(ptf.usePatchList_),
    reservoirMoles_(ptf.reservoirMoles_),
    reservoirX_(ptf.reservoirX_),
    reservoirY_(ptf.reservoirY_),
    relax_(ptf.relax_),
    minValue_(ptf.minValue_),
    maxValue_(ptf.maxValue_),
    minReservoirMoles_(ptf.minReservoirMoles_),
    minFlux_(ptf.minFlux_),
    log_(ptf.log_),
    lastReservoirUpdateIndex_(ptf.lastReservoirUpdateIndex_),
    lastA12_(ptf.lastA12_),
    lastLiquidToReservoir_(ptf.lastLiquidToReservoir_),
    lastGasToReservoir_(ptf.lastGasToReservoir_),
    lastLiquidFromReservoir_(ptf.lastLiquidFromReservoir_),
    lastGasFromReservoir_(ptf.lastGasFromReservoir_),
    lastSpeciesToReservoir_(ptf.lastSpeciesToReservoir_),
    lastSpeciesFromReservoir_(ptf.lastSpeciesFromReservoir_),
    lastMolesResidual_(ptf.lastMolesResidual_),
    lastSpeciesResidual_(ptf.lastSpeciesResidual_),
    lastLiquidXToReservoir_(ptf.lastLiquidXToReservoir_),
    lastGasXToReservoir_(ptf.lastGasXToReservoir_),
    writeBalance_(ptf.writeBalance_),
    useRawSpeciesFlux_(ptf.useRawSpeciesFlux_),
    UName_(ptf.UName_),
    alphaLiquidName_(ptf.alphaLiquidName_),
    rhoLiquidName_(ptf.rhoLiquidName_),
    rhoGasName_(ptf.rhoGasName_),
    species1LiquidName_(ptf.species1LiquidName_),
    species2LiquidName_(ptf.species2LiquidName_),
    species1GasName_(ptf.species1GasName_),
    species2GasName_(ptf.species2GasName_),
    W1_(ptf.W1_),
    W2_(ptf.W2_)
{}


inletOutletBoilerFvPatchScalarField::
inletOutletBoilerFvPatchScalarField
(
    const inletOutletBoilerFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    liquidFluxName_(ptf.liquidFluxName_),
    gasFluxName_(ptf.gasFluxName_),
    A12Name_(ptf.A12Name_),
    groupName_(ptf.groupName_),
    patchNames_(ptf.patchNames_),
    usePatchList_(ptf.usePatchList_),
    reservoirMoles_(ptf.reservoirMoles_),
    reservoirX_(ptf.reservoirX_),
    reservoirY_(ptf.reservoirY_),
    relax_(ptf.relax_),
    minValue_(ptf.minValue_),
    maxValue_(ptf.maxValue_),
    minReservoirMoles_(ptf.minReservoirMoles_),
    minFlux_(ptf.minFlux_),
    log_(ptf.log_),
    lastReservoirUpdateIndex_(ptf.lastReservoirUpdateIndex_),
    lastA12_(ptf.lastA12_),
    lastLiquidToReservoir_(ptf.lastLiquidToReservoir_),
    lastGasToReservoir_(ptf.lastGasToReservoir_),
    lastLiquidFromReservoir_(ptf.lastLiquidFromReservoir_),
    lastGasFromReservoir_(ptf.lastGasFromReservoir_),
    lastSpeciesToReservoir_(ptf.lastSpeciesToReservoir_),
    lastSpeciesFromReservoir_(ptf.lastSpeciesFromReservoir_),
    lastMolesResidual_(ptf.lastMolesResidual_),
    lastSpeciesResidual_(ptf.lastSpeciesResidual_),
    lastLiquidXToReservoir_(ptf.lastLiquidXToReservoir_),
    lastGasXToReservoir_(ptf.lastGasXToReservoir_),
    writeBalance_(ptf.writeBalance_),
    useRawSpeciesFlux_(ptf.useRawSpeciesFlux_),
    UName_(ptf.UName_),
    alphaLiquidName_(ptf.alphaLiquidName_),
    rhoLiquidName_(ptf.rhoLiquidName_),
    rhoGasName_(ptf.rhoGasName_),
    species1LiquidName_(ptf.species1LiquidName_),
    species2LiquidName_(ptf.species2LiquidName_),
    species1GasName_(ptf.species1GasName_),
    species2GasName_(ptf.species2GasName_),
    W1_(ptf.W1_),
    W2_(ptf.W2_)
{}


void inletOutletBoilerFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Reservoir is advanced by compositionPredictor AFTER the composition
    // solve via advanceReservoir(), using post-solve X1.  Here we only set
    // the BC coefficients from the current (previous-step) reservoir state.
    refValue() = reservoirInletValue();
    refGrad() = 0;

    // The BC behaves as inletOutlet for both phases:
    // phases entering the CFD domain receive the flux-weighted reservoir
    // composition; phases leaving the CFD domain use zeroGradient.
    valueFraction() = reservoirInletFraction();

    mixedFvPatchScalarField::updateCoeffs();
}


void inletOutletBoilerFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntryIfDifferent<word>
    (
        os,
        "liquidFlux",
        "alphaCLPhi1",
        liquidFluxName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        "gasFlux",
        "alphaCGPhi2",
        gasFluxName_
    );

    writeEntryIfDifferent<word>(os, "A12", "A12", A12Name_);

    if (!groupName_.empty())
    {
        writeEntry(os, "group", groupName_);
    }

    if (usePatchList_)
    {
        writeEntry(os, "patches", patchNames_);
    }

    writeEntry(os, "reservoirMoles", reservoirMoles_);
    writeEntry(os, "reservoirX", reservoirX_);
    writeEntry(os, "reservoirY", reservoirY_);

    writeEntryIfDifferent<scalar>(os, "relax", 1, relax_);
    writeEntryIfDifferent<scalar>(os, "minValue", 0, minValue_);
    writeEntryIfDifferent<scalar>(os, "maxValue", 1, maxValue_);

    writeEntryIfDifferent<scalar>
    (
        os,
        "minReservoirMoles",
        SMALL,
        minReservoirMoles_
    );

    writeEntryIfDifferent<scalar>(os, "minFlux", VSMALL, minFlux_);
    writeEntryIfDifferent<Switch>(os, "log", true, log_);
    writeEntryIfDifferent<Switch>(os, "writeBalance", true, writeBalance_);

    // Restart/diagnostic entries
    writeEntry(os, "lastA12", lastA12_);
    writeEntry(os, "lastLiquidToReservoir", lastLiquidToReservoir_);
    writeEntry(os, "lastGasToReservoir", lastGasToReservoir_);
    writeEntry(os, "lastLiquidFromReservoir", lastLiquidFromReservoir_);
    writeEntry(os, "lastGasFromReservoir", lastGasFromReservoir_);
    writeEntry(os, "lastSpeciesToReservoir", lastSpeciesToReservoir_);
    writeEntry(os, "lastSpeciesFromReservoir", lastSpeciesFromReservoir_);
    writeEntry(os, "lastMolesResidual", lastMolesResidual_);
    writeEntry(os, "lastSpeciesResidual", lastSpeciesResidual_);
    writeEntry(os, "lastLiquidXToReservoir", lastLiquidXToReservoir_);
    writeEntry(os, "lastGasXToReservoir", lastGasXToReservoir_);



    writeEntryIfDifferent<word>
    (
        os,
        "alphaLiquid",
        "alpha.liquid",
        alphaLiquidName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        "rhoLiquid",
        "rho.liquid",
        rhoLiquidName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        "rhoGas",
        "rho.gas",
        rhoGasName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        "species1Liquid",
        "cyclohexane.liquid",
        species1LiquidName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        "species2Liquid",
        "nHeptane.liquid",
        species2LiquidName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        "species1Gas",
        "cyclohexane.gas",
        species1GasName_
    );

    writeEntryIfDifferent<word>
    (
        os,
        "species2Gas",
        "nHeptane.gas",
        species2GasName_
    );
    writeEntryIfDifferent<Switch>
    (
        os,
        "useRawSpeciesFlux",
        true,
        useRawSpeciesFlux_
    );
    
    writeEntryIfDifferent<word>(os, "U", "U", UName_);

    writeEntryIfDifferent<scalar>(os, "W1", 84.16, W1_*1000.0);
    writeEntryIfDifferent<scalar>(os, "W2", 100.2, W2_*1000.0);
    writeEntry(os, "value", *this);

}


void inletOutletBoilerFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchScalarField::operator=
    (
        valueFraction()*refValue()
      + (scalar(1) - valueFraction())*ptf
    );
}


makePatchTypeField
(
    fvPatchScalarField,
    inletOutletBoilerFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //


