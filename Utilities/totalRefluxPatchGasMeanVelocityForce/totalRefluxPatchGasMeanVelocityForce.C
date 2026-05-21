/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
   \\  /    A nd           |
    \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "totalRefluxPatchGasMeanVelocityForce.H"

#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "IFstream.H"
#include "IOdictionary.H"
#include "timeIOdictionary.H"
#include "zero.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(totalRefluxPatchGasMeanVelocityForce, 0);

    addToRunTimeSelectionTable
    (
        fvConstraint,
        totalRefluxPatchGasMeanVelocityForce,
        dictionary
    );
}
}


Foam::scalar Foam::fv::totalRefluxPatchGasMeanVelocityForce::limitedAlpha
(
    const scalar alpha
)
{
    return max(scalar(0), min(scalar(1), alpha));
}


bool Foam::fv::totalRefluxPatchGasMeanVelocityForce::isGasLikeName
(
    const word& name
)
{
    return
        name == "air"
     || name == "gas"
     || name == "vapour"
     || name == "vapor"
     || name == "steam";
}


bool Foam::fv::totalRefluxPatchGasMeanVelocityForce::isAlphaFieldName
(
    const word& name
)
{
    if (name == "alpha1")
    {
        return true;
    }

    if (name.find("alpha.") != 0)
    {
        return false;
    }

    if (name.size() > 2 && name.substr(name.size() - 2) == "_0")
    {
        return false;
    }

    return true;
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::readCoeffs
(
    const dictionary& dict
)
{
    if (dict.found("patches"))
    {
        patches_ = dict.lookup<wordList>("patches");
    }
    else if (dict.found("patch"))
    {
        patches_.setSize(1);
        patches_[0] = dict.lookup<word>("patch");
    }
    else
    {
        FatalErrorInFunction
            << "Either 'patch' or 'patches' must be specified"
            << exit(FatalError);
    }

    forAll(patches_, patchNamei)
    {
        if (mesh().boundaryMesh().findIndex(patches_[patchNamei]) < 0)
        {
            FatalErrorInFunction
                << "Cannot find patch " << patches_[patchNamei]
                << exit(FatalError);
        }
    }

    UName_ = dict.lookupOrDefault<word>("U", "U");

    // Optional overrides. With no alpha specified, the code looks for
    // alpha.liquid, alpha.water, alpha1, then the first non-gas alpha.* field.
    alphaName_ = dict.lookupOrDefault<word>("alpha", word::null);
    alphaPhiName_ = dict.lookupOrDefault<word>("alphaPhi", word::null);
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    // Robust default: use alpha*U projected onto gasDirection instead of
    // solver-internal alphaPhi, which can contain very large intermediate
    // values on NCC/coupled patches during pressure/alpha corrections.
    liquidFluxMode_ =
        dict.lookupOrDefault<word>("liquidFluxMode", word("directionalAlphaU"));

    if (liquidFluxMode_ != "directionalAlphaU" && liquidFluxMode_ != "alphaPhi")
    {
        FatalErrorInFunction
            << "Unsupported liquidFluxMode " << liquidFluxMode_ << nl
            << "Use directionalAlphaU or alphaPhi."
            << exit(FatalError);
    }

    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    rhoLiquidFieldName_ =
        dict.lookupOrDefault<word>("rhoLiquidField", "rho.liquid");
    rhoGasFieldName_ =
        dict.lookupOrDefault<word>("rhoGasField", "rho.gas");

    molarTotalReflux_ =
        dict.lookupOrDefault<Switch>("molarTotalReflux", true);

    species1LiquidName_ =
        dict.lookupOrDefault<word>("species1Liquid", "cyclohexane.liquid");
    species2LiquidName_ =
        dict.lookupOrDefault<word>("species2Liquid", "nHeptane.liquid");
    species1GasName_ =
        dict.lookupOrDefault<word>("species1Gas", "cyclohexane.gas");
    species2GasName_ =
        dict.lookupOrDefault<word>("species2Gas", "nHeptane.gas");

    W1_ = dict.lookupOrDefault<scalar>("W1", 84.16)/1000.0;
    W2_ = dict.lookupOrDefault<scalar>("W2", 100.2)/1000.0;

    // Optional constants. If they are not given, rho.<phase> fields are tried;
    // if they are not available, the mixture field rho is used to estimate the
    // pure-phase densities from rho(alpha).
    rhoLiquid_ = dict.lookupOrDefault<scalar>("rhoLiquid", -1);
    rhoGas_ = dict.lookupOrDefault<scalar>("rhoGas", -1);

    densityLiquidAlphaCutoff_ =
        dict.lookupOrDefault<scalar>("densityLiquidAlphaCutoff", 0.99);
    densityGasAlphaCutoff_ =
        dict.lookupOrDefault<scalar>("densityGasAlphaCutoff", 0.01);

    gasDirection_ = dict.lookupOrDefault<vector>("gasDirection", vector(0, 1, 0));

    if (mag(gasDirection_) <= SMALL)
    {
        FatalErrorInFunction
            << "gasDirection must be non-zero"
            << exit(FatalError);
    }

    gasDirection_ /= mag(gasDirection_);

    // Formula 4.26/4.27 usually uses 0.01 s. A value <=0 gives instantaneous
    // B_L = Q_L/A_0 with no averaging window.
    averagingInterval_ = dict.lookupOrDefault<scalar>("averagingInterval", 0.01);

    relaxation_ = dict.lookupOrDefault<scalar>("relaxation", 0.05);

    // User requirement: cells with little gas shall be counted in the
    // measurement but receive no momentum source.  A smooth ramp avoids
    // source discontinuities at the interface.
    alphaCutoff_ = dict.lookupOrDefault<scalar>("alphaCutoff", 0.2);
    useCutoff_ = dict.lookupOrDefault<Switch>("useCutoff", true);
    alphaRamp_ = dict.lookupOrDefault<scalar>("alphaRamp", 0.1);

    // Total reflux: liquid load is positive when liquid moves opposite to
    // gasDirection.  For alphaPhi mode the default remains the old signed
    // patch-normal behaviour unless overridden.
    const bool defaultCounterCurrent = (liquidFluxMode_ != "alphaPhi");
    counterCurrent_ =
        dict.lookupOrDefault<Switch>("counterCurrent", defaultCounterCurrent);

    // Negative values mean unlimited. These are optional numerical safety caps.
    maxGasLoading_ = dict.lookupOrDefault<scalar>("maxGasLoading", -1);
    maxPressureGradient_ =
        dict.lookupOrDefault<scalar>("maxPressureGradient", -1);
    maxDeltaPressureGradient_ =
        dict.lookupOrDefault<scalar>("maxDeltaPressureGradient", -1);

    resetGradient_ = dict.lookupOrDefault<Switch>("resetGradient", false);

    // Window-integral balance controller.  The base target is still computed
    // from the instantaneous/filtered reflux, but this controller integrates
    // the mismatch
    //
    //     eDot = mDown - mGasUp
    //
    // over balanceWindow and adds a superficial gas-loading correction to
    // recover the accumulated mass error over the following windows.
    balanceControl_ =
        dict.lookupOrDefault<Switch>("balanceControl", true);

    balanceWindow_ =
        dict.lookupOrDefault<scalar>("balanceWindow", 0.01);

    balanceTolerance_ =
        dict.lookupOrDefault<scalar>("balanceTolerance", 1e-8);

    balanceRelativeTolerance_ =
        dict.lookupOrDefault<scalar>("balanceRelativeTolerance", 1e-4);

    balanceIntegralGain_ =
        dict.lookupOrDefault<scalar>("balanceIntegralGain", 0.25);

    maxBalanceCorrection_ =
        dict.lookupOrDefault<scalar>("maxBalanceCorrection", -1);

    maxBalanceCorrectionFraction_ =
        dict.lookupOrDefault<scalar>("maxBalanceCorrectionFraction", -1);
}


Foam::word Foam::fv::totalRefluxPatchGasMeanVelocityForce::actualAlphaName() const
{
    if (alphaName_ != word::null)
    {
        return alphaName_;
    }

    if (mesh().foundObject<volScalarField>("alpha.liquid"))
    {
        return "alpha.liquid";
    }

    if (mesh().foundObject<volScalarField>("alpha.water"))
    {
        return "alpha.water";
    }

    if (mesh().foundObject<volScalarField>("alpha1"))
    {
        return "alpha1";
    }

    const wordList alphaNames(mesh().toc<volScalarField>());

    // Prefer a non-gas alpha.* field. This keeps alpha.air/alpha.gas from being
    // selected as liquid when both alpha.water and alpha.air exist.
    forAll(alphaNames, i)
    {
        const word& fld = alphaNames[i];

        if (isAlphaFieldName(fld))
        {
            const word phase(IOobject::group(fld));

            if (!isGasLikeName(phase))
            {
                return fld;
            }
        }
    }

    forAll(alphaNames, i)
    {
        if (isAlphaFieldName(alphaNames[i]))
        {
            return alphaNames[i];
        }
    }

    FatalErrorInFunction
        << "Cannot auto-detect a liquid alpha field for " << name() << nl
        << "Set 'alpha' explicitly, e.g. alpha alpha.water;" << nl
        << "Available volScalarField objects are " << alphaNames
        << exit(FatalError);

    return word::null;
}


Foam::word Foam::fv::totalRefluxPatchGasMeanVelocityForce::actualAlphaPhiName
(
    const volScalarField& alpha1
) const
{
    if (alphaPhiName_ != word::null)
    {
        return alphaPhiName_;
    }

    if (alpha1.group() != word::null)
    {
        const word phaseAlphaPhiName
        (
            IOobject::groupName("alphaPhi", alpha1.group())
        );

        if (mesh().foundObject<surfaceScalarField>(phaseAlphaPhiName))
        {
            return phaseAlphaPhiName;
        }
    }

    if (mesh().foundObject<surfaceScalarField>("alphaPhi1"))
    {
        return "alphaPhi1";
    }

    if (mesh().foundObject<surfaceScalarField>("alphaPhi"))
    {
        return "alphaPhi";
    }

    return word::null;
}


Foam::word Foam::fv::totalRefluxPatchGasMeanVelocityForce::actualGasPhaseName
(
    const volScalarField& alpha1
) const
{
    if (mesh().foundObject<volScalarField>("alpha.gas"))
    {
        return "gas";
    }

    if (mesh().foundObject<volScalarField>("alpha.air"))
    {
        return "air";
    }

    const wordList alphaNames(mesh().toc<volScalarField>());

    forAll(alphaNames, i)
    {
        const word& fld = alphaNames[i];

        if (isAlphaFieldName(fld) && fld != alpha1.name())
        {
            const word phase(IOobject::group(fld));

            if (isGasLikeName(phase))
            {
                return phase;
            }
        }
    }

    forAll(alphaNames, i)
    {
        const word& fld = alphaNames[i];

        if (isAlphaFieldName(fld) && fld != alpha1.name())
        {
            return IOobject::group(fld);
        }
    }

    return word::null;
}


Foam::word Foam::fv::totalRefluxPatchGasMeanVelocityForce::actualGasAlphaPhiName
(
    const volScalarField& alpha1
) const
{
    const word gasPhase(actualGasPhaseName(alpha1));

    if (gasPhase != word::null)
    {
        const word phaseAlphaPhiName
        (
            IOobject::groupName("alphaPhi", gasPhase)
        );

        if (mesh().foundObject<surfaceScalarField>(phaseAlphaPhiName))
        {
            return phaseAlphaPhiName;
        }
    }

    if (mesh().foundObject<surfaceScalarField>("alphaPhi2"))
    {
        return "alphaPhi2";
    }

    return word::null;
}


Foam::scalar Foam::fv::totalRefluxPatchGasMeanVelocityForce::gasFluxWeight
(
    const scalar alpha1
) const
{
    // This is only for measuring superficial gas velocity.
    // Do not apply alphaCutoff_ here.
    return scalar(1) - limitedAlpha(alpha1);
}


Foam::scalar Foam::fv::totalRefluxPatchGasMeanVelocityForce::sourceWeight
(
    const scalar alpha1
) const
{
    const scalar aG = scalar(1) - limitedAlpha(alpha1);

    if (!useCutoff_)
    {
        return aG;
    }

    if (aG <= alphaCutoff_)
    {
        return 0;
    }

    if (alphaRamp_ <= SMALL)
    {
        return aG;
    }

    const scalar x = min
    (
        scalar(1),
        max(scalar(0), (aG - alphaCutoff_)/max(alphaRamp_, SMALL))
    );

    // Smoothstep ramp: zero slope at both ends.
    const scalar gate = x*x*(scalar(3) - scalar(2)*x);

    return aG*gate;
}


Foam::scalar Foam::fv::totalRefluxPatchGasMeanVelocityForce::sourceWeight
(
    const volScalarField& alpha1,
    const label celli
) const
{
    return sourceWeight(alpha1[celli]);
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::writeProps
(
    const scalar gradP
) const
{
    if (mesh().time().writeTime())
    {
        timeIOdictionary propsDict
        (
            IOobject
            (
                name() + "Properties",
                mesh().time().name(),
                "uniform",
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );

        propsDict.add("gradient", gradP);
        propsDict.add("targetUbar", targetUbar_);
        propsDict.add("gasDirection", gasDirection_);
        propsDict.add("liquidFlowRate", liquidFlowRate_);
        propsDict.add("patchArea", patchArea_);
        propsDict.add("liquidLoading", liquidLoading_);
        propsDict.add("accumulatedEquivalentGasVolume", accumulatedEquivalentGasVolume_);
        propsDict.add("accumulatedTime", accumulatedTime_);
        propsDict.add("rhoLiquid", rhoLiquidPatchMean_);
        propsDict.add("rhoGas", rhoGasPatchMean_);
        propsDict.add("densitySource", densitySource_);
        propsDict.add("totalMassFlowRateDown", totalMassFlowRateDown_);
        propsDict.add("liquidVolumeFlowRateDown", liquidVolumeFlowRateDown_);
        propsDict.add("gasVolumeFlowRateDown", gasVolumeFlowRateDown_);
        propsDict.add("instantGasLoadingTarget", instantGasLoadingTarget_);
        propsDict.add("gasLoadingTarget", gasLoadingTarget_);
        propsDict.add("gasMassFlowRateUp", gasMassFlowRateUp_);
        propsDict.add("gasVolumeFlowRateUp", gasVolumeFlowRateUp_);
        propsDict.add("totalMolarFlowRateDown", totalMolarFlowRateDown_);
        propsDict.add("gasMolarFlowRateUp", gasMolarFlowRateUp_);
        propsDict.add("gasMolarConcentrationPatchMean", gasMolarConcentrationPatchMean_);
        propsDict.add("molarTotalReflux", molarTotalReflux_);
        propsDict.add("balanceErrorRate", balanceErrorRate_);
        propsDict.add("balanceWindowErrorMass", balanceWindowErrorMass_);
        propsDict.add("balanceWindowTime", balanceWindowTime_);
        propsDict.add("balanceCorrectionLoading", balanceCorrectionLoading_);

        propsDict.regIOobject::write();
    }
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::measureDownwardTotalReflux
(
    const volVectorField& U,
    const volScalarField& alpha1,
    scalar& totalMassFlowDown,
    scalar& liquidVolumeFlowDown,
    scalar& gasVolumeFlowDown,
    scalar& gasVolumeFlowUp,
    scalar& patchArea
) const
{
    const word alphaPhiLName(actualAlphaPhiName(alpha1));
    const word alphaPhiGName(actualGasAlphaPhiName(alpha1));

    const surfaceScalarField* alphaPhiLPtr = nullptr;
    const surfaceScalarField* alphaPhiGPtr = nullptr;
    const surfaceScalarField* phiPtr = nullptr;

    if (liquidFluxMode_ == "alphaPhi")
    {
        if (alphaPhiLName != word::null)
        {
            alphaPhiLPtr =
                &mesh().lookupObject<surfaceScalarField>(alphaPhiLName);
        }

        if (alphaPhiGName != word::null)
        {
            alphaPhiGPtr =
                &mesh().lookupObject<surfaceScalarField>(alphaPhiGName);
        }

        if
        (
            (!alphaPhiLPtr || !alphaPhiGPtr)
         && mesh().foundObject<surfaceScalarField>(phiName_)
        )
        {
            phiPtr = &mesh().lookupObject<surfaceScalarField>(phiName_);
        }

        if (!alphaPhiLPtr && !phiPtr)
        {
            FatalErrorInFunction
                << "Cannot find liquid alphaPhi field for " << alpha1.name()
                << " and cannot find fallback phi field " << phiName_ << nl
                << "Use liquidFluxMode directionalAlphaU for NCC patches, "
                << "or provide alphaPhi/phi explicitly." << nl
                << "Available surfaceScalarField objects are "
                << mesh().toc<surfaceScalarField>()
                << exit(FatalError);
        }
    }

    scalar sumMDown = 0;
    scalar sumQLDown = 0;
    scalar sumQGDown = 0;
    scalar sumQGUp = 0;
    scalar sumA = 0;

    forAll(patches_, patchNamei)
    {
        const word& patchName = patches_[patchNamei];
        const label patchi = mesh().boundaryMesh().findIndex(patchName);

        const scalarField& magSf = mesh().boundary()[patchi].magSf();
        const fvPatchScalarField& alpha1p = alpha1.boundaryField()[patchi];
        const fvPatchVectorField& Up = U.boundaryField()[patchi];

        const fvsPatchScalarField* alphaPhiLpPtr = nullptr;
        const fvsPatchScalarField* alphaPhiGpPtr = nullptr;
        const fvsPatchScalarField* phipPtr = nullptr;

        if (liquidFluxMode_ == "alphaPhi")
        {
            if (alphaPhiLPtr)
            {
                alphaPhiLpPtr = &alphaPhiLPtr->boundaryField()[patchi];
            }

            if (alphaPhiGPtr)
            {
                alphaPhiGpPtr = &alphaPhiGPtr->boundaryField()[patchi];
            }

            if (phiPtr)
            {
                phipPtr = &phiPtr->boundaryField()[patchi];
            }
        }

        forAll(magSf, facei)
        {
            const scalar A = magSf[facei];

            const scalar aL = limitedAlpha(alpha1p[facei]);
            const scalar aG = scalar(1) - aL;

            scalar phiL = 0;
            scalar phiG = 0;

            if (liquidFluxMode_ == "alphaPhi")
            {
                if (alphaPhiLpPtr)
                {
                    phiL = (*alphaPhiLpPtr)[facei];
                }
                else
                {
                    phiL = aL*(*phipPtr)[facei];
                }

                if (alphaPhiGpPtr)
                {
                    phiG = (*alphaPhiGpPtr)[facei];
                }
                else if (phipPtr)
                {
                    phiG = (*phipPtr)[facei] - phiL;
                }
                else
                {
                    phiG = 0;
                }
            }
            else
            {
                // Use the same patch-normal volume flux as the apparatus BCs:
                // positive U&Sf leaves the CFD domain into the apparatus.
                const scalar phi = Up[facei] & mesh().Sf().boundaryField()[patchi][facei];

                phiL = aL*phi;
                phiG = aG*phi;
            }

            // Apparatus convention: positive = column -> apparatus.
            // Gas return from apparatus to column is negative gas flux.
            const scalar qLDown = max( phiL, scalar(0));
            const scalar qGDown = max( phiG, scalar(0));
            const scalar qGUp   = max(-phiG, scalar(0));

            sumQLDown += qLDown;
            sumQGDown += qGDown;
            sumQGUp += qGUp;

            sumMDown +=
                rhoLiquidPatchMean_*qLDown
              + rhoGasPatchMean_*qGDown;

            sumA += A;
        }
    }

    reduce(sumMDown, sumOp<scalar>());
    reduce(sumQLDown, sumOp<scalar>());
    reduce(sumQGDown, sumOp<scalar>());
    reduce(sumQGUp, sumOp<scalar>());
    reduce(sumA, sumOp<scalar>());

    totalMassFlowDown = sumMDown;
    liquidVolumeFlowDown = sumQLDown;
    gasVolumeFlowDown = sumQGDown;
    gasVolumeFlowUp = sumQGUp;
    patchArea = sumA;
}



bool Foam::fv::totalRefluxPatchGasMeanVelocityForce::molarFluxFieldsAvailable() const
{
    return
        mesh().foundObject<volScalarField>(rhoLiquidFieldName_)
     && mesh().foundObject<volScalarField>(rhoGasFieldName_)
     && mesh().foundObject<volScalarField>(species1LiquidName_)
     && mesh().foundObject<volScalarField>(species2LiquidName_)
     && mesh().foundObject<volScalarField>(species1GasName_)
     && mesh().foundObject<volScalarField>(species2GasName_);
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::measureDownwardTotalRefluxMolar
(
    const volVectorField& U,
    const volScalarField& alpha1,
    scalar& totalMolarFlowDown,
    scalar& gasMolarFlowUp,
    scalar& liquidVolumeFlowDown,
    scalar& gasVolumeFlowDown,
    scalar& gasVolumeFlowUp,
    scalar& gasMolarConcentration,
    scalar& patchArea
) const
{
    if (!molarFluxFieldsAvailable())
    {
        FatalErrorInFunction
            << "molarTotalReflux requires the following volScalarFields:" << nl
            << "  rhoLiquidField=" << rhoLiquidFieldName_ << " found="
            << mesh().foundObject<volScalarField>(rhoLiquidFieldName_) << nl
            << "  rhoGasField=" << rhoGasFieldName_ << " found="
            << mesh().foundObject<volScalarField>(rhoGasFieldName_) << nl
            << "  species1Liquid=" << species1LiquidName_ << " found="
            << mesh().foundObject<volScalarField>(species1LiquidName_) << nl
            << "  species2Liquid=" << species2LiquidName_ << " found="
            << mesh().foundObject<volScalarField>(species2LiquidName_) << nl
            << "  species1Gas=" << species1GasName_ << " found="
            << mesh().foundObject<volScalarField>(species1GasName_) << nl
            << "  species2Gas=" << species2GasName_ << " found="
            << mesh().foundObject<volScalarField>(species2GasName_) << nl
            << exit(FatalError);
    }

    if (W1_ <= SMALL || W2_ <= SMALL)
    {
        FatalErrorInFunction
            << "W1/W2 must be positive. W1=" << W1_ << ", W2=" << W2_
            << exit(FatalError);
    }

    const volScalarField& rhoL = mesh().lookupObject<volScalarField>(rhoLiquidFieldName_);
    const volScalarField& rhoG = mesh().lookupObject<volScalarField>(rhoGasFieldName_);
    const volScalarField& Y1L = mesh().lookupObject<volScalarField>(species1LiquidName_);
    const volScalarField& Y2L = mesh().lookupObject<volScalarField>(species2LiquidName_);
    const volScalarField& Y1G = mesh().lookupObject<volScalarField>(species1GasName_);
    const volScalarField& Y2G = mesh().lookupObject<volScalarField>(species2GasName_);

    scalar sumNDown = 0;
    scalar sumNGasUp = 0;
    scalar sumQLDown = 0;
    scalar sumQGDown = 0;
    scalar sumQGUp = 0;
    scalar sumA = 0;
    scalar sumCGAG = 0;
    scalar sumAG = 0;

    forAll(patches_, patchNamei)
    {
        const word& patchName = patches_[patchNamei];
        const label patchi = mesh().boundaryMesh().findIndex(patchName);

        const scalarField& magSf = mesh().boundary()[patchi].magSf();
        const vectorField& Sfp = mesh().Sf().boundaryField()[patchi];
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvPatchScalarField& alpha1p = alpha1.boundaryField()[patchi];
        const fvPatchScalarField& rhoLp = rhoL.boundaryField()[patchi];
        const fvPatchScalarField& rhoGp = rhoG.boundaryField()[patchi];
        const fvPatchScalarField& Y1Lp = Y1L.boundaryField()[patchi];
        const fvPatchScalarField& Y2Lp = Y2L.boundaryField()[patchi];
        const fvPatchScalarField& Y1Gp = Y1G.boundaryField()[patchi];
        const fvPatchScalarField& Y2Gp = Y2G.boundaryField()[patchi];

        forAll(magSf, facei)
        {
            const scalar A = magSf[facei];
            const scalar aL = limitedAlpha(alpha1p[facei]);
            const scalar aG = scalar(1) - aL;
            const scalar phi = Up[facei] & Sfp[facei];

            const scalar CL = rhoLp[facei]*(Y1Lp[facei]/W1_ + Y2Lp[facei]/W2_);
            const scalar CG = rhoGp[facei]*(Y1Gp[facei]/W1_ + Y2Gp[facei]/W2_);

            const scalar qL = aL*phi;
            const scalar qG = aG*phi;
            const scalar nL = CL*qL;
            const scalar nG = CG*qG;

            sumNDown += max(nL, scalar(0)) + max(nG, scalar(0));
            sumNGasUp += max(-nG, scalar(0));

            sumQLDown += max(qL, scalar(0));
            sumQGDown += max(qG, scalar(0));
            sumQGUp += max(-qG, scalar(0));

            sumCGAG += aG*CG*A;
            sumAG += aG*A;
            sumA += A;
        }
    }

    reduce(sumNDown, sumOp<scalar>());
    reduce(sumNGasUp, sumOp<scalar>());
    reduce(sumQLDown, sumOp<scalar>());
    reduce(sumQGDown, sumOp<scalar>());
    reduce(sumQGUp, sumOp<scalar>());
    reduce(sumCGAG, sumOp<scalar>());
    reduce(sumAG, sumOp<scalar>());
    reduce(sumA, sumOp<scalar>());

    totalMolarFlowDown = sumNDown;
    gasMolarFlowUp = sumNGasUp;
    liquidVolumeFlowDown = sumQLDown;
    gasVolumeFlowDown = sumQGDown;
    gasVolumeFlowUp = sumQGUp;
    gasMolarConcentration = sumAG > SMALL ? sumCGAG/sumAG : 0;
    patchArea = sumA;
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::updateBalanceController
(
    const volVectorField& U,
    const volScalarField& alpha1
) const
{
    if (!balanceControl_)
    {
        return;
    }

    scalar eDot = 0;
    scalar reference = SMALL;
    scalar A0 = 0;
    scalar denom = SMALL;

    if (molarTotalReflux_)
    {
        scalar nDown = 0;
        scalar nGasUp = 0;
        scalar qLDown = 0;
        scalar qGDown = 0;
        scalar qGUp = 0;
        scalar CG = 0;

        measureDownwardTotalRefluxMolar
        (
            U,
            alpha1,
            nDown,
            nGasUp,
            qLDown,
            qGDown,
            qGUp,
            CG,
            A0
        );

        totalMolarFlowRateDown_ = nDown;
        gasMolarFlowRateUp_ = nGasUp;
        gasMolarConcentrationPatchMean_ = CG;
        gasVolumeFlowRateUp_ = qGUp;
        gasMassFlowRateUp_ = rhoGasPatchMean_*qGUp;

        eDot = nDown - nGasUp;
        reference = max(max(nDown, nGasUp), SMALL);
        denom = max(CG*A0, SMALL);
    }
    else
    {
        scalar mDown = 0;
        scalar qLDown = 0;
        scalar qGDown = 0;
        scalar qGUp = 0;

        measureDownwardTotalReflux
        (
            U,
            alpha1,
            mDown,
            qLDown,
            qGDown,
            qGUp,
            A0
        );

        const scalar mGasUp = rhoGasPatchMean_*qGUp;
        gasVolumeFlowRateUp_ = qGUp;
        gasMassFlowRateUp_ = mGasUp;

        eDot = mDown - mGasUp;
        reference = max(max(mDown, mGasUp), SMALL);
        denom = max(rhoGasPatchMean_*A0, SMALL);
    }

    balanceErrorRate_ = eDot;

    if (A0 <= SMALL || denom <= SMALL)
    {
        return;
    }

    const scalar dt = mesh().time().deltaTValue();
    const label timeIndex = mesh().time().timeIndex();

    if (dt <= SMALL || timeIndex == lastBalanceUpdateTimeIndex_)
    {
        return;
    }

    balanceWindowErrorMass_ += eDot*dt;
    balanceWindowTime_ += dt;
    lastBalanceUpdateTimeIndex_ = timeIndex;

    const scalar window = max(balanceWindow_, dt);

    if (balanceWindowTime_ < window)
    {
        return;
    }

    const scalar eDotAve = balanceWindowErrorMass_/max(balanceWindowTime_, SMALL);

    const scalar tolerance = max
    (
        balanceTolerance_,
        balanceRelativeTolerance_*reference
    );

    if (mag(eDotAve) > tolerance)
    {
        // eDot > 0 means more molar flow leaves the column than enters as gas.
        // Increase gas loading so the following window returns the inventory.
        const scalar dJ = balanceIntegralGain_*eDotAve/denom;

        balanceCorrectionLoading_ += dJ;

        // Relative cap: fraction of current gas loading target.
        // Takes priority over absolute cap when target > 0.
        if
        (
            maxBalanceCorrectionFraction_ > 0
         && gasLoadingTarget_ > SMALL
        )
        {
            const scalar relMax =
                maxBalanceCorrectionFraction_ * gasLoadingTarget_;

            balanceCorrectionLoading_ =
                balanceCorrectionLoading_ < 0
              ? max(-relMax, balanceCorrectionLoading_)
              : min( relMax, balanceCorrectionLoading_);
        }
        else if
        (
            maxBalanceCorrection_ > 0
         && mag(balanceCorrectionLoading_) > maxBalanceCorrection_
        )
        {
            balanceCorrectionLoading_ =
                balanceCorrectionLoading_ < 0
              ? -maxBalanceCorrection_
              :  maxBalanceCorrection_;
        }
    }

    Info<< "Total-reflux window balance:"
        << (molarTotalReflux_ ? " basis=molar" : " basis=mass")
        << ", errorAve = " << eDotAve
        << ", tolerance = " << tolerance
        << ", window = " << balanceWindowTime_
        << ", balanceCorrection jG = " << balanceCorrectionLoading_
        << " (cap: "
        << (
               maxBalanceCorrectionFraction_ > 0 && gasLoadingTarget_ > SMALL
             ? maxBalanceCorrectionFraction_*gasLoadingTarget_
             : maxBalanceCorrection_
           )
        << " m/s)"
        << endl;

    balanceWindowErrorMass_ = 0;
    balanceWindowTime_ = 0;
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::updateTotalRefluxLoading
(
    const volVectorField& U,
    const volScalarField& alpha1
) const
{
    scalar qGTarget = 0;
    scalar A0 = 0;

    if (molarTotalReflux_)
    {
        scalar nDown = 0;
        scalar nGasUp = 0;
        scalar qLDown = 0;
        scalar qGDown = 0;
        scalar qGUp = 0;
        scalar CG = 0;

        measureDownwardTotalRefluxMolar
        (
            U,
            alpha1,
            nDown,
            nGasUp,
            qLDown,
            qGDown,
            qGUp,
            CG,
            A0
        );

        totalMolarFlowRateDown_ = nDown;
        gasMolarFlowRateUp_ = nGasUp;
        gasMolarConcentrationPatchMean_ = CG;
        liquidVolumeFlowRateDown_ = qLDown;
        gasVolumeFlowRateDown_ = qGDown;
        gasVolumeFlowRateUp_ = qGUp;
        gasMassFlowRateUp_ = rhoGasPatchMean_*qGUp;
        balanceErrorRate_ = nDown - nGasUp;
        liquidFlowRate_ = qLDown;
        patchArea_ = A0;

        if (CG <= SMALL)
        {
            FatalErrorInFunction
                << "Gas total molar concentration is too small for molar total reflux."
                << exit(FatalError);
        }

        // Molar total reflux:
        //     NGasUp,target = NDown
        //     QGasUp,target = NDown / CG
        qGTarget = nDown/CG;
    }
    else
    {
        scalar mDown = 0;
        scalar qLDown = 0;
        scalar qGDown = 0;
        scalar qGUp = 0;

        measureDownwardTotalReflux(U, alpha1, mDown, qLDown, qGDown, qGUp, A0);

        totalMassFlowRateDown_ = mDown;
        liquidVolumeFlowRateDown_ = qLDown;
        gasVolumeFlowRateDown_ = qGDown;
        gasVolumeFlowRateUp_ = qGUp;
        gasMassFlowRateUp_ = rhoGasPatchMean_*qGUp;
        balanceErrorRate_ = mDown - gasMassFlowRateUp_;
        liquidFlowRate_ = qLDown;
        patchArea_ = A0;

        if (rhoGasPatchMean_ <= SMALL)
        {
            FatalErrorInFunction
                << "Gas density is too small for mass-based total reflux."
                << exit(FatalError);
        }

        qGTarget = mDown/rhoGasPatchMean_;
    }

    if (A0 <= SMALL)
    {
        WarningInFunction
            << "Selected patch area is too small for " << name()
            << ". Setting total-reflux target to zero."
            << endl;

        instantGasLoadingTarget_ = 0;
        gasLoadingTarget_ = 0;
        liquidLoading_ = 0;
        liquidLoadingInitialized_ = true;
        return;
    }

    const scalar jGbase = max(scalar(0), qGTarget/A0);
    const scalar jGinst = max(scalar(0), jGbase + balanceCorrectionLoading_);

    instantGasLoadingTarget_ = jGinst;

    const scalar dt = mesh().time().deltaTValue();
    const label timeIndex = mesh().time().timeIndex();

    if (!liquidLoadingInitialized_ || averagingInterval_ <= SMALL)
    {
        gasLoadingTarget_ = jGinst;
        liquidLoading_ = gasLoadingTarget_;

        liquidLoadingInitialized_ = true;
        lastIntegratedTimeIndex_ = timeIndex;
        accumulatedEquivalentGasVolume_ = 0;
        accumulatedTime_ = 0;

        return;
    }

    if (timeIndex != lastIntegratedTimeIndex_)
    {
        if (dt > SMALL)
        {
            const scalar beta =
                min(scalar(1), dt/max(averagingInterval_, SMALL));

            gasLoadingTarget_ =
                (scalar(1) - beta)*gasLoadingTarget_
              + beta*jGinst;

            liquidLoading_ = gasLoadingTarget_;

            accumulatedEquivalentGasVolume_ += qGTarget*dt;
            accumulatedTime_ += dt;
        }

        lastIntegratedTimeIndex_ = timeIndex;
    }
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::updateDensities
(
    const volScalarField& alpha1
) const
{
    if (rhoLiquid_ > 0 && rhoGas_ > 0)
    {
        rhoLiquidPatchMean_ = rhoLiquid_;
        rhoGasPatchMean_ = rhoGas_;
        densitySource_ = "constants";
        return;
    }

    word rhoLName(rhoLiquidFieldName_);
    word rhoGName(rhoGasFieldName_);

    if (rhoLName == word::null && alpha1.group() != word::null)
    {
        const word candidate(IOobject::groupName("rho", alpha1.group()));

        if (mesh().foundObject<volScalarField>(candidate))
        {
            rhoLName = candidate;
        }
    }

    if (rhoGName == word::null)
    {
        const word gasPhase(actualGasPhaseName(alpha1));

        if (gasPhase != word::null)
        {
            const word candidate(IOobject::groupName("rho", gasPhase));

            if (mesh().foundObject<volScalarField>(candidate))
            {
                rhoGName = candidate;
            }
        }
    }

    if
    (
        rhoLName != word::null
     && rhoGName != word::null
     && mesh().foundObject<volScalarField>(rhoLName)
     && mesh().foundObject<volScalarField>(rhoGName)
    )
    {
        const volScalarField& rhoL =
            mesh().lookupObject<volScalarField>(rhoLName);

        const volScalarField& rhoG =
            mesh().lookupObject<volScalarField>(rhoGName);

        scalar sumA = 0;
        scalar sumRhoLA = 0;
        scalar sumRhoGA = 0;

        forAll(patches_, patchNamei)
        {
            const word& patchName = patches_[patchNamei];
            const label patchi = mesh().boundaryMesh().findIndex(patchName);

            const scalarField& magSf = mesh().boundary()[patchi].magSf();
            const fvPatchScalarField& rhoLp = rhoL.boundaryField()[patchi];
            const fvPatchScalarField& rhoGp = rhoG.boundaryField()[patchi];

            forAll(magSf, facei)
            {
                const scalar A = magSf[facei];

                sumA += A;
                sumRhoLA += rhoLp[facei]*A;
                sumRhoGA += rhoGp[facei]*A;
            }
        }

        reduce(sumA, sumOp<scalar>());
        reduce(sumRhoLA, sumOp<scalar>());
        reduce(sumRhoGA, sumOp<scalar>());

        if (sumA > SMALL)
        {
            rhoLiquidPatchMean_ = sumRhoLA/sumA;
            rhoGasPatchMean_ = sumRhoGA/sumA;
            densitySource_ = "phaseRhoFields";
            return;
        }
    }

    if (!mesh().foundObject<volScalarField>(rhoName_))
    {
        FatalErrorInFunction
            << "Cannot determine phase densities for " << name() << nl
            << "Either set rhoLiquid and rhoGas, provide rhoLiquidField/rhoGasField, "
            << "or make sure the mixture density field '" << rhoName_
            << "' exists." << nl
            << "Available volScalarField objects are "
            << mesh().toc<volScalarField>()
            << exit(FatalError);
    }

    const volScalarField& rho = mesh().lookupObject<volScalarField>(rhoName_);
    const scalarField& V = mesh().V();

    scalar sumVL = 0;
    scalar sumRhoL = 0;
    scalar sumVG = 0;
    scalar sumRhoG = 0;

    scalar S = 0;
    scalar Sa = 0;
    scalar Sr = 0;
    scalar Saa = 0;
    scalar Sar = 0;

    forAll(alpha1, celli)
    {
        const scalar a = limitedAlpha(alpha1[celli]);
        const scalar r = rho[celli];
        const scalar w = V[celli];

        if (a >= densityLiquidAlphaCutoff_)
        {
            sumVL += w;
            sumRhoL += r*w;
        }

        if (a <= densityGasAlphaCutoff_)
        {
            sumVG += w;
            sumRhoG += r*w;
        }

        S += w;
        Sa += a*w;
        Sr += r*w;
        Saa += a*a*w;
        Sar += a*r*w;
    }

    reduce(sumVL, sumOp<scalar>());
    reduce(sumRhoL, sumOp<scalar>());
    reduce(sumVG, sumOp<scalar>());
    reduce(sumRhoG, sumOp<scalar>());
    reduce(S, sumOp<scalar>());
    reduce(Sa, sumOp<scalar>());
    reduce(Sr, sumOp<scalar>());
    reduce(Saa, sumOp<scalar>());
    reduce(Sar, sumOp<scalar>());

    if (sumVL > SMALL && sumVG > SMALL)
    {
        rhoLiquidPatchMean_ = sumRhoL/sumVL;
        rhoGasPatchMean_ = sumRhoG/sumVG;
        densitySource_ = "mixtureRhoPureCells";
        return;
    }

    const scalar denom = S*Saa - Sa*Sa;

    if (S > SMALL && mag(denom) > SMALL)
    {
        const scalar slope = (S*Sar - Sa*Sr)/denom;
        const scalar intercept = (Sr - slope*Sa)/S;

        rhoGasPatchMean_ = intercept;
        rhoLiquidPatchMean_ = intercept + slope;
        densitySource_ = "mixtureRhoLinearFit";

        if (rhoLiquidPatchMean_ > 0 && rhoGasPatchMean_ > 0)
        {
            return;
        }
    }

    FatalErrorInFunction
        << "Could not infer valid phase densities from mixture field "
        << rhoName_ << ". Set rhoLiquid and rhoGas explicitly."
        << nl << "rhoLiquid = " << rhoLiquidPatchMean_
        << ", rhoGas = " << rhoGasPatchMean_
        << exit(FatalError);
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::updateTarget
(
    const volVectorField& U,
    const volScalarField& alpha1
) const
{
    updateDensities(alpha1);

    updateTotalRefluxLoading(U, alpha1);

    targetUbar_ = gasDirection_*gasLoadingTarget_;

    const scalar targetMag = mag(targetUbar_);

    if (maxGasLoading_ > 0 && targetMag > maxGasLoading_)
    {
        WarningInFunction
            << "Capping total-reflux target gas loading from "
            << targetMag << " to " << maxGasLoading_ << nl
            << "Increase maxGasLoading or set it negative to disable this cap."
            << endl;

        targetUbar_ *= maxGasLoading_/targetMag;
    }
}


Foam::scalar Foam::fv::totalRefluxPatchGasMeanVelocityForce::patchGasUbarAve
(
    const volVectorField& U,
    const volScalarField& alpha1
) const
{
    scalar sumA = 0;
    scalar sumQGUp = 0;

    forAll(patches_, patchNamei)
    {
        const word& patchName = patches_[patchNamei];
        const label patchi = mesh().boundaryMesh().findIndex(patchName);

        const scalarField& magSf = mesh().boundary()[patchi].magSf();
        const vectorField& Sfp = mesh().Sf().boundaryField()[patchi];
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvPatchScalarField& alpha1p = alpha1.boundaryField()[patchi];

        forAll(magSf, facei)
        {
            const scalar A = magSf[facei];
            const scalar aG = scalar(1) - limitedAlpha(alpha1p[facei]);
            const scalar phi = Up[facei] & Sfp[facei];
            const scalar qG = aG*phi;

            sumA += A;
            sumQGUp += max(-qG, scalar(0));
        }
    }

    reduce(sumA, sumOp<scalar>());
    reduce(sumQGUp, sumOp<scalar>());

    if (sumA <= SMALL)
    {
        WarningInFunction
            << "Selected patch area is too small for " << name()
            << ". Returning zero patch-averaged gas velocity." << endl;

        return 0;
    }

    return sumQGUp/sumA;
}


Foam::scalar Foam::fv::totalRefluxPatchGasMeanVelocityForce::patchResponseAve
(
    const volScalarField& rA,
    const volScalarField& alpha1
) const
{
    scalar sumA = 0;
    scalar sumResponse = 0;

    forAll(patches_, patchNamei)
    {
        const word& patchName = patches_[patchNamei];
        const label patchi = mesh().boundaryMesh().findIndex(patchName);

        const scalarField& magSf = mesh().boundary()[patchi].magSf();
        const fvPatchScalarField& alpha1p = alpha1.boundaryField()[patchi];
        const fvPatchScalarField& rAp = rA.boundaryField()[patchi];

        forAll(magSf, facei)
        {
            const scalar A = magSf[facei];

            sumA += A;
            sumResponse += sourceWeight(alpha1p[facei])*rAp[facei]*A;
        }
    }

    reduce(sumA, sumOp<scalar>());
    reduce(sumResponse, sumOp<scalar>());

    if (sumA <= SMALL)
    {
        return 0;
    }

    return sumResponse/sumA;
}


Foam::scalar Foam::fv::totalRefluxPatchGasMeanVelocityForce::zoneResponseAve
(
    const volScalarField& rA,
    const volScalarField& alpha1
) const
{
    const labelList& cells = zone_.zone();
    const scalarField& cv = mesh().V();
    const scalarField& rAI = rA;

    scalar sumW = 0;
    scalar sumWR = 0;

    forAll(cells, i)
    {
        const label celli = cells[i];
        const scalar w = sourceWeight(alpha1, celli);
        const scalar volW = w*cv[celli];

        sumW += volW;
        sumWR += volW*rAI[celli];
    }

    reduce(sumW, sumOp<scalar>());
    reduce(sumWR, sumOp<scalar>());

    if (sumW <= SMALL)
    {
        return 0;
    }

    return sumWR/sumW;
}


Foam::fv::totalRefluxPatchGasMeanVelocityForce::
totalRefluxPatchGasMeanVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(sourceName, modelType, mesh, dict),
    zone_(mesh, coeffs(dict)),
    patches_(),
    UName_(word::null),
    alphaName_(word::null),
    alphaPhiName_(word::null),
    liquidFluxMode_("directionalAlphaU"),
    phiName_("phi"),
    rhoName_("rho"),
    rhoLiquidFieldName_("rho.liquid"),
    rhoGasFieldName_("rho.gas"),
    molarTotalReflux_(true),
    species1LiquidName_("cyclohexane.liquid"),
    species2LiquidName_("nHeptane.liquid"),
    species1GasName_("cyclohexane.gas"),
    species2GasName_("nHeptane.gas"),
    W1_(84.16/1000.0),
    W2_(100.2/1000.0),
    rhoLiquid_(-1),
    rhoGas_(-1),
    densityLiquidAlphaCutoff_(0.99),
    densityGasAlphaCutoff_(0.01),
    gasDirection_(0, 1, 0),
    averagingInterval_(0.01),
    relaxation_(0.05),
    alphaCutoff_(0.2),
    useCutoff_(true),
    alphaRamp_(0.1),
    counterCurrent_(true),
    maxGasLoading_(-1),
    maxPressureGradient_(-1),
    maxDeltaPressureGradient_(-1),
    lastIntegratedTimeIndex_(-1),
    accumulatedEquivalentGasVolume_(0),
    accumulatedTime_(0),
    liquidFlowRate_(0),
    patchArea_(0),
    instantLiquidLoading_(0),
    liquidLoading_(0),
    liquidLoadingInitialized_(false),
    totalMassFlowRateDown_(0),
    liquidVolumeFlowRateDown_(0),
    gasVolumeFlowRateDown_(0),
    instantGasLoadingTarget_(0),
    gasLoadingTarget_(0),
    balanceControl_(true),
    balanceWindow_(0.01),
    balanceTolerance_(1e-8),
    balanceRelativeTolerance_(1e-4),
    balanceIntegralGain_(0.25),
    maxBalanceCorrection_(-1),
    maxBalanceCorrectionFraction_(-1),
    gasMassFlowRateUp_(0),
    gasVolumeFlowRateUp_(0),
    totalMolarFlowRateDown_(0),
    gasMolarFlowRateUp_(0),
    gasMolarConcentrationPatchMean_(0),
    balanceErrorRate_(0),
    balanceWindowErrorMass_(0),
    balanceWindowTime_(0),
    balanceCorrectionLoading_(0),
    lastBalanceUpdateTimeIndex_(-1),
    targetUbar_(Zero),
    rhoLiquidPatchMean_(0),
    rhoGasPatchMean_(0),
    densitySource_(word::null),
    resetGradient_(false),
    gradP0_(0),
    dGradP_(0),
    rAPtr_(nullptr)
{
    readCoeffs(coeffs(dict));

    IFstream propsFile
    (
        mesh.time().timePath()/"uniform"/(this->name() + "Properties")
    );

    if (propsFile.good() && !resetGradient_)
    {
        Info<< "Reading pressure gradient and loading state from file" << endl;

        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("gradient") >> gradP0_;

        if (propsDict.readIfPresent("gasLoadingTarget", gasLoadingTarget_))
        {
            liquidLoading_ = gasLoadingTarget_;
            liquidLoadingInitialized_ = true;
        }
        else if (propsDict.readIfPresent("liquidLoading", liquidLoading_))
        {
            gasLoadingTarget_ = liquidLoading_;
            liquidLoadingInitialized_ = true;
        }

        propsDict.readIfPresent
        (
            "accumulatedEquivalentGasVolume",
            accumulatedEquivalentGasVolume_
        );
        propsDict.readIfPresent("accumulatedTime", accumulatedTime_);
    }

    Info<< "Initial pressure gradient = " << gradP0_ << nl
        << "Total reflux source is field-based and solver-independent." << nl
        << "Target definition: "
        << (molarTotalReflux_ ? word("NGasUp,target = NDown") : word("rhoG*QGup,target = rhoL*QLdown + rhoG*QGdown")) << nl
        << "Window balance correction: "
        << (balanceControl_ ? word("on") : word("off"))
        << ", window = " << balanceWindow_
        << " s, tolerance = " << balanceTolerance_
        << (molarTotalReflux_ ? " mol/s" : " kg/s")
        << ", relative tolerance = " << balanceRelativeTolerance_
        << ", gain = " << balanceIntegralGain_ << nl
        << "Liquid alpha field = "
        << (alphaName_ == word::null ? word("auto") : alphaName_)
        << ", liquid alphaPhi field = "
        << (alphaPhiName_ == word::null ? word("auto") : alphaPhiName_)
        << ", fallback phi = " << phiName_ << nl
        << "Flux measurement mode = " << liquidFluxMode_
        << ", counterCurrent = " << counterCurrent_ << nl
        << "Target gas loading is low-pass filtered with time scale "
        << averagingInterval_ << " s" << nl
        << "Gas target direction = " << gasDirection_
        << endl;
}


Foam::wordList Foam::fv::totalRefluxPatchGasMeanVelocityForce::
constrainedFields() const
{
    return wordList(1, UName_);
}


bool Foam::fv::totalRefluxPatchGasMeanVelocityForce::constrain
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (fieldName != UName_)
    {
        return false;
    }

    const word alphaFieldName(actualAlphaName());

    if (!mesh().foundObject<volScalarField>(alphaFieldName))
    {
        FatalErrorInFunction
            << "Cannot find liquid alpha field " << alphaFieldName << nl
            << "Available volScalarField objects are "
            << mesh().toc<volScalarField>()
            << exit(FatalError);
    }

    const volScalarField& alpha1 =
        mesh().lookupObject<volScalarField>(alphaFieldName);

    const volVectorField& U =
        mesh().lookupObject<volVectorField>(UName_);

    updateTarget(U, alpha1);

    if (mag(targetUbar_) <= SMALL)
    {
        dGradP_ = -gradP0_;
    }

    const labelList& cells = zone_.zone();
    const scalar gradP = gradP0_ + dGradP_;

    volVectorField::Internal Su
    (
        IOobject
        (
            name() + fieldName + "Sup",
            mesh().time().name(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    forAll(cells, i)
    {
        const label celli = cells[i];
        Su[celli] = sourceWeight(alpha1, celli)*gasDirection_*gradP;
    }

    eqn -= Su;

    if (rAPtr_.empty())
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name() + ":rA",
                    mesh().time().name(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                1/eqn.A()
            )
        );
    }
    else
    {
        rAPtr_() = 1/eqn.A();
    }

    gradP0_ += dGradP_;
    dGradP_ = 0;

    return true;
}


bool Foam::fv::totalRefluxPatchGasMeanVelocityForce::constrain
(
    volVectorField& U
) const
{
    if (U.name() != UName_)
    {
        return false;
    }

    if (rAPtr_.empty())
    {
        FatalErrorInFunction
            << "rA field is not available. The matrix constrain stage must be "
            << "called before the field correction stage." << nl
            << exit(FatalError);
    }

    const word alphaFieldName(actualAlphaName());
    const volScalarField& alpha1 =
        mesh().lookupObject<volScalarField>(alphaFieldName);

    const volScalarField& rA = rAPtr_();
    const scalarField& rAU = rA;
    const labelList& cells = zone_.zone();

    updateTarget(U, alpha1);

    scalar jGbarAve = this->patchGasUbarAve(U, alpha1);
    const scalar jGtarget = gasDirection_ & targetUbar_;

    if (mag(targetUbar_) <= SMALL)
    {
        Info<< "Total-reflux gas pressure gradient source:"
            << " mDown = " << totalMassFlowRateDown_
            << ", QLdown = " << liquidVolumeFlowRateDown_
            << ", QGdown = " << gasVolumeFlowRateDown_
            << ", QGup = " << gasVolumeFlowRateUp_
            << ", balanceErr = " << balanceErrorRate_
            << ", target jG = " << jGtarget
            << ", actual jG = " << jGbarAve
            << ", pressure gradient = " << gradP0_ + dGradP_
            << endl;

        writeProps(0);
        return true;
    }

    scalar rAUeff = patchResponseAve(rA, alpha1);

    if (rAUeff <= SMALL)
    {
        rAUeff = zoneResponseAve(rA, alpha1);
    }

    if (rAUeff <= SMALL)
    {
        WarningInFunction
            << "Effective gas response is too small. No forcing update "
            << "applied for " << name() << endl;

        return true;
    }

    dGradP_ = relaxation_*(jGtarget - jGbarAve)/max(rAUeff, SMALL);

    if (maxDeltaPressureGradient_ > 0 && mag(dGradP_) > maxDeltaPressureGradient_)
    {
        dGradP_ =
            dGradP_ < 0 ? -maxDeltaPressureGradient_ : maxDeltaPressureGradient_;
    }

    if (maxPressureGradient_ > 0 && mag(gradP0_ + dGradP_) > maxPressureGradient_)
    {
        const scalar limitedGradP =
            (gradP0_ + dGradP_) < 0
          ? -maxPressureGradient_
          :  maxPressureGradient_;

        dGradP_ = limitedGradP - gradP0_;
    }

    forAll(cells, i)
    {
        const label celli = cells[i];
        const scalar w = sourceWeight(alpha1, celli);

        U[celli] += w*gasDirection_*rAU[celli]*dGradP_;
    }

    U.correctBoundaryConditions();

    // Update the integral balance controller with the corrected velocity.
    // Its correction affects the target in the following calls/time windows.
    updateBalanceController(U, alpha1);

    jGbarAve = this->patchGasUbarAve(U, alpha1);

    const scalar gradP = gradP0_ + dGradP_;

    Info<< "Total-reflux gas pressure gradient source:"
        << " mDown = " << totalMassFlowRateDown_
        << ", QLdown = " << liquidVolumeFlowRateDown_
        << ", QGdown = " << gasVolumeFlowRateDown_
        << ", QGup = " << gasVolumeFlowRateUp_
        << ", mGasUp = " << gasMassFlowRateUp_
        << ", balanceErr = " << balanceErrorRate_
        << ", balanceCorrJ = " << balanceCorrectionLoading_
        << ", filtered target jG = " << gasLoadingTarget_
        << ", filter time = " << averagingInterval_
        << ", jG target/actual= " << jGtarget << "/" << jGbarAve
        << ", pressure gradient = " << gradP
        << endl;

    writeProps(gradP);

    return true;
}


bool Foam::fv::totalRefluxPatchGasMeanVelocityForce::movePoints()
{
    zone_.movePoints();
    return true;
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::topoChange
(
    const polyTopoChangeMap& map
)
{
    zone_.topoChange(map);
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::mapMesh
(
    const polyMeshMap& map
)
{
    zone_.mapMesh(map);
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::distribute
(
    const polyDistributionMap& map
)
{
    zone_.distribute(map);
}


bool Foam::fv::totalRefluxPatchGasMeanVelocityForce::read
(
    const dictionary& dict
)
{
    if (fvConstraint::read(dict))
    {
        zone_.read(coeffs(dict));
        readCoeffs(coeffs(dict));
        return true;
    }

    return false;
}

