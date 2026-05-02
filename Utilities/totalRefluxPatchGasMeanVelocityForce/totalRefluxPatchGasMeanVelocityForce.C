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
        dict.lookupOrDefault<word>("rhoLiquidField", word::null);
    rhoGasFieldName_ =
        dict.lookupOrDefault<word>("rhoGasField", word::null);

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


Foam::scalar Foam::fv::totalRefluxPatchGasMeanVelocityForce::gasFluxWeight
(
    const scalar alpha1
) const
{
    // This is used for measuring the superficial gas velocity.  It must not
    // include the alphaCutoff_, because low-gas-fraction faces still belong to
    // the physical flux measurement.
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
        propsDict.add("accumulatedLiquidVolume", accumulatedLiquidVolume_);
        propsDict.add("accumulatedTime", accumulatedTime_);
        propsDict.add("rhoLiquid", rhoLiquidPatchMean_);
        propsDict.add("rhoGas", rhoGasPatchMean_);
        propsDict.add("densitySource", densitySource_);
        propsDict.regIOobject::write();
    }
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::measureLiquidFlowRate
(
    const volVectorField& U,
    const volScalarField& alpha1,
    scalar& liquidFlowRate,
    scalar& patchArea
) const
{
    const word alphaPhiFieldName(actualAlphaPhiName(alpha1));

    const surfaceScalarField* alphaPhi1Ptr = nullptr;
    const surfaceScalarField* phiPtr = nullptr;

    if (liquidFluxMode_ == "alphaPhi")
    {
        if (alphaPhiFieldName != word::null)
        {
            alphaPhi1Ptr =
                &mesh().lookupObject<surfaceScalarField>(alphaPhiFieldName);
        }
        else
        {
            if (!mesh().foundObject<surfaceScalarField>(phiName_))
            {
                FatalErrorInFunction
                    << "Cannot find liquid phase flux field alphaPhi for "
                    << alpha1.name() << " and cannot find fallback flux field "
                    << phiName_ << nl
                    << "Set 'alphaPhi' explicitly or switch to "
                    << "liquidFluxMode directionalAlphaU."
                    << nl << "Available surfaceScalarField objects are "
                    << mesh().toc<surfaceScalarField>()
                    << exit(FatalError);
            }

            phiPtr = &mesh().lookupObject<surfaceScalarField>(phiName_);
        }
    }

    scalar sumQL = 0;
    scalar sumA = 0;

    forAll(patches_, patchNamei)
    {
        const word& patchName = patches_[patchNamei];
        const label patchi = mesh().boundaryMesh().findIndex(patchName);

        const scalarField& magSf = mesh().boundary()[patchi].magSf();
        const fvPatchScalarField& alpha1p = alpha1.boundaryField()[patchi];
        const fvPatchVectorField& Up = U.boundaryField()[patchi];

        const fvsPatchScalarField* alphaPhi1pPtr = nullptr;
        const fvsPatchScalarField* phipPtr = nullptr;

        if (liquidFluxMode_ == "alphaPhi")
        {
            if (alphaPhi1Ptr)
            {
                alphaPhi1pPtr = &alphaPhi1Ptr->boundaryField()[patchi];
            }
            else
            {
                phipPtr = &phiPtr->boundaryField()[patchi];
            }
        }

        forAll(magSf, facei)
        {
            scalar phiL = 0;

            if (liquidFluxMode_ == "alphaPhi")
            {
                if (alphaPhi1pPtr)
                {
                    phiL = (*alphaPhi1pPtr)[facei];
                }
                else
                {
                    phiL = limitedAlpha(alpha1p[facei])*(*phipPtr)[facei];
                }
            }
            else
            {
                // Robust default: physical directional liquid flow through the
                // selected cross-section.  This avoids NCC/coupled-patch
                // alphaPhi spikes and is independent of patch-normal direction.
                phiL = limitedAlpha(alpha1p[facei])
                    *(Up[facei] & gasDirection_)*magSf[facei];
            }

            sumQL += phiL;
            sumA += magSf[facei];
        }
    }

    reduce(sumQL, sumOp<scalar>());
    reduce(sumA, sumOp<scalar>());

    liquidFlowRate = sumQL;
    patchArea = sumA;
}


void Foam::fv::totalRefluxPatchGasMeanVelocityForce::updateLiquidLoading
(
    const volVectorField& U,
    const volScalarField& alpha1
) const
{
    scalar QL = 0;
    scalar A0 = 0;
    measureLiquidFlowRate(U, alpha1, QL, A0);

    liquidFlowRate_ = QL;
    patchArea_ = A0;

    if (A0 <= SMALL)
    {
        WarningInFunction
            << "Selected patch area is too small for " << name()
            << ". Setting liquid loading and total-reflux gas target to zero."
            << endl;

        instantLiquidLoading_ = 0;
        liquidLoading_ = 0;
        liquidLoadingInitialized_ = true;
        return;
    }

    // directionalAlphaU gives QL signed along gasDirection.  Total reflux
    // needs the liquid load opposing gasDirection, hence -QL/A0.
    // alphaPhi mode keeps the old signed patch-normal definition unless
    // counterCurrent is enabled explicitly.
    scalar jLinst = QL/A0;

    if (counterCurrent_)
    {
        jLinst = max(scalar(0), -jLinst);
    }

    instantLiquidLoading_ = jLinst;

    const label timeIndex = mesh().time().timeIndex();

    if (!liquidLoadingInitialized_ || averagingInterval_ <= SMALL)
    {
        liquidLoading_ = jLinst;
        liquidLoadingInitialized_ = true;
        lastIntegratedTimeIndex_ = timeIndex;
        accumulatedLiquidVolume_ = 0;
        accumulatedTime_ = 0;
        return;
    }

    if (timeIndex != lastIntegratedTimeIndex_)
    {
        const scalar dt = mesh().time().deltaTValue();

        if (dt > SMALL)
        {
            // Exponential physical-time filter.  It is smoother than the old
            // rectangular reset window and avoids step changes every interval.
            const scalar beta = min(scalar(1), dt/max(averagingInterval_, SMALL));
            liquidLoading_ = (scalar(1) - beta)*liquidLoading_ + beta*jLinst;

            accumulatedLiquidVolume_ += QL*dt;
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
        const volScalarField& rhoL = mesh().lookupObject<volScalarField>(rhoLName);
        const volScalarField& rhoG = mesh().lookupObject<volScalarField>(rhoGName);

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
    updateLiquidLoading(U, alpha1);
    updateDensities(alpha1);

    if (rhoGasPatchMean_ <= SMALL)
    {
        FatalErrorInFunction
            << "Gas density is too small for " << name()
            << ". Add rhoGas or check the mixture rho field."
            << exit(FatalError);
    }

    const scalar uSoll = (rhoLiquidPatchMean_/rhoGasPatchMean_)*liquidLoading_;

    targetUbar_ = gasDirection_*uSoll;

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
    scalar sumAjG = 0;

    forAll(patches_, patchNamei)
    {
        const word& patchName = patches_[patchNamei];
        const label patchi = mesh().boundaryMesh().findIndex(patchName);

        const scalarField& magSf = mesh().boundary()[patchi].magSf();
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvPatchScalarField& alpha1p = alpha1.boundaryField()[patchi];

        forAll(magSf, facei)
        {
            const scalar A = magSf[facei];
            const scalar w = gasFluxWeight(alpha1p[facei]);

            sumA += A;
            sumAjG += w*(gasDirection_ & Up[facei])*A;
        }
    }

    reduce(sumA, sumOp<scalar>());
    reduce(sumAjG, sumOp<scalar>());

    if (sumA <= SMALL)
    {
        WarningInFunction
            << "Selected patch area is too small for " << name()
            << ". Returning zero patch-averaged gas velocity." << endl;

        return 0;
    }

    return sumAjG/sumA;
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
    rhoLiquidFieldName_(word::null),
    rhoGasFieldName_(word::null),
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
    accumulatedLiquidVolume_(0),
    accumulatedTime_(0),
    liquidFlowRate_(0),
    patchArea_(0),
    instantLiquidLoading_(0),
    liquidLoading_(0),
    liquidLoadingInitialized_(false),
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

        if (propsDict.readIfPresent("liquidLoading", liquidLoading_))
        {
            liquidLoadingInitialized_ = true;
        }
        propsDict.readIfPresent("accumulatedLiquidVolume", accumulatedLiquidVolume_);
        propsDict.readIfPresent("accumulatedTime", accumulatedTime_);
    }

    Info<< "Initial pressure gradient = " << gradP0_ << nl
        << "Total reflux source is field-based and solver-independent." << nl
        << "Liquid alpha field = "
        << (alphaName_ == word::null ? word("auto") : alphaName_)
        << ", liquid alphaPhi field = "
        << (alphaPhiName_ == word::null ? word("auto") : alphaPhiName_)
        << ", fallback phi = " << phiName_ << nl
        << "Liquid loading mode = " << liquidFluxMode_
        << ", counterCurrent = " << counterCurrent_ << nl
        << "Liquid loading is low-pass filtered with time scale "
        << averagingInterval_ << " s" << nl
        << "Gas target direction = " << gasDirection_ << nl
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

    const volScalarField& alpha1 = mesh().lookupObject<volScalarField>(alphaFieldName);

    const volVectorField& U = mesh().lookupObject<volVectorField>(UName_);

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
    const volScalarField& alpha1 = mesh().lookupObject<volScalarField>(alphaFieldName);
    const volScalarField& rA = rAPtr_();
    const scalarField& rAU = rA;
    const labelList& cells = zone_.zone();

    if (mag(targetUbar_) <= SMALL)
    {
        Info<< "Total reflux gas target is zero for " << name()
            << ": alpha = " << alpha1.name()
            << ", liquid QL = " << liquidFlowRate_
            << ", patch area = " << patchArea_
            << ", liquid loading B_L = " << liquidLoading_
            << ", target jG = " << targetUbar_ << endl;

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

    const scalar jGbarAve = this->patchGasUbarAve(U, alpha1);
    const scalar jGtarget = gasDirection_ & targetUbar_;

    dGradP_ = relaxation_*(jGtarget - jGbarAve)/max(rAUeff, SMALL);

    if (maxDeltaPressureGradient_ > 0 && mag(dGradP_) > maxDeltaPressureGradient_)
    {
        dGradP_ =
            dGradP_ < 0 ? -maxDeltaPressureGradient_ : maxDeltaPressureGradient_;
    }

    if (maxPressureGradient_ > 0 && mag(gradP0_ + dGradP_) > maxPressureGradient_)
    {
        const scalar limitedGradP =
            (gradP0_ + dGradP_) < 0 ? -maxPressureGradient_ : maxPressureGradient_;

        dGradP_ = limitedGradP - gradP0_;
    }

    forAll(cells, i)
    {
        const label celli = cells[i];
        const scalar w = sourceWeight(alpha1, celli);

        U[celli] += w*gasDirection_*rAU[celli]*dGradP_;
    }

    const scalar gradP = gradP0_ + dGradP_;

    Info<< "Total-reflux gas pressure gradient source: alpha = "
        << alpha1.name()
        << ", alphaPhi = " << actualAlphaPhiName(alpha1)
        << ", liquid QL = " << liquidFlowRate_
        << ", patch area A0 = " << patchArea_
        << ", instantaneous liquid loading B_L_inst = " << instantLiquidLoading_
        << ", filtered liquid loading B_L = " << liquidLoading_
        << ", filter time = " << averagingInterval_
        << ", response = " << rAUeff
        << ", accumulated time = " << accumulatedTime_
        << ", rhoL = " << rhoLiquidPatchMean_
        << ", rhoG = " << rhoGasPatchMean_
        << ", density source = " << densitySource_
        << ", gas direction = " << gasDirection_
        << ", target jG = " << jGtarget
        << ", uncorrected jG = " << jGbarAve
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
    else
    {
        return false;
    }
}

// ************************************************************************* //
