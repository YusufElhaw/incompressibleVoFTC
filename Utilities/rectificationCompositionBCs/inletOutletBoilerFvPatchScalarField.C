/*---------------------------------------------------------------------------*\
  Dynamic boiler inlet/outlet composition boundary condition
\*---------------------------------------------------------------------------*/

#include "inletOutletBoilerFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "typeInfo.H"

namespace Foam
{

scalar inletOutletBoilerFvPatchScalarField::clamp(const scalar x) const
{
    return min(max(x, minValue_), maxValue_);
}


wordList inletOutletBoilerFvPatchScalarField::phaseNames() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    IOdictionary phaseProperties
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    wordList phases;
    phaseProperties.lookup("phases") >> phases;

    if (phases.size() != 2)
    {
        FatalErrorInFunction
            << "Expected exactly two phases in constant/phaseProperties, "
            << "but found " << phases << nl
            << "This boundary condition assumes one liquid phase and one gas phase."
            << exit(FatalError);
    }

    return phases;
}


word inletOutletBoilerFvPatchScalarField::alphaFieldName
(
    const word& phaseName
) const
{
    return IOobject::groupName("alpha", phaseName);
}


word inletOutletBoilerFvPatchScalarField::alphaPhiName
(
    const word& phaseName
) const
{
    return IOobject::groupName("alphaPhi", phaseName);
}


word inletOutletBoilerFvPatchScalarField::liquidPhaseName() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const wordList phases = phaseNames();

    wordList candidates(phases.size());
    label n = 0;

    // Primary criterion: the alpha.<phase> volScalarField registered by the
    // solver. In incompressibleVoFTC there is only one such field, and it is
    // the liquid phase fraction.
    forAll(phases, i)
    {
        const word aName = alphaFieldName(phases[i]);

        if (mesh.foundObject<volScalarField>(aName))
        {
            candidates[n++] = phases[i];
        }
    }

    // Fallback for early construction/restart edge cases: check whether an
    // alpha.<phase> file exists in the current time directory.
    if (n == 0)
    {
        forAll(phases, i)
        {
            IOobject aHeader
            (
                alphaFieldName(phases[i]),
                mesh.time().name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );

            if (aHeader.headerOk())
            {
                candidates[n++] = phases[i];
            }
        }
    }

    if (n == 1)
    {
        return candidates[0];
    }

    candidates.setSize(n);

    FatalErrorInFunction
        << "Could not determine the liquid phase from the alpha.<phase> field." << nl
        << "The phase order in constant/phaseProperties is not used to choose "
        << "the liquid phase." << nl
        << "Exactly one alpha.<phase> field must exist and it is treated as "
        << "the liquid phase." << nl
        << "phases = " << phases << nl
        << "matching alpha fields = " << candidates << nl
        << exit(FatalError);

    return word();
}


word inletOutletBoilerFvPatchScalarField::gasPhaseName
(
    const word& liquidPhaseName
) const
{
    const wordList phases = phaseNames();

    forAll(phases, i)
    {
        if (phases[i] != liquidPhaseName)
        {
            return phases[i];
        }
    }

    FatalErrorInFunction
        << "Could not determine the gas phase. Liquid phase is "
        << liquidPhaseName << ", phases are " << phases
        << exit(FatalError);

    return word();
}


void inletOutletBoilerFvPatchScalarField::setAutomaticFluxNames()
{
    // incompressibleVoFTC convention for this BC library:
    // The single alpha.<phase> field present in the case/registry identifies
    // the liquid phase. The other phase in phaseProperties is treated as gas.
    const word liquid = liquidPhaseName();
    const word gas = gasPhaseName(liquid);

    // Boiler/reboiler:
    // - measure outgoing liquid composition
    // - impose incoming gas composition
    samplePhiName_ = alphaPhiName(liquid);
    phiName_       = alphaPhiName(gas);
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
    else if (dict.found("otherBoilerPatches"))
    {
        wordList others;
        dict.lookup("otherBoilerPatches") >> others;

        patchNames_.setSize(others.size() + 1);
        patchNames_[0] = patch().name();
        forAll(others, i)
        {
            patchNames_[i + 1] = others[i];
        }
        usePatchList_ = true;
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
            const label n = patchNames_.size();
            patchNames_.setSize(n + 1);
            patchNames_[n] = patch().name();
        }
    }
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
    const volScalarField& x =
        mesh.lookupObject<volScalarField>(internalField().name());

    const label nPatches = x.boundaryField().size();
    wordList names(nPatches);
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


label inletOutletBoilerFvPatchScalarField::patchID(const word& name) const
{
    const fvBoundaryMesh& bm = patch().boundaryMesh();

    forAll(bm, patchi)
    {
        if (bm[patchi].name() == name)
        {
            return patchi;
        }
    }

    FatalErrorInFunction
        << "Patch " << name << " was requested by boundary condition "
        << type() << " on patch " << patch().name()
        << ", but it does not exist in the mesh boundary."
        << exit(FatalError);

    return -1;
}


scalar inletOutletBoilerFvPatchScalarField::outgoingLiquidMean
(
    const wordList& patchNames
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const volScalarField& x =
        mesh.lookupObject<volScalarField>(internalField().name());

    const surfaceScalarField& samplePhi =
        mesh.lookupObject<surfaceScalarField>(samplePhiName_);

    scalar num = 0;
    scalar den = 0;

    forAll(patchNames, i)
    {
        const label patchi = patchID(patchNames[i]);

        const fvsPatchField<scalar>& samplePhip =
            samplePhi.boundaryField()[patchi];

        const scalarField xPatchInternal
        (
            x.boundaryField()[patchi].patchInternalField()
        );

        forAll(samplePhip, facei)
        {
            const scalar w = max(samplePhip[facei], scalar(0));
            num += w*xPatchInternal[facei];
            den += w;
        }
    }

    reduce(num, sumOp<scalar>());
    reduce(den, sumOp<scalar>());

    if (den > minFlux_)
    {
        return clamp(num/den);
    }

    return measuredLiquidValue_;
}


scalar inletOutletBoilerFvPatchScalarField::outgoingA12Mean
(
    const wordList& patchNames
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (!mesh.foundObject<volScalarField>(A12Name_))
    {
        FatalErrorInFunction
            << "No volScalarField named " << A12Name_
            << " was found for " << type()
            << " on patch " << patch().name() << nl
            << "This version of inletOutletBoiler always uses the field A12 "
            << "generated by compositionPredictor. Dictionary entries "
            << "A12Field, useA12Field and constant A12 are not supported."
            << exit(FatalError);
    }

    const volScalarField& A12Field =
        mesh.lookupObject<volScalarField>(A12Name_);

    const surfaceScalarField& samplePhi =
        mesh.lookupObject<surfaceScalarField>(samplePhiName_);

    scalar num = 0;
    scalar den = 0;

    forAll(patchNames, i)
    {
        const label patchi = patchID(patchNames[i]);

        const fvsPatchField<scalar>& samplePhip =
            samplePhi.boundaryField()[patchi];

        const scalarField A12PatchInternal
        (
            A12Field.boundaryField()[patchi].patchInternalField()
        );

        forAll(samplePhip, facei)
        {
            const scalar w = max(samplePhip[facei], scalar(0));
            num += w*A12PatchInternal[facei];
            den += w;
        }
    }

    reduce(num, sumOp<scalar>());
    reduce(den, sumOp<scalar>());

    if (den > minFlux_)
    {
        return max(num/den, SMALL);
    }

    return max(measuredA12_, SMALL);
}


scalar inletOutletBoilerFvPatchScalarField::oldInletMean
(
    const wordList& patchNames
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const volScalarField& x =
        mesh.lookupObject<volScalarField>(internalField().name());

    scalar sumValue = 0;
    scalar nValue = 0;

    forAll(patchNames, i)
    {
        const label patchi = patchID(patchNames[i]);
        const fvPatchScalarField& pf = x.boundaryField()[patchi];

        if (isA<inletOutletBoilerFvPatchScalarField>(pf))
        {
            const inletOutletBoilerFvPatchScalarField& bpf =
                refCast<const inletOutletBoilerFvPatchScalarField>(pf);

            sumValue += bpf.inletValue_;
            nValue += 1;
        }
    }

    reduce(sumValue, sumOp<scalar>());
    reduce(nValue, sumOp<scalar>());

    if (nValue > SMALL)
    {
        return clamp(sumValue/nValue);
    }

    return inletValue_;
}


scalar inletOutletBoilerFvPatchScalarField::vapourEquilibrium
(
    const scalar x,
    const scalar A12
) const
{
    const scalar xb = clamp(x);
    const scalar Ab = max(A12, SMALL);

    return clamp(Ab*xb/(1 + (Ab - 1)*xb));
}


inletOutletBoilerFvPatchScalarField::inletOutletBoilerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_(word()),
    samplePhiName_(word()),
    groupName_(word()),
    patchNames_(),
    usePatchList_(false),
    A12Name_("A12"),
    relax_(1),
    inletValue_(0),
    measuredLiquidValue_(0),
    measuredA12_(1),
    minValue_(0),
    maxValue_(1),
    minFlux_(VSMALL)
{
    setAutomaticFluxNames();

    refValue() = inletValue_;
    refGrad() = 0;
    valueFraction() = 0;
    fvPatchScalarField::operator=(inletValue_);
}


inletOutletBoilerFvPatchScalarField::inletOutletBoilerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    phiName_(word()),
    samplePhiName_(word()),
    groupName_(dict.lookupOrDefault<word>("group", word())),
    patchNames_(),
    usePatchList_(false),
    A12Name_("A12"),
    relax_(dict.lookupOrDefault<scalar>("relax", 1)),
    inletValue_(dict.lookupOrDefault<scalar>("inletValue", 0)),
    measuredLiquidValue_
    (
        dict.lookupOrDefault<scalar>("measuredLiquidValue", inletValue_)
    ),
    measuredA12_(dict.lookupOrDefault<scalar>("measuredA12", 1)),
    minValue_(dict.lookupOrDefault<scalar>("minValue", 0)),
    maxValue_(dict.lookupOrDefault<scalar>("maxValue", 1)),
    minFlux_(dict.lookupOrDefault<scalar>("minFlux", VSMALL))
{
    setAutomaticFluxNames();
    readPatchSelection(dict);

    inletValue_ = clamp(inletValue_);
    measuredLiquidValue_ = clamp(measuredLiquidValue_);
    measuredA12_ = max(measuredA12_, SMALL);
    relax_ = min(max(relax_, scalar(0)), scalar(1));

    refValue() = inletValue_;
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
        fvPatchScalarField::operator=(inletValue_);
    }
}


inletOutletBoilerFvPatchScalarField::inletOutletBoilerFvPatchScalarField
(
    const inletOutletBoilerFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    samplePhiName_(ptf.samplePhiName_),
    groupName_(ptf.groupName_),
    patchNames_(ptf.patchNames_),
    usePatchList_(ptf.usePatchList_),
    A12Name_(ptf.A12Name_),
    relax_(ptf.relax_),
    inletValue_(ptf.inletValue_),
    measuredLiquidValue_(ptf.measuredLiquidValue_),
    measuredA12_(ptf.measuredA12_),
    minValue_(ptf.minValue_),
    maxValue_(ptf.maxValue_),
    minFlux_(ptf.minFlux_)
{}


inletOutletBoilerFvPatchScalarField::inletOutletBoilerFvPatchScalarField
(
    const inletOutletBoilerFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    samplePhiName_(ptf.samplePhiName_),
    groupName_(ptf.groupName_),
    patchNames_(ptf.patchNames_),
    usePatchList_(ptf.usePatchList_),
    A12Name_(ptf.A12Name_),
    relax_(ptf.relax_),
    inletValue_(ptf.inletValue_),
    measuredLiquidValue_(ptf.measuredLiquidValue_),
    measuredA12_(ptf.measuredA12_),
    minValue_(ptf.minValue_),
    maxValue_(ptf.maxValue_),
    minFlux_(ptf.minFlux_)
{}


void inletOutletBoilerFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const wordList patchNames = selectedPatchNames();

    measuredLiquidValue_ = outgoingLiquidMean(patchNames);
    measuredA12_ = outgoingA12Mean(patchNames);

    const scalar target = vapourEquilibrium(measuredLiquidValue_, measuredA12_);
    const scalar oldValue = oldInletMean(patchNames);

    inletValue_ = clamp((1 - relax_)*oldValue + relax_*target);

    refValue() = inletValue_;
    refGrad() = 0;
    valueFraction() = neg(phip);

    mixedFvPatchScalarField::updateCoeffs();
}


void inletOutletBoilerFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    if (!groupName_.empty())
    {
        writeEntry(os, "group", groupName_);
    }

    if (usePatchList_)
    {
        writeEntry(os, "patches", patchNames_);
    }

    writeEntryIfDifferent<scalar>(os, "relax", 1, relax_);
    writeEntry(os, "inletValue", inletValue_);
    writeEntry(os, "measuredLiquidValue", measuredLiquidValue_);
    writeEntry(os, "measuredA12", measuredA12_);
    writeEntryIfDifferent<scalar>(os, "minValue", 0, minValue_);
    writeEntryIfDifferent<scalar>(os, "maxValue", 1, maxValue_);
    writeEntryIfDifferent<scalar>(os, "minFlux", VSMALL, minFlux_);
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
      + (1 - valueFraction())*ptf
    );
}


makePatchTypeField
(
    fvPatchScalarField,
    inletOutletBoilerFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
