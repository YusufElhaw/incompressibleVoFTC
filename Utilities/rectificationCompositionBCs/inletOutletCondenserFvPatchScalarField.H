/*---------------------------------------------------------------------------*\
  Dynamic condenser inlet/outlet composition boundary condition
\*---------------------------------------------------------------------------*/

#ifndef inletOutletCondenserFvPatchScalarField_H
#define inletOutletCondenserFvPatchScalarField_H

#include "mixedFvPatchFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class inletOutletCondenserFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class inletOutletCondenserFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    //- Flux used for inlet/outlet switching of this field.
    //  Negative flux means inlet and therefore fixed refValue.
    word phiName_;

    //- Flux used to measure the outgoing vapour composition.
    //  Positive flux is treated as outflow from the CFD domain.
    word samplePhiName_;

    //- Optional group name. If set and no explicit patches list is given,
    //  all inletOutletCondensator patches with the same group are coupled.
    word groupName_;

    //- Explicit list of coupled patches. If empty, use groupName_ or only self.
    wordList patchNames_;

    //- True if patchNames_ was explicitly supplied by the user.
    bool usePatchList_;

    //- Relaxation factor for the dynamic inlet value: 1 = no smoothing
    scalar relax_;

    //- Current dynamic liquid inlet value
    scalar inletValue_;

    //- Last measured outgoing vapour value before relaxation
    scalar measuredValue_;

    //- Bounds for the dynamic inlet value
    scalar minValue_;
    scalar maxValue_;

    //- Small flux below which the previous inlet value is kept
    scalar minFlux_;

    //- Return phase names from constant/phaseProperties.
    wordList phaseNames() const;

    //- Return alpha.<phaseName>.
    word alphaFieldName(const word& phaseName) const;

    //- Return alphaPhi.<phaseName>.
    word alphaPhiName(const word& phaseName) const;

    //- Detect the liquid phase from the single registered/current-time
    //  alpha.<phase> field. The phase order in phaseProperties is not used
    //  to choose the liquid phase.
    word liquidPhaseName() const;

    //- Return the other phase from phaseProperties as gas.
    word gasPhaseName(const word& liquidPhaseName) const;

    //- Set automatic flux names.
    //  Condenser: phi = alphaPhi.liquidPhase, samplePhi = alphaPhi.gasPhase.
    void setAutomaticFluxNames();

    //- Set explicit patch list from dictionary, including legacy keyword
    void readPatchSelection(const dictionary&);

    //- Return the coupled patch names for this update
    wordList selectedPatchNames() const;

    //- Return patch ID or abort with a clear error
    label patchID(const word&) const;

    //- Compute an outgoing-flux-weighted vapour mean over all selected patches
    scalar outgoingMean(const wordList& patchNames) const;

    //- Return group-mean old inlet value for consistent relaxation
    scalar oldInletMean(const wordList& patchNames) const;

    //- Clamp a value to [minValue_, maxValue_]
    scalar clamp(const scalar x) const;

public:

    TypeName("inletOutletCondensator");

    // Constructors

    inletOutletCondenserFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    inletOutletCondenserFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    inletOutletCondenserFvPatchScalarField
    (
        const inletOutletCondenserFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fieldMapper&
    );

    inletOutletCondenserFvPatchScalarField
    (
        const inletOutletCondenserFvPatchScalarField&
    ) = delete;

    inletOutletCondenserFvPatchScalarField
    (
        const inletOutletCondenserFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new inletOutletCondenserFvPatchScalarField(*this, iF)
        );
    }

    //- Return true: this patch field is altered by assignment
    virtual bool assignable() const
    {
        return true;
    }

    // Member Functions

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;

    // Member Operators

    virtual void operator=(const fvPatchField<scalar>&);
};

} // End namespace Foam

#endif

// ************************************************************************* //
