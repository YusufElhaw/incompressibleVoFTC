#include "prghNonConformalCyclicPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "fvcSnGrad.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "fvMesh.H"
#include "fvSchemes.H"
#include "snGradScheme.H"

namespace Foam
{

defineTypeNameAndDebug(prghNonConformalCyclicPressureFvPatchScalarField, 0);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    prghNonConformalCyclicPressureFvPatchScalarField,
    dictionary
);


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

prghNonConformalCyclicPressureFvPatchScalarField::
prghNonConformalCyclicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nonConformalCyclicFvPatchField<scalar>(p, iF),
    rhoName_("rho"),
    rhoInf_(NaN),
    jump_(p.size(), Zero)
{}


prghNonConformalCyclicPressureFvPatchScalarField::
prghNonConformalCyclicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nonConformalCyclicFvPatchField<scalar>(p, iF, dict),
    rhoName_
    (
        this->cyclicPatch().owner()
      ? dict.lookupOrDefault<word>("rho", "rho")
      : word::null
    ),
    rhoInf_
    (
        this->cyclicPatch().owner()
      ? dict.lookup<scalar>("rhoInf", dimDensity)
      : NaN
    ),
    jump_(p.size(), Zero)
{
    if (dict.found("jump"))
    {
        jump_ = scalarField("jump", iF.dimensions(), dict, p.size());
    }
    else
    {
        jump_ = scalarField(p.size(), Zero);
}

   // updateCoeffs();
}


prghNonConformalCyclicPressureFvPatchScalarField::
prghNonConformalCyclicPressureFvPatchScalarField
(
    const prghNonConformalCyclicPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    nonConformalCyclicFvPatchField<scalar>(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_),
    rhoInf_(ptf.rhoInf_),
    jump_(mapper(ptf.jump_))
{}


prghNonConformalCyclicPressureFvPatchScalarField::
prghNonConformalCyclicPressureFvPatchScalarField
(
    const prghNonConformalCyclicPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nonConformalCyclicFvPatchField<scalar>(ptf, iF),
    rhoName_(ptf.rhoName_),
    rhoInf_(ptf.rhoInf_),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //

tmp<Field<scalar>>
prghNonConformalCyclicPressureFvPatchScalarField::patchNeighbourField
(
    const Pstream::commsTypes
) const
{
    const VolField<scalar>& vf =
        static_cast<const VolField<scalar>&>(this->internalField());

    const labelUList& nbrFaceCells =
        this->cyclicPatch().neighbFvPatch().faceCells();

    return
        this->transform().transform(Field<scalar>(vf, nbrFaceCells))
      + jump_;
}


void prghNonConformalCyclicPressureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    nonConformalCyclicFvPatchField<scalar>::map(ptf, mapper);

    const prghNonConformalCyclicPressureFvPatchScalarField& rhs =
        refCast<const prghNonConformalCyclicPressureFvPatchScalarField>(ptf);

    mapper(jump_, rhs.jump_);
}


void prghNonConformalCyclicPressureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    nonConformalCyclicFvPatchField<scalar>::reset(ptf);

    const prghNonConformalCyclicPressureFvPatchScalarField& rhs =
        refCast<const prghNonConformalCyclicPressureFvPatchScalarField>(ptf);

    jump_.reset(rhs.jump_);
}


void prghNonConformalCyclicPressureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const volScalarField& vf =
        static_cast<const volScalarField&>(this->internalField());

    const prghNonConformalCyclicPressureFvPatchScalarField& nbrPf =
        refCast<const prghNonConformalCyclicPressureFvPatchScalarField>
        (
            this->nbrPatchField()
        );

    const label patchi = this->patch().index();
    const label nbrPatchi = nbrPf.patch().index();

    // Buoyancy fields
    const volScalarField& rhoVf =
        this->db().lookupObject<volScalarField>
        (
            this->cyclicPatch().owner() ? rhoName_ : nbrPf.rhoName_
        );

    const volScalarField::Boundary& rhoBf = rhoVf.boundaryField();

    const surfaceScalarField::Boundary& ghfBf =
        this->db().lookupObject<surfaceScalarField>("ghf").boundaryField();

    // Pressure solution fields
    const surfaceScalarField::Boundary& rAUfBf =
        this->db().lookupObject<surfaceScalarField>("rAUf").boundaryField();

    const surfaceScalarField::Boundary& phiHbyABf =
        this->db().lookupObject<surfaceScalarField>("phiHbyA").boundaryField();

    // Delta coeffs
    const tmp<surfaceScalarField> deltaCoeffsSf =
        fv::snGradScheme<scalar>::New
        (
            vf.mesh(),
            vf.mesh().schemes().snGrad(this->internalField().name())
        )->deltaCoeffs(vf);

    const scalarField& deltaCoeffsPf =
        deltaCoeffsSf->boundaryField()[patchi];

    // Same formula as prghCyclicPressure
    jump_ =
        (rhoBf[patchi] - (this->cyclicPatch().owner() ? rhoInf_ : nbrPf.rhoInf_))
       *(ghfBf[nbrPatchi] - ghfBf[patchi])
      + (
            phiHbyABf[patchi]
           /rAUfBf[patchi]
           /this->patch().magSf()
          + phiHbyABf[nbrPatchi]
           /rAUfBf[nbrPatchi]
           /nbrPf.patch().magSf()
        )
       *this->patch().weights()
       /deltaCoeffsPf;

    nonConformalCyclicFvPatchField<scalar>::updateCoeffs();
}


void prghNonConformalCyclicPressureFvPatchScalarField::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes
) const
{
    const labelUList& faceCells = this->patch().faceCells();
    const labelUList& nbrFaceCells = this->cyclicPatch().neighbFvPatch().faceCells();

    scalarField nbrPf(psiInternal, nbrFaceCells);

    // Only apply jump to the original field
    if (&psiInternal == &this->primitiveField())
    {
        nbrPf += jump_;
    }

    forAll(faceCells, facei)
    {
        result[faceCells[facei]] -= coeffs[facei]*nbrPf[facei];
    }
}


void prghNonConformalCyclicPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    if (this->cyclicPatch().owner())
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntry(os, "rhoInf", rhoInf_);
    }

    writeEntry(os, "value", *this);
}


} // End namespace Foam