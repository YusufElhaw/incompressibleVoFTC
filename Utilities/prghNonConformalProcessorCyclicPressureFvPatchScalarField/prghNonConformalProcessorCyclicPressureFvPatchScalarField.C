/*#include "prghNonConformalProcessorCyclicPressureFvPatchScalarField.H"
#include "prghNonConformalCyclicPressureFvPatchScalarField.H"

#include "fvMeshStitcherTools.H"
#include "nonConformalFvPatch.H"

#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "fvMesh.H"
#include "fvSchemes.H"
#include "snGradScheme.H"
#include "UIPstream.H"
#include "UOPstream.H"

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    prghNonConformalProcessorCyclicPressureFvPatchScalarField
);


// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(p, iF),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p)),
    rhoName_("rho"),
    rhoInf_(NaN),
    rhoInfSpecified_(false),
    owner_(procPatch_.owner()),
    jump_(p.size(), Zero),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(p, iF, dict),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p)),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    rhoInf_
    (
        dict.found("rhoInf")
      ? dict.lookup<scalar>("rhoInf", dimDensity)
      : NaN
    ),
    rhoInfSpecified_(dict.found("rhoInf")),
    owner_(procPatch_.owner()),
    jump_(p.size(), Zero),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{
    if (dict.found("jump"))
    {
        jump_ = scalarField("jump", iF.dimensions(), dict, p.size());
    }

    evaluateNoUpdateCoeffs();
}


prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(ptf, p, iF, mapper),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p)),
    rhoName_(ptf.rhoName_),
    rhoInf_(ptf.rhoInf_),
    rhoInfSpecified_(ptf.rhoInfSpecified_),
    owner_(ptf.owner_),
    jump_(mapper(ptf.jump_)),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(ptf, iF),
    procPatch_
    (
        refCast<const nonConformalProcessorCyclicFvPatch>(this->patch())
    ),
    rhoName_(ptf.rhoName_),
    rhoInf_(ptf.rhoInf_),
    rhoInfSpecified_(ptf.rhoInfSpecified_),
    owner_(ptf.owner_),
    jump_(ptf.jump_),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


// * * * * * * * * * * * * * * Private member functions * * * * * * * * * //

void
prghNonConformalProcessorCyclicPressureFvPatchScalarField::updateJumpFromReferPatch()
{
    const fvBoundaryMesh& fvbm = this->patch().boundaryMesh();

    const label referPatchi = procPatch_.referPatchIndex();

    if (referPatchi < 0 || referPatchi >= fvbm.size())
    {
        FatalErrorInFunction
            << "Invalid referPatchIndex " << referPatchi
            << " for patch " << this->patch().name()
            << abort(FatalError);
    }

    const fvPatch& referFvp = fvbm[referPatchi];

    if (!isA<nonConformalFvPatch>(referFvp))
    {
        FatalErrorInFunction
            << "Referred patch " << referFvp.name()
            << " of patch " << this->patch().name()
            << " is not a nonConformalFvPatch."
            << abort(FatalError);
    }

    const nonConformalFvPatch& thisNcFvp =
        refCast<const nonConformalFvPatch>(this->patch());

    const nonConformalFvPatch& referNcFvp =
        refCast<const nonConformalFvPatch>(referFvp);

    if (thisNcFvp.origPatchIndex() != referNcFvp.origPatchIndex())
    {
        FatalErrorInFunction
            << "Patch " << this->patch().name()
            << " and its referred patch " << referFvp.name()
            << " do not refer to the same original patch."
            << abort(FatalError);
    }

    const volScalarField& pRgh =
        static_cast<const volScalarField&>(this->internalField());

    volScalarField::Boundary& pRghBf =
        const_cast<volScalarField&>(pRgh).boundaryFieldRef();

    fvPatchScalarField& referBasePf = pRghBf[referPatchi];

    if (!isA<prghNonConformalCyclicPressureFvPatchScalarField>(referBasePf))
    {
        FatalErrorInFunction
            << "Referred patch field on patch " << referFvp.name()
            << " is type " << referBasePf.type()
            << ", but prghNonConformalCyclicPressure was expected."
            << abort(FatalError);
    }

    // Make sure the referred NCC patch has already computed its jump
    referBasePf.updateCoeffs();

    const prghNonConformalCyclicPressureFvPatchScalarField& referPrghPf =
        refCast<const prghNonConformalCyclicPressureFvPatchScalarField>
        (
            referBasePf
        );

    const fvPatch& origFvp = referNcFvp.origPatch();

    const labelList referOrigPatchFace
    (
        referNcFvp.polyFaces() - origFvp.start()
    );

    const labelList thisOrigPatchFace
    (
        thisNcFvp.polyFaces() - origFvp.start()
    );

    const scalarField& referArea = referNcFvp.patch().magSf();

    const scalarField origJumpSum
    (
        fvMeshStitcherTools::fieldRMapSum
        (
            referPrghPf.jump()*referArea,
            origFvp.size(),
            referOrigPatchFace
        )
    );

    const scalarField origAreaSum
    (
        fvMeshStitcherTools::fieldRMapSum
        (
            referArea,
            origFvp.size(),
            referOrigPatchFace
        )
    );

    scalarField origJump(origFvp.size(), Zero);

    forAll(origJump, origFacei)
    {
        if (origAreaSum[origFacei] > rootVSmall)
        {
            origJump[origFacei] =
                origJumpSum[origFacei]/origAreaSum[origFacei];
        }
    }

    jump_ = fvMeshStitcherTools::fieldMap(origJump, thisOrigPatchFace);
}


// * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * //

tmp<Field<scalar>>
prghNonConformalProcessorCyclicPressureFvPatchScalarField::patchNeighbourField
(
    const Pstream::commsTypes
) const
{
    if (debug && !this->ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name()
            << " outstanding request."
            << abort(FatalError);
    }

    tmp<Field<scalar>> tNbr(new scalarField(*this));
    tNbr.ref() += jump_;
    return tNbr;
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    nonConformalProcessorCyclicFvPatchField<scalar>::map(ptf, mapper);

    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& rhs =
        refCast<const prghNonConformalProcessorCyclicPressureFvPatchScalarField>
        (
            ptf
        );

    mapper(jump_, rhs.jump_);
    owner_ = rhs.owner_;
    rhoName_ = rhs.rhoName_;
    rhoInf_ = rhs.rhoInf_;
    rhoInfSpecified_ = rhs.rhoInfSpecified_;
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    nonConformalProcessorCyclicFvPatchField<scalar>::reset(ptf);

    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& rhs =
        refCast<const prghNonConformalProcessorCyclicPressureFvPatchScalarField>
        (
            ptf
        );

    jump_.reset(rhs.jump_);
    owner_ = rhs.owner_;
    rhoName_ = rhs.rhoName_;
    rhoInf_ = rhs.rhoInf_;
    rhoInfSpecified_ = rhs.rhoInfSpecified_;
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    nonConformalProcessorCyclicFvPatchField<scalar>::initEvaluate(commsType);
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    updateCoeffs();
    nonConformalProcessorCyclicFvPatchField<scalar>::evaluate(commsType);
}


tmp<Field<scalar>>
prghNonConformalProcessorCyclicPressureFvPatchScalarField::snGrad() const
{
    return snGrad(this->patch().deltaCoeffs());
}


tmp<Field<scalar>>
prghNonConformalProcessorCyclicPressureFvPatchScalarField::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    const scalarField pif(this->patchInternalField());

    scalarField nbrPf(*this);
    nbrPf += jump_;

    return deltaCoeffs*(nbrPf - pif);
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    updateJumpFromReferPatch();

    nonConformalProcessorCyclicFvPatchField<scalar>::updateCoeffs();
}


void
prghNonConformalProcessorCyclicPressureFvPatchScalarField::initInterfaceMatrixUpdate
(
    scalarField&,
    const scalarField& psiInternal,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    this->patch().patchInternalField(psiInternal, scalarSendBuf_);

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        if (debug && !this->ready())
        {
            FatalErrorInFunction
                << "On patch " << procPatch_.name()
                << " outstanding request."
                << abort(FatalError);
        }

        scalarReceiveBuf_.setSize(scalarSendBuf_.size());

        outstandingRecvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<char*>(scalarReceiveBuf_.begin()),
            scalarReceiveBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(scalarSendBuf_.begin()),
            scalarSendBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.compressedSend(commsType, scalarSendBuf_);
    }

    const_cast
    <prghNonConformalProcessorCyclicPressureFvPatchScalarField&>(*this)
        .updatedMatrix() = false;
}


void
prghNonConformalProcessorCyclicPressureFvPatchScalarField::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    const labelUList& faceCells = this->patch().faceCells();

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            UPstream::waitRequest(outstandingRecvRequest_);
        }

        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;

        transformCoupleField(scalarReceiveBuf_, cmpt);

        if (&psiInternal == &this->primitiveField())
        {
            scalarReceiveBuf_ += jump_;
        }

        forAll(faceCells, facei)
        {
            result[faceCells[facei]] -= coeffs[facei]*scalarReceiveBuf_[facei];
        }
    }
    else
    {
        scalarField nbrPsi
        (
            procPatch_.compressedReceive<scalar>(commsType, this->size())()
        );

        transformCoupleField(nbrPsi, cmpt);

        if (&psiInternal == &this->primitiveField())
        {
            nbrPsi += jump_;
        }

        forAll(faceCells, facei)
        {
            result[faceCells[facei]] -= coeffs[facei]*nbrPsi[facei];
        }
    }

    const_cast
    <prghNonConformalProcessorCyclicPressureFvPatchScalarField&>(*this)
        .updatedMatrix() = true;
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "patchType", this->patch().type());
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    if (owner_ && rhoInfSpecified_)
    {
        writeEntry(os, "rhoInf", rhoInf_);
    }

    writeEntry(os, "jump", jump_);
    writeEntry(os, "value", *this);
}

} */// End namespace Foam
/* //Gemini
 #include "prghNonConformalProcessorCyclicPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "fvMesh.H"
#include "fvSchemes.H"
#include "snGradScheme.H"
#include "UIPstream.H"
#include "UOPstream.H"
#include "volFields.H" // Erforderlich für den Zugriff auf volScalarField

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    prghNonConformalProcessorCyclicPressureFvPatchScalarField
);


// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(p, iF),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p)),
    rhoName_("rho"),
    rhoInf_(NaN),
    rhoInfSpecified_(false),
    owner_(procPatch_.owner()),
    jump_(p.size(), Zero),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(p, iF, dict),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p)),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    rhoInf_
    (
        dict.found("rhoInf")
      ? dict.lookup<scalar>("rhoInf", dimDensity)
      : NaN
    ),
    rhoInfSpecified_(dict.found("rhoInf")),
    owner_(procPatch_.owner()),
    jump_(p.size(), Zero),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{
    if (dict.found("jump"))
    {
        jump_ = scalarField("jump", iF.dimensions(), dict, p.size());
    }

    evaluateNoUpdateCoeffs();
}


prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(ptf, p, iF, mapper),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p)),
    rhoName_(ptf.rhoName_),
    rhoInf_(ptf.rhoInf_),
    rhoInfSpecified_(ptf.rhoInfSpecified_),
    owner_(ptf.owner_),
    jump_(mapper(ptf.jump_)),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(ptf, iF),
    procPatch_
    (
        refCast<const nonConformalProcessorCyclicFvPatch>(this->patch())
    ),
    rhoName_(ptf.rhoName_),
    rhoInf_(ptf.rhoInf_),
    rhoInfSpecified_(ptf.rhoInfSpecified_),
    owner_(ptf.owner_),
    jump_(ptf.jump_),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


// * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * //

tmp<Field<scalar>>
prghNonConformalProcessorCyclicPressureFvPatchScalarField::patchNeighbourField
(
    const Pstream::commsTypes commsType
) const
{
    if (debug && !this->ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name()
            << " outstanding request."
            << abort(FatalError);
    }

    // 1. Hole die Werte vom anderen Prozessor (ist evtl. nur eine const Referenz)
    tmp<Field<scalar>> tNbrBase =
        nonConformalProcessorCyclicFvPatchField<scalar>::patchNeighbourField(commsType);

    // 2. Erzeuge eine explizite Kopie des Feldes im Speicher
    tmp<Field<scalar>> tNbr(new Field<scalar>(tNbrBase()));

    // 3. Jetzt können wir sicher unseren hydrostatischen Sprung addieren
    tNbr.ref() += jump_;
    
    return tNbr;
}

void prghNonConformalProcessorCyclicPressureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    nonConformalProcessorCyclicFvPatchField<scalar>::map(ptf, mapper);

    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& rhs =
        refCast<const prghNonConformalProcessorCyclicPressureFvPatchScalarField>
        (
            ptf
        );

    mapper(jump_, rhs.jump_);
    owner_ = rhs.owner_;
    rhoName_ = rhs.rhoName_;
    rhoInf_ = rhs.rhoInf_;
    rhoInfSpecified_ = rhs.rhoInfSpecified_;
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    nonConformalProcessorCyclicFvPatchField<scalar>::reset(ptf);

    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& rhs =
        refCast<const prghNonConformalProcessorCyclicPressureFvPatchScalarField>
        (
            ptf
        );

    jump_.reset(rhs.jump_);
    owner_ = rhs.owner_;
    rhoName_ = rhs.rhoName_;
    rhoInf_ = rhs.rhoInf_;
    rhoInfSpecified_ = rhs.rhoInfSpecified_;
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    nonConformalProcessorCyclicFvPatchField<scalar>::initEvaluate(commsType);
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    updateCoeffs();
    nonConformalProcessorCyclicFvPatchField<scalar>::evaluate(commsType);
}


tmp<Field<scalar>>
prghNonConformalProcessorCyclicPressureFvPatchScalarField::snGrad() const
{
    return snGrad(this->patch().deltaCoeffs());
}


tmp<Field<scalar>>
prghNonConformalProcessorCyclicPressureFvPatchScalarField::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    const scalarField pif(this->patchInternalField());

    scalarField nbrPf(*this);
    nbrPf += jump_;

    return deltaCoeffs*(nbrPf - pif);
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatchScalarField& rhop =
        this->patch().template lookupPatchField<volScalarField, scalar>(rhoName_);

    const fvPatchScalarField& ghp =
        this->patch().template lookupPatchField<volScalarField, scalar>("gh");

    scalar usedRhoInf = rhoInfSpecified_ ? rhoInf_ : 0.0;

    // KORREKTUR: Wir entfernen die Mittelung der Dichte mit dem Nachbar-Patch!
    // Nutze ausschließlich die exakte, lokale Dichte (rhop) dieser Zelle, 
    // um das hydrostatische Gleichgewicht für diese spezifische Face zu erhalten.
    jump_ =
        (rhop - usedRhoInf)
       *(ghp - ghp.patchNeighbourField());

    nonConformalProcessorCyclicFvPatchField<scalar>::updateCoeffs();
}


void
prghNonConformalProcessorCyclicPressureFvPatchScalarField::initInterfaceMatrixUpdate
(
    scalarField&,
    const scalarField& psiInternal,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    this->patch().patchInternalField(psiInternal, scalarSendBuf_);

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        if (debug && !this->ready())
        {
            FatalErrorInFunction
                << "On patch " << procPatch_.name()
                << " outstanding request."
                << abort(FatalError);
        }

        scalarReceiveBuf_.setSize(scalarSendBuf_.size());

        outstandingRecvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<char*>(scalarReceiveBuf_.begin()),
            scalarReceiveBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            reinterpret_cast<const char*>(scalarSendBuf_.begin()),
            scalarSendBuf_.byteSize(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.compressedSend(commsType, scalarSendBuf_);
    }

    const_cast
    <prghNonConformalProcessorCyclicPressureFvPatchScalarField&>(*this)
        .updatedMatrix() = false;
}


void
prghNonConformalProcessorCyclicPressureFvPatchScalarField::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    const labelUList& faceCells = this->patch().faceCells();

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            UPstream::waitRequest(outstandingRecvRequest_);
        }

        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;

        transformCoupleField(scalarReceiveBuf_, cmpt);

        if (&psiInternal == &this->primitiveField())
        {
            scalarReceiveBuf_ += jump_;
        }

        forAll(faceCells, facei)
        {
            result[faceCells[facei]] -= coeffs[facei]*scalarReceiveBuf_[facei];
        }
    }
    else
    {
        scalarField nbrPsi
        (
            procPatch_.compressedReceive<scalar>(commsType, this->size())()
        );

        transformCoupleField(nbrPsi, cmpt);

        if (&psiInternal == &this->primitiveField())
        {
            nbrPsi += jump_;
        }

        forAll(faceCells, facei)
        {
            result[faceCells[facei]] -= coeffs[facei]*nbrPsi[facei];
        }
    }

    const_cast
    <prghNonConformalProcessorCyclicPressureFvPatchScalarField&>(*this)
        .updatedMatrix() = true;
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "patchType", this->patch().type());
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    if (owner_ && rhoInfSpecified_)
    {
        writeEntry(os, "rhoInf", rhoInf_);
    }

    writeEntry(os, "jump", jump_);
    writeEntry(os, "value", *this);
}

} // End namespace Foam
*/













/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  13
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
#include "prghNonConformalProcessorCyclicPressureFvPatchScalarField.H"
#include "prghNonConformalCyclicPressureFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "nonConformalProcessorCyclicFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "volFields.H"
#include "surfaceFields.H"

#include <cmath>
#include "fvMesh.H"
#include "fvSchemes.H"
#include "snGradScheme.H"


namespace Foam
{

defineTypeNameAndDebug
(
    prghNonConformalProcessorCyclicPressureFvPatchScalarField,
    0
);

addToRunTimeSelectionTable
(
    fvPatchScalarField,
    prghNonConformalProcessorCyclicPressureFvPatchScalarField,
    dictionary
);


// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(p, iF),
    rhoName_("rho"),
    rhoInf_(NaN),
    jump_(p.size(), Zero)
{}


prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(p, iF, dict),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    rhoInf_
    (
        dict.found("rhoInf")
      ? dict.lookup<scalar>("rhoInf", dimDensity)
      : NaN
    ),
    jump_(p.size(), Zero)
{}


prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_),
    rhoInf_(ptf.rhoInf_),
    jump_(mapper(ptf.jump_))
{}


prghNonConformalProcessorCyclicPressureFvPatchScalarField::
prghNonConformalProcessorCyclicPressureFvPatchScalarField
(
    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nonConformalProcessorCyclicFvPatchField<scalar>(ptf, iF),
    rhoName_(ptf.rhoName_),
    rhoInf_(ptf.rhoInf_),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * //

tmp<Field<scalar>>
prghNonConformalProcessorCyclicPressureFvPatchScalarField::patchNeighbourField
(
    const Pstream::commsTypes commsType
) const
{
    return
        nonConformalProcessorCyclicFvPatchField<scalar>::patchNeighbourField
        (
            commsType
        );
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    nonConformalProcessorCyclicFvPatchField<scalar>::map(ptf, mapper);

    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& rhs =
        refCast
        <
            const prghNonConformalProcessorCyclicPressureFvPatchScalarField
        >(ptf);

    mapper(jump_, rhs.jump_);
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    nonConformalProcessorCyclicFvPatchField<scalar>::reset(ptf);

    const prghNonConformalProcessorCyclicPressureFvPatchScalarField& rhs =
        refCast
        <
            const prghNonConformalProcessorCyclicPressureFvPatchScalarField
        >(ptf);

    jump_.reset(rhs.jump_);
}



void prghNonConformalProcessorCyclicPressureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if
    (
        !this->db().foundObject<surfaceScalarField>("ghf")
     || !this->db().foundObject<surfaceScalarField>("rAUf")
     || !this->db().foundObject<surfaceScalarField>("phiHbyA")
    )
    {
        jump_ = scalarField(this->patch().size(), Zero);
        nonConformalProcessorCyclicFvPatchField<scalar>::updateCoeffs();
        return;
    }

    const volScalarField& vf =
        static_cast<const volScalarField&>(this->internalField());

    const label patchi = this->patch().index();

    if (!this->db().foundObject<volScalarField>(rhoName_))
    {
        jump_ = scalarField(this->patch().size(), Zero);
        nonConformalProcessorCyclicFvPatchField<scalar>::updateCoeffs();
        return;
    }

    const volScalarField& rhoVf =
        this->db().lookupObject<volScalarField>(rhoName_);

    const volScalarField::Boundary& rhoBf = rhoVf.boundaryField();

    const surfaceScalarField& ghf =
        this->db().lookupObject<surfaceScalarField>("ghf");

    const surfaceScalarField::Boundary& ghfBf = ghf.boundaryField();

    const tmp<surfaceScalarField::Boundary> tGhfNbrBf =
        ghf.boundaryField().boundaryNeighbourField();

    const scalarField ghfNbr(tGhfNbrBf()[patchi]);

    const scalarField hydroTerm =
        (rhoBf[patchi] - rhoInf_)*(ghfNbr - ghfBf[patchi]);

    jump_ = hydroTerm;

    nonConformalProcessorCyclicFvPatchField<scalar>::updateCoeffs();
if (this->patch().name() == "procBoundary0to1throughNCC_R4S1_on_PatchBR4S1"
 || this->patch().name() == "procBoundary1to0throughNCC_R4S1_on_PatchTR4S1")
{
Pout<< " patch " << this->patch().name()
    << " hydro min/max = " << gMin(hydroTerm) << " " << gMax(hydroTerm)
   // << " flux min/max = " << gMin(fluxTerm)  << " " << gMax(fluxTerm)
    << " jump min/max = " << gMin(jump_)     << " " << gMax(jump_)
    << endl;
}
}



void prghNonConformalProcessorCyclicPressureFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        updateCoeffs();
    }

    nonConformalProcessorCyclicFvPatchField<scalar>::evaluate(commsType);

    if (this->patch().name() == "procBoundary0to1throughNCC_R4S1_on_PatchBR4S1"
 || this->patch().name() == "procBoundary1to0throughNCC_R4S1_on_PatchTR4S1")
{
    Pout<< " evaluate patch " << this->patch().name()
    << " before add min/max = " << gMin(*this) << " " << gMax(*this)
    << endl;

}
    operator+=(jump_);


        if (this->patch().name() == "procBoundary0to1throughNCC_R4S1_on_PatchBR4S1"
 || this->patch().name() == "procBoundary1to0throughNCC_R4S1_on_PatchTR4S1")
{
    Pout<< " evaluate patch " << this->patch().name()
    << " after add min/max = " << gMin(*this) << " " << gMax(*this)
    << endl;

}
}

void prghNonConformalProcessorCyclicPressureFvPatchScalarField::
updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    // Normal processor NCC contribution first
    nonConformalProcessorCyclicFvPatchField<scalar>::updateInterfaceMatrix
    (
        result,
        psiInternal,
        coeffs,
        cmpt,
        commsType
    );

    // Add only the missing jump term for the original field
    if (&psiInternal == &this->primitiveField())
    {
        const labelUList& faceCells = this->patch().faceCells();

        forAll(faceCells, facei)
        {
            result[faceCells[facei]] -= coeffs[facei]*jump_[facei];
        }
    }
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "patchType", word("nonConformalProcessorCyclic"));
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    if (!std::isnan(rhoInf_))
    {
        writeEntry(os, "rhoInf", rhoInf_);
    }

    writeEntry(os, "value", *this);
}


} // End namespace Foam