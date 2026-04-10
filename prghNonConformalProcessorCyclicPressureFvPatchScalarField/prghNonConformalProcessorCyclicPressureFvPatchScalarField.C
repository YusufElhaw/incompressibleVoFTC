/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "prghNonConformalProcessorCyclicPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "fvcSnGrad.H"
#include "fvMesh.H"
#include "fvSchemes.H"
#include "snGradScheme.H"
#include "UIPstream.H"
#include "UOPstream.H"

#include <cmath>

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


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
    rhoInf_(0),
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
    rhoInf_(dict.lookupOrDefault<scalar>("rhoInf", 0)),
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


// * * * * * * * * * * * * * * Private Helpers * * * * * * * * * * * * * * //

tmp<scalarField>
prghNonConformalProcessorCyclicPressureFvPatchScalarField::swapWithNeighbour
(
    const scalarField& fld
) const
{
    tmp<scalarField> tfld(new scalarField(fld.size(), Zero));
    scalarField& nbrFld = tfld.ref();

    if (!Pstream::parRun())
    {
        nbrFld = fld;
        return tfld;
    }

    label recvRequest = UPstream::nRequests();
    UIPstream::read
    (
        Pstream::commsTypes::nonBlocking,
        procPatch_.neighbProcNo(),
        reinterpret_cast<char*>(nbrFld.begin()),
        nbrFld.byteSize(),
        procPatch_.tag(),
        procPatch_.comm()
    );

    label sendRequest = UPstream::nRequests();
    UOPstream::write
    (
        Pstream::commsTypes::nonBlocking,
        procPatch_.neighbProcNo(),
        reinterpret_cast<const char*>(fld.begin()),
        fld.byteSize(),
        procPatch_.tag(),
        procPatch_.comm()
    );

    if (recvRequest >= 0 && recvRequest < Pstream::nRequests())
    {
        UPstream::waitRequest(recvRequest);
    }

    if (sendRequest >= 0 && sendRequest < Pstream::nRequests())
    {
        UPstream::waitRequest(sendRequest);
    }

    return tfld;
}


scalar
prghNonConformalProcessorCyclicPressureFvPatchScalarField::swapWithNeighbour
(
    const scalar s
) const
{
    scalarField local(1, s);
    return swapWithNeighbour(local)()[0];
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

tmp<Field<scalar>>
prghNonConformalProcessorCyclicPressureFvPatchScalarField::patchNeighbourField
(
    const Pstream::commsTypes
) const
{
    tmp<Field<scalar>> tLocal(this->patchInternalField());
    const scalarField nbrInternal(swapWithNeighbour(tLocal())());

    tmp<Field<scalar>> tNbr(new Field<scalar>(nbrInternal + jump_));
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
    const Pstream::commsTypes
)
{
    // no-op:
    // for this custom patch we exchange the needed neighbour values directly
    // inside patchNeighbourField()/updateInterfaceMatrix().
}



void prghNonConformalProcessorCyclicPressureFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    updateCoeffs();

    const tmp<Field<scalar>> tNbr = patchNeighbourField(commsType);

    // Set the boundary values directly from neighbour-internal + jump
    this->operator==(tNbr());
}


void prghNonConformalProcessorCyclicPressureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    Pout<< "prghNCPCC updateCoeffs patch=" << this->patch().name()
        << " owner_=" << (owner_ ? "true" : "false")
        << " rhoInfSpecified_=" << (rhoInfSpecified_ ? "true" : "false")
        << " neighbProcNo=" << procPatch_.neighbProcNo()
        << nl;

    const volScalarField& vf =
        static_cast<const volScalarField&>(this->internalField());

    const label patchi = this->patch().index();
    const label nFaces = this->patch().size();

    const volScalarField& rho =
        this->db().lookupObject<volScalarField>(rhoName_);
    const volScalarField::Boundary& rhoBf = rho.boundaryField();

    const surfaceScalarField& ghf =
        this->db().lookupObject<surfaceScalarField>("ghf");
    const surfaceScalarField::Boundary& ghfBf = ghf.boundaryField();

    const surfaceScalarField& rAUf =
        this->db().lookupObject<surfaceScalarField>("rAUf");
    const surfaceScalarField::Boundary& rAUfBf = rAUf.boundaryField();

    const surfaceScalarField& phiHbyA =
        this->db().lookupObject<surfaceScalarField>("phiHbyA");
    const surfaceScalarField::Boundary& phiHbyABf = phiHbyA.boundaryField();

    const tmp<surfaceScalarField> tDeltaCoeffs =
        fv::snGradScheme<scalar>::New
        (
            vf.mesh(),
            vf.mesh().schemes().snGrad(vf.name())
        )->deltaCoeffs(vf);

    const scalarField& deltaCoeffs = tDeltaCoeffs->boundaryField()[patchi];
    const scalarField& magSf = this->patch().magSf();
    const scalarField& w = this->patch().weights();

    const scalarField& rhoLocal = rhoBf[patchi];
    const scalarField& ghfLocal = ghfBf[patchi];
    const scalarField& rAUfLocal = rAUfBf[patchi];
    const scalarField& phiHbyALocal = phiHbyABf[patchi];

    if
    (
        rhoLocal.size() != nFaces
     || ghfLocal.size() != nFaces
     || rAUfLocal.size() != nFaces
     || phiHbyALocal.size() != nFaces
     || deltaCoeffs.size() != nFaces
     || magSf.size() != nFaces
     || w.size() != nFaces
    )
    {
        FatalErrorInFunction
            << "Size mismatch on patch " << this->patch().name() << nl
            << "patch size = " << nFaces << nl
            << "rhoLocal size = " << rhoLocal.size() << nl
            << "ghfLocal size = " << ghfLocal.size() << nl
            << "rAUfLocal size = " << rAUfLocal.size() << nl
            << "phiHbyALocal size = " << phiHbyALocal.size() << nl
            << "deltaCoeffs size = " << deltaCoeffs.size() << nl
            << "magSf size = " << magSf.size() << nl
            << "weights size = " << w.size()
            << abort(FatalError);
    }

    scalarField fluxLocal(nFaces, Zero);

    forAll(fluxLocal, facei)
    {
        if (!std::isfinite(rhoLocal[facei]))
        {
            FatalErrorInFunction
                << "rhoLocal non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (!std::isfinite(ghfLocal[facei]))
        {
            FatalErrorInFunction
                << "ghfLocal non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (!std::isfinite(rAUfLocal[facei]))
        {
            FatalErrorInFunction
                << "rAUfLocal non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (!std::isfinite(phiHbyALocal[facei]))
        {
            FatalErrorInFunction
                << "phiHbyALocal non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (!std::isfinite(magSf[facei]))
        {
            FatalErrorInFunction
                << "magSf non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (!std::isfinite(deltaCoeffs[facei]))
        {
            FatalErrorInFunction
                << "deltaCoeffs non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (!std::isfinite(w[facei]))
        {
            FatalErrorInFunction
                << "weights non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (mag(rAUfLocal[facei]) <= VSMALL)
        {
            FatalErrorInFunction
                << "rAUfLocal near zero on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (mag(magSf[facei]) <= VSMALL)
        {
            FatalErrorInFunction
                << "magSf near zero on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (mag(deltaCoeffs[facei]) <= VSMALL)
        {
            FatalErrorInFunction
                << "deltaCoeffs near zero on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        fluxLocal[facei] =
            phiHbyALocal[facei]/rAUfLocal[facei]/magSf[facei];

        if (!std::isfinite(fluxLocal[facei]))
        {
            FatalErrorInFunction
                << "fluxLocal non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }
    }

    const scalarField rhoNbr(swapWithNeighbour(rhoLocal)());
    const scalarField ghfNbr(swapWithNeighbour(ghfLocal)());
    const scalarField fluxNbr(swapWithNeighbour(fluxLocal)());

    if
    (
        rhoNbr.size() != nFaces
     || ghfNbr.size() != nFaces
     || fluxNbr.size() != nFaces
    )
    {
        FatalErrorInFunction
            << "Neighbour size mismatch on patch " << this->patch().name() << nl
            << "patch size = " << nFaces << nl
            << "rhoNbr size = " << rhoNbr.size() << nl
            << "ghfNbr size = " << ghfNbr.size() << nl
            << "fluxNbr size = " << fluxNbr.size()
            << abort(FatalError);
    }

    forAll(rhoNbr, facei)
    {
        if (!std::isfinite(rhoNbr[facei]))
        {
            FatalErrorInFunction
                << "rhoNbr non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (!std::isfinite(ghfNbr[facei]))
        {
            FatalErrorInFunction
                << "ghfNbr non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }

        if (!std::isfinite(fluxNbr[facei]))
        {
            FatalErrorInFunction
                << "fluxNbr non-finite on patch " << this->patch().name()
                << ", face " << facei
                << abort(FatalError);
        }
    }

    if (owner_ && !rhoInfSpecified_)
    {
        FatalErrorInFunction
            << "Patch " << this->patch().name()
            << " is owner-side, but rhoInf is missing."
            << abort(FatalError);
    }

    scalarField rhoInfSendField(1, owner_ ? rhoInf_ : 0.0);
    scalarField rhoInfRecvField(swapWithNeighbour(rhoInfSendField)());

    const scalar rhoInfUsed = owner_ ? rhoInf_ : rhoInfRecvField[0];

    if (!std::isfinite(rhoInfUsed))
    {
        FatalErrorInFunction
            << "rhoInfUsed non-finite on patch " << this->patch().name()
            << abort(FatalError);
    }

    jump_.setSize(nFaces);

    forAll(jump_, facei)
    {
        const bool rhoLocalFinite  = std::isfinite(rhoLocal[facei]);
        const bool ghfLocalFinite  = std::isfinite(ghfLocal[facei]);
        const bool ghfNbrFinite    = std::isfinite(ghfNbr[facei]);
        const bool fluxLocalFinite = std::isfinite(fluxLocal[facei]);
        const bool fluxNbrFinite   = std::isfinite(fluxNbr[facei]);
        const bool wFinite         = std::isfinite(w[facei]);
        const bool dFinite         = std::isfinite(deltaCoeffs[facei]);
        const bool rhoInfFinite    = std::isfinite(rhoInfUsed);

        if
        (
            !rhoLocalFinite
         || !ghfLocalFinite
         || !ghfNbrFinite
         || !fluxLocalFinite
         || !fluxNbrFinite
         || !wFinite
         || !dFinite
         || !rhoInfFinite
        )
        {
            FatalErrorInFunction
                << "Non-finite input before jump evaluation on patch "
                << this->patch().name() << ", face " << facei << nl
                << "rhoLocal finite  = " << rhoLocalFinite << nl
                << "ghfLocal finite  = " << ghfLocalFinite << nl
                << "ghfNbr finite    = " << ghfNbrFinite << nl
                << "fluxLocal finite = " << fluxLocalFinite << nl
                << "fluxNbr finite   = " << fluxNbrFinite << nl
                << "w finite         = " << wFinite << nl
                << "deltaCoeff finite= " << dFinite << nl
                << "rhoInfUsed finite= " << rhoInfFinite
                << abort(FatalError);
        }

        jump_[facei] =
            (rhoLocal[facei] - rhoInfUsed)
           *(ghfNbr[facei] - ghfLocal[facei])
          + (fluxLocal[facei] + fluxNbr[facei])
           *w[facei]/deltaCoeffs[facei];

        if (!std::isfinite(jump_[facei]))
        {
            FatalErrorInFunction
                << "Non-finite jump on patch " << this->patch().name()
                << ", face " << facei << nl
                << "rhoLocal finite  = " << rhoLocalFinite << nl
                << "ghfLocal finite  = " << ghfLocalFinite << nl
                << "ghfNbr finite    = " << ghfNbrFinite << nl
                << "fluxLocal finite = " << fluxLocalFinite << nl
                << "fluxNbr finite   = " << fluxNbrFinite << nl
                << "w finite         = " << wFinite << nl
                << "deltaCoeff finite= " << dFinite << nl
                << "rhoInfUsed finite= " << rhoInfFinite
                << abort(FatalError);
        }
    }

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

        //if (&psiInternal == &this->primitiveField())
        //{
            scalarReceiveBuf_ += jump_;
        //}

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

        //if (&psiInternal == &this->primitiveField())
        //{
            nbrPsi += jump_;
        //}

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

    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);

    if (owner_ && rhoInfSpecified_)
    {
        writeEntry(os, "rhoInf", rhoInf_);
    }

    writeEntry(os, "jump", jump_);
    writeEntry(os, "value", *this);
}


} // End namespace Foam

// ************************************************************************* //