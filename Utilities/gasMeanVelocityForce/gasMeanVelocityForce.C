/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "gasMeanVelocityForce.H"

#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IFstream.H"
#include "IOdictionary.H"
#include "timeIOdictionary.H"
#include "volFields.H"
#include "zero.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(gasMeanVelocityForce, 0);
    addToRunTimeSelectionTable(fvConstraint, gasMeanVelocityForce, dictionary);
}
}


void Foam::fv::gasMeanVelocityForce::readCoeffs(const dictionary& dict)
{
    UName_ = dict.lookupOrDefault<word>("U", "U");
    alphaName_ = dict.lookup<word>("alpha");
    Ubar_ = dict.lookup<vector>("Ubar");
    relaxation_ = dict.lookupOrDefault<scalar>("relaxation", 1);
    alphaCutoff_ = dict.lookupOrDefault<scalar>("alphaCutoff", 0);
    useCutoff_ = dict.lookupOrDefault<Switch>("useCutoff", false);
}

Foam::scalar Foam::fv::gasMeanVelocityForce::weight(const scalar alphaL) const
{
    const scalar aG = max(scalar(0), min(scalar(1), scalar(1) - alphaL));

    if (useCutoff_ && aG <= alphaCutoff_)
    {
        return 0;
    }

    return aG;
}


Foam::scalar Foam::fv::gasMeanVelocityForce::weight
(
    const volScalarField& alpha,
    const label celli
) const
{
    return weight(alpha[celli]);
}


Foam::scalar Foam::fv::gasMeanVelocityForce::magUbarGasAve
(
    const volVectorField& U,
    const volScalarField& alpha
) const
{
    const labelList& cells = zone_.zone(); 
    const scalarField& cv = mesh().V();  // field of cell volumes (fvMesh::V())
    const vector direction = normalised(Ubar_);

    scalar jG = 0;

    forAll(cells, i)
    {
        const label celli = cells[i];
        jG += weight(alpha, celli)*(direction & U[celli])*cv[celli];
    }

    reduce(jG, sumOp<scalar>());
    jG /= zone_.V();    // devision by total volume of the cell zone (fvCellZone::V())

    return jG;
}


void Foam::fv::gasMeanVelocityForce::writeProps(const scalar gradP) const
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
        propsDict.regIOobject::write();
    }
}


Foam::fv::gasMeanVelocityForce::gasMeanVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(sourceName, modelType, mesh, dict),
    zone_(mesh, coeffs(dict)),
    UName_(word::null),
    alphaName_(word::null),
    Ubar_(vector::uniform(NaN)),
    relaxation_(NaN),
    alphaCutoff_(0),
    useCutoff_(false),
    gradP0_(0),
    dGradP_(0),
    rAPtr_(nullptr)
{
    readCoeffs(coeffs(dict));

    IFstream propsFile
    (
        mesh.time().timePath()/"uniform"/(this->name() + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "Reading pressure gradient from file" << endl;
        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("gradient") >> gradP0_;
    }

    Info<< "Initial pressure gradient = " << gradP0_ << nl << endl;
}


Foam::wordList Foam::fv::gasMeanVelocityForce::constrainedFields() const
{
    return wordList(1, UName_);
}


bool Foam::fv::gasMeanVelocityForce::constrain
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (fieldName != UName_)
    {
        return false;
    }

    const volScalarField& alpha = mesh().lookupObject<volScalarField>(alphaName_);
    const labelList& cells = zone_.zone();
    const scalar gradP = gradP0_ + dGradP_;
    const vector direction = normalised(Ubar_);

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
        Su[celli] = weight(alpha, celli)*direction*gradP;
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


bool Foam::fv::gasMeanVelocityForce::constrain(volVectorField& U) const
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

    const volScalarField& alpha = mesh().lookupObject<volScalarField>(alphaName_);
    const scalarField& rAU = rAPtr_();
    const labelList& cells = zone_.zone();
    const scalarField& cv = mesh().V();
    const vector direction = normalised(Ubar_);

    scalar rAUeff = 0;
    scalar weightVol = 0;

    forAll(cells, i)
    {
        const label celli = cells[i];
        const scalar w = weight(alpha, celli);
        rAUeff += w*rAU[celli]*cv[celli];
        weightVol += w*cv[celli];
    }

    reduce(rAUeff, sumOp<scalar>());
    reduce(weightVol, sumOp<scalar>());

    if (weightVol <= SMALL)
    {
        WarningInFunction
            << "Effective gas volume in zone is too small. No forcing update "
            << "applied for " << name() << endl;
        return true;
    }

    rAUeff /= weightVol;

    const scalar jGbarAve = this->magUbarGasAve(U, alpha);
    dGradP_ = relaxation_*(mag(Ubar_) - jGbarAve)/max(rAUeff, SMALL);

    forAll(cells, i)
    {
        const label celli = cells[i];
        const scalar w = weight(alpha, celli);
        U[celli] += w*direction*rAU[celli]*dGradP_;
    }

    const scalar gradP = gradP0_ + dGradP_;

    Info<< "Gas pressure gradient source: uncorrected jG = " << jGbarAve
        << ", target jG = " << mag(Ubar_)
        << ", pressure gradient = " << gradP << endl;

    writeProps(gradP);

    return true;
}


bool Foam::fv::gasMeanVelocityForce::movePoints()
{
    zone_.movePoints();
    return true;
}


void Foam::fv::gasMeanVelocityForce::topoChange(const polyTopoChangeMap& map)
{
    zone_.topoChange(map);
}


void Foam::fv::gasMeanVelocityForce::mapMesh(const polyMeshMap& map)
{
    zone_.mapMesh(map);
}


void Foam::fv::gasMeanVelocityForce::distribute(const polyDistributionMap& map)
{
    zone_.distribute(map);
}


bool Foam::fv::gasMeanVelocityForce::read(const dictionary& dict)
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
