/*---------------------------------------------------------------------------*\
  Reservoir Fenske function object for incompressibleVoFTC.
\*---------------------------------------------------------------------------*/

#include "reservoirFenske.H"
#include "addToRunTimeSelectionTable.H"
#include "OSspecific.H"
#include <cmath>

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(reservoirFenske, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        reservoirFenske,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private member functions  * * * * * * * * * * //

Foam::label
Foam::functionObjects::reservoirFenske::patchID
(
    const word& patchName
) const
{
    const label patchi = mesh().boundaryMesh().findIndex(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch " << patchName
            << " not found in mesh boundary."
            << exit(FatalError);
    }

    return patchi;
}


bool
Foam::functionObjects::reservoirFenske::volumeIntegral
(
    scalar& result,
    const word& fieldName
) const
{
    if (!mesh().foundObject<volScalarField>(fieldName))
    {
        Info<< type() << " " << name()
            << ": field " << fieldName << " not found yet, skipping."
            << nl;

        return false;
    }

    const volScalarField& f =
        mesh().lookupObject<volScalarField>(fieldName);

    const scalarField& V  = mesh().V();
    const scalarField& fi = f.internalField();

    scalar sum = 0;

    forAll(fi, celli)
    {
        sum += V[celli]*fi[celli];
    }

    reduce(sum, sumOp<scalar>());

    result = sum;
    return true;
}


bool
Foam::functionObjects::reservoirFenske::volumeGeoMean
(
    scalar& result,
    const word& fieldName
) const
{
    if (!mesh().foundObject<volScalarField>(fieldName))
    {
        Info<< type() << " " << name()
            << ": field " << fieldName << " not found yet, skipping."
            << nl;

        return false;
    }

    const volScalarField& f =
        mesh().lookupObject<volScalarField>(fieldName);

    const scalarField& V  = mesh().V();
    const scalarField& fi = f.internalField();

    scalar sumLog = 0;
    scalar sumV   = 0;

    forAll(fi, celli)
    {
        if (fi[celli] > SMALL)
        {
            sumLog += V[celli]*std::log(fi[celli]);
            sumV   += V[celli];
        }
    }

    reduce(sumLog, sumOp<scalar>());
    reduce(sumV,   sumOp<scalar>());

    if (sumV > VSMALL)
    {
        result = std::exp(sumLog/sumV);
        return true;
    }

    return false;
}


bool
Foam::functionObjects::reservoirFenske::readCondenserReservoir
(
    scalar& ND,
    scalar& xD,
    const volScalarField& X,
    const wordList& patchNames
) const
{
    forAll(patchNames, i)
    {
        const label patchi = patchID(patchNames[i]);

        const fvPatchScalarField& pf = X.boundaryField()[patchi];

        if (isA<inletOutletCondenserFvPatchScalarField>(pf))
        {
            const inletOutletCondenserFvPatchScalarField& bc =
                refCast<const inletOutletCondenserFvPatchScalarField>(pf);

            ND = bc.reservoirMoles();
            xD = bc.reservoirX();

            return true;
        }
    }

    Info<< type() << " " << name()
        << ": no inletOutletCondenser BC found on topReservoirPatches "
        << patchNames << nl;

    return false;
}


bool
Foam::functionObjects::reservoirFenske::readBoilerReservoir
(
    scalar& NB,
    scalar& xB,
    scalar& yB,
    const volScalarField& X,
    const wordList& patchNames
) const
{
    forAll(patchNames, i)
    {
        const label patchi = patchID(patchNames[i]);

        const fvPatchScalarField& pf = X.boundaryField()[patchi];

        if (isA<inletOutletBoilerFvPatchScalarField>(pf))
        {
            const inletOutletBoilerFvPatchScalarField& bc =
                refCast<const inletOutletBoilerFvPatchScalarField>(pf);

            NB = bc.reservoirMoles();
            xB = bc.reservoirX();
            yB = bc.reservoirY();

            return true;
        }
    }

    Info<< type() << " " << name()
        << ": no inletOutletBoiler BC found on bottomReservoirPatches "
        << patchNames << nl;

    return false;
}


Foam::scalar
Foam::functionObjects::reservoirFenske::fenskeNmin
(
    const scalar xD,
    const scalar xB,
    const scalar A12
)
{
    return
        std::log((xD/(1 - xD))*((1 - xB)/xB))
       /std::log(A12);
}


bool
Foam::functionObjects::reservoirFenske::validFenskeInputs
(
    const scalar xD,
    const scalar xB,
    const scalar A12
)
{
    return
        xD > 0
     && xD < 1
     && xB > 0
     && xB < 1
     && xD > xB
     && A12 > 1;
}


void
Foam::functionObjects::reservoirFenske::writeValueOrNa
(
    std::ofstream& os,
    const bool valid,
    const scalar value
)
{
    if (valid)
    {
        os << value;
    }
    else
    {
        os << "nan";
    }
}


void
Foam::functionObjects::reservoirFenske::readCoeffs
(
    const dictionary& dict
)
{
    speciesField_ =
        dict.lookupOrDefault<word>("speciesField", "cyclohexane.X");

    A12Field_ =
        dict.lookupOrDefault<word>("A12Field", "A12");

    columnTotalMolesField_ =
        dict.lookupOrDefault<word>("columnTotalMolesField", "nMolesTotal");

    columnSpecies1MolesField_ =
        dict.lookupOrDefault<word>("columnSpecies1MolesField", "nMoles1");

    columnSpecies2MolesField_ =
        dict.lookupOrDefault<word>("columnSpecies2MolesField", "nMoles2");

    dict.lookup("topReservoirPatches")    >> topReservoirPatches_;
    dict.lookup("bottomReservoirPatches") >> bottomReservoirPatches_;

    packingLength_ =
        dict.lookupOrDefault<scalar>("packingLength", 0.05);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reservoirFenske::
reservoirFenske
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    speciesField_("cyclohexane.X"),
    A12Field_("A12"),
    columnTotalMolesField_("nMolesTotal"),
    columnSpecies1MolesField_("nMoles1"),
    columnSpecies2MolesField_("nMoles2"),
    topReservoirPatches_(),
    bottomReservoirPatches_(),
    packingLength_(0.05),
    initialized_(false),
    NSystem0_(0),
    N1System0_(0),
    N2System0_(0)
{
    readCoeffs(dict);
}


Foam::functionObjects::reservoirFenske::
~reservoirFenske()
{}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * //

Foam::wordList
Foam::functionObjects::reservoirFenske::fields() const
{
    return wordList
    {
        speciesField_,
        A12Field_,
        columnTotalMolesField_,
        columnSpecies1MolesField_,
        columnSpecies2MolesField_
    };
}


bool Foam::functionObjects::reservoirFenske::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    readCoeffs(dict);
    return true;
}


bool Foam::functionObjects::reservoirFenske::execute()
{
    // Species mole-fraction field must exist
    if (!mesh().foundObject<volScalarField>(speciesField_))
    {
        Info<< type() << " " << name()
            << ": field " << speciesField_ << " not found yet, skipping."
            << nl;

        return true;
    }

    const volScalarField& X =
        mesh().lookupObject<volScalarField>(speciesField_);

    // --- Reservoir states --------------------------------------------------

    scalar ND = 0;
    scalar xD = 0;
    const bool validTop =
        readCondenserReservoir(ND, xD, X, topReservoirPatches_);

    scalar NB = 0;
    scalar xB = 0;
    scalar yB = 0;
    const bool validBottom =
        readBoilerReservoir(NB, xB, yB, X, bottomReservoirPatches_);

    if (!validTop || !validBottom)
    {
        return true;
    }

    // --- Column molar inventory (volume integrals) -------------------------

    scalar NColumn  = 0;
    scalar N1Column = 0;
    scalar N2Column = 0;

    const bool validColumn =
        volumeIntegral(NColumn,  columnTotalMolesField_)
     && volumeIntegral(N1Column, columnSpecies1MolesField_)
     && volumeIntegral(N2Column, columnSpecies2MolesField_);

    if (!validColumn)
    {
        return true;
    }

    // --- Geometric mean A12 ------------------------------------------------

    scalar A12ColumnGeo = 0;
    const bool validA12 = volumeGeoMean(A12ColumnGeo, A12Field_);

    // --- System totals -----------------------------------------------------

    const scalar N1Top    = ND*xD;
    const scalar N2Top    = ND*(1 - xD);
    const scalar N1Bottom = NB*xB;
    const scalar N2Bottom = NB*(1 - xB);

    const scalar NSystem  = NColumn  + ND + NB;
    const scalar N1System = N1Column + N1Top + N1Bottom;
    const scalar N2System = N2Column + N2Top + N2Bottom;

    // Capture baseline at the first valid time step
    if (!initialized_ && NSystem > VSMALL)
    {
        initialized_ = true;
        NSystem0_    = NSystem;
        N1System0_   = N1System;
        N2System0_   = N2System;
    }

    const scalar dN  = initialized_ ? NSystem  - NSystem0_  : 0;
    const scalar dN1 = initialized_ ? N1System - N1System0_ : 0;
    const scalar dN2 = initialized_ ? N2System - N2System0_ : 0;

    const scalar relN  =
        initialized_ && mag(NSystem0_)  > VSMALL ? dN /NSystem0_  : 0;
    const scalar relN1 =
        initialized_ && mag(N1System0_) > VSMALL ? dN1/N1System0_ : 0;
    const scalar relN2 =
        initialized_ && mag(N2System0_) > VSMALL ? dN2/N2System0_ : 0;

    // --- Fenske / HETP -----------------------------------------------------
    // Nmin (Fenske) includes the boiler as one equilibrium stage.
    // Npacking = Nmin - 1 counts only the column stages.
    // HETP = packingLength / Npacking  (valid only when Nmin > 1).

    scalar NminColumnGeo  = 0;
    scalar HETPColumnGeo  = 0;

    const bool validNmin =
        validA12
     && validFenskeInputs(xD, xB, A12ColumnGeo);

    if (validNmin)
    {
        NminColumnGeo = fenskeNmin(xD, xB, A12ColumnGeo);
    }

    // HETP requires at least one column stage above the boiler (Nmin > 1).
    const bool validHETP =
        validNmin
     && NminColumnGeo > scalar(1) + VSMALL
     && packingLength_ > VSMALL;

    if (validHETP)
    {
        HETPColumnGeo = packingLength_ / (NminColumnGeo - scalar(1));
    }

    // --- Write (master only) -----------------------------------------------

    if (Pstream::master())
    {
        const fileName dir =
            mesh().time().path()
           /"postProcessing"
           /name()
           /"0";

        mkDir(dir);

        const fileName file = dir/"values.dat";
        const bool newFile  = !isFile(file);

        std::ofstream os(file.c_str(), std::ios::app);

        if (newFile)
        {
            os  << "# Time"
                << " ND xD"
                << " NB xB yB"
                << " NColumn N1Column N2Column"
                << " NSystem N1System N2System"
                << " dN dN1 dN2"
                << " relN relN1 relN2"
                << " A12ColumnGeo"
                << " NminColumnGeo"
                << " HETPColumnGeo"
                << std::endl;
        }

        os  << mesh().time().value() << " "
            << ND << " " << xD << " "
            << NB << " " << xB << " " << yB << " "
            << NColumn  << " " << N1Column  << " " << N2Column  << " "
            << NSystem  << " " << N1System  << " " << N2System  << " "
            << dN       << " " << dN1       << " " << dN2       << " "
            << relN     << " " << relN1     << " " << relN2     << " ";

        writeValueOrNa(os, validA12,  A12ColumnGeo);
        os << " ";
        writeValueOrNa(os, validNmin, NminColumnGeo);
        os << " ";
        writeValueOrNa(os, validHETP, HETPColumnGeo);
        os << std::endl;
    }

    Info<< type() << " " << name() << ": "
        << "xD=" << xD
        << ", xB=" << xB
        << ", yB=" << yB
        << ", Nsys=" << NSystem
        << ", dN=" << dN
        << ", dN1=" << dN1
        << ", dN2=" << dN2
        << ", relN=" << relN
        << ", A12geo=";

    if (validA12)  { Info << A12ColumnGeo;  }
    else           { Info << "nan";         }

    Info << ", Nmin=";

    if (validNmin) { Info << NminColumnGeo; }
    else           { Info << "nan";         }

    Info << nl;

    return true;
}


bool Foam::functionObjects::reservoirFenske::write()
{
    return true;
}

// ************************************************************************* //
