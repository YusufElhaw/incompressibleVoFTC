/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Derived solver module
     \\/     M anipulation  |
-------------------------------------------------------------------------------
This Document is part of the Master thesis of Yusuf Elhaw.

License
    This file is derived from OpenFOAM and released under the GNU General
    Public License, version 3 or later.
\*---------------------------------------------------------------------------*/

#include "incompressibleVoFCyclic.H"
#include "IOdictionary.H"
#include "localEulerDdtScheme.H"
#include "fvCorrectPhi.H"
#include "geometricZeroField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleVoFCyclic, 0);
    addToRunTimeSelectionTable(solver, incompressibleVoFCyclic, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleVoFCyclic::incompressibleVoFCyclic(fvMesh& mesh)
:
    incompressibleVoF(mesh),
    pPressureReferencePtr_(nullptr)
{
Info<< "Reading static pressure p for incompressibleVoFCyclic" << nl;

volScalarField pFromFile
(
    IOobject
    (
        "p",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh
);

p.primitiveFieldRef() = pFromFile.primitiveField();

forAll(p.boundaryField(), patchi)
{
    p.boundaryFieldRef().set
    (
        patchi,
        pFromFile.boundaryField()[patchi].clone(p).ptr()
    );
}

Info<< "p patch types:" << nl
    << p.boundaryField().types() << nl << endl;
    // The pressure equation now solves for p and uses pEqn.flux().
    mesh.schemes().setFluxRequired(p.name());

    // Direct pressure reference: p is the solved pressure, not p_rgh.
    pPressureReferencePtr_.reset(new Foam::pressureReference(p, pimple.dict()));

    // Keep p_rgh consistent for writing/post-processing and for any inherited
    // code path which only expects the field to exist.
    p_rgh_ = p - rho*buoyancy.gh;

    // The base incompressibleVoF constructor performs this correction with
    // p_rgh.  Re-apply it here with p so that moving/topology-changing cases
    // and corrected initial fluxes are consistent with the direct-p form.
    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }

    Info<< "Using direct static pressure p in momentum and pressure corrector"
        << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleVoFCyclic::~incompressibleVoFCyclic()
{}


// ************************************************************************* //
