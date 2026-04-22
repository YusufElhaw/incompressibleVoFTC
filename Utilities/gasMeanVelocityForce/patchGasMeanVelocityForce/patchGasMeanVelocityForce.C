/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/
#include "patchGasMeanVelocityForce.H"

#include "addToRunTimeSelectionTable.H"
#include "processorCyclicPolyPatch.H"
#include "volFields.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(patchGasMeanVelocityForce, 0);
    addToRunTimeSelectionTable(fvConstraint, patchGasMeanVelocityForce, dictionary);
}
}


void Foam::fv::patchGasMeanVelocityForce::readCoeffs(const dictionary& dict)
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

    forAll(patches_, i)
    {
        if (mesh().boundaryMesh().findIndex(patches_[i]) < 0)
        {
            FatalErrorInFunction
                << "Cannot find patch " << patches_[i]
                << exit(FatalError);
        }
    }
}


Foam::fv::patchGasMeanVelocityForce::patchGasMeanVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    gasMeanVelocityForce(sourceName, modelType, mesh, dict),
    patches_()
{
    readCoeffs(coeffs(dict));
}


Foam::scalar Foam::fv::patchGasMeanVelocityForce::magUbarGasAve
(
    const volVectorField& U,
    const volScalarField& alpha
) const
{
    /*const label patchi = mesh().boundaryMesh().findIndex(patch_);
    const vector direction = normalised(Ubar());

    scalar sumA = 0;
    scalar sumAjG = 0;

    // Main patch
    {
        const scalarField& magSf = mesh().boundary()[patchi].magSf();
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvPatchScalarField& alphaLp = alpha.boundaryField()[patchi];

        forAll(magSf, facei)
        {
            const scalar w = weight(alphaLp[facei]);   // alphaG on patch face
            sumA += magSf[facei];
            sumAjG += w*(direction & Up[facei])*magSf[facei];
        }
    }

    // If this is a cyclic patch in parallel, add processorCyclic pieces
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    if (Pstream::parRun() && isA<cyclicPolyPatch>(patches[patchi]))
    {
        const labelList processorCyclicPatches
        (
            processorCyclicPolyPatch::patchIDs(patch_, patches)
        );

        forAll(processorCyclicPatches, pcpi)
        {
            const label pPatchi = processorCyclicPatches[pcpi];

            const scalarField& magSf = mesh().boundary()[pPatchi].magSf();
            const fvPatchVectorField& Up = U.boundaryField()[pPatchi];
            const fvPatchScalarField& alphaLp = alpha.boundaryField()[pPatchi];

            forAll(magSf, facei)
            {
                const scalar w = weight(alphaLp[facei]);
                sumA += magSf[facei];
                sumAjG += w*(direction & Up[facei])*magSf[facei];
            }
        }
    }

    reduce(sumA, sumOp<scalar>());
    reduce(sumAjG, sumOp<scalar>());

    if (sumA <= SMALL)
    {
        WarningInFunction
            << "Patch area is too small for " << name()
            << ". Returning zero patch-averaged gas velocity." << endl;
        return 0;
    }

    return sumAjG/sumA;*/
    const vector direction = normalised(Ubar());

    scalar sumA = 0;
    scalar sumAjG = 0;

    forAll(patches_, patchNamei)
    {
        const word& patchName = patches_[patchNamei];
        const label patchi = mesh().boundaryMesh().findIndex(patchName);

        const scalarField& magSf = mesh().boundary()[patchi].magSf();
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvPatchScalarField& alphaLp = alpha.boundaryField()[patchi];

        forAll(magSf, facei)
        {
            const scalar w = weight(alphaLp[facei]);   // alphaG = 1 - alphaL
            sumA += magSf[facei];
            sumAjG += w*(direction & Up[facei])*magSf[facei];
        }
    }

    reduce(sumA, sumOp<scalar>());
    reduce(sumAjG, sumOp<scalar>());

    if (sumA <= SMALL)
    {
        WarningInFunction
            << "Total selected patch area is too small for " << name()
            << ". Returning zero patch-averaged gas velocity." << endl;
        return 0;
    }

    return sumAjG/sumA;
}


bool Foam::fv::patchGasMeanVelocityForce::read(const dictionary& dict)
{
    if (gasMeanVelocityForce::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
