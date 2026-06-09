/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
This Document is part of the Master thesis of Yusuf Elhaw.

\*---------------------------------------------------------------------------*/
#include "patchGasMeanVelocityForce.H"

#include "addToRunTimeSelectionTable.H"
#include "processorCyclicPolyPatch.H"
#include "volFields.H"
#include "DynamicList.H"

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
    const vector direction = normalised(Ubar());

    scalar sumA = 0;
    scalar sumAjG = 0;

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    DynamicList<label> measurementPatchIDs(patches_.size());

    auto appendPatchID =
    [&measurementPatchIDs](const label patchi)
    {
        if (patchi < 0)
        {
            return;
        }

        forAll(measurementPatchIDs, i)
        {
            if (measurementPatchIDs[i] == patchi)
            {
                return;
            }
        }

        measurementPatchIDs.append(patchi);
    };

    // 1. Collect the user-selected patches
    // 2. In parallel, also collect processorCyclic pieces belonging to them
    forAll(patches_, patchNamei)
    {
        const word& patchName = patches_[patchNamei];

        const label patchi = pbm.findIndex(patchName);

        if (patchi < 0)
        {
            FatalErrorInFunction
                << "Cannot find patch " << patchName
                << exit(FatalError);
        }

        appendPatchID(patchi);

        if (Pstream::parRun())
        {
            // Standard OpenFOAM helper for processor-cyclic pieces
            const labelList processorCyclicPatchIDs
            (
                processorCyclicPolyPatch::patchIDs(patchName, pbm)
            );

            forAll(processorCyclicPatchIDs, pcpi)
            {
                appendPatchID(processorCyclicPatchIDs[pcpi]);
            }

            // Fallback for NCC-generated names such as:
            // procBoundary0to10throughNCC_R1S1_on_PatchBR1S1
            forAll(pbm, pPatchi)
            {
                if
                (
                    isA<processorCyclicPolyPatch>(pbm[pPatchi])
                 && pbm[pPatchi].name().find("through") != string::npos
                 && pbm[pPatchi].name().find(patchName) != string::npos
                )
                {
                    appendPatchID(pPatchi);
                }
            }
        }
    }

    // Sum contributions over all selected physical/NCC/processorCyclic pieces
    forAll(measurementPatchIDs, i)
    {
        const label patchi = measurementPatchIDs[i];

        const scalarField& magSf = mesh().boundary()[patchi].magSf();
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvPatchScalarField& alphaLp = alpha.boundaryField()[patchi];

        forAll(magSf, facei)
        {
            const scalar w = weight(alphaLp[facei]); // gas weighting from alphaLiquid

            // Superficial gas velocity:
            // denominator remains total selected area, not gas-weighted area
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
