/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Utility: adjustContinuity
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "wordList.H"
#include "OFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

static bool containsName(const wordList& names, const word& name)
{
    forAll(names, i)
    {
        if (names[i] == name)
        {
            return true;
        }
    }
    return false;
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Symmetrically adjusts selected U patch normal velocities "
        "to satisfy global continuity."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary adjustDict
    (
        IOobject
        (
            "adjustContinuityDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word fieldName(adjustDict.lookupOrDefault<word>("field", "U"));
    word timeName(adjustDict.lookupOrDefault<word>("time", runTime.name()));
    wordList patchNames(adjustDict.lookup("patches"));

    const scalar tolerance = adjustDict.lookupOrDefault<scalar>("tolerance", 1e-8);
    const scalar split = adjustDict.lookupOrDefault<scalar>("split", 0.5);
    const Switch dryRun = adjustDict.lookupOrDefault<Switch>("dryRun", false);
    const Switch writeReport = adjustDict.lookupOrDefault<Switch>("writeReport", true);

    if (split < 0 || split > 1)
    {
        FatalErrorInFunction
            << "split must be between 0 and 1. Current value: " << split << nl
            << "split = 0.5 means half correction on outflow and half on inflow." << nl
            << exit(FatalError);
    }

    if (patchNames.size() == 0)
    {
        FatalErrorInFunction
            << "No patches listed in system/adjustContinuityDict" << nl
            << exit(FatalError);
    }

    Info<< nl << "Reading field " << fieldName << " from time " << timeName << nl << endl;

    volVectorField U
    (
        IOobject
        (
            fieldName,
            timeName,
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Selected correction patches:" << nl;
    forAll(patchNames, i)
    {
        const label patchi = mesh.boundaryMesh().findIndex(patchNames[i]);
        if (patchi < 0)
        {
            FatalErrorInFunction
                << "Patch " << patchNames[i] << " not found in mesh." << nl
                << "Available patches are:" << nl << mesh.boundaryMesh().names() << nl
                << exit(FatalError);
        }
        Info<< "  " << patchNames[i] << nl;
    }
    Info<< endl;

    scalar selectedIn = 0;
    scalar selectedOut = 0;
    scalar fixedIn = 0;
    scalar fixedOut = 0;
    scalar totalFluxAbs = 0;
    scalar netFlux = 0;

    const surfaceVectorField& Sf = mesh.Sf();

    forAll(mesh.boundaryMesh(), patchi)
    {
        const bool selected = containsName(patchNames, mesh.boundaryMesh()[patchi].name());
        const vectorField& Up = U.boundaryField()[patchi];
        const vectorField& Sfp = Sf.boundaryField()[patchi];

        forAll(Up, facei)
        {
            const scalar phi = (Up[facei] & Sfp[facei]);
            netFlux += phi;
            totalFluxAbs += mag(phi);

            if (phi < 0)
            {
                if (selected)
                {
                    selectedIn += -phi;
                }
                else
                {
                    fixedIn += -phi;
                }
            }
            else
            {
                if (selected)
                {
                    selectedOut += phi;
                }
                else
                {
                    fixedOut += phi;
                }
            }
        }
    }

    const scalar initialRelativeError =
        totalFluxAbs > vSmall ? mag(netFlux)/totalFluxAbs : 0;

    Info<< "Initial global boundary flux balance" << nl
        << "  selected inflow       : " << selectedIn << nl
        << "  selected outflow      : " << selectedOut << nl
        << "  non-selected inflow   : " << fixedIn << nl
        << "  non-selected outflow  : " << fixedOut << nl
        << "  net flux sum(phi)     : " << netFlux << nl
        << "  total flux sum(|phi|) : " << totalFluxAbs << nl
        << "  relative error        : " << initialRelativeError << nl
        << endl;

    if (initialRelativeError <= tolerance)
    {
        Info<< "Continuity already satisfies tolerance " << tolerance << ". No change." << nl << endl;
        return 0;
    }

    const scalar selectedFluxAbs = selectedIn + selectedOut;
    if (selectedFluxAbs <= vSmall)
    {
        FatalErrorInFunction
            << "Selected patches have no non-zero normal flux to correct." << nl
            << exit(FatalError);
    }

    scalar outFraction = split;
    scalar inFraction = 1 - split;

    if (selectedOut <= vSmall)
    {
        outFraction = 0;
        inFraction = 1;
        Info<< "Selected patches contain no outflow. Applying correction only to inflow." << nl;
    }
    else if (selectedIn <= vSmall)
    {
        outFraction = 1;
        inFraction = 0;
        Info<< "Selected patches contain no inflow. Applying correction only to outflow." << nl;
    }

    scalar outScale = 1;
    scalar inScale = 1;

    // Positive phi is outflow. Negative phi is inflow.
    // Need: netFlux + selectedOut*(outScale - 1) - selectedIn*(inScale - 1) = 0
    if (selectedOut > vSmall)
    {
        outScale = 1 - outFraction*netFlux/selectedOut;
    }

    if (selectedIn > vSmall)
    {
        inScale = 1 + inFraction*netFlux/selectedIn;
    }

    if (outScale < 0 || inScale < 0)
    {
        FatalErrorInFunction
            << "Correction would reverse the direction of at least one selected face group." << nl
            << "  outScale = " << outScale << nl
            << "  inScale  = " << inScale << nl
            << "Choose more correction patches or change split." << nl
            << exit(FatalError);
    }

    Info<< "Correction factors" << nl
        << "  split applied to outflow : " << outFraction << nl
        << "  split applied to inflow  : " << inFraction << nl
        << "  outflow scale            : " << outScale << nl
        << "  inflow scale             : " << inScale << nl
        << endl;

    label changedFaces = 0;

    forAll(patchNames, i)
    {
        const label patchi = mesh.boundaryMesh().findIndex(patchNames[i]);

        vectorField& Up = U.boundaryFieldRef()[patchi];
        const vectorField& Sfp = Sf.boundaryField()[patchi];

        forAll(Up, facei)
        {
            const scalar phi = (Up[facei] & Sfp[facei]);

            scalar newPhi = phi;
            if (phi > 0)
            {
                newPhi = outScale*phi;
            }
            else if (phi < 0)
            {
                newPhi = inScale*phi;
            }

            const scalar dPhi = newPhi - phi;
            if (mag(dPhi) > small)
            {
                const scalar magSf2 = magSqr(Sfp[facei]);
                if (magSf2 > vSmall)
                {
                    // Only the normal component is changed.
                    // Tangential velocity is preserved.
                    Up[facei] += (dPhi/magSf2)*Sfp[facei];
                    ++changedFaces;
                }
            }
        }
    }

    scalar finalNetFlux = 0;
    scalar finalTotalFluxAbs = 0;
    scalar finalSelectedIn = 0;
    scalar finalSelectedOut = 0;
    scalar finalFixedIn = 0;
    scalar finalFixedOut = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const bool selected = containsName(patchNames, mesh.boundaryMesh()[patchi].name());
        const vectorField& Up = U.boundaryField()[patchi];
        const vectorField& Sfp = Sf.boundaryField()[patchi];

        forAll(Up, facei)
        {
            const scalar phi = (Up[facei] & Sfp[facei]);
            finalNetFlux += phi;
            finalTotalFluxAbs += mag(phi);

            if (phi < 0)
            {
                if (selected)
                {
                    finalSelectedIn += -phi;
                }
                else
                {
                    finalFixedIn += -phi;
                }
            }
            else
            {
                if (selected)
                {
                    finalSelectedOut += phi;
                }
                else
                {
                    finalFixedOut += phi;
                }
            }
        }
    }

    const scalar finalRelativeError =
        finalTotalFluxAbs > vSmall ? mag(finalNetFlux)/finalTotalFluxAbs : 0;

    Info<< "Final global boundary flux balance" << nl
        << "  selected inflow       : " << finalSelectedIn << nl
        << "  selected outflow      : " << finalSelectedOut << nl
        << "  non-selected inflow   : " << finalFixedIn << nl
        << "  non-selected outflow  : " << finalFixedOut << nl
        << "  net flux sum(phi)     : " << finalNetFlux << nl
        << "  total flux sum(|phi|) : " << finalTotalFluxAbs << nl
        << "  relative error        : " << finalRelativeError << nl
        << "  changed faces         : " << changedFaces << nl
        << endl;

    if (finalRelativeError > tolerance)
    {
        FatalErrorInFunction
            << "Continuity still above tolerance after correction." << nl
            << "  tolerance      : " << tolerance << nl
            << "  relative error : " << finalRelativeError << nl
            << exit(FatalError);
    }

    if (writeReport)
    {
        mkDir(runTime.path()/"postProcessing"/"adjustContinuity"/timeName);
        OFstream report(runTime.path()/"postProcessing"/"adjustContinuity"/timeName/"report.dat");

        report
            << "# adjustContinuity report\n"
            << "field " << fieldName << "\n"
            << "time " << timeName << "\n"
            << "patches " << patchNames << "\n"
            << "tolerance " << tolerance << "\n"
            << "initialNetFlux " << netFlux << "\n"
            << "initialTotalFluxAbs " << totalFluxAbs << "\n"
            << "initialRelativeError " << initialRelativeError << "\n"
            << "selectedIn " << selectedIn << "\n"
            << "selectedOut " << selectedOut << "\n"
            << "nonSelectedIn " << fixedIn << "\n"
            << "nonSelectedOut " << fixedOut << "\n"
            << "outScale " << outScale << "\n"
            << "inScale " << inScale << "\n"
            << "finalNetFlux " << finalNetFlux << "\n"
            << "finalTotalFluxAbs " << finalTotalFluxAbs << "\n"
            << "finalRelativeError " << finalRelativeError << "\n"
            << "changedFaces " << changedFaces << "\n";
    }

    if (dryRun)
    {
        Info<< "dryRun = true: not writing " << fieldName << "." << nl << endl;
    }
    else
    {
        Info<< "Writing corrected " << fieldName << " to time " << timeName << nl << endl;
        U.write();
    }

    Info<< "End" << nl << endl;
    return 0;
}

// ************************************************************************* //
