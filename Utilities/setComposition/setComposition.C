/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     | Website:  https://openfoam.org
   \\  /    A nd           | Copyright (C) 2026
    \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    setComposition

Description
    Initialise two species mole-fraction fields and four species/phase mass-
    fraction fields for a two-phase/two-species incompressibleVoFTC case.

    The utility reads system/setCompositionDict, reads phase names from
    constant/phaseProperties, detects the liquid/phase1 field from the existing
    alpha.<phase> file, reads the two species, defaultSpecie entries and
    molWeight values from constant/physicalProperties.<phase>, applies the
    dummy-value logic used by compositionPredictor, converts mole fractions to
    mass fractions and writes:

        species1.X
        species2.X
        species1.phase1
        species2.phase1
        species1.phase2
        species2.phase2

    In interface cells, species1.X/species2.X are weighted with alpha*C,
    where C is the phase molar concentration. For incompressible phases C is
    computed from rho and mixture molWeight. For perfectGas phases C is
    computed cell-wise from p/(R*T) when p/T fields are available.

\*---------------------------------------------------------------------------*/

#include "setComposition.H"

#include "timeSelector.H"
#include "systemDict.H"
#include "Pstream.H"
#include "OSspecific.H"
#include "DynamicList.H"

#include <cctype>
#include <cstddef>

using namespace Foam;
using namespace Foam::setCompositionTools;

namespace Foam
{
namespace setCompositionTools
{

word readWordEntry(const dictionary& dict, const word& key)
{
    ITstream& is = dict.lookup(key);
    token t(is);

    if (t.isWord())
    {
        return t.wordToken();
    }

    if (t.isString())
    {
        return word(t.stringToken());
    }

    FatalIOErrorInFunction(is)
        << "Entry '" << key << "' must be a word or a quoted string."
        << nl << exit(FatalIOError);

    return word::null;
}

bool sameWordNoCase(const word& a, const word& b)
{
    if (a.size() != b.size())
    {
        return false;
    }

    for (std::size_t i = 0; i < a.size(); ++i)
    {
        const unsigned char ca = static_cast<unsigned char>(a[i]);
        const unsigned char cb = static_cast<unsigned char>(b[i]);

        if (std::tolower(ca) != std::tolower(cb))
        {
            return false;
        }
    }

    return true;
}

bool findEntryNameNoCase(const dictionary& dict, const word& wanted, word& actual)
{
    if (dict.found(wanted))
    {
        actual = wanted;
        return true;
    }

    const wordList toc(dict.toc());

    forAll(toc, i)
    {
        if (sameWordNoCase(toc[i], wanted))
        {
            actual = toc[i];
            return true;
        }
    }

    return false;
}

bool containsWord(const wordList& words, const word& w)
{
    forAll(words, i)
    {
        if (words[i] == w)
        {
            return true;
        }
    }

    return false;
}
label patchIDByName
(
    const fvMesh& mesh,
    const word& patchName
)
{
    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].name() == patchName)
        {
            return patchi;
        }
    }

    return -1;
}
void checkPatchListConflicts
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    wordList inletOutletPatches;
    wordList fixedValuePatches;
    wordList zeroGradientPatches;

    if (dict.found("inletOutletPatches"))
    {
        dict.lookup("inletOutletPatches") >> inletOutletPatches;
    }

    if (dict.found("fixedValuePatches"))
    {
        dict.lookup("fixedValuePatches") >> fixedValuePatches;
    }

    if (dict.found("zeroGradientPatches"))
    {
        dict.lookup("zeroGradientPatches") >> zeroGradientPatches;
    }

    // Check if all listed patches exist
    forAll(inletOutletPatches, i)
    {
        if (patchIDByName(mesh, inletOutletPatches[i]) < 0)
        {
            FatalErrorInFunction
                << "Patch '" << inletOutletPatches[i]
                << "' listed in inletOutletPatches does not exist in the mesh."
                << nl << exit(FatalError);
        }
    }

    forAll(fixedValuePatches, i)
    {
        if (patchIDByName(mesh, fixedValuePatches[i]) < 0)
        {
            FatalErrorInFunction
                << "Patch '" << fixedValuePatches[i]
                << "' listed in fixedValuePatches does not exist in the mesh."
                << nl << exit(FatalError);
        }
    }

    forAll(zeroGradientPatches, i)
    {
        if (patchIDByName(mesh, zeroGradientPatches[i]) < 0)
        {
            FatalErrorInFunction
                << "Patch '" << zeroGradientPatches[i]
                << "' listed in zeroGradientPatches does not exist in the mesh."
                << nl << exit(FatalError);
        }
    }

    // Check conflicts between inletOutletPatches and fixedValuePatches
    forAll(inletOutletPatches, i)
    {
        const word& p = inletOutletPatches[i];

        if (containsWord(fixedValuePatches, p))
        {
            FatalErrorInFunction
                << "Patch '" << p << "' is listed in both "
                << "inletOutletPatches and fixedValuePatches." << nl
                << "This is ambiguous. A patch may only appear in one "
                << "composition boundary-condition list." << nl
                << "Remove the patch from one of the lists in "
                << "system/setCompositionDict."
                << nl << exit(FatalError);
        }

        if (containsWord(zeroGradientPatches, p))
        {
            FatalErrorInFunction
                << "Patch '" << p << "' is listed in both "
                << "inletOutletPatches and zeroGradientPatches." << nl
                << "This is ambiguous. Remove the patch from one of the lists in "
                << "system/setCompositionDict."
                << nl << exit(FatalError);
        }
    }

    // Check conflicts between fixedValuePatches and zeroGradientPatches
    forAll(fixedValuePatches, i)
    {
        const word& p = fixedValuePatches[i];

        if (containsWord(zeroGradientPatches, p))
        {
            FatalErrorInFunction
                << "Patch '" << p << "' is listed in both "
                << "fixedValuePatches and zeroGradientPatches." << nl
                << "This is ambiguous. Remove the patch from one of the lists in "
                << "system/setCompositionDict."
                << nl << exit(FatalError);
        }
    }
}

bool containsDynamicWord(const DynamicList<word>& words, const word& w)
{
    forAll(words, i)
    {
        if (words[i] == w)
        {
            return true;
        }
    }

    return false;
}

void appendUnique(DynamicList<word>& words, const word& w)
{
    if (!containsDynamicWord(words, w))
    {
        words.append(w);
    }
}

wordList dynamicToWordList(const DynamicList<word>& words)
{
    wordList result(words.size());

    forAll(words, i)
    {
        result[i] = words[i];
    }

    return result;
}

bool containsWordNoCase(const wordList& words, const word& w, word& actual)
{
    forAll(words, i)
    {
        if (sameWordNoCase(words[i], w))
        {
            actual = words[i];
            return true;
        }
    }

    return false;
}

wordList readPhases(const Time& runTime, const fvMesh& mesh, const dictionary& dict)
{
    if (dict.found("phases"))
    {
        wordList phases;
        dict.lookup("phases") >> phases;

        if (phases.size() != 2)
        {
            FatalErrorInFunction
                << "Entry 'phases' in setCompositionDict must contain exactly two phases."
                << nl << exit(FatalError);
        }

        Info<< "Phases read from setCompositionDict/phases: " << phases << nl;
        return phases;
    }

    IOobject phaseHeader
    (
        "phaseProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    IOdictionary phaseProperties(phaseHeader);

    wordList phases;
    phaseProperties.lookup("phases") >> phases;

    if (phases.size() != 2)
    {
        FatalErrorInFunction
            << "constant/phaseProperties must contain exactly two phases, e.g.\n"
            << "    phases (liquid gas);"
            << nl << exit(FatalError);
    }

    Info<< "Phases read from constant/phaseProperties: " << phases << nl;
    return phases;
}

word findPhase1FromAlphaField(const Time& runTime, const fvMesh& mesh, const wordList& phases)
{
    label found = -1;
    label nFound = 0;

    forAll(phases, i)
    {
        const word alphaName = word("alpha.") + phases[i];

        IOobject alphaHeader
        (
            alphaName,
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );

        if (alphaHeader.headerOk())
        {
            if (found == -1)
            {
                found = i;
            }

            ++nFound;
        }
    }

    if (found == -1)
    {
        FatalErrorInFunction
            << "Could not find any alpha.<phase> field in time directory '"
            << runTime.name() << "' for phases " << phases << ".\n"
            << "Expected for example: " << runTime.name() << "/alpha."
            << phases[0] << nl << exit(FatalError);
    }

    if (nFound > 1)
    {
        Info<< "Warning: more than one alpha.<phase> field was found. "
            << "Using alpha." << phases[found] << " as phase1." << nl;
    }

    Info<< "phase1 detected from existing field alpha." << phases[found]
        << ": " << phases[found] << nl;

    return phases[found];
}

word physicalPropertiesObject
(
    const Time& runTime,
    const fvMesh& mesh,
    const word& phase,
    const bool required
)
{
    wordList candidates(4);
    candidates[0] = word("physicalProperties.") + phase;
    candidates[1] = word("thermophysicalProperties.") + phase;
    candidates[2] = "physicalProperties";
    candidates[3] = "thermophysicalProperties";

    forAll(candidates, i)
    {
        IOobject io
        (
            candidates[i],
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );

        if (io.headerOk())
        {
            return candidates[i];
        }
    }

    if (required)
    {
        FatalErrorInFunction
            << "Could not find constant/physicalProperties." << phase
            << " or constant/thermophysicalProperties." << phase << "."
            << nl << exit(FatalError);
    }

    return word::null;
}

wordList readSpeciesFromPhysicalProperties
(
    const Time& runTime,
    const fvMesh& mesh,
    const word& phase
)
{
    const word objectName = physicalPropertiesObject(runTime, mesh, phase, true);

    IOdictionary props
    (
        IOobject
        (
            objectName,
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (!props.found("species"))
    {
        FatalErrorInFunction
            << "Could not find entry 'species' in constant/" << objectName << "."
            << nl << exit(FatalError);
    }

    wordList species;
    props.lookup("species") >> species;

    if (species.size() != 2)
    {
        FatalErrorInFunction
            << "Entry 'species' in constant/" << objectName
            << " must contain exactly two species. Found: " << species
            << nl << exit(FatalError);
    }

    Info<< "Species read from constant/" << objectName << ": " << species << nl;
    return species;
}

word readDefaultSpecie
(
    const Time& runTime,
    const fvMesh& mesh,
    const word& phase,
    const wordList& species
)
{
    const word objectName = physicalPropertiesObject(runTime, mesh, phase, true);

    IOdictionary props
    (
        IOobject
        (
            objectName,
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word defaultSpecie = species[0];

    if (props.found("defaultSpecie"))
    {
        defaultSpecie = readWordEntry(props, "defaultSpecie");
    }
    else
    {
        Info<< "Warning: constant/" << objectName
            << " has no defaultSpecie. Using " << defaultSpecie
            << " for dummy values in phase " << phase << "." << nl;
    }

    word actualDefault;

    if (!containsWordNoCase(species, defaultSpecie, actualDefault))
    {
        FatalErrorInFunction
            << "defaultSpecie '" << defaultSpecie << "' in constant/"
            << objectName << " is not contained in species list " << species << "."
            << nl << exit(FatalError);
    }

    if (actualDefault != defaultSpecie)
    {
        Info<< "defaultSpecie '" << defaultSpecie << "' in constant/"
            << objectName << " matched species name '" << actualDefault
            << "' ignoring case." << nl;
    }

    Info<< "defaultSpecie for phase " << phase << " read from constant/"
        << objectName << ": " << actualDefault << nl;

    return actualDefault;
}

bool readMolWeightFromFile
(
    const Time& runTime,
    const fvMesh& mesh,
    const word& objectName,
    const word& species,
    scalar& W
)
{
    IOobject io
    (
        objectName,
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (!io.headerOk())
    {
        return false;
    }

    IOdictionary thermoDict(io);

    word actualSpecies;
    if (!findEntryNameNoCase(thermoDict, species, actualSpecies))
    {
        return false;
    }

    if (!thermoDict.isDict(actualSpecies))
    {
        return false;
    }

    const dictionary& spDict = thermoDict.subDict(actualSpecies);

    if (!spDict.isDict("specie"))
    {
        return false;
    }

    const dictionary& specieDict = spDict.subDict("specie");

    if (!specieDict.found("molWeight"))
    {
        return false;
    }

    W = readScalar(specieDict.lookup("molWeight"));
    return true;
}

bool readMolWeightFromDict
(
    const dictionary& dict,
    const word& species,
    scalar& W
)
{
    if (!dict.found("molecularWeights"))
    {
        return false;
    }

    const dictionary& mw = dict.subDict("molecularWeights");

    word actualSpecies;
    if (!findEntryNameNoCase(mw, species, actualSpecies))
    {
        return false;
    }

    W = readScalar(mw.lookup(actualSpecies));
    return true;
}

scalar readMolWeight
(
    const Time& runTime,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& species,
    const word& phase1,
    const word& phase2,
    fileName& source
)
{
    scalar W = -1;

    const word object1 = physicalPropertiesObject(runTime, mesh, phase1, false);
    const word object2 = physicalPropertiesObject(runTime, mesh, phase2, false);

    if (object1 != word::null && readMolWeightFromFile(runTime, mesh, object1, species, W))
    {
        source = fileName("constant")/object1;
        Info<< "molWeight(" << species << ") = " << W
            << " read from " << source << nl;
        return W;
    }

    if (object2 != word::null && readMolWeightFromFile(runTime, mesh, object2, species, W))
    {
        source = fileName("constant")/object2;
        Info<< "molWeight(" << species << ") = " << W
            << " read from " << source << nl;
        return W;
    }

    if (readMolWeightFromDict(dict, species, W))
    {
        source = fileName("system")/"setCompositionDict";
        Info<< "molWeight(" << species << ") = " << W
            << " read from " << source << nl;
        return W;
    }

    FatalErrorInFunction
        << "Could not find molWeight for species '" << species << "'.\n"
        << "The preferred location is:\n"
        << "    constant/physicalProperties." << phase1 << "\n"
        << "or\n"
        << "    constant/physicalProperties." << phase2 << "\n"
        << "inside:\n"
        << "    " << species << "\n"
        << "    {\n"
        << "        specie\n"
        << "        {\n"
        << "            molWeight <value>;\n"
        << "        }\n"
        << "    }\n"
        << "As fallback, you may define molecularWeights in system/setCompositionDict."
        << nl << exit(FatalError);

    return W;
}

scalar readOptionalSpeciesPhaseValue
(
    const dictionary& dict,
    const word& species,
    const word& phase,
    const bool required,
    bool& found
)
{
    found = false;

    word actualSpecies;
    if (!findEntryNameNoCase(dict, species, actualSpecies) || !dict.isDict(actualSpecies))
    {
        if (required)
        {
            FatalErrorInFunction
                << "Missing species sub-dictionary '" << species << "' in setCompositionDict.\n"
                << "The species names must correspond to the names read from "
                << "constant/physicalProperties.<phase>."
                << nl << exit(FatalError);
        }

        return 0;
    }

    const dictionary& sd = dict.subDict(actualSpecies);

    word actualPhase;
    if (!findEntryNameNoCase(sd, phase, actualPhase))
    {
        return 0;
    }

    found = true;

    if (actualSpecies != species || actualPhase != phase)
    {
        Info<< "Using setCompositionDict entry '" << actualSpecies << "."
            << actualPhase << "' for '" << species << "." << phase
            << "' ignoring case." << nl;
    }

    return readScalar(sd.lookup(actualPhase));
}

void normalisePhaseMoles
(
    const dictionary& dict,
    const word& species1,
    const word& species2,
    const word& phase,
    scalar& x1,
    scalar& x2
)
{
    bool has1 = false;
    bool has2 = false;

    x1 = readOptionalSpeciesPhaseValue(dict, species1, phase, false, has1);
    x2 = readOptionalSpeciesPhaseValue(dict, species2, phase, false, has2);

    if (!has1 && !has2)
    {
        FatalErrorInFunction
            << "Neither '" << species1 << "." << phase << "' nor '"
            << species2 << "." << phase << "' is defined in system/setCompositionDict.\n"
            << "Add for example:\n"
            << "    " << species1 << "\n"
            << "    {\n"
            << "        " << phase << " 0.5;\n"
            << "    }\n"
            << nl << exit(FatalError);
    }

    if (has1 && !has2)
    {
        x2 = 1 - x1;
        has2 = true;
    }
    else if (!has1 && has2)
    {
        x1 = 1 - x2;
        has1 = true;
    }

    if (x1 < -SMALL || x2 < -SMALL)
    {
        FatalErrorInFunction
            << "Negative mole fraction for phase '" << phase << "'.\n"
            << species1 << " = " << x1 << ", "
            << species2 << " = " << x2
            << nl << exit(FatalError);
    }

    x1 = max(x1, scalar(0));
    x2 = max(x2, scalar(0));

    const scalar sumX = x1 + x2;

    if (sumX <= VSMALL)
    {
        FatalErrorInFunction
            << "Mole-fraction sum for phase '" << phase << "' is zero."
            << nl << exit(FatalError);
    }

    if (mag(sumX - 1) > 1e-12)
    {
        Info<< "Normalising mole fractions for phase '" << phase << "': "
            << species1 << "=" << x1 << ", "
            << species2 << "=" << x2 << ", sum=" << sumX << nl;

        x1 /= sumX;
        x2 /= sumX;
    }
}

scalar moleToMassFraction
(
    const scalar x1,
    const scalar x2,
    const scalar W1,
    const scalar W2,
    const label speciesIndex
)
{
    const scalar denom = x1*W1 + x2*W2;

    if (denom <= VSMALL)
    {
        FatalErrorInFunction
            << "Invalid mole-to-mass conversion denominator."
            << nl << exit(FatalError);
    }

    if (speciesIndex == 0)
    {
        return (x1*W1)/denom;
    }

    return (x2*W2)/denom;
}

scalar mixtureMolWeight
(
    const scalar x1,
    const scalar x2,
    const scalar W1,
    const scalar W2
)
{
    const scalar Wmix = x1*W1 + x2*W2;

    if (Wmix <= VSMALL)
    {
        FatalErrorInFunction
            << "Invalid mixture molWeight."
            << nl << exit(FatalError);
    }

    return Wmix;
}

word readEquationOfState
(
    const Time& runTime,
    const fvMesh& mesh,
    const word& phase
)
{
    const word objectName = physicalPropertiesObject(runTime, mesh, phase, true);

    IOdictionary props
    (
        IOobject
        (
            objectName,
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word eos = word::null;

    if (props.isDict("thermoType"))
    {
        const dictionary& thermoType = props.subDict("thermoType");

        if (thermoType.found("equationOfState"))
        {
            eos = readWordEntry(thermoType, "equationOfState");
        }
    }

    if (eos == word::null && props.found("equationOfState"))
    {
        eos = readWordEntry(props, "equationOfState");
    }

    if (eos == word::null)
    {
        Info<< "No equationOfState entry found for phase " << phase
            << " in constant/" << objectName
            << ". Using rho/molWeight for molar concentration." << nl;

        return word("none");
    }

    Info<< "equationOfState(" << phase << ") = " << eos
        << " read from constant/" << objectName << nl;

    return eos;
}

bool readRhoFromPhaseProperties
(
    const Time& runTime,
    const fvMesh& mesh,
    const word& phase,
    scalar& rho,
    fileName& source
)
{
    const word objectName = physicalPropertiesObject(runTime, mesh, phase, false);

    if (objectName == word::null)
    {
        return false;
    }

    IOobject io
    (
        objectName,
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (!io.headerOk())
    {
        return false;
    }

    IOdictionary props(io);

    if (!props.found("rho"))
    {
        return false;
    }

    rho = readScalar(props.lookup("rho"));
    source = fileName("constant")/objectName;

    return true;
}

word findScalarFieldObject
(
    const Time& runTime,
    const fvMesh& mesh,
    const wordList& candidates,
    const bool required,
    const word& description
)
{
    forAll(candidates, i)
    {
        IOobject io
        (
            candidates[i],
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );

        if (io.headerOk())
        {
            Info<< description << " field read from "
                << runTime.name() << "/" << candidates[i] << nl;

            return candidates[i];
        }
    }

    if (required)
    {
        FatalErrorInFunction
            << "Could not find " << description << " field in time directory '"
            << runTime.name() << "'. Tried: " << candidates << nl
            << exit(FatalError);
    }

    return word::null;
}

wordList pressureFieldCandidates(const dictionary& dict)
{
    const word base = dict.lookupOrDefault<word>("pressureField", word("p_rgh"));

    DynamicList<word> names;
    appendUnique(names, base);
    appendUnique(names, word(base + ".orig"));
    appendUnique(names, word("p_rgh"));
    appendUnique(names, word("p_rgh.orig"));
    appendUnique(names, word("p"));
    appendUnique(names, word("p.orig"));

    return dynamicToWordList(names);
}

wordList temperatureFieldCandidates(const dictionary& dict)
{
    const word base = dict.lookupOrDefault<word>("temperatureField", word("T"));

    DynamicList<word> names;
    appendUnique(names, base);
    appendUnique(names, word(base + ".orig"));
    appendUnique(names, word("T"));
    appendUnique(names, word("T.orig"));

    return dynamicToWordList(names);
}

scalar pressureOffsetForPhase
(
    const dictionary& dict,
    const word& phase
)
{
    if (dict.isDict("pressureOffset"))
    {
        const dictionary& p0 = dict.subDict("pressureOffset");

        word actualPhase;
        if (findEntryNameNoCase(p0, phase, actualPhase))
        {
            return readScalar(p0.lookup(actualPhase));
        }

        return 0.0;
    }

    if (dict.found("pressureOffset"))
    {
        return readScalar(dict.lookup("pressureOffset"));
    }

    return 0.0;
}

void phaseMolarConcentration
(
    const Time& runTime,
    const fvMesh& mesh,
    const dictionary& dict,
    const compositionInput& c,
    const word& phase,
    const label nCells,
    scalarField& CInternal,
    PtrList<scalarField>& CBoundary
)
{
    CInternal.setSize(nCells);
    CBoundary.setSize(mesh.boundary().size());

    const bool isPhase1 = phase == c.phase1;

    const scalar x1 = isPhase1 ? c.x1Phase1 : c.x1Phase2;
    const scalar x2 = isPhase1 ? c.x2Phase1 : c.x2Phase2;
    const scalar Wmix = mixtureMolWeight(x1, x2, c.W1, c.W2);

    const word eos = readEquationOfState(runTime, mesh, phase);
    const bool perfectGas = sameWordNoCase(eos, word("perfectGas"));

    if (perfectGas)
    {
        const word pName = findScalarFieldObject
        (
            runTime,
            mesh,
            pressureFieldCandidates(dict),
            false,
            word("pressure")
        );

        const word TName = findScalarFieldObject
        (
            runTime,
            mesh,
            temperatureFieldCandidates(dict),
            false,
            word("temperature")
        );

        if (pName != word::null && TName != word::null)
        {
            volScalarField p
            (
                IOobject
                (
                    pName,
                    runTime.name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            volScalarField T
            (
                IOobject
                (
                    TName,
                    runTime.name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            const scalar pOffset = pressureOffsetForPhase(dict, phase);
            const scalar R = dict.lookupOrDefault<scalar>
            (
                "universalGasConstant",
                8314.46261815324
            );

            const scalarField& pIf = p.primitiveField();
            const scalarField& TIf = T.primitiveField();

            forAll(CInternal, celli)
            {
                const scalar pAbs = pIf[celli] + pOffset;
                const scalar Ti = TIf[celli];

                if (Ti <= VSMALL)
                {
                    FatalErrorInFunction
                        << "Temperature for perfectGas phase '" << phase
                        << "' is not positive in cell " << celli
                        << ": T=" << Ti << nl << exit(FatalError);
                }

                if (pAbs <= VSMALL)
                {
                    FatalErrorInFunction
                        << "Pressure for perfectGas phase '" << phase
                        << "' is not positive in cell " << celli
                        << ": pField=" << pIf[celli]
                        << ", pressureOffset=" << pOffset
                        << ", pAbs=" << pAbs << nl
                        << "If your p_rgh is relative/gauge pressure, add e.g." << nl
                        << "    pressureOffset 101325;" << nl
                        << "to system/setCompositionDict."
                        << nl << exit(FatalError);
                }

                CInternal[celli] = pAbs/(R*Ti);
            }

            forAll(mesh.boundary(), patchi)
            {
                const fvPatchScalarField& pPf = p.boundaryField()[patchi];
                const fvPatchScalarField& TPf = T.boundaryField()[patchi];

                tmp<scalarField> tpInternal = pPf.patchInternalField();
                tmp<scalarField> tTInternal = TPf.patchInternalField();

                const scalarField& pInternalPatch = tpInternal();
                const scalarField& TInternalPatch = tTInternal();

                const word pPatchType = pPf.type();
                const word TPatchType = TPf.type();

                const bool useInternalP =
                (
                    pPatchType == "fixedFluxPressure"
                || pPatchType == "fixedFluxExtrapolatedPressure"
                || pPatchType == "zeroGradient"
                || pPatchType == "calculated"
                );

                const bool useInternalT =
                (
                    TPatchType == "zeroGradient"
                || TPatchType == "calculated"
                );

                CBoundary.set(patchi, new scalarField(pPf.size(), 0.0));
                scalarField& CPf = CBoundary[patchi];

                forAll(CPf, facei)
                {
                    scalar pRaw = useInternalP ? pInternalPatch[facei] : pPf[facei];
                    scalar Ti   = useInternalT ? TInternalPatch[facei] : TPf[facei];

                    // Safety fallback:
                    // Some pressure boundary conditions store a dummy value such as 0.
                    // If the patch value is unusable but the adjacent cell value is valid,
                    // use the adjacent cell value.
                    if ((pRaw + pOffset) <= VSMALL && (pInternalPatch[facei] + pOffset) > VSMALL)
                    {
                        pRaw = pInternalPatch[facei];
                    }

                    if (Ti <= VSMALL && TInternalPatch[facei] > VSMALL)
                    {
                        Ti = TInternalPatch[facei];
                    }

                    const scalar pAbs = pRaw + pOffset;

                    if (Ti <= VSMALL)
                    {
                        FatalErrorInFunction
                            << "Patch temperature for perfectGas phase '" << phase
                            << "' is not positive on patch '"
                            << mesh.boundary()[patchi].name()
                            << "' face " << facei
                            << ": T=" << Ti
                            << ", TInternal=" << TInternalPatch[facei]
                            << ", TPatchType=" << TPatchType
                            << nl << exit(FatalError);
                    }

                    if (pAbs <= VSMALL)
                    {
                        FatalErrorInFunction
                            << "Patch pressure for perfectGas phase '" << phase
                            << "' is not positive on patch '"
                            << mesh.boundary()[patchi].name()
                            << "' face " << facei
                            << ": pPatch=" << pPf[facei]
                            << ", pInternal=" << pInternalPatch[facei]
                            << ", pUsed=" << pRaw
                            << ", pressureOffset=" << pOffset
                            << ", pAbs=" << pAbs
                            << ", pPatchType=" << pPatchType << nl
                            << "If both patch and internal pressure are relative/gauge pressure, "
                            << "add for example:" << nl
                            << "    pressureOffset 101325;" << nl
                            << "to system/setCompositionDict."
                            << nl << exit(FatalError);
                    }

                    CPf[facei] = pAbs/(R*Ti);
                }
            }

            Info<< "Molar concentration C(" << phase
                << ") calculated cell-wise and patch-wise from perfectGas: C=p/(R*T), "
                << "pField=" << pName << ", TField=" << TName
                << ", pressureOffset=" << pOffset
                << ", R=" << R << nl;

            return;
        }

        Info<< "Warning: phase " << phase
            << " uses equationOfState perfectGas, but p/T fields were not both found. "
            << "Falling back to rho/molWeight if rho is available." << nl;
    }

    scalar rho = -1;
    fileName rhoSource;

    if (!readRhoFromPhaseProperties(runTime, mesh, phase, rho, rhoSource))
    {
        FatalErrorInFunction
            << "Could not calculate molar concentration C for phase '" << phase << "'." << nl
            << "For incompressible phases, define rho in constant/physicalProperties."
            << phase << "." << nl
            << "For perfectGas phases, provide p/T fields, e.g. 0/p_rgh and 0/T "
            << "or 0/p_rgh.orig and 0/T.orig."
            << nl << exit(FatalError);
    }

    if (rho <= VSMALL)
    {
        FatalErrorInFunction
            << "rho for phase '" << phase << "' must be positive. Found rho="
            << rho << " in " << rhoSource << nl << exit(FatalError);
    }

    const scalar C0 = rho/Wmix;

    forAll(CInternal, celli)
    {
        CInternal[celli] = C0;
    }

    forAll(mesh.boundary(), patchi)
    {
        CBoundary.set
        (
            patchi,
            new scalarField(mesh.boundary()[patchi].size(), C0)
        );
    }

    Info<< "rho(" << phase << ") = " << rho
        << " read from " << rhoSource << nl
        << "Mmix(" << phase << ") = " << Wmix
        << " from mole fractions and molWeights" << nl
        << "Molar concentration C(" << phase << ") = rho/Mmix = "
        << C0 << " used for cells and patches" << nl;
}

void calculateXMoleFractionBoundaryValues
(
    const fvMesh& mesh,
    const volScalarField& alpha,
    const scalar alphaTol,
    const compositionInput& c,
    const PtrList<scalarField>& CPhase1Boundary,
    const PtrList<scalarField>& CPhase2Boundary,
    PtrList<scalarField>& X1Boundary,
    PtrList<scalarField>& X2Boundary
)
{
    X1Boundary.setSize(mesh.boundary().size());
    X2Boundary.setSize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        const fvPatchScalarField& aPatch = alpha.boundaryField()[patchi];
        const scalarField& C1Patch = CPhase1Boundary[patchi];
        const scalarField& C2Patch = CPhase2Boundary[patchi];

        X1Boundary.set(patchi, new scalarField(aPatch.size(), 0.0));
        X2Boundary.set(patchi, new scalarField(aPatch.size(), 0.0));

        scalarField& X1Patch = X1Boundary[patchi];
        scalarField& X2Patch = X2Boundary[patchi];

        forAll(aPatch, facei)
        {
            const scalar ac = min(max(aPatch[facei], scalar(0)), scalar(1));

            if (ac >= 1 - alphaTol)
            {
                X1Patch[facei] = c.x1Phase1;
                X2Patch[facei] = c.x2Phase1;
            }
            else if (ac <= alphaTol)
            {
                X1Patch[facei] = c.x1Phase2;
                X2Patch[facei] = c.x2Phase2;
            }
            else
            {
                const scalar w1 = ac*C1Patch[facei];
                const scalar w2 = (1 - ac)*C2Patch[facei];
                const scalar wSum = w1 + w2;

                if (wSum <= VSMALL)
                {
                    FatalErrorInFunction
                        << "Invalid patch alpha*C denominator on patch '"
                        << mesh.boundary()[patchi].name()
                        << "' face " << facei
                        << ": alpha=" << ac
                        << ", C(" << c.phase1 << ")=" << C1Patch[facei]
                        << ", C(" << c.phase2 << ")=" << C2Patch[facei]
                        << nl << exit(FatalError);
                }

                X1Patch[facei] = (w1*c.x1Phase1 + w2*c.x1Phase2)/wSum;
                X2Patch[facei] = 1 - X1Patch[facei];
            }
        }
    }

    Info<< "Patch values for " << c.species1 << ".X and " << c.species2
        << ".X calculated with (alpha*C*x)_sum/(alpha*C)_sum. "
        << "For perfectGas phases, patch p/T values are used in C=p/(R*T)." << nl;
}

compositionInput readCompositionInput
(
    const Time& runTime,
    const fvMesh& mesh,
    const dictionary& dict
)
{
    compositionInput c;

    const wordList allPhases = readPhases(runTime, mesh, dict);

    c.phase1 = findPhase1FromAlphaField(runTime, mesh, allPhases);

    if (c.phase1 == allPhases[0])
    {
        c.phase2 = allPhases[1];
    }
    else if (c.phase1 == allPhases[1])
    {
        c.phase2 = allPhases[0];
    }
    else
    {
        FatalErrorInFunction
            << "Internal error while ordering phases."
            << nl << exit(FatalError);
    }

    Info<< "phase2 selected from phaseProperties as the other phase: "
        << c.phase2 << nl;

    const wordList species = readSpeciesFromPhysicalProperties(runTime, mesh, c.phase1);

    c.species1 = species[0];
    c.species2 = species[1];

    c.defaultSpeciePhase1 = readDefaultSpecie(runTime, mesh, c.phase1, species);
    c.defaultSpeciePhase2 = readDefaultSpecie(runTime, mesh, c.phase2, species);

    c.W1 = readMolWeight(runTime, mesh, dict, c.species1, c.phase1, c.phase2, c.W1Source);
    c.W2 = readMolWeight(runTime, mesh, dict, c.species2, c.phase1, c.phase2, c.W2Source);

    normalisePhaseMoles(dict, c.species1, c.species2, c.phase1, c.x1Phase1, c.x2Phase1);
    normalisePhaseMoles(dict, c.species1, c.species2, c.phase2, c.x1Phase2, c.x2Phase2);

    c.y1Phase1 = moleToMassFraction(c.x1Phase1, c.x2Phase1, c.W1, c.W2, 0);
    c.y2Phase1 = moleToMassFraction(c.x1Phase1, c.x2Phase1, c.W1, c.W2, 1);
    c.y1Phase2 = moleToMassFraction(c.x1Phase2, c.x2Phase2, c.W1, c.W2, 0);
    c.y2Phase2 = moleToMassFraction(c.x1Phase2, c.x2Phase2, c.W1, c.W2, 1);

    return c;
}

bool isConstraintPatchType(const word& meshPatchType, const word& fieldPatchType)
{
    return
    (
        meshPatchType == "empty"
     || meshPatchType == "cyclic"
     || meshPatchType == "cyclicAMI"
     || meshPatchType == "processor"
     || meshPatchType == "processorCyclic"
     || meshPatchType == "symmetryPlane"
     || meshPatchType == "symmetry"
     || meshPatchType == "wedge"
     || fieldPatchType == "empty"
     || fieldPatchType == "cyclic"
     || fieldPatchType == "cyclicAMI"
     || fieldPatchType == "processor"
     || fieldPatchType == "processorCyclic"
     || fieldPatchType == "symmetryPlane"
     || fieldPatchType == "symmetry"
     || fieldPatchType == "wedge"
    );
}

scalar phaseValue
(
    const compositionInput& c,
    const word& species,
    const word& phase,
    const bool massFraction
)
{
    if (massFraction)
    {
        if (species == c.species1 && phase == c.phase1) return c.y1Phase1;
        if (species == c.species2 && phase == c.phase1) return c.y2Phase1;
        if (species == c.species1 && phase == c.phase2) return c.y1Phase2;
        if (species == c.species2 && phase == c.phase2) return c.y2Phase2;
    }
    else
    {
        if (species == c.species1 && phase == c.phase1) return c.x1Phase1;
        if (species == c.species2 && phase == c.phase1) return c.x2Phase1;
        if (species == c.species1 && phase == c.phase2) return c.x1Phase2;
        if (species == c.species2 && phase == c.phase2) return c.x2Phase2;
    }

    FatalErrorInFunction
        << "Unknown species/phase combination: species='" << species
        << "', phase='" << phase << "'."
        << nl << exit(FatalError);

    return 0;
}

scalar dummyValue(const compositionInput& c, const word& species, const word& phase)
{
    if (phase == c.phase1)
    {
        return species == c.defaultSpeciePhase1 ? 1.0 : 0.0;
    }

    if (phase == c.phase2)
    {
        return species == c.defaultSpeciePhase2 ? 1.0 : 0.0;
    }

    FatalErrorInFunction
        << "Unknown phase '" << phase << "' for dummy value."
        << nl << exit(FatalError);

    return 0;
}

word inletPhaseForPatch
(
    const compositionInput& c,
    const dictionary& dict,
    const word& patchName
)
{
    word defaultPhase = c.phase1;

    if (dict.found("defaultInletPhase"))
    {
        defaultPhase = readWordEntry(dict, "defaultInletPhase");
    }

    word inletPhase = defaultPhase;

    if (dict.isDict("patchInletPhases"))
    {
        const dictionary& p = dict.subDict("patchInletPhases");

        word actualPatch;
        if (findEntryNameNoCase(p, patchName, actualPatch))
        {
            ITstream& is = p.lookup(actualPatch);
            token t(is);

            if (t.isWord())
            {
                inletPhase = t.wordToken();
            }
            else if (t.isString())
            {
                inletPhase = word(t.stringToken());
            }
            else
            {
                FatalIOErrorInFunction(is)
                    << "patchInletPhases entry for patch '" << patchName
                    << "' must be a word."
                    << nl << exit(FatalIOError);
            }
        }
    }

    if (sameWordNoCase(inletPhase, c.phase1))
    {
        return c.phase1;
    }

    if (sameWordNoCase(inletPhase, c.phase2))
    {
        return c.phase2;
    }

    FatalErrorInFunction
        << "Patch '" << patchName << "' uses inlet phase '" << inletPhase
        << "', but valid phases are '" << c.phase1 << "' and '"
        << c.phase2 << "'."
        << nl << exit(FatalError);

    return c.phase1;
}

scalar inletValueForPatch
(
    const compositionInput& c,
    const dictionary& dict,
    const word& patchName,
    const word& fieldSpecies,
    const word& fieldPhase,
    const bool massFraction
)
{
    if (massFraction)
    {
        return phaseValue(c, fieldSpecies, fieldPhase, true);
    }

    const word inletPhase = inletPhaseForPatch(c, dict, patchName);
    return phaseValue(c, fieldSpecies, inletPhase, false);
}

void writeFieldHeader(OFstream& os, const word& fieldName)
{
    os  << "/*--------------------------------*- C++ -*----------------------------------*\\" << nl
        << "| =========                 |                                                 |" << nl
        << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |" << nl
        << "|  \\    /   O peration     | Generated by setComposition                    |" << nl
        << "|   \\  /    A nd           |                                                 |" << nl
        << "|    \\/     M anipulation  |                                                 |" << nl
        << "\\*---------------------------------------------------------------------------*/" << nl
        << "FoamFile" << nl
        << "{" << nl
        << "    version     2.0;" << nl
        << "    format      ascii;" << nl
        << "    class       volScalarField;" << nl
        << "    object      " << fieldName << ";" << nl
        << "}" << nl
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
        << nl << nl;
}

word patchScalarType
(
    const fvMesh& mesh,
    const volScalarField& alpha,
    const label patchi,
    const dictionary& dict
)
{
    const word& patchName = mesh.boundary()[patchi].name();
    const word meshPatchType = mesh.boundaryMesh()[patchi].type();
    const word fieldPatchType = alpha.boundaryField()[patchi].type();

    // Constraint patches must keep their type, e.g. empty, cyclic, symmetry
    if (isConstraintPatchType(meshPatchType, fieldPatchType))
    {
        return fieldPatchType;
    }

    wordList inletOutletPatches;
    if (dict.found("inletOutletPatches"))
    {
        dict.lookup("inletOutletPatches") >> inletOutletPatches;

        if (containsWord(inletOutletPatches, patchName))
        {
            return "inletOutlet";
        }
    }

    wordList fixedValuePatches;
    if (dict.found("fixedValuePatches"))
    {
        dict.lookup("fixedValuePatches") >> fixedValuePatches;

        if (containsWord(fixedValuePatches, patchName))
        {
            return "fixedValue";
        }
    }

    wordList zeroGradientPatches;
    if (dict.found("zeroGradientPatches"))
    {
        dict.lookup("zeroGradientPatches") >> zeroGradientPatches;

        if (containsWord(zeroGradientPatches, patchName))
        {
            return "zeroGradient";
        }
    }

    if (dict.found("defaultPatchType"))
    {
        return readWordEntry(dict, "defaultPatchType");
    }

    return "zeroGradient";
}

void writeScalarFieldFile
(
    const fvMesh& mesh,
    const Time& runTime,
    const volScalarField& alpha,
    const dictionary& dict,
    const word& fieldName,
    const scalarField& internalValues,
    const word& fieldSpecies,
    const word& fieldPhase,
    const bool massFraction,
    const compositionInput& c,
    const PtrList<scalarField>* boundaryValuesPtr
)
{
    IOobject io
    (
        fieldName,
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    mkDir(io.objectPath(false).path());

    OFstream os(io.objectPath(false));

    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open field file for writing: " << io.objectPath(false)
            << nl << exit(FatalError);
    }

    os.precision(12);

    writeFieldHeader(os, fieldName);

    os << "dimensions      [0 0 0 0 0 0 0];" << nl << nl;

    os  << "internalField   nonuniform List<scalar>" << nl
        << internalValues.size() << nl
        << "(" << nl;

    forAll(internalValues, i)
    {
        os << internalValues[i] << nl;
    }

    os  << ")" << nl
        << ";" << nl << nl;

    os  << "boundaryField" << nl
        << "{" << nl;

    const word phiName = dict.lookupOrDefault<word>("phi", word("phi"));

    forAll(mesh.boundary(), patchi)
    {
        const word& patchName = mesh.boundary()[patchi].name();
        const word patchType = patchScalarType(mesh, alpha, patchi, dict);

        os  << "    " << patchName << nl
            << "    {" << nl
            << "        type        " << patchType << ";" << nl;

        if (patchType == "inletOutlet")
        {
            const scalar inletVal = inletValueForPatch
            (
                c,
                dict,
                patchName,
                fieldSpecies,
                fieldPhase,
                massFraction
            );

            os  << "        phi         " << phiName << ";" << nl
                << "        inletValue  uniform " << inletVal << ";" << nl;

            if (boundaryValuesPtr)
            {
                const scalarField& patchValues = (*boundaryValuesPtr)[patchi];

                os  << "        value       nonuniform List<scalar>" << nl
                    << patchValues.size() << nl
                    << "(" << nl;

                forAll(patchValues, facei)
                {
                    os << patchValues[facei] << nl;
                }

                os  << ")" << nl
                    << ";" << nl;
            }
            else
            {
                os  << "        value       uniform " << inletVal << ";" << nl;
            }
        }
        else if (patchType == "fixedValue")
        {
            const scalar val = inletValueForPatch
            (
                c,
                dict,
                patchName,
                fieldSpecies,
                fieldPhase,
                massFraction
            );

            os << "        value       uniform " << val << ";" << nl;
        }

        os  << "    }" << nl;
    }

    os  << "}" << nl << nl
        << "// ************************************************************************* //" << nl;

    Info<< "Written field: " << io.objectPath(false) << nl;
}

} // End namespace setCompositionTools
} // End namespace Foam

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    timeSelector::select0(runTime, args);

    #include "createRegionMeshNoChangers.H"

    const dictionary setCompositionDict(systemDict("setCompositionDict", args, mesh));
    checkPatchListConflicts(mesh, setCompositionDict);

    const compositionInput c = readCompositionInput(runTime, mesh, setCompositionDict);

    const scalar alphaTol = setCompositionDict.lookupOrDefault<scalar>("alphaTol", 1e-4);
    const word alphaName = "alpha." + c.phase1;

    volScalarField alpha
    (
        IOobject
        (
            alphaName,
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< nl
        << "setComposition" << nl
        << "    phase1 from alpha field     : " << c.phase1 << nl
        << "    phase2                      : " << c.phase2 << nl
        << "    alpha field                 : " << alphaName << nl
        << "    species1                    : " << c.species1 << nl
        << "    species2                    : " << c.species2 << nl
        << "    defaultSpecie(" << c.phase1 << ")        : " << c.defaultSpeciePhase1 << nl
        << "    defaultSpecie(" << c.phase2 << ")        : " << c.defaultSpeciePhase2 << nl
        << "    molWeight(" << c.species1 << ")          : " << c.W1
        << "  [" << c.W1Source << "]" << nl
        << "    molWeight(" << c.species2 << ")          : " << c.W2
        << "  [" << c.W2Source << "]" << nl
        << "    x(" << c.species1 << "," << c.phase1 << ")       : " << c.x1Phase1 << nl
        << "    x(" << c.species2 << "," << c.phase1 << ")       : " << c.x2Phase1 << nl
        << "    x(" << c.species1 << "," << c.phase2 << ")       : " << c.x1Phase2 << nl
        << "    x(" << c.species2 << "," << c.phase2 << ")       : " << c.x2Phase2 << nl
        << "    Y(" << c.species1 << "," << c.phase1 << ")       : " << c.y1Phase1 << nl
        << "    Y(" << c.species2 << "," << c.phase1 << ")       : " << c.y2Phase1 << nl
        << "    Y(" << c.species1 << "," << c.phase2 << ")       : " << c.y1Phase2 << nl
        << "    Y(" << c.species2 << "," << c.phase2 << ")       : " << c.y2Phase2 << nl
        << "    dummy in missing " << c.phase1 << "        : "
        << c.defaultSpeciePhase1 << "=1, other=0" << nl
        << "    dummy in missing " << c.phase2 << "        : "
        << c.defaultSpeciePhase2 << "=1, other=0" << nl
        << "    alphaTol                   : " << alphaTol << nl
        << "    mixed-cell X weighting     : (alpha*C*x)_sum/(alpha*C)_sum" << nl
        << endl;

    const scalarField& a = alpha.primitiveField();

    scalarField CPhase1(a.size(), 0.0);
    scalarField CPhase2(a.size(), 0.0);
    PtrList<scalarField> CPhase1Boundary(mesh.boundary().size());
    PtrList<scalarField> CPhase2Boundary(mesh.boundary().size());

    phaseMolarConcentration
    (
        runTime,
        mesh,
        setCompositionDict,
        c,
        c.phase1,
        a.size(),
        CPhase1,
        CPhase1Boundary
    );

    phaseMolarConcentration
    (
        runTime,
        mesh,
        setCompositionDict,
        c,
        c.phase2,
        a.size(),
        CPhase2,
        CPhase2Boundary
    );

    scalarField X1(a.size(), 0.0);
    scalarField X2(a.size(), 0.0);
    PtrList<scalarField> X1Boundary(mesh.boundary().size());
    PtrList<scalarField> X2Boundary(mesh.boundary().size());

    scalarField Y1Phase1(a.size(), 0.0);
    scalarField Y2Phase1(a.size(), 0.0);
    scalarField Y1Phase2(a.size(), 0.0);
    scalarField Y2Phase2(a.size(), 0.0);

    forAll(a, celli)
    {
        const scalar ac = min(max(a[celli], scalar(0)), scalar(1));

        if (ac >= 1 - alphaTol)
        {
            // Pure phase1: phase2 fields are dummy values.
            X1[celli] = c.x1Phase1;
            X2[celli] = c.x2Phase1;

            Y1Phase1[celli] = c.y1Phase1;
            Y2Phase1[celli] = c.y2Phase1;
            Y1Phase2[celli] = dummyValue(c, c.species1, c.phase2);
            Y2Phase2[celli] = dummyValue(c, c.species2, c.phase2);
        }
        else if (ac <= alphaTol)
        {
            // Pure phase2: phase1 fields are dummy values.
            X1[celli] = c.x1Phase2;
            X2[celli] = c.x2Phase2;

            Y1Phase1[celli] = dummyValue(c, c.species1, c.phase1);
            Y2Phase1[celli] = dummyValue(c, c.species2, c.phase1);
            Y1Phase2[celli] = c.y1Phase2;
            Y2Phase2[celli] = c.y2Phase2;
        }
        else
        {
            // Mixed interface cell: both phases exist, so both sets are real.
            // X is the total mixture mole fraction and must be weighted with
            // alpha*C, not with alpha alone.
            const scalar w1 = ac*CPhase1[celli];
            const scalar w2 = (1 - ac)*CPhase2[celli];
            const scalar wSum = w1 + w2;

            if (wSum <= VSMALL)
            {
                FatalErrorInFunction
                    << "Invalid alpha*C denominator in mixed cell " << celli
                    << ": alpha=" << ac
                    << ", C(" << c.phase1 << ")=" << CPhase1[celli]
                    << ", C(" << c.phase2 << ")=" << CPhase2[celli]
                    << nl << exit(FatalError);
            }

            X1[celli] = (w1*c.x1Phase1 + w2*c.x1Phase2)/wSum;
            X2[celli] = 1 - X1[celli];

            Y1Phase1[celli] = c.y1Phase1;
            Y2Phase1[celli] = c.y2Phase1;
            Y1Phase2[celli] = c.y1Phase2;
            Y2Phase2[celli] = c.y2Phase2;
        }
    }

    calculateXMoleFractionBoundaryValues
    (
        mesh,
        alpha,
        alphaTol,
        c,
        CPhase1Boundary,
        CPhase2Boundary,
        X1Boundary,
        X2Boundary
    );

    writeScalarFieldFile
    (
        mesh,
        runTime,
        alpha,
        setCompositionDict,
        c.species1 + ".X",
        X1,
        c.species1,
        word::null,
        false,
        c,
        &X1Boundary
    );

    writeScalarFieldFile
    (
        mesh,
        runTime,
        alpha,
        setCompositionDict,
        c.species2 + ".X",
        X2,
        c.species2,
        word::null,
        false,
        c,
        &X2Boundary
    );

    writeScalarFieldFile
    (
        mesh,
        runTime,
        alpha,
        setCompositionDict,
        c.species1 + "." + c.phase1,
        Y1Phase1,
        c.species1,
        c.phase1,
        true,
        c,
        nullptr
    );

    writeScalarFieldFile
    (
        mesh,
        runTime,
        alpha,
        setCompositionDict,
        c.species2 + "." + c.phase1,
        Y2Phase1,
        c.species2,
        c.phase1,
        true,
        c,
        nullptr
    );

    writeScalarFieldFile
    (
        mesh,
        runTime,
        alpha,
        setCompositionDict,
        c.species1 + "." + c.phase2,
        Y1Phase2,
        c.species1,
        c.phase2,
        true,
        c,
        nullptr
    );

    writeScalarFieldFile
    (
        mesh,
        runTime,
        alpha,
        setCompositionDict,
        c.species2 + "." + c.phase2,
        Y2Phase2,
        c.species2,
        c.phase2,
        true,
        c,
        nullptr
    );

    Info<< nl << "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
