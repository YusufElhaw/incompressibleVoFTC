#ifndef setComposition_H
#define setComposition_H

#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "wordList.H"
#include "dictionary.H"
#include "error.H"
#include "token.H"
#include "fileName.H"
#include "PtrList.H"

namespace Foam
{
namespace setCompositionTools
{

struct compositionInput
{
    word species1;
    word species2;

    word phase1;
    word phase2;

    word defaultSpeciePhase1;
    word defaultSpeciePhase2;

    scalar W1;
    scalar W2;

    fileName W1Source;
    fileName W2Source;

    scalar x1Phase1;
    scalar x2Phase1;
    scalar x1Phase2;
    scalar x2Phase2;

    scalar y1Phase1;
    scalar y2Phase1;
    scalar y1Phase2;
    scalar y2Phase2;
};

word readWordEntry(const dictionary& dict, const word& key);

bool sameWordNoCase(const word& a, const word& b);

bool findEntryNameNoCase(const dictionary& dict, const word& wanted, word& actual);

wordList readPhases(const Time& runTime, const fvMesh& mesh, const dictionary& dict);

word findPhase1FromAlphaField(const Time& runTime, const fvMesh& mesh, const wordList& phases);

wordList readSpeciesFromPhysicalProperties
(
    const Time& runTime,
    const fvMesh& mesh,
    const word& phase
);

scalar moleToMassFraction
(
    const scalar x1,
    const scalar x2,
    const scalar W1,
    const scalar W2,
    const label speciesIndex
);

compositionInput readCompositionInput
(
    const Time& runTime,
    const fvMesh& mesh,
    const dictionary& dict
);

bool containsWord(const wordList& words, const word& w);

bool isConstraintPatchType(const word& meshPatchType, const word& fieldPatchType);

scalar phaseValue
(
    const compositionInput& c,
    const word& species,
    const word& phase,
    const bool massFraction
);

scalar inletValueForPatch
(
    const compositionInput& c,
    const dictionary& dict,
    const word& patchName,
    const word& fieldSpecies,
    const word& fieldPhase,
    const bool massFraction
);

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
);

} // End namespace setCompositionTools
} // End namespace Foam

#endif
