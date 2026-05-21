/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Function object for total-reflux rectification simulations.

    Reads reservoir states (ND, xD, NB, xB, yB) from the active
    inletOutletCondenser / inletOutletBoiler BC objects and integrates the
    molar density fields (nMolesTotal, nMoles1, nMoles2) over the column
    volume.  From these it computes:

      - System-wide mole totals:  NSystem = NColumn + ND + NB
      - Conservation drift:       dN, dN1, dN2 relative to the first time step
      - Volume-weighted geometric mean A12 over the column
      - Fenske minimum stages:    Nmin = ln(xD(1-xB)/xB(1-xD)) / ln(A12_geo)
                                  (includes boiler as one equilibrium stage)
      - Packing stages:           Npacking = Nmin - 1  (column stages only)
      - HETP:                     HETPColumnGeo = packingLength / Npacking
                                  (valid only when Nmin > 1)

    Results are appended to
      postProcessing/<functionObjectName>/0/values.dat

    Usage in system/functions:

        myFenske
        {
            type            reservoirFenske;
            libs            ("libreservoirFenske.so");

            executeControl  timeStep;
            executeInterval 1;
            writeControl    none;           // file write happens in execute()

            speciesField    cyclohexane.X;
            A12Field        A12;

            columnTotalMolesField    nMolesTotal;
            columnSpecies1MolesField nMoles1;
            columnSpecies2MolesField nMoles2;

            // One representative patch per apparatus is sufficient.
            // All patches in a group share the same reservoir state.
            topReservoirPatches    ( bottom_Top_on_top    );
            bottomReservoirPatches ( bottom_Top_on_bottom );

            packingLength   0.05;   // [m] total packing height — must match mesh
        }

\*---------------------------------------------------------------------------*/

#ifndef reservoirFenske_H
#define reservoirFenske_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "inletOutletCondenserFvPatchScalarField.H"
#include "inletOutletBoilerFvPatchScalarField.H"

#include <fstream>

namespace Foam
{
namespace functionObjects
{

class reservoirFenske
:
    public fvMeshFunctionObject
{
    // Configuration
    word        speciesField_;
    word        A12Field_;
    word        columnTotalMolesField_;
    word        columnSpecies1MolesField_;
    word        columnSpecies2MolesField_;
    wordList    topReservoirPatches_;
    wordList    bottomReservoirPatches_;
    scalar      packingLength_;

    // Baseline captured at the first valid time step for drift reporting.
    // Reset each time foamRun starts (relative to the current run, not t=0
    // of the original simulation if restarting).
    bool   initialized_;
    scalar NSystem0_;
    scalar N1System0_;
    scalar N2System0_;

    // Private helpers

    //- Return patch index or FatalError.
    label patchID(const word& patchName) const;

    //- Sum fieldName * cellVolume over all cells (parallel-safe via reduce).
    bool volumeIntegral(scalar& result, const word& fieldName) const;

    //- Volume-weighted geometric mean: exp(sum V*ln(f) / sum V).
    //  Cells with f <= SMALL are excluded.
    bool volumeGeoMean(scalar& result, const word& fieldName) const;

    //- Read condenser state from the first inletOutletCondenser BC found on
    //  patchNames.  With group coupling all patches share the same state, so
    //  one representative patch is sufficient.
    bool readCondenserReservoir
    (
        scalar& ND,
        scalar& xD,
        const volScalarField& X,
        const wordList& patchNames
    ) const;

    //- Read boiler state from the first inletOutletBoiler BC found on
    //  patchNames.
    bool readBoilerReservoir
    (
        scalar& NB,
        scalar& xB,
        scalar& yB,
        const volScalarField& X,
        const wordList& patchNames
    ) const;

    //- Fenske equation (result includes boiler as one equilibrium stage).
    static scalar fenskeNmin(scalar xD, scalar xB, scalar A12);

    //- Guard: all inputs in valid range for Fenske (0 < xB < xD < 1, A12 > 1).
    static bool validFenskeInputs(scalar xD, scalar xB, scalar A12);

    //- Write scalar or "nan" to an already-open stream.
    static void writeValueOrNa(std::ofstream& os, bool valid, scalar value);

    //- Read configuration from dict (shared by constructor and read()).
    void readCoeffs(const dictionary& dict);

public:

    TypeName("reservoirFenske");

    reservoirFenske
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual ~reservoirFenske();

    virtual bool read(const dictionary& dict);

    //- Fields required by this function object.
    virtual wordList fields() const;

    //- Compute and append one row to values.dat.
    virtual bool execute();

    //- No-op: all output is produced in execute().
    virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif
// ************************************************************************* //
