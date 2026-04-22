/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    patchGasMeanVelocityForce

    A user fvConstraint derived from the OpenFOAM Foundation v13
    meanVelocityForce model. It regulates a superficial gas velocity / gas
    loading by weighting the driving force with alpha_G = 1 - alphaLiquid.

    Intended for incompressibleVoF-style single-U simulations where the user
    wants to preserve periodic/NCC-cyclic boundaries and drive a prescribed
    upward gas loading through a representative periodic packing segment.

    The controlled quantity is

        jG = (1/V_zone) * sum(alpha.Gas * (e.U) * V_cells)

    which corresponds to the superficial gas velocity in the forcing direction
    for a periodic zone of uniform cross-section.
    Example usage:
    \verbatim
    gasDrive
    {
        type            patchGasMeanVelocityForce;
        cellZone        all;            // the cellZone where the force is applied
        patch           inlet;          // the patch where jG controlled
        patches ();

        U               U;
        alpha           alpha.liquid;   // always the liquid VoF field from 0/
        Ubar            (0 1 0);      // targeted gas loading / mean gas velocity [m/s]
        relaxation      0.1;

        useCutoff       false;  // optional
        alphaCutoff     1e-3;   // only used if useCutoff=true. Cells with alphaG <= alphaCutoff contribute zero weight.
    }
    \endverbatim

    \verbatim
    momentumSource
    {
        type            patchGasMeanVelocityForce;

        cellZone        all;
        U               U;
        alpha           alpha.water;    // always the liquid VoF field from 0/

        
        Ubar            (0 1.5 0);      // targeted gas loading / mean gas velocity [m/s]

        // multiple patches
        // by NCC-Coupling write the names of the NCC patches, not the original patches
        patches 
        (
            bottom_01_on_top
            bottom_02_on_top
            bottom_03_on_top
            bottom_04_on_top
        );

        relaxation      0.02;

        useCutoff       false;
        alphaCutoff     1e-3;
    }

    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef patchGasMeanVelocityForce_H
#define patchGasMeanVelocityForce_H

#include "gasMeanVelocityForce.H"

namespace Foam
{
namespace fv
{

class patchGasMeanVelocityForce
:
    public gasMeanVelocityForce
{
    // Private data

        //- Patch names used for the control measurement
        wordList patches_;


    // Private Member Functions

        void readCoeffs(const dictionary& dict);


protected:

        virtual scalar magUbarGasAve
        (
            const volVectorField& U,
            const volScalarField& alpha
        ) const;

public:

    TypeName("patchGasMeanVelocityForce");

    patchGasMeanVelocityForce
    (
        const word& sourceName,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict
    );

    patchGasMeanVelocityForce(const patchGasMeanVelocityForce&) = delete;

    virtual bool read(const dictionary& dict);

    void operator=(const patchGasMeanVelocityForce&) = delete;
};

} // End namespace fv
} // End namespace Foam

#endif
// ************************************************************************* //
