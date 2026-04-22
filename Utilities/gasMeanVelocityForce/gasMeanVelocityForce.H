/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    gasMeanVelocityForce

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
        type            gasMeanVelocityForce;
        cellZone        all;

        U               U;
        alpha           alpha.liquid;   // always the liquid VoF field from 0/
        Ubar            (0 1 0);      // targeted gas loading / mean gas velocity [m/s]
        relaxation      0.1;

        useCutoff       false;  // optional
        alphaCutoff     1e-3;   // only used if useCutoff=true. Cells with alphaG <= alphaCutoff contribute zero weight.
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef gasMeanVelocityForce_H
#define gasMeanVelocityForce_H

#include "fvConstraint.H"
#include "fvCellZone.H"
#include "volFields.H"

namespace Foam
{
namespace fv
{

class gasMeanVelocityForce
:
    public fvConstraint
{
    // Private data

        //- The cellZone the fvConstraint applies to
        fvCellZone zone_;

        //- Name of the velocity field
        word UName_;

        //- Name of the liquid alpha field; gas fraction is alphaG = 1 - alpha
        word alphaName_;

        //- Desired superficial gas velocity vector [m/s]
        vector Ubar_;

        //- Relaxation factor for the forcing update
        scalar relaxation_;

        //- Optional lower cutoff for alphaG. Cells with alphaG <= alphaCutoff_
        //  contribute zero weight.
        scalar alphaCutoff_;

        //- Whether to use a hard cutoff mask in addition to alphaG weighting
        bool useCutoff_;

        //- Pressure gradient before correction
        mutable scalar gradP0_;

        //- Change in pressure gradient
        mutable scalar dGradP_;

        //- Matrix 1/A coefficients field pointer
        mutable autoPtr<volScalarField> rAPtr_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const scalar gradP) const;
   
protected:

    //- Gas weighting factor from liquid volume fraction alphaL
    scalar weight(const scalar alphaL) const;

    //- Gas weighting factor for a cell
    scalar weight(const volScalarField& alpha, const label celli) const;

    //- Calculate the current superficial gas velocity projected on Ubar_
    virtual scalar magUbarGasAve
    (
        const volVectorField& U,
        const volScalarField& alpha
    ) const;

public:

    //- Runtime type information
    TypeName("gasMeanVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        gasMeanVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        gasMeanVelocityForce(const gasMeanVelocityForce&) = delete;


    // Member Functions

        //- Disallow default bitwise copy construction
        const vector& Ubar() const
        {
            return Ubar_;
        }

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Add the momentum source and set the 1/A coefficient field
        virtual bool constrain
        (
            fvMatrix<vector>& eqn,
            const word& fieldName
        ) const;

        //- Correct the pressure gradient and the velocity field
        virtual bool constrain(volVectorField& U) const;

        //- Update for mesh motion
        virtual bool movePoints();

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Read source dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gasMeanVelocityForce&) = delete;
};

} // End namespace fv
} // End namespace Foam

#endif

// ************************************************************************* //
