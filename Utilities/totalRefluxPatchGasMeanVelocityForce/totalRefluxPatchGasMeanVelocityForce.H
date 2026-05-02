/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
   \\  /    A nd           |
    \\/     M anipulation  |
-------------------------------------------------------------------------------

    totalRefluxPatchGasMeanVelocityForce

    Solver-independent fvConstraint for total reflux operation in two-phase VOF
    cases.  It uses only fields from the objectRegistry and is intended to run
    with incompressible/compressible VOF-style two-phase solvers.

    Robust default liquid-load measurement:

        Q_L = sum_patch alphaLiquid*(U & gasDirection)*magSf
        B_L = max(0, -Q_L/A_0)       // counterCurrent total reflux mode
        U_G,soll = (rho_L/rho_G) B_L gasDirection

    Compared to an alphaPhi-only measurement this avoids extremely large
    intermediate alphaPhi values that can occur on NCC/coupled patches during
    pressure/alpha corrections.  alphaPhi mode is still available for cases
    where exact surfaceFieldValue-like patch-normal fluxes are desired.

\*---------------------------------------------------------------------------*/

#ifndef totalRefluxPatchGasMeanVelocityForce_H
#define totalRefluxPatchGasMeanVelocityForce_H

#include "fvConstraint.H"
#include "fvCellZone.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
namespace fv
{

class totalRefluxPatchGasMeanVelocityForce
:
    public fvConstraint
{
    // Private Data

        //- Cell zone where the force is applied
        fvCellZone zone_;

        //- Patch names used for liquid loading measurement and gas control
        wordList patches_;

        //- Name of the velocity field
        word UName_;

        //- Name of liquid volume fraction field. Empty means auto-detect
        word alphaName_;

        //- Optional name of liquid phase flux field. Empty means auto-detect
        word alphaPhiName_;

        //- Liquid loading measurement mode: directionalAlphaU or alphaPhi
        word liquidFluxMode_;

        //- Name of mixture volume flux field for alpha*phi fallback
        word phiName_;

        //- Name of mixture density field for automatic density estimation
        word rhoName_;

        //- Optional phase density field names
        word rhoLiquidFieldName_;
        word rhoGasFieldName_;

        //- Optional explicit phase densities. Negative means not specified
        scalar rhoLiquid_;
        scalar rhoGas_;

        //- Cutoffs used for estimating rhoL/rhoG from mixture rho(alpha)
        scalar densityLiquidAlphaCutoff_;
        scalar densityGasAlphaCutoff_;

        //- Fixed gas target direction. Default is upward y: (0 1 0)
        vector gasDirection_;

        //- Time interval for liquid loading averaging [s]. <=0 means instantaneous
        scalar averagingInterval_;

        //- Relaxation factor for the pressure-gradient update
        scalar relaxation_;

        //- Optional lower cutoff for alphaG. Cells with alphaG <= alphaCutoff_
        //  contribute zero force.
        scalar alphaCutoff_;

        //- Whether to suppress the source below alphaCutoff_
        bool useCutoff_;

        //- Smooth ramp width above alphaCutoff_ for source activation
        scalar alphaRamp_;

        //- For total reflux, use only liquid motion opposite to gasDirection
        bool counterCurrent_;

        //- Optional cap for |targetUbar_|. A negative value disables the cap.
        scalar maxGasLoading_;

        //- Optional cap for |gradP|. A negative value disables the cap.
        scalar maxPressureGradient_;

        //- Optional cap for |dGradP| per velocity correction. Negative disables.
        scalar maxDeltaPressureGradient_;

        //- Last time-index used for integrating the liquid flux
        mutable label lastIntegratedTimeIndex_;

        //- Accumulated liquid volume over the current averaging window [m3]
        mutable scalar accumulatedLiquidVolume_;

        //- Accumulated time over the current averaging window [s]
        mutable scalar accumulatedTime_;

        //- Last measured liquid volume flow rate Q_L [m3/s]
        mutable scalar liquidFlowRate_;

        //- Last selected patch/segment area A_0 [m2]
        mutable scalar patchArea_;

        //- Last instantaneous liquid loading B_L,inst [m/s]
        mutable scalar instantLiquidLoading_;

        //- Last filtered liquid loading B_L [m/s]
        mutable scalar liquidLoading_;

        //- True after the liquid loading filter has been initialized
        mutable bool liquidLoadingInitialized_;

        //- Last computed target superficial gas velocity vector
        mutable vector targetUbar_;

        //- Last phase densities used for the reflux ratio
        mutable scalar rhoLiquidPatchMean_;
        mutable scalar rhoGasPatchMean_;

        //- Text describing the last density source used
        mutable word densitySource_;

        //- Ignore restart file and start pressure-gradient controller from zero
        bool resetGradient_;

        //- Pressure gradient before correction
        mutable scalar gradP0_;

        //- Change in pressure gradient
        mutable scalar dGradP_;

        //- Matrix 1/A coefficients field pointer
        mutable autoPtr<volScalarField> rAPtr_;


    // Private Member Functions

        //- Read coefficients from dictionary
        void readCoeffs(const dictionary& dict);

        //- Return true if a word looks like a liquid alpha field name
        static bool isAlphaFieldName(const word& name);

        //- Return true if the word looks like a gas/vapour phase name
        static bool isGasLikeName(const word& name);

        //- Return selected/autodetected liquid alpha field name
        word actualAlphaName() const;

        //- Return the liquid alpha-flux field name actually used, or word::null
        word actualAlphaPhiName(const volScalarField& alpha1) const;

        //- Return a likely gas phase/group name for diagnostics and density fields
        word actualGasPhaseName(const volScalarField& alpha1) const;

        //- Write the pressure gradient and loading state to file for restarts
        void writeProps(const scalar gradP) const;

        //- Gas superficial-flux weighting factor from liquid volume fraction alpha1
        scalar gasFluxWeight(const scalar alpha1) const;

        //- Smooth source weighting factor from liquid volume fraction alpha1
        scalar sourceWeight(const scalar alpha1) const;

        //- Smooth source weighting factor for a cell
        scalar sourceWeight(const volScalarField& alpha1, const label celli) const;

        //- Clamp a volume fraction into [0, 1]
        static scalar limitedAlpha(const scalar alpha);

        //- Measure selected-patch signed liquid flow rate and area
        void measureLiquidFlowRate
        (
            const volVectorField& U,
            const volScalarField& alpha1,
            scalar& liquidFlowRate,
            scalar& patchArea
        ) const;

        //- Update the liquid-loading average according to formula 4.26/4.27
        void updateLiquidLoading
        (
            const volVectorField& U,
            const volScalarField& alpha1
        ) const;

        //- Compute phase densities from constants, phase fields or mixture rho
        void updateDensities(const volScalarField& alpha1) const;

        //- Compute the total-reflux gas target from averaged liquid loading
        void updateTarget
        (
            const volVectorField& U,
            const volScalarField& alpha1
        ) const;

        //- Calculate the current patch superficial gas velocity in target direction
        scalar patchGasUbarAve
        (
            const volVectorField& U,
            const volScalarField& alpha1
        ) const;

        //- Patch-based response d(jG_patch)/d(gradP), using rA on patches
        scalar patchResponseAve
        (
            const volScalarField& rA,
            const volScalarField& alpha1
        ) const;

        //- Zone fallback response d(jG)/d(gradP), used if patch response is zero
        scalar zoneResponseAve
        (
            const volScalarField& rA,
            const volScalarField& alpha1
        ) const;


public:

    //- Runtime type information
    TypeName("totalRefluxPatchGasMeanVelocityForce");

    // Constructors

        //- Construct from explicit source name and mesh
        totalRefluxPatchGasMeanVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        totalRefluxPatchGasMeanVelocityForce
        (
            const totalRefluxPatchGasMeanVelocityForce&
        ) = delete;


    // Member Functions

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Add the momentum source and set the 1/A coefficient field
        virtual bool constrain
        (
            fvMatrix<vector>& eqn,
            const word& fieldName
        ) const;

        //- Correct the velocity field and update the pressure gradient
        virtual bool constrain(volVectorField& U) const;

        // Mesh changes
        virtual bool movePoints();
        virtual void topoChange(const polyTopoChangeMap& map);
        virtual void mapMesh(const polyMeshMap& map);
        virtual void distribute(const polyDistributionMap& map);

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};

} // End namespace fv
} // End namespace Foam

#endif

// ************************************************************************* //
