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

    Two-phase control:

    Phase 1 (pre-SS): feed-forward
        jGbase = nDown / (CG * A0)      [molar mode]
               = mDown / (rhoG * A0)    [mass mode]
        smoothed by EMA filter (averagingInterval)

    I-controller (always active once N_ref is known):
        e_rel = (N_boiler(t) - N_ref) / N_ref
        jGcorrection_ += K * e_rel * jGeff / tau * dt
        jGcorrection_ = clamp(jGcorrection_, ±antiWindupFraction * jGeff)
        jGtarget = max(0, jGbase + jGcorrection_)
    K is dimensionless and scale-invariant. jGeff = gasLoadingTarget_ (smoothed).

    Always: rate limiter on |d(jGtarget)/dt| <= maxDeltaGasLoading.
    No velocity cap, no window balance controller.

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

        //- Patch names used for measurement and gas control
        wordList patches_;

        //- Resolved patch indices: includes original NCC patches from patches_
        //  plus nonConformalProcessorCyclic patches that refer to them (parallel)
        labelList patchIDs_;

        //- Name of the velocity field
        word UName_;

        //- Name of liquid volume fraction field. Empty means auto-detect
        word alphaName_;

        //- Optional name of liquid phase alpha flux field. Empty means auto-detect
        word alphaPhiName_;

        //- Liquid/gas flux measurement mode: directionalAlphaU or alphaPhi
        word liquidFluxMode_;

        //- Name of mixture volume flux field for alpha*phi fallback
        word phiName_;

        //- Name of mixture density field for automatic density estimation
        word rhoName_;

        //- Optional phase density field names
        word rhoLiquidFieldName_;
        word rhoGasFieldName_;

        //- If true, target total reflux by total molar flux, not mass flux
        Switch molarTotalReflux_;

        //- Phase/species fields for molar total reflux closure
        word species1LiquidName_;
        word species2LiquidName_;
        word species1GasName_;
        word species2GasName_;

        //- Molecular weights [kg/mol] internally. Dict input is kg/kmol.
        scalar W1_;
        scalar W2_;

        //- Optional explicit phase densities. Negative means not specified
        scalar rhoLiquid_;
        scalar rhoGas_;

        //- Cutoffs used for estimating rhoL/rhoG from mixture rho(alpha)
        scalar densityLiquidAlphaCutoff_;
        scalar densityGasAlphaCutoff_;

        //- Fixed gas target direction. Default is upward y: (0 1 0)
        vector gasDirection_;

        //- Time interval for EMA loading averaging [s]. <=0 means instantaneous
        scalar averagingInterval_;

        //- Relaxation factor for the pressure-gradient update
        scalar relaxation_;

        //- Optional lower cutoff for alphaG source activation
        scalar alphaCutoff_;

        //- Whether to suppress the source below alphaCutoff_
        bool useCutoff_;

        //- Smooth ramp width above alphaCutoff_ for source activation
        scalar alphaRamp_;

        //- For total reflux, downward/outgoing means opposite to gasDirection
        bool counterCurrent_;

        //- Max allowed |d(jGtarget)/dt| [m/s/s]. <=0 disables rate limiter.
        scalar maxDeltaGasLoading_;

        //- Last time-index used for integrating the measured reflux / I-controller
        mutable label lastIntegratedTimeIndex_;

        //- Accumulated equivalent target gas volume over averaging window [m3]
        mutable scalar accumulatedEquivalentGasVolume_;

        //- Accumulated time over averaging window [s]
        mutable scalar accumulatedTime_;

        //- Last downward liquid volume flow rate [m3/s]
        mutable scalar liquidFlowRate_;

        //- Last selected patch/segment area [m2]
        mutable scalar patchArea_;

        //- Last filtered equivalent gas loading [m/s]
        mutable scalar liquidLoading_;

        //- True after filter initialization
        mutable bool liquidLoadingInitialized_;

        //- Total downward mass flow rate [kg/s]
        mutable scalar totalMassFlowRateDown_;

        //- Downward liquid volume flow rate [m3/s]
        mutable scalar liquidVolumeFlowRateDown_;

        //- Downward gas volume flow rate [m3/s]
        mutable scalar gasVolumeFlowRateDown_;

        //- Instantaneous target gas loading [m/s]
        mutable scalar instantGasLoadingTarget_;

        //- Filtered target gas loading [m/s]
        mutable scalar gasLoadingTarget_;

        //- Previous timestep's final jG target, used by rate limiter [m/s]
        mutable scalar jGtargetPrev_;

        //- True after steady-state detection
        mutable bool steadyStateReached_;

        //- Optional patch names measured ONLY for SS detection.
        //  Typically the top/condenser patch. Empty means use control patches.
        wordList steadyStatePatchNames_;

        //- Resolved patch IDs for SS detection (incl. NCC processor patches).
        labelList steadyStatePatchIDs_;

        //- Relative threshold: |eDot_avg|/flow_ref < relTol triggers SS.
        //  Negative (default) disables auto-detection.
        scalar balanceSteadyStateRelTol_;

        //- Minimum absolute flow [mol/s or kg/s] before SS windows are counted.
        scalar balanceSteadyStateMinFlow_;

        //- Consecutive SS detection windows below relTol required to confirm SS.
        label balanceSteadyStateConfirmCount_;

        //- Penalty subtracted from consecutiveSteadyWindows_ per failed window.
        label balanceSteadyStateViolationPenalty_;

        //- Counter of consecutive qualifying windows.
        mutable label consecutiveSteadyWindows_;

        //- Boiler moles [mol] at the moment SS was detected.
        //  Negative when not yet captured (no SS or boiler BC not found).
        mutable scalar reservoirMolesAtSS_;

        //- I-controller: accumulated jG correction [m/s] since SS detection.
        //  d(jGcorrection_)/dt = K * e / (CG * A0 * tau), e = N_boiler - N_ref
        mutable scalar jGcorrection_;

        //- Gain for the I-controller on boiler moles [-]
        scalar steadyStateIntegralGain_;

        //- Timescale [s] normalising the I gain: d(jGcorr)/dt = K*e_rel*jGeff/tau
        scalar steadyStateIntegralTimescale_;

        //- Anti-windup fraction: |jGcorrection| <= antiWindupFraction * jGeff
        scalar antiWindupFraction_;

        //- Width of SS detection window [s] (dict key: balanceWindow)
        scalar ssDetectWindow_;

        //- Accumulated eDot*dt in the current SS detection window [mol or kg]
        mutable scalar ssDetectWindowError_;

        //- Elapsed time in current SS detection window [s]
        mutable scalar ssDetectWindowTime_;

        //- Last time index used for SS detection window accumulation
        mutable label lastSSDetectTimeIndex_;

        //- Last upward gas mass flow through the measured patch [kg/s]
        mutable scalar gasMassFlowRateUp_;

        //- Last upward gas volume flow through the measured patch [m3/s]
        mutable scalar gasVolumeFlowRateUp_;

        //- Last downward total molar flow, from column into apparatus [mol/s]
        mutable scalar totalMolarFlowRateDown_;

        //- Last upward gas molar flow, from apparatus into column [mol/s]
        mutable scalar gasMolarFlowRateUp_;

        //- Gas total molar concentration used at the patch [mol/m3]
        mutable scalar gasMolarConcentrationPatchMean_;

        //- Last balance error rate. If molarTotalReflux=true: mol/s; else kg/s
        mutable scalar balanceErrorRate_;

        //- Last computed target superficial gas velocity vector
        mutable vector targetUbar_;

        //- Last phase densities used
        mutable scalar rhoLiquidPatchMean_;
        mutable scalar rhoGasPatchMean_;

        //- Text describing density source
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

        //- Build patchIDs_ from patches_ (original NCC) plus
        //  nonConformalProcessorCyclic patches that refer to them (parallel)
        void buildPatchIDs();

        //- Build steadyStatePatchIDs_ from steadyStatePatchNames_ (same logic).
        //  Falls back to patchIDs_ when steadyStatePatchNames_ is empty.
        void buildSteadyStatePatchIDs();

        //- Read reservoir moles directly from inletOutletBoilerFvPatchScalarField
        //  found at the control patches. Returns -1 if not available.
        //  In parallel: broadcasts from owning processor via maxOp reduce.
        scalar readReservoirMoles() const;

        //- Measure instantaneous eDot = nDown - nGasUp (molar) or mDown - mGasUp
        //  (mass) at the given patch indices. Used for SS detection and integral.
        scalar measureInstantEDot
        (
            const volVectorField& U,
            const volScalarField& alpha1,
            const labelList& pids
        ) const;

        //- Accumulate SS detection window; set steadyStateReached_ when confirmed
        void detectSteadyState
        (
            const volVectorField& U,
            const volScalarField& alpha1
        ) const;

        //- Return true if a word looks like a liquid alpha field name
        static bool isAlphaFieldName(const word& name);

        //- Return true if the word looks like a gas/vapour phase name
        static bool isGasLikeName(const word& name);

        //- Return selected/autodetected liquid alpha field name
        word actualAlphaName() const;

        //- Return liquid alpha-flux field name actually used, or word::null
        word actualAlphaPhiName(const volScalarField& alpha1) const;

        //- Return gas alpha-flux field name if available, otherwise word::null
        word actualGasAlphaPhiName(const volScalarField& alpha1) const;

        //- Return a likely gas phase/group name
        word actualGasPhaseName(const volScalarField& alpha1) const;

        //- Write pressure gradient and measured state to uniform dictionary
        void writeProps(const scalar gradP) const;

        //- Gas superficial-flux weighting factor from liquid alpha
        scalar gasFluxWeight(const scalar alpha1) const;

        //- Smooth source weighting factor from liquid alpha
        scalar sourceWeight(const scalar alpha1) const;

        //- Smooth source weighting factor for a cell
        scalar sourceWeight(const volScalarField& alpha1, const label celli) const;

        //- Clamp alpha into [0,1]
        static scalar limitedAlpha(const scalar alpha);

        //- Measure total downward/outgoing mass reflux:
        //      mDown = rhoL*QLdown + rhoG*QGdown
        void measureDownwardTotalReflux
        (
            const volVectorField& U,
            const volScalarField& alpha1,
            scalar& totalMassFlowDown,
            scalar& liquidVolumeFlowDown,
            scalar& gasVolumeFlowDown,
            scalar& gasVolumeFlowUp,
            scalar& patchArea
        ) const;

        //- True if all fields needed for molar reflux closure exist
        bool molarFluxFieldsAvailable() const;

        //- Measure patch-normal molar total reflux using exactly U&Sf, alpha, rho, Y, W.
        //  Positive U&Sf leaves the CFD domain into the apparatus.
        //  The upward gas return from apparatus into the CFD domain is negative gas flux.
        void measureDownwardTotalRefluxMolar
        (
            const volVectorField& U,
            const volScalarField& alpha1,
            scalar& totalMolarFlowDown,
            scalar& gasMolarFlowUp,
            scalar& liquidVolumeFlowDown,
            scalar& gasVolumeFlowDown,
            scalar& gasVolumeFlowUp,
            scalar& gasMolarConcentration,
            scalar& patchArea
        ) const;

        //- Update SS detection (called after each velocity correction)
        void updateBalanceController
        (
            const volVectorField& U,
            const volScalarField& alpha1
        ) const;

        //- Update target gas loading from measured total downward mass flow
        void updateTotalRefluxLoading
        (
            const volVectorField& U,
            const volScalarField& alpha1
        ) const;

        //- Compute phase densities from constants, phase fields or mixture rho
        void updateDensities(const volScalarField& alpha1) const;

        //- Compute target gas velocity
        void updateTarget
        (
            const volVectorField& U,
            const volScalarField& alpha1
        ) const;

        //- Current patch superficial gas velocity in target direction
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

        //- Zone fallback response d(jG)/d(gradP)
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
