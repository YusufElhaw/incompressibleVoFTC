/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
This Document is part of the Master thesis of Yusuf Elhaw.

Description
    Finite-holdup total-condenser/reflux-drum boundary condition for the
    light-key mole-fraction field X.

    Sign convention:
      - positive phase molar flux  = phase leaves CFD domain into reservoir
      - negative phase molar flux  = phase enters CFD domain from reservoir

    Reservoir balances:
        dN_D/dt      = sum(phi_L^+ + phi_G^+) - sum(phi_L^- + phi_G^-)

        d(N_D*x_D)/dt =
            sum((phi_L^+ + phi_G^+) * X_patchInternal)
          - sum((phi_L^- + phi_G^-) * x_D)

    Patch behaviour:
      - incoming reservoir flow gets refValue = reservoirX
      - outgoing flow uses zeroGradient behaviour
      - by default, incoming liquid and incoming gas are both treated as
        reservoir outflow for species conservation; gas inlet treatment can be
        disabled with applyToGasInlet false.

\*---------------------------------------------------------------------------*/
#ifndef inletOutletCondenserFvPatchScalarField_H
#define inletOutletCondenserFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

class inletOutletCondenserFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    //- Liquid total molar flux [mol/s]. Positive leaves CFD domain.
    word liquidFluxName_;

    //- Gas total molar flux [mol/s]. Positive leaves CFD domain.
    word gasFluxName_;

    //- Optional group name. If set and no explicit patches list is given,
    //  all patches of this boundary-condition type on the same field
    //  with the same group are coupled.
    word groupName_;

    //- Coupled top/condenser patches. Empty means only self.
    wordList patchNames_;

    //- True if patchNames_ was explicitly supplied.
    bool usePatchList_;

    //- Condenser/reflux-drum total holdup [mol].
    scalar reservoirMoles_;

    //- Reflux-drum mole fraction of species 1.
    scalar reservoirX_;

    //- Under-relaxation for reservoir update.
    //  1 = physical explicit Euler; <1 = numerical damping.
    scalar relax_;

    //- Bounds for composition.
    scalar minValue_;
    scalar maxValue_;

    //- Minimum total holdup [mol] to avoid division by zero.
    scalar minReservoirMoles_;

    //- Small flux below which faces are treated as no-flow for switching.
    scalar minFlux_;

    //- If true, gas entering the CFD domain also receives reservoirX.
    //  This makes the species balance consistent for all incoming phases.
    Switch applyToGasInlet_;

    //- Print reservoir diagnostics.
    Switch log_;

    //- Last time index at which the reservoir was advanced.
    label lastReservoirUpdateIndex_;

    // Diagnostics
    scalar lastLiquidToReservoir_;
    scalar lastGasToReservoir_;
    scalar lastLiquidFromReservoir_;
    scalar lastGasFromReservoir_;
    scalar lastSpeciesToReservoir_;
    scalar lastSpeciesFromReservoir_;
    scalar lastMolesResidual_;
    scalar lastSpeciesResidual_;

    //- Write detailed time-resolved apparatus balance
    Switch writeBalance_;
    
    //- If true, reservoir N/N1 fluxes are measured from rho/Y/W/U/Sf
    //  instead of alphaCLPhi/alphaCGPhi times the transported X field.
    Switch useRawSpeciesFlux_;

    //- Raw flux field names
    word UName_;
    word alphaLiquidName_;
    word rhoLiquidName_;
    word rhoGasName_;
    word species1LiquidName_;
    word species2LiquidName_;
    word species1GasName_;
    word species2GasName_;

    //- Molecular weights [kg/mol] internally.
    //  Dictionary input is kg/kmol and converted by /1000.
    scalar W1_;
    scalar W2_;

    //- Name of VLE relative volatility field in the mesh registry.
    word A12Name_;

    //- Clamp scalar to [minValue_, maxValue_].
    scalar clamp(const scalar x) const;

    //- VLE: y = A12*x / (1 + (A12-1)*x), clamped.
    scalar yVLE(scalar x, scalar A12) const;

    //- Solve for liquid mole fraction xL given overall X, Lto, Gto, A12.
    //  Uses VLE equilibrium: X*(Lto+Gto) = Lto*xL + Gto*yVLE(xL,A12).
    scalar xLiquidFromXtotal(scalar X, scalar Lto, scalar Gto, scalar A12) const;

    //- Compute area-weighted mean A12 over selected patches.
    scalar meanA12(const wordList& patchNames) const;

    //- Return selected coupled patches.
    wordList selectedPatchNames() const;

    //- Return patch ID or abort.
    label patchID(const word& patchName) const;

    //- Lookup a surfaceScalarField or abort with available names.
    const surfaceScalarField& lookupSurfaceFlux(const word& name) const;

    //- Read selected patches from dictionary.
    void readPatchSelection(const dictionary& dict);

    //- Write detailed apparatus balance diagnostics into postProcessing
    void writeApparatusBalance
    (
        const scalar dt,
        const scalar NOld,
        const scalar xOld,
        const scalar NNew,
        const scalar xNew,
        const scalar nLtoReservoir,
        const scalar nGtoReservoir,
        const scalar nLfromReservoir,
        const scalar nGfromReservoir,
        const scalar n1ToReservoir,
        const scalar n1FromReservoir
    ) const;

        //- Return true if all raw species-flux fields exist.
    bool compositionFluxFieldsAvailable() const;

    //- Measure raw molar N and N1 fluxes using rho/Y/W/U/Sf.
    //  Positive flux leaves CFD domain into the reservoir.
    void measureCompositionPredictorFluxes
    (
        const wordList& patchNames,
        scalar& nLtoReservoir,
        scalar& nGtoReservoir,
        scalar& nLfromReservoir,
        scalar& nGfromReservoir,
        scalar& n1LtoReservoir,
        scalar& n1GtoReservoir,
        scalar& n1LfromReservoir,
        scalar& n1GfromReservoir
    ) const;

    //- Inlet/outlet fraction for incoming reservoir flow.
    tmp<scalarField> reservoirInletFraction() const;

public:

    TypeName("inletOutletCondenser");

    // Constructors

    inletOutletCondenserFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    inletOutletCondenserFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    inletOutletCondenserFvPatchScalarField
    (
        const inletOutletCondenserFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fieldMapper&
    );

    inletOutletCondenserFvPatchScalarField
    (
        const inletOutletCondenserFvPatchScalarField&
    ) = delete;

    inletOutletCondenserFvPatchScalarField
    (
        const inletOutletCondenserFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new inletOutletCondenserFvPatchScalarField(*this, iF)
        );
    }
    
    scalar reservoirMoles() const { return reservoirMoles_; }

    scalar reservoirX() const { return reservoirX_; }

    //- Advance the condenser reservoir using post-solve composition.
    //  Called explicitly from compositionPredictor after the solve.
    void advanceReservoir();

    //- Correct reservoir state for NCPC molar non-conservation.
    //  deltaN  [mol]: moles to add (negative removes moles from reservoir).
    //  deltaN1 [mol]: species-1 moles to add correspondingly.
    void correctReservoir(const scalar deltaN, const scalar deltaN1);

    //- In parallel runs, directly set reservoir inlet composition on
    //  nonConformalProcessorCyclic faces for incoming flow.
    //  Must be called AFTER X.correctBoundaryConditions() (MPI exchange)
    //  and BEFORE fvMatrix assembly on each time step / Picard iteration.
    void applyInletValuesOnNCPCPatches(volScalarField& X) const;

    //- Update coefficients.
    virtual void updateCoeffs();

    // Zero-gradient for diffusion operators only.
    // Convective inletOutlet behaviour (outflow=zeroGradient, inflow=fixedValue
    // via refValue) is unaffected because fvc::div uses patchFaceValues(), not
    // snGrad(). Only fvm::laplacian uses these methods, so zeroing them
    // eliminates the untracked diffusive flux at the apparatus boundary without
    // touching species transport by convection.
    virtual tmp<scalarField> snGrad() const
    {
        return tmp<scalarField>(new scalarField(this->size(), scalar(0)));
    }

    virtual tmp<scalarField> gradientInternalCoeffs() const
    {
        return tmp<scalarField>(new scalarField(this->size(), scalar(0)));
    }

    virtual tmp<scalarField> gradientBoundaryCoeffs() const
    {
        return tmp<scalarField>(new scalarField(this->size(), scalar(0)));
    }

    //- Write.
    virtual void write(Ostream&) const;

    //- Assignment.
    virtual void operator=(const fvPatchField<scalar>&);
};

} // End namespace Foam

#endif

// ************************************************************************* //


