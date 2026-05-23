/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Dynamic reboiler finite-holdup inletOutlet boundary condition for the
    light-key mole-fraction field X.

    Sign convention:
      - positive phase molar flux = phase leaves CFD domain into boiler
      - negative phase molar flux = phase enters CFD domain from boiler

    The boiler contains a liquid reservoir with total holdup N_B and liquid
    composition x_B. Vapour leaving the boiler and entering the CFD domain is
    in VLE with the reservoir liquid:

        y_B = A12*x_B/(1 + (A12 - 1)*x_B)

    Complete reservoir balances:
        dN_B/dt =
            L_toB + G_toB
          - L_fromB - G_fromB

        d(N_B*x_B)/dt =
            L_toB*x_L,out + G_toB*y_G,out
          - L_fromB*x_B   - G_fromB*y_B

    Patch behaviour:
      - incoming liquid from boiler receives x_B
      - incoming gas from boiler receives y_B
      - if both enter through the same face, refValue is flux-weighted
      - outgoing faces use zeroGradient behaviour

\*---------------------------------------------------------------------------*/
#ifndef inletOutletBoilerFvPatchScalarField_H
#define inletOutletBoilerFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

class inletOutletBoilerFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    //- Liquid total molar flux [mol/s]. Positive leaves CFD domain.
    word liquidFluxName_;

    //- Gas total molar flux [mol/s]. Positive leaves CFD domain.
    word gasFluxName_;

    //- Relative volatility field created by compositionPredictor.
    word A12Name_;

    //- Optional group name. If set and no explicit patches list is given,
    //  all patches of this boundary-condition type on the same field
    //  with the same group are coupled.
    word groupName_;

    //- Coupled boiler patches. Empty means only self.
    wordList patchNames_;

    //- True if patchNames_ was explicitly supplied.
    bool usePatchList_;

    //- Boiler liquid holdup [mol].
    scalar reservoirMoles_;

    //- Boiler liquid mole fraction of species 1.
    scalar reservoirX_;

    //- Current vapour mole fraction in VLE with reservoirX.
    scalar reservoirY_;

    //- Under-relaxation for reservoir update.
    //  1 = physical explicit Euler; <1 = numerical damping.
    scalar relax_;

    //- Bounds for composition.
    scalar minValue_;
    scalar maxValue_;

    //- Minimum boiler holdup [mol] to avoid division by zero.
    scalar minReservoirMoles_;

    //- Small flux below which faces are treated as no-flow.
    scalar minFlux_;

    //- Print reservoir diagnostics.
    Switch log_;

    //- Last time index at which the reservoir was advanced.
    label lastReservoirUpdateIndex_;

    // Diagnostics
    scalar lastA12_;
    scalar lastLiquidToReservoir_;
    scalar lastGasToReservoir_;
    scalar lastLiquidFromReservoir_;
    scalar lastGasFromReservoir_;
    scalar lastSpeciesToReservoir_;
    scalar lastSpeciesFromReservoir_;
    scalar lastMolesResidual_;
    scalar lastSpeciesResidual_;
    scalar lastLiquidXToReservoir_;
    scalar lastGasXToReservoir_;
    
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

    //- Clamp scalar to [minValue_, maxValue_].
    scalar clamp(const scalar x) const;

    //- VLE map y = A*x/(1 + (A - 1)*x).
    scalar yVLE(const scalar x, const scalar A) const;

    //- Return selected coupled patches.
    wordList selectedPatchNames() const;

    //- Return patch ID or abort.
    label patchID(const word& patchName) const;

    //- Lookup a surfaceScalarField or abort with available names.
    const surfaceScalarField& lookupSurfaceFlux(const word& name) const;

    //- Read selected patches from dictionary.
    void readPatchSelection(const dictionary& dict);

    //- Compute patch/group mean A12.
    scalar meanA12(const wordList& patchNames) const;


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
    
    //- Inlet/outlet fraction for all phases entering from the boiler.
    tmp<scalarField> reservoirInletFraction() const;

    //- Flux-weighted inlet composition for phases entering from the boiler.
    tmp<scalarField> reservoirInletValue() const;

public:

    TypeName("inletOutletBoiler");

    // Constructors

    inletOutletBoilerFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    inletOutletBoilerFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    inletOutletBoilerFvPatchScalarField
    (
        const inletOutletBoilerFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fieldMapper&
    );

    inletOutletBoilerFvPatchScalarField
    (
        const inletOutletBoilerFvPatchScalarField&
    ) = delete;

    inletOutletBoilerFvPatchScalarField
    (
        const inletOutletBoilerFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new inletOutletBoilerFvPatchScalarField(*this, iF)
        );
    }
    
    scalar reservoirMoles() const { return reservoirMoles_; }

    scalar reservoirX() const { return reservoirX_; }

    scalar reservoirY() const { return reservoirY_; }

    //- Advance the boiler reservoir using post-solve composition.
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


