/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Parallel NCC variant of prghNonConformalCyclicPressure for p_rgh on
    nonConformalProcessorCyclic patches.

    This patch field mirrors the logic of prghCyclicPressure and of the custom
    prghNonConformalCyclicPressure, but it uses processor communication for the
    remote side and therefore also works when the two coupled NCC halves live on
    different MPI ranks.
\*---------------------------------------------------------------------------*/

#ifndef prghNonConformalProcessorCyclicPressureFvPatchScalarField_H
#define prghNonConformalProcessorCyclicPressureFvPatchScalarField_H

#include "nonConformalProcessorCyclicFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

class prghNonConformalProcessorCyclicPressureFvPatchScalarField
:
    public nonConformalProcessorCyclicFvPatchField<scalar>
{
    // Private data

        //- Local reference cast into the processor patch
        const nonConformalProcessorCyclicFvPatch& procPatch_;

        //- Name of the density field. Defaults to rho on both sides.
        word rhoName_;

        //- Owner-side reference density
        scalar rhoInf_;
        bool rhoInfSpecified_;

        //- True on the owner side of the non-conformal cyclic pair
        bool owner_;

        //- Jump in value from the neighbouring side to this side
        scalarField jump_;

        //- Scalar send buffer for custom matrix updates
        mutable scalarField scalarSendBuf_;

        //- Scalar receive buffer for custom matrix updates
        mutable scalarField scalarReceiveBuf_;

        //- Outstanding request ids for custom matrix updates
        mutable label outstandingSendRequest_;
        mutable label outstandingRecvRequest_;


    // Private Member Functions

        //- Exchange a scalar field with the neighbour processor
        tmp<scalarField> swapWithNeighbour(const scalarField& fld) const;

        //- Exchange a single scalar with the neighbour processor
        scalar swapWithNeighbour(const scalar s) const;


public:

    TypeName("prghNonConformalProcessorCyclicPressure");


    // Constructors

        prghNonConformalProcessorCyclicPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        prghNonConformalProcessorCyclicPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        prghNonConformalProcessorCyclicPressureFvPatchScalarField
        (
            const prghNonConformalProcessorCyclicPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        prghNonConformalProcessorCyclicPressureFvPatchScalarField
        (
            const prghNonConformalProcessorCyclicPressureFvPatchScalarField&
        ) = delete;

        prghNonConformalProcessorCyclicPressureFvPatchScalarField
        (
            const prghNonConformalProcessorCyclicPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new prghNonConformalProcessorCyclicPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Return neighbour coupled field including the hydrostatic jump
        virtual tmp<Field<scalar>> patchNeighbourField
        (
            const Pstream::commsTypes
        ) const;

        //- Map the given fvPatchField onto this fvPatchField
        virtual void map(const fvPatchScalarField&, const fieldMapper&);

        //- Reset the fvPatchField to the given fvPatchField
        virtual void reset(const fvPatchScalarField&);

        //- Start processor exchange for the boundary values
        virtual void initEvaluate(const Pstream::commsTypes commsType);

        //- Complete processor exchange and then add the p_rgh jump
        virtual void evaluate(const Pstream::commsTypes commsType);

        //- Update the patch coefficients and recompute the jump
        virtual void updateCoeffs();

        //- Start matrix-interface communication for scalar solves
        virtual void initInterfaceMatrixUpdate
        (
            scalarField& result,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;

        //- Complete matrix-interface communication for scalar solves
        virtual void updateInterfaceMatrix
        (
            scalarField& result,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;

        //- Write dictionary entries
        virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif

// ************************************************************************* //
