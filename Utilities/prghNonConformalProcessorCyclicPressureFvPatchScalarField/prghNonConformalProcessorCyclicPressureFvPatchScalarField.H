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



sed -i '/nonConformalProcessorCyclic/c\
type            prghNonConformalProcessorCyclicPressure;\
patchType       nonConformalProcessorCyclic;\
rho             rho;\
rhoInf          1;\
value           uniform 0;' processor* /0/p_rgh


\*---------------------------------------------------------------------------*/
/*#ifndef prghNonConformalProcessorCyclicPressureFvPatchScalarField_H
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

        //- Reference to processor patch
        const nonConformalProcessorCyclicFvPatch& procPatch_;

        //- Kept for dictionary compatibility
        word rhoName_;
        scalar rhoInf_;
        bool rhoInfSpecified_;
        bool owner_;

        //- Pressure jump transported from the referred NCC patch
        scalarField jump_;

        //- Buffers for matrix communication
        mutable scalarField scalarSendBuf_;
        mutable scalarField scalarReceiveBuf_;
        mutable label outstandingSendRequest_;
        mutable label outstandingRecvRequest_;


    // Private member functions

        //- Pull jump from referred prghNonConformalCyclic patch,
        //  reverse-map to orig patch faces, then map onto this processor patch
        void updateJumpFromReferPatch();


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
            const prghNonConformalProcessorCyclicPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member functions

        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new prghNonConformalProcessorCyclicPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

        const scalarField& jump() const
        {
            return jump_;
        }

        virtual tmp<Field<scalar>> patchNeighbourField
        (
            const Pstream::commsTypes
        ) const;

        virtual void map
        (
            const fvPatchScalarField&,
            const fieldMapper&
        );

        virtual void reset
        (
            const fvPatchScalarField&
        );

        virtual void initEvaluate
        (
            const Pstream::commsTypes commsType
        );

        virtual void evaluate
        (
            const Pstream::commsTypes commsType
        );

        virtual tmp<Field<scalar>> snGrad() const;

        virtual tmp<Field<scalar>> snGrad
        (
            const scalarField& deltaCoeffs
        ) const;

        virtual void updateCoeffs();

        virtual void initInterfaceMatrixUpdate
        (
            scalarField& result,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;

        virtual void updateInterfaceMatrix
        (
            scalarField& result,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;

        virtual void write(Ostream&) const;
};

} // End namespace Foam 

#endif*/

//Gemini 
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Parallel NCC variant of prghNonConformalCyclicPressure for p_rgh on
    nonConformalProcessorCyclic patches.
\*---------------------------------------------------------------------------*/
/*
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

        //- Reference to processor patch
        const nonConformalProcessorCyclicFvPatch& procPatch_;

        //- Kept for dictionary compatibility
        word rhoName_;
        scalar rhoInf_;
        bool rhoInfSpecified_;
        bool owner_;

        //- Pressure jump calculated from local and neighbour fields
        scalarField jump_;

        //- Buffers for matrix communication
        mutable scalarField scalarSendBuf_;
        mutable scalarField scalarReceiveBuf_;
        mutable label outstandingSendRequest_;
        mutable label outstandingRecvRequest_;

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
            const prghNonConformalProcessorCyclicPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member functions

        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new prghNonConformalProcessorCyclicPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

        const scalarField& jump() const
        {
            return jump_;
        }

        virtual tmp<Field<scalar>> patchNeighbourField
        (
            const Pstream::commsTypes
        ) const;

        virtual void map
        (
            const fvPatchScalarField&,
            const fieldMapper&
        );

        virtual void reset
        (
            const fvPatchScalarField&
        );

        virtual void initEvaluate
        (
            const Pstream::commsTypes commsType
        );

        virtual void evaluate
        (
            const Pstream::commsTypes commsType
        );

        virtual tmp<Field<scalar>> snGrad() const;

        virtual tmp<Field<scalar>> snGrad
        (
            const scalarField& deltaCoeffs
        ) const;

        virtual void updateCoeffs();

        virtual void initInterfaceMatrixUpdate
        (
            scalarField& result,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;

        virtual void updateInterfaceMatrix
        (
            scalarField& result,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;

        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif*/










/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  13
     \\/     M anipulation  |
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

        //- Kept for dictionary compatibility / write()
        word rhoName_;

        //- Kept for dictionary compatibility / write()
        scalar rhoInf_;

        //- Jump copied from the referred local NCC patch
        scalarField jump_;


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


    // Member functions

        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new prghNonConformalProcessorCyclicPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

        virtual tmp<Field<scalar>> patchNeighbourField
        (
            const Pstream::commsTypes
        ) const;

        virtual void map
        (
            const fvPatchScalarField&,
            const fieldMapper&
        );

        virtual void reset
        (
            const fvPatchScalarField&
        );

        virtual void updateCoeffs();

        virtual void evaluate
        (
            const Pstream::commsTypes commsType =
                Pstream::commsTypes::blocking
        );

        virtual void updateInterfaceMatrix
        (
            scalarField& result,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;

        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif