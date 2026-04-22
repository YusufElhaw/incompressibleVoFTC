/*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     | Website:  https://openfoam.org
     \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
      \\/     M anipulation  |
 -------------------------------------------------------------------------------
 License
     This file is part of OpenFOAM.
  
     OpenFOAM is free software: you can redistribute it and/or modify it
     under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.
  
     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
     for more details.
  
     You should have received a copy of the GNU General Public License
     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
  
 Class
     Foam::prghNonConformalCyclicPressureFvPatchScalarField
  
 Description
     This boundary condition provides a cyclic condition for p_rgh in NCC
     "non-conformal cyclic" situations. It applies corrections to the value 
     and gradient on both sides of the cyclic to account for the non-cylicity 
     of the gravitational force.
  
     This condition is only needed when the cyclic patches have a transformation
     and a normal component in the direction of gravity. If the cyclic patches
     are orthogonal to the direction gravity, then a normal cyclic boundary
     condition can be used instead.
  
     Care must be taken when using this boundary condition that the simulation
     is actually cyclic. The following constraints apply:
  
     - Both cyclic patches must be oriented in the same way with respect to
       gravity. In practice this means that applicability is limited to cyclics
       with translational transformations.
  
     - The model cannot have any dependence on the absolute value of the
       pressure field. The absolute value of the pressure, in reality, varies
       between each repetition of the geometry; it is not actually formally
       cyclic. Only the gradient of the pressure field can be truly cyclic. This
       model is therefore only valid if the absolute value of the pressure is
       arbitrary, and only the gradient has an effect on the solution. This is
       the case for incompressible multiphase solutions or incompressible
       Boussinesq-like models of density variation. It is not true if (for
       example) a compressible thermodynamic model is being used.
  
 Example of the boundary condition specification usage:
    patches top and bottom are walls in blockMeshDict, but are coupled together 
    in createNonConformalCouples utility with the following specification:
    1- in 0/p_rgh:
        bottom
        {
            type           fixedFluxPressure;
            gradient        uniform 0;
            value           uniform 0;
        }
        top
        {
            type            fixedFluxPressure;
            gradient        uniform 0;
            value           uniform 0;
        }    
    
    2- in 0/U:
        bottom
        {
            type            slip; 
        }
        top
        {
            type            slip;
        }
    
    3- in 0/alpha.* or any other fields:
        bottom
        {
            type            zeroGradient;
        }
        top
        {
            type            zeroGradient;
        }
       
    4- in controlDict:  //add the library and couple definition
     libs
        (
            "libprghNonConformalCyclicPressure.so"
        );

    5- in createNonConformalCouplesDict:
    fields  yes; //    only for one processor. for more then one processor, should be "no" and the patch fields should be added manually to the patches in the couple definition as shown below.

    nonConformalCouples
    {
        bottom_Top // or any other name for the couple
        {
            owner
                {
                    patch       bottom;

                    patchFields
                    {
                        p_rgh
                        {
                            type        prghNonConformalCyclicPressure;
                            patchType   nonConformalCyclic;
                            rho         rho;
                            rhoInf      1;
                            value       uniform 0;
                        }
                    }
                }
            neighbour
                {
                    patch       top;

                    patchFields
                    {
                        p_rgh
                        {
                            type        prghNonConformalCyclicPressure;
                            patchType   nonConformalCyclic;
                            value       uniform 0;
                        }
                    }
                }
            transform       translational;
            separation      (0 -0.025 0); // or what ever the separation between the cyclic patches is
        }    
    }
6-  for more than one processor, do the decomposPar first. do the cut only 
    parallel to the gravity direction. so that the coupeld Patchs for example
    ("top" and "bottom") are in the same processor.
    example for 48 processors: (8 1 6). where the gravity in the y direction and
    the normal to the coupled patches "top" and "bottom" are in the y direction.  
    Then the patch field for p_rgh should be added manually to the patches in the
    couple definition in the processors folders as shown below, and "fields no" 
    should be used in createNonConformalCouplesDict to avoid adding the patch 
    fields to all patches in the couple: nonConformalCouples.

    Then in the processors folders under /0/p_rgh
   you get something like this:
            bottom_Top_on_bottom
            {
                type nonConformalCyclic;
            }
            bottom_Top_on_top
            {
                type nonConformalCyclic;
            }
        change it to this:
            bottom_Top_on_bottom
            {
                type        prghNonConformalCyclicPressure;
                patchType   nonConformalCyclic;
                rho         rho;
                rhoInf      1;
                value       uniform 0;
            }
            bottom_Top_on_top
            {
                type        prghNonConformalCyclicPressure;
                patchType   nonConformalCyclic;
                value       uniform 0;
            }
        you can  just use this line in the terminal to do it automatically for all processors:
        \verbatim    
sed -i '/nonConformalCyclic/c\
type      prghNonConformalCyclicPressure;\
patchType   nonConformalCyclic;\
rho         rho;\
rhoInf      1;\
value       uniform 0;' processor* /0/p_rgh
        \endverbatim //           ^delete this character
        then start the simulation in parallel as usual. 
        before running reconstructPar you have to change those back in every processor folder in every time to 
        "type nonConformalCyclic;". you can use this line in the terminal to do it automatically for all time steps and all processors:
        \verbatim
perl -0pi -e 's/(\{\s*)[^{}]*?type\s+prghNonConformalCyclicPressure;[^{}]*?(\s*\})/$1    type            nonConformalCyclic;\n$2/gs' processor* / * /p_rgh
        \endverbatim  //                                                                                                                                   ^ ^ ^ delete these characters
        
        then do the reconstructPar as usual.
  
 SourceFiles
     prghNonConformalCyclicPressureFvPatchScalarField.C
  
 \*---------------------------------------------------------------------------*/

#ifndef prghNonConformalCyclicPressureFvPatchScalarField_H
#define prghNonConformalCyclicPressureFvPatchScalarField_H

#include "nonConformalCyclicFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

class prghNonConformalCyclicPressureFvPatchScalarField
:
    public nonConformalCyclicFvPatchField<scalar>
{
    // Private data

        //- Density field name
        word rhoName_;

        //- Reference density
        scalar rhoInf_;

        //- Pressure jump (neighbour minus owner)
        scalarField jump_;


public:

    TypeName("prghNonConformalCyclicPressure");


    // Constructors

        //- Construct from patch and internal field
        prghNonConformalCyclicPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        prghNonConformalCyclicPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        prghNonConformalCyclicPressureFvPatchScalarField
        (
            const prghNonConformalCyclicPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Copy constructor setting internal field reference
        prghNonConformalCyclicPressureFvPatchScalarField
        (
            const prghNonConformalCyclicPressureFvPatchScalarField&,
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
                new prghNonConformalCyclicPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

        //- Return neighbour field including jump
        virtual tmp<Field<scalar>> patchNeighbourField
        (
            const Pstream::commsTypes
        ) const;

        //- Map from self
        virtual void map
        (
            const fvPatchScalarField&,
            const fieldMapper&
        );

        //- Reset from self
        virtual void reset
        (
            const fvPatchScalarField&
        );

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Update matrix contribution for scalar solves
        virtual void updateInterfaceMatrix
        (
            scalarField& result,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;
        const scalarField& jump() const
        {
            return jump_;
        }
        
        //- Write
        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
