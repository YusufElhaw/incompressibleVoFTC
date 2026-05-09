/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Derived solver module
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derived from OpenFOAM. It is free software: you can
    redistribute it and/or modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

Class
    Foam::solvers::incompressibleVoFCyclic

Description
    incompressibleVoF variant which solves the momentum and pressure correction
    directly with the static pressure p instead of the buoyant pressure p_rgh.

    The OpenFOAM VoF base classes still construct/read p_rgh; this module keeps
    p_rgh as a compatibility/output field, but the overridden momentum predictor
    and pressure corrector use p consistently.

\*---------------------------------------------------------------------------*/

#ifndef incompressibleVoFCyclic_H
#define incompressibleVoFCyclic_H

#include "incompressibleVoF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                    Class incompressibleVoFCyclic Declaration
\*---------------------------------------------------------------------------*/

class incompressibleVoFCyclic
:
    public incompressibleVoF
{
protected:

    //- Pressure reference based on the static pressure p
    autoPtr<Foam::pressureReference> pPressureReferencePtr_;


    // Protected Member Functions

        //- Return the pressure reference based on p, not p_rgh
        virtual const Foam::pressureReference& pressureReference() const override
        {
            return pPressureReferencePtr_();
        }

public:

    //- Runtime type information
    TypeName("incompressibleVoFCyclic");


    // Constructors

        //- Construct from region mesh
        incompressibleVoFCyclic(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        incompressibleVoFCyclic(const incompressibleVoFCyclic&) = delete;


    //- Destructor
    virtual ~incompressibleVoFCyclic();


    // Member Functions

        //- Construct and solve the momentum predictor with p
        virtual void momentumPredictor() override;

        //- Construct and solve the pressure equation directly for p
        virtual void pressureCorrector() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const incompressibleVoFCyclic&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
