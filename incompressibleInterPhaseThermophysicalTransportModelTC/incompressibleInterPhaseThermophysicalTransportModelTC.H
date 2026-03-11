/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::incompressibleInterPhaseThermophysicalTransportModelTC

Description
    incompressible two-phase VoF thermophysicalTransportModel

    Provides the mixture thermal conductivity and diffusivity for
    Fourier's temperature gradient heat flux model.

SourceFiles
    incompressibleInterPhaseThermophysicalTransportModelTC.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleInterPhaseThermophysicalTransportModelTC_H
#define incompressibleInterPhaseThermophysicalTransportModelTC_H

#include "thermophysicalTransportModel.H"
#include "massDiffusivity.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleInterPhaseTransportModelTC;
class massDiffusivity;  

/*---------------------------------------------------------------------------*\
    Class incompressibleInterPhaseThermophysicalTransportModelTC Declaration
\*---------------------------------------------------------------------------*/

class incompressibleInterPhaseThermophysicalTransportModelTC
:
    public thermophysicalTransportModel
{
    // Private Data

        //- incompressibleInterPhaseTransportModel
        const incompressibleInterPhaseTransportModelTC& turbulence_;
        massDiffusivity massDiffusivity_;



public:

    // Constructors

        //- Construct from components
        incompressibleInterPhaseThermophysicalTransportModelTC
        (
            const incompressibleInterPhaseTransportModelTC& turbulence
        );

        //- Disallow default bitwise copy construction
        incompressibleInterPhaseThermophysicalTransportModelTC
        (
            const incompressibleInterPhaseThermophysicalTransportModelTC&
        ) = delete;


    // Member Functions
         
        //- Read model coefficients if they have changed
        virtual bool read();
        
        // Heat Transfer

            //- Effective thermal turbulent conductivity
            //  of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff() const ; // german lambda = kappa

            //- Effective thermal turbulent conductivity
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff(const label patchi) const;// german lambda = kappa

            //- Return the effective temperature transport coefficient
            //  derived from the phase internal energy equations i.e. <kappa/Cv>
            tmp<volScalarField> alphaEff() const;   // german a= alpha Temperaturleitfähigkeit

            //- Return the heat flux [W/m^2]
            virtual tmp<surfaceScalarField> q() const;

            //- Return the patch heat flux [W/m^2]
            virtual tmp<scalarField> q(const label patchi) const;

            //- Return the patch heat flux correction [W/m^2]
            //  For isotropic or patch-aligned thermal conductivity qCorr is null
            virtual tmp<scalarField> qCorr(const label patchi) const;

            //- Return the source term for the energy equation
            virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;
        // Mass Transfer
        
            //- Effective Mass turbulent conductivity
            //  of mixture [m^2/s]
            virtual tmp<volScalarField> DEff() const ;      // analog to a= alpha Temperaturleitfähigkeit 
            virtual tmp<volScalarField> D1Eff() const ;      // analog to a= alpha Temperaturleitfähigkeit 
            virtual tmp<volScalarField> D2Eff() const ;      // analog to a= alpha Temperaturleitfähigkeit 

           //- Effective turbulent Diffusion coeffcient
           //  of mixture for patch [m^2/s]
           virtual tmp<scalarField> DEff(const label patchi) const;
           virtual tmp<scalarField> D1Eff(const label patchi) const;
           virtual tmp<scalarField> D2Eff(const label patchi) const;

           
//            //- Return the Mass flux [Mol/(m^2*s)]
//            virtual tmp<surfaceScalarField> j() const;
//
//            //- Return the patch Mass flux [Mol/(m^2*s)]
//            virtual tmp<scalarField> j(const label patchi) const;
//
//            //- Return the patch Mass flux correction [Mol/(m^2*s)]
//            //  For isotropic or patch-aligned Mass conductivity jCorr is null
//            virtual tmp<scalarField> jCorr(const label patchi) const;
//            
            //- Predict the transport coefficients if possible
            //  without solving thermophysical transport model equations
            virtual void predict();

            //- Solve the thermophysical transport model equations
            //  and correct the transport coefficients
            virtual void correct();


        // Member Operators

            //- Disallow default bitwise assignment
            void operator=
            (
                const incompressibleInterPhaseThermophysicalTransportModelTC&
            ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
