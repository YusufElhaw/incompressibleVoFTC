/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::incompressibleTwoPhaseVoFMixtureT

Description
    Class to represent a mixture of two constant density phases

SourceFiles
    incompressibleTwoPhaseVoFMixtureT.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoPhaseVoFMixtureTC_H
#define incompressibleTwoPhaseVoFMixtureTC_H

#include "twoPhaseVoFMixture.H"
#include "incompressibleTwoPhases.H"
#include "viscosityModel.H"
#include "rhoFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class incompressibleTwoPhaseVoFMixtureT Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhaseVoFMixtureTC
:
    public twoPhaseVoFMixture,
    virtual public incompressibleTwoPhases,
    public viscosity
{
// Private Data

        Switch totalInternalEnergy_;

        //- Pressure
        volScalarField p_;

        //- Mixture temperature
        volScalarField T_;
        ////
        // --- Species (always 2) ---
        word species1Name_;
        word species2Name_;

        //- Species mass-fractions per phase
        //  File names in 0/ and time/:
        //    <species1>.<phase1>, <species2>.<phase1>
        //    <species1>.<phase2>, <species2>.<phase2>
        autoPtr<volScalarField> s1Phase1_;
        autoPtr<volScalarField> s2Phase1_;
        autoPtr<volScalarField> s1Phase2_;
        autoPtr<volScalarField> s2Phase2_;

        //- Indices in thermo species tables (per phase thermo)
        label s1IndexThermo1_;
        label s2IndexThermo1_;
        label s1IndexThermo2_;
        label s2IndexThermo2_;
        
      //dimensionedScalar residualAlpha1_;
      //dimensionedScalar residualAlpha2_;

        ////////
        //- Thermo-package of phase 1
        autoPtr<rhoFluidMulticomponentThermo> thermo1_;

        //- Thermo-package of phase 2
        autoPtr<rhoFluidMulticomponentThermo> thermo2_;

    //- Viscosity model for phase 1
    autoPtr<viscosityModel> nuModel1_;

    //- Viscosity model for phase 2
    autoPtr<viscosityModel> nuModel2_;

    //- Constant density of phase 1
    dimensionedScalar rho1_;

    //- Constant density of phase 2
    dimensionedScalar rho2_;

    //- Mixture density
    volScalarField rho_;

    //- Mixture viscosity
    volScalarField nu_;
 
public:

    TypeName("incompressibleTwoPhaseVoFMixtureTC");


    // Constructors

        //- Construct from a mesh
        incompressibleTwoPhaseVoFMixtureTC(const fvMesh& mesh);


    //- Destructor
    virtual ~incompressibleTwoPhaseVoFMixtureTC()
    {}


    // Member Functions
        //- Return true to solve for total internal energy
        //  return false to solve for internal energy
        bool totalInternalEnergy() const
        {
            return totalInternalEnergy_;
        }
                //- Return pressure [Pa]
        volScalarField& p()
        {
            return p_;
        }

        //- Return mixture temperature [K]
        volScalarField& T()
        {
            return T_;
        }

        //- Return mixture temperature [K]
        const volScalarField& T() const
        {
            return T_;
        }

        //- Return the thermo for phase 1
        const rhoFluidMulticomponentThermo& thermo1() const
        {
            return thermo1_();
        }

        //- Return the thermo for phase 2
        const rhoFluidMulticomponentThermo& thermo2() const
        {
            return thermo2_();
        }

        //- Return the thermo for phase 1
        rhoFluidMulticomponentThermo& thermo1()
        {
            return thermo1_();
        }

        //- Return the thermo for phase 2
        rhoFluidMulticomponentThermo& thermo2()
        {
            return thermo2_();
        }
        //- Correct the thermodynamics of each phase
        virtual void correctThermo();
        
        //- Correct the composition of each phase
        virtual void correctComposition();
        
            /////
            // Species accessors (mass fractions)
            volScalarField& s1Phase1() { return s1Phase1_(); }
            volScalarField& s2Phase1() { return s2Phase1_(); }
            volScalarField& s1Phase2() { return s1Phase2_(); }
            volScalarField& s2Phase2() { return s2Phase2_(); }

            const volScalarField& s1Phase1() const { return s1Phase1_(); }
            const volScalarField& s2Phase1() const { return s2Phase1_(); }
            const volScalarField& s1Phase2() const { return s1Phase2_(); }
            const volScalarField& s2Phase2() const { return s2Phase2_(); }

            const word& species1Name() const { return species1Name_; }
            const word& species2Name() const { return species2Name_; }
            
          //  const dimensionedScalar& residualAlpha1() const;
           // const dimensionedScalar& residualAlpha2() const;


            //////

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Return the mixture density
        virtual const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }
        
        //- Correct the mixture density and laminar viscosity
        virtual void correct();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
