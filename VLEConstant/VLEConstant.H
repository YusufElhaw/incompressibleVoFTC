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
    Foam::VLEConstant

Description
    class to manage the value of Vapur Liquid Equilibrium (VLE) constant based
    on the Rauolt's law. 
    y 𝜑_L p  =  x 𝛾 𝜑_G p_sat exp{1/(RT)*  ∫( v_0 )dp }
    y and x are the molfraction in the gas phase and the liquid phase 
    respectivly. 
    𝛾 is the activity coeffecient. 
    𝜑_L and 𝜑_G are the fugacity factors.
    p_sat is the saturation pressure of the Species.
    exp{1/(RT)*  ∫( v_0 )dp|s,0-P } is the poynting correction.

    The VLE constant is in the formation of K=y/x. 
    For non-polar and chemical similar mixtures and low pressures can the 
    fugacity factors and the Poynting factor and the activity coeffecient be 
    assumed to be 1. 
    
    For more details see:
    Themis Matsoukas: Fundamentals of Chemical Engineering Thermodynamics.


SourceFiles
    VLEConstant.C

\*---------------------------------------------------------------------------*/

#ifndef VLEConstant_H
#define VLEConstant_H

#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "incompressibleInterPhaseTransportModelTC.H"
#include "incompressibleTwoPhaseVoFMixtureTC.H"

namespace Foam
{

class VLEConstant
{
    const incompressibleInterPhaseTransportModelTC& turbulence_;
    const incompressibleTwoPhaseVoFMixtureTC& mixture_;

    // parameters per-
    dimensionedScalar Dm1_, Dm2_;  // [m2/s]

    //  turbulence per-phase flags (read once)
    bool turb1_, turb2_;

    // dictionaries (cached in registry)
    const IOdictionary& readThermophysicalTransportDict(const word& phaseName) const;
    
    // checking from momentumTransport both Phases for turbulence
    bool checkPhaseTurbulence(const word& phaseName) const;

public:
    VLEConstant(const incompressibleInterPhaseTransportModelTC& turb);

    // saturation pressure
    virtual void pSat(); 

    // activity coeffecient 
    virtual void gamma();

    // fugacity coeffecient of the liquid phase
    virtual void phiL();

    // fugacity coeffecient of the Gas phase
    virtual void phiG();

    // poynting correction
    virtual void Poy();


};

} // End namespace Foam

#endif

