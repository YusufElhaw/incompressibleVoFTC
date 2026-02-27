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
    class to manage the value of Vapour Liquid Equilibrium (VLE) constant based
    on the Raoult's law. 
    y 𝜑_L p  =  x 𝛾 𝜑_G p_sat exp{1/(RT)*  ∫( v_0 )dp }
    y and x are the molfraction in the gas phase and the liquid phase 
    respectively. 
    𝛾 is the activity coefficient. 
    𝜑_L and 𝜑_G are the fugacity factors.
    p_sat is the saturation pressure of the Species.
    exp{1/(RT)*  ∫( v_0 )dp|s,0-P } is the poynting correction.

    The VLE constant is in the formation of K=y/x. 
    For non-polar and chemical similar mixtures and low pressures can the 
    fugacity factors and the Poynting factor and the activity coefficient be 
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
#include "dictionary.H"
#include "autoPtr.H"
#include "word.H"

#include "saturationPressureModel.H"
#include "incompressibleInterPhaseTransportModelTC.H"
#include "incompressibleTwoPhaseVoFMixtureTC.H"

namespace Foam
{

class VLEConstant
{
    const incompressibleInterPhaseTransportModelTC& turbulence_;
    const incompressibleTwoPhaseVoFMixtureTC& mixture;

    /// Names
    word phaseL_;     // e.g. "phase1" / "liquid"
    word phaseG_;     // e.g. "phase2" / "gas"
    word specie1_;    // optional label/name for user readability
    word specie2_;
    //scalar KMin_;
    //scalar KMax_;

    // Saturation pressure models (one per species)
    autoPtr<saturationPressureModel> pSat1Model_;
    autoPtr<saturationPressureModel> pSat2Model_;

    const IOdictionary& readDictOrDie(const word& dictName) const;
    const dictionary& coeffsDict(const word& dictName) const;

public:
    // Construct from transport + mixture.
    // dictName: file in constant/, "VLEProperties"
    VLEConstant
    (
        const incompressibleInterPhaseTransportModelTC& turb,
        const incompressibleTwoPhaseVoFMixtureTC& mix,
        const word& dictName = "VLEProperties"
    );

    // pSat(T) fields
    tmp<volScalarField> pSat1(const volScalarField& T) const;
    tmp<volScalarField> pSat2(const volScalarField& T) const;

    // K = pSat/p (Raoult, idealised: gamma=phiL=phiG=Poy=1)
    tmp<volScalarField> K1(const volScalarField& p, const volScalarField& T) const;
    tmp<volScalarField> K2(const volScalarField& p, const volScalarField& T) const;

    // Convenience: use phase temperatures (if desired)
    tmp<volScalarField> K1(const volScalarField& p) const;
    tmp<volScalarField> K2(const volScalarField& p) const;

    // Accessors
    const word& phaseL() const { return phaseL_; }
    const word& phaseG() const { return phaseG_; }
    // Destructor 
    virtual ~VLEConstant();

    // Member Functions

    // activity coefficient 
    tmp<scalarField> gamma(const label patchi) const;
    tmp<volScalarField> gamma(const volScalarField& T) const;

    // fugacity coefficient of the liquid phase
    tmp<scalarField> phiL (const label patchi) const;
    tmp<volScalarField> phiL (const volScalarField& p, const volScalarField& T) const;

    // fugacity coefficient of the Gas phase
    tmp<scalarField> phiG (const label patchi) const;
    tmp<volScalarField> phiG (const volScalarField& p, const volScalarField& T) const;

    // poynting correction
    tmp<scalarField> Poy  (const label patchi) const;
    tmp<volScalarField> Poy (const volScalarField& p, const volScalarField& T) const;

};
} // End namespace Foam

#endif
