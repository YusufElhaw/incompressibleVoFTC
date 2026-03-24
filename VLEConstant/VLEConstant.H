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
    y 𝜑_G p  =  x 𝛾 𝜑_L p_sat exp{1/(RT) * ∫( v_0 )dp }
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

    Usage:
    in constant make a File with the name: "VLEProperties" 
    VLEConstantCoeffs
    {
    <Species1> // same name in physicalProperties.<phase1>
    {
        pSat1
        {
            type constant;  // or Antoine or ArdenBuck
            value 1;        // for constant write the value Always in bar.
            // ... Antoine-Koeffizienten wie in OpenFOAM saturationModels
        }
    }

    <Species2> // same name in physicalProperties.<phase1>
    { 
        pSat2 
        {
            type Antoine;   // pSat= exp(A + B/(C + T))
            A    ;
            B    ;
            C    ; 
        }
    }
    }

    Attention:
    1-  for pSat type constant, write the saturation pressure in bar.
    2-  look carefully on the units of the Antoine equation. in this solver it is 
        assumed, that the saturation pressure is in bar. 
    3-  Antoine equation calculates with the natural logarithm in the form: 
        ln(pSat) = A + B/(C + T) or 
        pSat = e^{A + B/(C + T)}
        If the constants provided are for Antoine equation in the decadic logarithm:
        log10(pSat) = A + B/(C + T) or 
        pSat= 10^{A_d + B_d/(C_d + T)} 
        you have to recalculate them as follows: 
        A = ln(10)*A_d
        B = ln(10)*B_d
        C = C_d
        and feed the A, B and C in the pSat subdictionary. 
    4-  be careful of the signs of the coeffecients A, B and C. They are all positive
        in this solver. If your equation looks like this: 
        ln(P) = A − (B / (T + C)), make sure to feed the solver with a negative B.
} 


SourceFiles
    VLEConstant.C

\*---------------------------------------------------------------------------*/

#ifndef VLEConstant_H
#define VLEConstant_H

//#include "volFields.H"
//#include "surfaceFields.H"
//#include "IOdictionary.H"
//#include "dictionary.H"
//#include "autoPtr.H"
//#include "word.H"

//#include "wordList.H"
//#include "HashTable.H"
//#include "fvMesh.H"

#include "saturationPressureModel.H"
#include "incompressibleInterPhaseTransportModelTC.H"
#include "incompressibleTwoPhaseVoFMixtureTC.H"

namespace Foam
{

class VLEConstant
{
     // Private data
    const incompressibleInterPhaseTransportModelTC& turbulence_;
    const incompressibleTwoPhaseVoFMixtureTC& mixture_;
    const fvMesh& mesh_;

    const word dictName_;
    const word speciesName_;

    autoPtr<IOdictionary> vleDictPtr_;
    autoPtr<saturationPressureModel> pSatModel_;

    // Helpers
    const IOdictionary& readDictOrDie(const word& dictName) const;
    const dictionary& speciesDictOrDie(const IOdictionary& d, const word& speciesName) const;

public:
    // Construct for a single species
    VLEConstant
    (
        const incompressibleInterPhaseTransportModelTC& turb,
        const incompressibleTwoPhaseVoFMixtureTC& mix,
        const word& speciesName,
        const word& dictName = "VLEProperties"
    );
    //- Disallow default bitwise copy construction
    VLEConstant (const VLEConstant&) = delete;
    void operator=(const VLEConstant&) = delete;

    // Destructor 
    virtual ~VLEConstant();

    // pSat(T) fields
    tmp<volScalarField> pSat(const volScalarField& T) const;

    // future work: Methode K to be changed, may be by overloading, after implemnting the fugacity, Poynting or activity
    //              so that a choose between the thermodynamic models be enabled (at best in the constant/VLEProperties document) 
    // K = pSat/p (Raoult, idealised: gamma=phiL=phiG=Poy=1)
    tmp<volScalarField> K(const volScalarField& p, const volScalarField& T) const;


    // Accessors
    //const word& phaseL() const { return mixture.thermo1().phaseName(); }
    //const word& phaseG() const { return mixture.thermo2().phaseName(); }
    
    const word& speciesName() const { return speciesName_; }

    // Member Functions

    // activity coefficient 
 //   tmp<scalarField> gamma(const label patchi) const;
    virtual tmp<scalarField>    gamma(const label patchi) const;
    virtual tmp<volScalarField> gamma(const volScalarField& Y1, const volScalarField& Y2) const;

    // fugacity coefficient of the liquid phase
    virtual tmp<scalarField>    phiL (const label patchi) const;
    virtual tmp<volScalarField> phiL (const volScalarField& p, const volScalarField& T) const;

    // fugacity coefficient of the Gas phase
    virtual tmp<scalarField>    phiG (const label patchi) const;
    virtual tmp<volScalarField> phiG (const volScalarField& p, const volScalarField& T) const;

    // poynting correction
    virtual tmp<scalarField>    Poy (const label patchi) const;
    virtual tmp<volScalarField> Poy (const volScalarField& p, const volScalarField& T) const;

};
} // End namespace Foam

#endif
