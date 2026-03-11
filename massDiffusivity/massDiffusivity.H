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
    Foam::MassDiffusivity

Description
    in constant use 2 files thermophysicalTransport.phase1 and 
    thermophysicalTransport.phase2. 
    exmaple of usage: 

    Dm		9e-6;   // molecular diffusion coeffcient: always obilgatory 
	Prt     1;      // turbulent Prandtl number 
                    // obilgatory only if the Phase was LES/RAS in momentumTransport
	Let		1;      // turbulent Prandtl number
                    // obilgatory only if the Phase was LES/RAS in momentumTransport

SourceFiles
    MassDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef massDiffusivity_H
#define massDiffusivity_H

#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "incompressibleInterPhaseTransportModelTC.H"
#include "incompressibleTwoPhaseVoFMixtureTC.H"

namespace Foam
{

class massDiffusivity
{
    const incompressibleInterPhaseTransportModelTC& turbulence_;
    const incompressibleTwoPhaseVoFMixtureTC& mixture_;

    // parameters per-phase (read once)
    dimensionedScalar Dm1_, Dm2_;  // [m2/s]
    scalar Prt1_, Prt2_;
    scalar Let1_, Let2_;

    //  turbulence per-phase flags (read once)
    bool turb1_, turb2_;

    // dictionaries (cached in registry)
    const IOdictionary& readThermophysicalTransportDict(const word& phaseName) const;
    
    // checking from momentumTransport both Phases for turbulence
    bool checkPhaseTurbulence(const word& phaseName) const;

public:
    massDiffusivity(const incompressibleInterPhaseTransportModelTC& turb);

    tmp<volScalarField> D1Eff() const;
    tmp<scalarField> D1Eff(const label patchi) const;

    tmp<volScalarField> D2Eff() const;
    tmp<scalarField> D2Eff(const label patchi) const;


    tmp<volScalarField> DEff() const;
    tmp<scalarField> DEff(const label patchi) const;
};

} // End namespace Foam

#endif

