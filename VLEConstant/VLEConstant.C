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

\*---------------------------------------------------------------------------*/

#include "VLEConstant.H"
#include "fvMesh.H"
#include "IOobject.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * constructor * * * * * * * * * * * * * * * * //

namespace Foam
{

const IOdictionary& VLEConstant::readDictOrDie(const word& dictName) const
{
    // Not used directly (we store a ptr), but kept for clarity if you extend later.
    FatalErrorInFunction
        << "readDictOrDie() should not be called without creating vleDictPtr_"
        << abort(FatalError);

    // dummy return to silence warnings (never reached)
    return *vleDictPtr_;
}


const dictionary& VLEConstant::speciesDictOrDie
(
    const IOdictionary& d,
    const word& speciesName
) const
{
    if (!d.found("VLEConstantCoeffs"))
    {
        FatalErrorInFunction
            << "Dictionary " << d.name()
            << " does not contain subDict 'VLEConstantCoeffs'."
            << nl << "File: " << d.objectPath()
            << abort(FatalError);
    }

    const dictionary& coeffs = d.subDict("VLEConstantCoeffs");

    if (!coeffs.found(speciesName))
    {
        FatalErrorInFunction
            << "VLEConstantCoeffs does not contain species subDict '"
            << speciesName << "'."
            << nl << "Available entries: " << coeffs.toc()
            << abort(FatalError);
    }

    return coeffs.subDict(speciesName);
}

VLEConstant::VLEConstant
(
    const incompressibleInterPhaseTransportModelTC& turb,
    const incompressibleTwoPhaseVoFMixtureTC& mix,
    const word& speciesName,
    const word& dictName

)
:
    turbulence_(turb),
    mixture_(mix),
    mesh_(mix.mesh()),
    dictName_(dictName),
    speciesName_(speciesName),
    vleDictPtr_(nullptr),
    pSatModel_(nullptr)
{
    vleDictPtr_.reset
    (
        new IOdictionary
        (
            IOobject
            (
                dictName_,
                mesh_.time().constant(),   // "constant"
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    const dictionary& spDict = speciesDictOrDie(*vleDictPtr_, speciesName_);

    if (!spDict.found("pSat"))
    {
        FatalErrorInFunction
            << "Species subDict '" << speciesName_
            << "' in " << vleDictPtr_->name()
            << " has no 'pSat' subDict."
            << abort(FatalError);
    }

    const dictionary& pSatDict = spDict.subDict("pSat");

    // ---- Factory: adjust here if your OF signature differs ----
    // Common variants:
    //   pSatModel_ = saturationPressureModel::New(pSatDict);
    //   pSatModel_ = saturationPressureModel::New("pSat", pSatDict);
    pSatModel_ = saturationPressureModel::New(pSatDict);
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::VLEConstant::~VLEConstant()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// pSat 
tmp<volScalarField> VLEConstant::pSat(const volScalarField& T) const
{
    tmp<volScalarField> tpSat = pSatModel_->pSat(T); // Antoine Equation or Constant [bar]
    return tpSat*1e5; // converts pSat from bar to pascal 
}

tmp<volScalarField> VLEConstant::K(const volScalarField& p, const volScalarField& T) const
{
    tmp<volScalarField> tpSat = this->pSat(T);

    tmp<volScalarField> tK
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("K", speciesName_),
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tpSat()/p
        )
    );

    // K is dimensionless if pSat and p are pressure
    return tK;
}

// gamma 
tmp<scalarField> VLEConstant::gamma(const label patchi) const
{
    NotImplemented;

    return tmp<scalarField>(nullptr);
}


tmp<volScalarField> VLEConstant::gamma(const volScalarField& Y1,const volScalarField& Y2) const
{
    NotImplemented;

    return tmp<volScalarField>(nullptr);
}

// fugacity liquid
tmp<scalarField> VLEConstant::phiL(const label patchi) const
{
    NotImplemented;

    return tmp<scalarField>(nullptr);
}

tmp<volScalarField> VLEConstant::phiL(const volScalarField& p, const volScalarField& T) const
{
     NotImplemented;

    return tmp<volScalarField>(nullptr);
}

// fugacity gas
tmp<scalarField> VLEConstant::phiG(const label patchi) const
{
    NotImplemented;

    return tmp<scalarField>(nullptr);
}

tmp<volScalarField> VLEConstant::phiG(const volScalarField& p, const volScalarField& T) const
{
    NotImplemented;

    return tmp<volScalarField>(nullptr);
}

// Poyinting correction
tmp<scalarField> VLEConstant::Poy(const label patchi) const
{
    NotImplemented;

    return tmp<scalarField>(nullptr);
}
tmp<volScalarField> VLEConstant::Poy(const volScalarField& p, const volScalarField& T) const
{
    NotImplemented;

    return tmp<volScalarField>(nullptr);
}

} // End namespace Foam
