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

namespace Foam
{

const IOdictionary& VLEConstant::readThermophysicalTransportDict
(
    const word& phaseName
) const
{
    const fvMesh& mesh = mixture_.mesh();
    const word dictName("thermophysicalTransport." + phaseName);

    if (mesh.foundObject<IOdictionary>(dictName))
    {
        return mesh.lookupObject<IOdictionary>(dictName);
    }

    IOdictionary* dictionaryPtr = new IOdictionary   // dictionary Pointer 
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    mesh.objectRegistry::store(dictionaryPtr);
    return *dictionaryPtr;
}


// Determine turbulence per phase from momentumTransport.<phaseName>
bool VLEConstant::checkPhaseTurbulence(const word& phaseName) const
{
    const fvMesh& mesh = mixture_.mesh();

    // helper to read simulationType from a given dict name
    auto readSimulationType = [&](const word& dictName, word& simType) -> bool
    {
        IOdictionary Dictionary
        (
            IOobject
            (
                dictName,
                mesh.time().constant(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );

        if (!Dictionary.headerOk())
        {
            return false;
        }

        simType = word(Dictionary.lookup("simulationType")); // must exist if file exists
        return true;
    };

    word simType("laminar");

    if (turbulence_.twoPhaseTransport_)
    {
        // try per-phase first
        if (readSimulationType("momentumTransport." + phaseName, simType))
        {
            return simType != "laminar";
        }

        // fallback to common momentumTransport
        if (readSimulationType("momentumTransport", simType))
        {
            return simType != "laminar";
        }

        return false;
    }
    else
    {
        // mixture-mode: only common momentumTransport
        if (readSimulationType("momentumTransport", simType))
        {
            return simType != "laminar";
        }

        return false;
    }
    
}

// * * * * * * * * * * * * * * * constructor * * * * * * * * * * * * * * * * //

VLEConstant::VLEConstant
(
)
  
{   
   
tmp<volScalarField> VLEConstant::pSat() const
{
   


}


tmp<scalarField> VLEConstant::gamma(const label patchi) const
{
    return 1;
}


tmp<scalarField> VLEConstant::phiL(const label patchi) const
{
    return 1;
}

tmp<scalarField> VLEConstant::phiG(const label patchi) const
{
    return 1;
}

tmp<scalarField> VLEConstant::Poy(const label patchi) const
{
    return 1;
}
} // End namespace Foam
