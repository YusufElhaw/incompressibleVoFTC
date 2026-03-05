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
#include "massDiffusivity.H"
#include "fvMesh.H"

namespace Foam
{

const IOdictionary& massDiffusivity::readThermophysicalTransportDict
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
bool massDiffusivity::checkPhaseTurbulence(const word& phaseName) const
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

massDiffusivity::massDiffusivity
(
    const incompressibleInterPhaseTransportModelTC& turbulence
)
:
    turbulence_(turbulence),
    mixture_(turbulence.mixture_),

    Dm1_(dimensionedScalar("Dm1", dimKinematicViscosity, 0)),
    Dm2_(dimensionedScalar("Dm2", dimKinematicViscosity, 0)),

    Prt1_(0.85),
    Prt2_(0.85),

    Let1_(1.0),
    Let2_(1.0),

    turb1_(false),
    turb2_(false)
{   
    // phase names from mixture
    const word& phase1Name = mixture_.phase1Name();
    const word& phase2Name = mixture_.phase2Name();

    // read thermo transport dicts (must exist)
    const dictionary& dictPhase1 = readThermophysicalTransportDict(phase1Name);  // file thermophysicalTransport.<phase1>
    const dictionary& dictPhase2 = readThermophysicalTransportDict(phase2Name);  // file thermophysicalTransport.<phase2>


    Dm1_ = dimensionedScalar("Dm1", dimKinematicViscosity, dictPhase1.lookup("Dm"));
    Dm2_ = dimensionedScalar("Dm2", dimKinematicViscosity, dictPhase2.lookup("Dm"));

    Info<<"Read molecular Diffusivty of "<< phase1Name << " phase Dm: "
        << Dm1_.value() <<" "<<Dm1_.dimensions()<<endl;
    Info<<"Read molecular Diffusivty of "<< phase2Name << " phase Dm: "
        << Dm2_.value() <<" "<<Dm2_.dimensions()<<endl;
    // turbulence flags from momentumTransport.<phaseName>
    turb1_ = checkPhaseTurbulence(phase1Name);
    turb2_ = checkPhaseTurbulence(phase2Name);

    // Prt/Let MUST exist if turbulent
    if (turb1_)
    {
        if (!dictPhase1.found("Prt"))
        {
            FatalErrorInFunction
                << "Phase " << phase1Name
                << " is turbulent, but turbulent Prandtl number 'Prt' is missing in constant/thermophysicalTransport."
                << phase1Name << nl << exit(FatalError);
        }
        if (!dictPhase1.found("Let"))
        {
            FatalErrorInFunction
                << "Phase " << phase1Name
                << " is turbulent, but turbulent Lewis number 'Let' is missing in constant/thermophysicalTransport."
                << phase1Name << nl << exit(FatalError);
        }

        dictPhase1.lookup("Prt") >> Prt1_;
        dictPhase1.lookup("Let") >> Let1_;

        Info<< "Read the turbulent Prandtl number of phase " << phase1Name
            << " Prt: " << Prt1_ << nl;
        Info<< "Read the turbulent Lewis number of phase " << phase1Name
            << " Let: " << Let1_ << nl;
    }

    if (turb2_)
    {
        if (!dictPhase2.found("Prt"))
        {
            FatalErrorInFunction
                << "Phase " << phase2Name
                << " is turbulent, but turbulent Prandtl number 'Prt' is missing in constant/thermophysicalTransport."
                << phase2Name << nl << exit(FatalError);
        }
        if (!dictPhase2.found("Let"))
        {
            FatalErrorInFunction
                << "Phase " << phase2Name
                << " is turbulent, but turbulent Lewis number 'Let' is missing in constant/thermophysicalTransport."
                << phase2Name << nl << exit(FatalError);
        }

        dictPhase2.lookup("Prt") >> Prt2_;
        dictPhase2.lookup("Let") >> Let2_;

        Info<< "Read the turbulent Prandtl number of phase " << phase2Name
            << " Prt: " << Prt2_ << nl;
        Info<< "Read the turbulent Lewis number of phase " << phase2Name
            << " Let: " << Let2_ << nl;

    }

}
tmp<volScalarField> massDiffusivity::DEff() const
{
//    const fvMesh& mesh = mixture_.mesh();

    volScalarField D1("D1", mixture_.alpha1()* Dm1_);
    volScalarField D2("D2", mixture_.alpha2()* Dm2_);

    if (turb1_)
    {
        D1 += mixture_.alpha1()*turbulence_.turbulence1_->nut()/(Prt1_*Let1_);
    }
    if (turb2_)
    {
        D2 += mixture_.alpha2()*turbulence_.turbulence2_->nut()/(Prt2_*Let2_);
    }

    return D1 + D2;
}



tmp<scalarField> massDiffusivity::DEff(const label patchi) const
{
    const scalarField a1 = mixture_.alpha1().boundaryField()[patchi];
    const scalarField a2 = mixture_.alpha2().boundaryField()[patchi];

    scalarField D1(a1.size(), Dm1_.value());
    scalarField D2(a1.size(), Dm2_.value());

    if (turb1_)
    {
        D1 += turbulence_.turbulence1_->nut(patchi)/(Prt1_*Let1_);
    }
    if (turb2_)
    {
        D2 += turbulence_.turbulence2_->nut(patchi)/(Prt2_*Let2_);
    }

    tmp<scalarField> t(new scalarField(a1*D1 + a2*D2));
    return t;
}

} // End namespace Foam
