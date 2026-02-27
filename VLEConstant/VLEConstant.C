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
    const fvMesh& mesh = mixture.mesh();

    // cache in registry (same pattern wie dein Stub)
    if (mesh.foundObject<IOdictionary>(dictName))
    {
        return mesh.lookupObject<IOdictionary>(dictName);
    }

    IOdictionary* dictPtr = new IOdictionary
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

    mesh.objectRegistry::store(dictPtr);
    return *dictPtr;
}

const dictionary& VLEConstant::coeffsDict(const word& dictName) const
{
    const IOdictionary& dict = readDictOrDie(dictName);

    if (!dict.found("VLEConstantCoeffs"))
    {
        FatalIOErrorInFunction(dict)
            << "Missing sub-dictionary 'VLEConstantCoeffs' in " << dictName
            << exit(FatalIOError);
    }

    return dict.subDict("VLEConstantCoeffs");
}

VLEConstant::VLEConstant
(
    const incompressibleInterPhaseTransportModelTC& turb,
    const incompressibleTwoPhaseVoFMixtureTC& mix,
    const word& dictName
)
:
    turbulence_(turb),
    mixture(mix),
    phaseL_("phase1"),
    phaseG_("phase2"),
    specie1_("specie1"),
    specie2_("specie2"),
   //KMin_(0.0),
   //KMax_(GREAT),
    pSat1Model_(nullptr),
    pSat2Model_(nullptr)
{
    const dictionary& coeffs = coeffsDict(dictName);

    coeffs.readIfPresent("phaseL", phaseL_);
    coeffs.readIfPresent("phaseG", phaseG_);
    coeffs.readIfPresent("specie1", specie1_);
    coeffs.readIfPresent("specie2", specie2_);
    //coeffs.readIfPresent("KMin", KMin_);
    //coeffs.readIfPresent("KMax", KMax_);

    if (!coeffs.found("pSat1") || !coeffs.found("pSat2"))
    {
        FatalIOErrorInFunction(readDictOrDie(dictName))
            << "VLEConstantCoeffs requires sub-dicts 'pSat1' and 'pSat2'\n"
            << "Each must define a saturationPressureModel, e.g. type Antoine/constant.\n"
            << exit(FatalIOError);
    }

    pSat1Model_ = saturationPressureModel::New("pSat1", coeffs);
    pSat2Model_ = saturationPressureModel::New("pSat2", coeffs);
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::VLEConstant::~VLEConstant()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// pSat 
tmp<volScalarField> VLEConstant::pSat1(const volScalarField& T) const
{
    return pSat1Model_->pSat(T);
}

tmp<volScalarField> VLEConstant::pSat2(const volScalarField& T) const
{
    return pSat2Model_->pSat(T);
}

// gamma 
// temporary as 1 

tmp<scalarField> VLEConstant::gamma(const label patchi) const
{
    const fvMesh& mesh = mixture.mesh();
    tmp<scalarField> t(new scalarField(mesh.boundary()[patchi].size(), 1.0));
    return t;
}


tmp<volScalarField> VLEConstant::gamma(const volScalarField& T) const
{
    tmp<volScalarField> t
    (
        new volScalarField
        (
            IOobject
            (
                "gamma",
                T.time().name(),
                T.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T.mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );
    return t;
}

// fugacity liquid
// temporary as 1 
tmp<scalarField> VLEConstant::phiL(const label patchi) const
{
    const fvMesh& mesh = mixture.mesh();
    tmp<scalarField> t(new scalarField(mesh.boundary()[patchi].size(), 1.0));
    return t;
}

tmp<volScalarField> VLEConstant::phiL(const volScalarField& p, const volScalarField& T) const
{
    tmp<volScalarField> t
    (
        new volScalarField
        (
            IOobject
            (
                "phiL",
                p.time().name(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p.mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );
    return t;
}

// fugacity gas
// temporary as 1 
tmp<scalarField> VLEConstant::phiG(const label patchi) const
{
    const fvMesh& mesh = mixture.mesh();
    tmp<scalarField> t(new scalarField(mesh.boundary()[patchi].size(), 1.0));
    return t;
}

tmp<volScalarField> VLEConstant::phiG(const volScalarField& p, const volScalarField& T) const
{
    tmp<volScalarField> t
    (
        new volScalarField
        (
            IOobject
            (
                "phiG",
                p.time().name(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p.mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );
    return t;
}

// Poyinting correction
// temporary as 1 
tmp<scalarField> VLEConstant::Poy(const label patchi) const
{
    const fvMesh& mesh = mixture.mesh();
    tmp<scalarField> t(new scalarField(mesh.boundary()[patchi].size(), 1.0));
    return t;
}
tmp<volScalarField> VLEConstant::Poy(const volScalarField& p, const volScalarField& T) const
{
    tmp<volScalarField> t
    (
        new volScalarField
        (
            IOobject
            (
                "Poy",
                p.time().name(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p.mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );
    return t;
}
tmp<volScalarField> VLEConstant::K1(const volScalarField& p, const volScalarField& T) const
{
    const tmp<volScalarField> tpSat = pSat1(T);
    const tmp<volScalarField> tGam  = gamma(T);
    const tmp<volScalarField> tPoy  = Poy(p, T);
    const tmp<volScalarField> tPhiL = phiL(p, T);
    const tmp<volScalarField> tPhiG = phiG(p, T);

    const dimensionedScalar pSmall(p.dimensions(), SMALL);

    tmp<volScalarField> tK
    (
        new volScalarField
        (
            IOobject
            (
                "K1",
                p.time().name(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (tpSat()*tGam()*tPoy()*(tPhiL()/max(tPhiG(), dimensionedScalar(dimless, SMALL))))
           /max(p, pSmall)
        )
    );

   // if (KMin_ > 0.0 || KMax_ < GREAT)
   // {
   //     tK.ref() = min(max(tK(), KMin_), KMax_);
   // }

    return tK;
}

tmp<volScalarField> VLEConstant::K2(const volScalarField& p, const volScalarField& T) const
{
    const tmp<volScalarField> tpSat = pSat2(T);
    const tmp<volScalarField> tGam  = gamma(T);
    const tmp<volScalarField> tPoy  = Poy(p, T);
    const tmp<volScalarField> tPhiL = phiL(p, T);
    const tmp<volScalarField> tPhiG = phiG(p, T);

    const dimensionedScalar pSmall(p.dimensions(), SMALL);

    tmp<volScalarField> tK
    (
        new volScalarField
        (
            IOobject
            (
                "K2",
                p.time().name(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (tpSat()*tGam()*tPoy()*(tPhiL()/max(tPhiG(), dimensionedScalar(dimless, SMALL))))
           /max(p, pSmall)
        )
    );

    //if (KMin_ > 0.0 || KMax_ < GREAT)
    //{
    //    tK.ref() = min(max(tK(), KMin_), KMax_);
    //}

    return tK;
}

tmp<volScalarField> VLEConstant::K1(const volScalarField& p) const
{
    return K1(p, mixture.thermo1().T());
}

tmp<volScalarField> VLEConstant::K2(const volScalarField& p) const
{
    return K2(p, mixture.thermo2().T());
}

} // End namespace Foam
