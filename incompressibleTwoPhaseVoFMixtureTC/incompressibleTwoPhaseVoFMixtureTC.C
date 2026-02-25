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

\*---------------------------------------------------------------------------*/

#include "incompressibleTwoPhaseVoFMixtureTC.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseVoFMixtureTC, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseVoFMixtureTC::incompressibleTwoPhaseVoFMixtureTC
(
    const fvMesh& mesh
)
:
    twoPhaseVoFMixture(mesh),
    totalInternalEnergy_
    (
        lookupOrDefault<Switch>("totalInternalEnergy", true)
    ),

    p_
    (
        IOobject
        (
            "p",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    T_
    (
        IOobject
        (
            "T",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    species1Name_(""),
    species2Name_(""),

    s1Phase1_(nullptr),
    s2Phase1_(nullptr),
    s1Phase2_(nullptr),
    s2Phase2_(nullptr),

    s1IndexThermo1_(-1),
    s2IndexThermo1_(-1),
    s1IndexThermo2_(-1),
    s2IndexThermo2_(-1),
    
    //residualAlpha1_("residualAlpha1", dimless, 0.0),
    //residualAlpha2_("residualAlpha2", dimless, 0.0),

    thermo1_(nullptr),
    thermo2_(nullptr),
    ////
    nuModel1_(viscosityModel::New(mesh, phase1Name())),
    nuModel2_(viscosityModel::New(mesh, phase2Name())),

    rho1_("rho", dimDensity, nuModel1_()),
    rho2_("rho", dimDensity, nuModel2_()),

    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 0)
    ),

    nu_
    (
        IOobject
        (
            "nu",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimKinematicViscosity, 0),
        calculatedFvPatchScalarField::typeName
    )
{
    
    {
        volScalarField T1
        (
            IOobject
            (
                IOobject::groupName("T", phase1Name()),
                mesh.time().name(),
                mesh
            ),
            T_,
            calculatedFvPatchScalarField::typeName
        );
        T1.write();
    }

    {
        volScalarField T2
        (
            IOobject
            (
                IOobject::groupName("T", phase2Name()),
                mesh.time().name(),
                mesh
            ),
            T_,
            calculatedFvPatchScalarField::typeName
        );
        T2.write();
    }

/////////
    
    const word timeName(mesh.time().name());

    // Note: we're writing files to be read in immediately afterwards.
    //       Avoid any thread-writing problems.
    // fileHandler().flush();

    thermo1_ = rhoFluidMulticomponentThermo::New(mesh, phase1Name());
    
    // Print the mole Weigts
    forAll(thermo1_->species(), i) 
    {
    Info<< "Reading "<<phase1Name()<<" phase species "<< i+1 <<" (" 
    << thermo1_->species()[i] << "): mole mass : "
    << thermo1_->Wi(i).value() <<" kg/kmol " << nl << endl;
    }
  
    thermo2_ = rhoFluidMulticomponentThermo::New(mesh, phase2Name());
    
    // Print the mole Weigts
    forAll(thermo2_->species(), i) 
    {
    Info<< "Reading "<<phase2Name()<<" phase species "<< i+1 <<" (" 
    << thermo2_->species()[i] << "): mole mass : "
    << thermo2_->Wi(i).value() << " kg/kmol " << nl << endl;
    }   
    
    species1Name_= thermo1_->species()[0];
    species2Name_= thermo1_->species()[1];
    
    //thermo1_->validate(phase1Name(), "e");     // making sure that physicalProperties of the Phase is on internalEnergy or Enthalpy
    //thermo2_->validate(phase2Name(), "e");     // making sure that physicalProperties of the Phase is on internalEnergy or Enthalpy
    
  
    // --- Map species indices in each thermo ---
    const speciesTable& sp1=thermo1_->species();
    s1IndexThermo1_ = sp1[species1Name_];
    s2IndexThermo1_ = sp1[species2Name_];

    const speciesTable& sp2=thermo2_->species();
    s1IndexThermo2_ = sp2[species1Name_];
    s2IndexThermo2_ = sp2[species2Name_];

    if (sp1.size() != 2 || sp2.size() != 2)
    {
        FatalErrorInFunction
            << "Expected exactly 2 species in both phases, got "
            << sp1.size() << " in  physicalProperties." << phase1Name()<< "and " 
            << sp2.size() << " in  physicalProperties." << phase2Name()<< nl
            << exit(FatalError);
    }
    
    s1Phase1_.reset
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName(species1Name_, phase1Name()),
                timeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    s2Phase1_.reset
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName(species2Name_, phase1Name()),
                timeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    s1Phase2_.reset
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName(species1Name_, phase2Name()),
                timeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    s2Phase2_.reset
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName(species2Name_, phase2Name()),
                timeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
        
    correct();  
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::incompressibleTwoPhaseVoFMixtureTC::correctThermo()
{
    thermo1_->T() = T_;
    thermo1_->he() = thermo1_->he(p_, T_);
    thermo1_->correct();

    thermo2_->T() = T_;
    thermo2_->he() = thermo2_->he(p_, T_);
    thermo2_->correct();
}

void Foam::incompressibleTwoPhaseVoFMixtureTC::correctComposition()
{
    // push mass fractions into each phase-thermo
    thermo1_->Y(s1IndexThermo1_) = s1Phase1_();
    thermo1_->Y(s2IndexThermo1_) = s2Phase1_();

    thermo2_->Y(s1IndexThermo2_) = s1Phase2_();
    thermo2_->Y(s2IndexThermo2_) = s2Phase2_();

    thermo1_->normaliseY();
    thermo2_->normaliseY();

    thermo1_->correct();
    thermo2_->correct();

}


bool Foam::incompressibleTwoPhaseVoFMixtureTC::read()
{
    if (twoPhaseVoFMixture::read())
    {
       // nuModel1_->lookup("rho") >> rho1_;
       // nuModel2_->lookup("rho") >> rho2_;

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::incompressibleTwoPhaseVoFMixtureTC::correct()
{

   // nuModel1_->correct();
   // nuModel2_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1(), scalar(0)), scalar(1))
    );

    const volScalarField limitedAlpha2 ("limitedAlpha2", scalar(1) - limitedAlpha1);

    rho_ = limitedAlpha1*thermo1().rho() + limitedAlpha2*thermo2().rho();

    // Average kinematic viscosity calculated from dynamic viscosity

    nu_ =
    ( limitedAlpha1*thermo1().rho()*thermo1().nu()      //was in incompressibleVoF "nuModel1_->nu()"
    + limitedAlpha2*thermo2().rho()*thermo2().nu()      //was in incompressibleVoF "nuModel2_->nu()"
    )/rho_;
    
}


// ************************************************************************* //
