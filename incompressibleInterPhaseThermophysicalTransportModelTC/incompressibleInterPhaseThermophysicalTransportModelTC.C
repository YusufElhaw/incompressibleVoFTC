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

#include "incompressibleInterPhaseThermophysicalTransportModelTC.H"
#include "incompressibleInterPhaseTransportModelTC.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleInterPhaseThermophysicalTransportModelTC::
incompressibleInterPhaseThermophysicalTransportModelTC
(
    const incompressibleInterPhaseTransportModelTC& turbulence
)
:
    thermophysicalTransportModel(turbulence.mixture_.mesh(), word::null),
    turbulence_(turbulence),
    massDiffusivity_(turbulence)

{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
        // * * * * * * * * * * * Heat Transfer * * * * * * * * * * * //
bool Foam::incompressibleInterPhaseThermophysicalTransportModelTC::read()
{
    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::kappaEff() const
{

    const incompressibleTwoPhaseVoFMixtureTC& mixture_ =
        turbulence_.mixture_;

    if (turbulence_.twoPhaseTransport_)
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *turbulence_.turbulence1_->nut()
            )
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *turbulence_.turbulence2_->nut()
            );
    }
    else
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *turbulence_.turbulence_->nut()
            )
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *turbulence_.turbulence_->nut()
            );
    }
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::kappaEff
(
    const label patchi
) const
{

    const incompressibleTwoPhaseVoFMixtureTC& mixture_ =
        turbulence_.mixture_;

    if (turbulence_.twoPhaseTransport_)
    {
        return
            mixture_.alpha1().boundaryField()[patchi]
           *(
                mixture_.thermo1().kappa().boundaryField()[patchi]
              + mixture_.thermo1().rho(patchi)
               *mixture_.thermo1().Cp().boundaryField()[patchi]
               *turbulence_.turbulence1_->nut(patchi)
            )
          + mixture_.alpha2().boundaryField()[patchi]
           *(
                mixture_.thermo2().kappa().boundaryField()[patchi]
              + mixture_.thermo2().rho(patchi)
               *mixture_.thermo2().Cp().boundaryField()[patchi]
               *turbulence_.turbulence2_->nut(patchi)
            );
    }
    else
    {
        return
            mixture_.alpha1().boundaryField()[patchi]
           *(
                mixture_.thermo1().kappa().boundaryField()[patchi]
              + mixture_.thermo1().rho(patchi)
               *mixture_.thermo1().Cp().boundaryField()[patchi]
               *turbulence_.turbulence_->nut(patchi)
            )
          + mixture_.alpha2().boundaryField()[patchi]
           *(
                mixture_.thermo2().kappa().boundaryField()[patchi]
              + mixture_.thermo2().rho(patchi)
               *mixture_.thermo2().Cp().boundaryField()[patchi]
               *turbulence_.turbulence_->nut(patchi)
            );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::alphaEff() const
{

    const incompressibleTwoPhaseVoFMixtureTC& mixture_ =
        turbulence_.mixture_;

    if (turbulence_.twoPhaseTransport_)
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *turbulence_.turbulence1_->nut()
            )/mixture_.thermo1().Cv()
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *turbulence_.turbulence2_->nut()
            )/mixture_.thermo2().Cv();
    }
    else
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *turbulence_.turbulence_->nut()
            )/mixture_.thermo1().Cv()
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *turbulence_.turbulence_->nut()
            )/mixture_.thermo2().Cv();
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::q() const
{
    return surfaceScalarField::New
    (
        "q",
        -fvc::interpolate(kappaEff())
        *fvc::snGrad(turbulence_.mixture_.T())
    );
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::q
(
    const label patchi
) const
{
    return
      - kappaEff(patchi)
       *turbulence_.mixture_.T().boundaryField()[patchi].snGrad();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::qCorr
(
    const label patchi
) const
{
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::divq
(
    volScalarField& he
) const
{
    NotImplemented;

    return tmp<fvScalarMatrix>(nullptr);
}


        // * * * * * * * * * * * Mass Transfer * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::DEff() const
{
    return massDiffusivity_.DEff();
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::D1Eff() const
{
    return massDiffusivity_.D1Eff();
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::D2Eff() const
{
    return massDiffusivity_.D2Eff();
}

Foam::tmp<Foam::scalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::DEff
(
    const label patchi
) const
{
   return massDiffusivity_.DEff(patchi);
}

Foam::tmp<Foam::scalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::D1Eff
(
    const label patchi
) const
{
   return massDiffusivity_.D1Eff(patchi);
}

Foam::tmp<Foam::scalarField>
Foam::incompressibleInterPhaseThermophysicalTransportModelTC::D2Eff
(
    const label patchi
) const
{
   return massDiffusivity_.D2Eff(patchi);
}
void Foam::incompressibleInterPhaseThermophysicalTransportModelTC::predict()
{}


void Foam::incompressibleInterPhaseThermophysicalTransportModelTC::correct()
{}


// ************************************************************************* //
