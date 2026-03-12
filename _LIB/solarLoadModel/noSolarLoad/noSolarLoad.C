/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "noSolarLoad.H"
#include "physicoChemicalConstants.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace solarLoad
    {
        defineTypeNameAndDebug(noSolarLoad, 0);
        addToSolarLoadRunTimeSelectionTables(noSolarLoad);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarLoad::noSolarLoad::noSolarLoad(const volScalarField& T)
:
    solarLoadModel(T)
{}


Foam::solarLoad::noSolarLoad::noSolarLoad
(
    const dictionary& dict,
    const volScalarField& T
)
:
    solarLoadModel(T)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solarLoad::noSolarLoad::~noSolarLoad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solarLoad::noSolarLoad::read()
{
    return solarLoadModel::read();
}


void Foam::solarLoad::noSolarLoad::calculate()
{
    // Do nothing
}


Foam::tmp<Foam::volScalarField> Foam::solarLoad::noSolarLoad::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "Rp",
                constant::physicoChemical::sigma.dimensions()/dimLength,
                0.0
            )
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::solarLoad::noSolarLoad::Ru() const
{
    return tmp<volScalarField::Internal>
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "Ru",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "Ru", dimMass/dimLength/pow3(dimTime), 0.0
            )
        )
    );
}


// ************************************************************************* //
