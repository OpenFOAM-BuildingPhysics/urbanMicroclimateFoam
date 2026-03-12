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

#include "noGrass.H"
#include "physicoChemicalConstants.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace grass
    {
        defineTypeNameAndDebug(noGrass, 0);
        addToGrassRunTimeSelectionTables(noGrass);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::grass::noGrass::noGrass(const volScalarField& T)
:
    grassModel(T)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::grass::noGrass::~noGrass()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::grass::noGrass::read()
{
    return grassModel::read();
}


void Foam::grass::noGrass::calculate
(
    const volScalarField& T_, 
    const volScalarField& w_,
    const volVectorField& U_
)
{
    // Do nothing
}

Foam::tmp<Foam::volScalarField> Foam::grass::noGrass::Sh() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimensionSet(1,-1,-3,0,0,0,0),
                0.0
            )
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::grass::noGrass::Cf() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Cf",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimensionSet(0,-1,0,0,0,0,0),
                0.0
            )
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::grass::noGrass::Sw() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Sw",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimensionSet(1,-3,-1,0,0,0,0),
                0.0
            )
        )
    );
}


// ************************************************************************* //
