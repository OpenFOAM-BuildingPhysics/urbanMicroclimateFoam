/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "solarLoadAbsorptionEmissionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace solarLoad
    {
        defineTypeNameAndDebug(solarLoadAbsorptionEmissionModel, 0);
        defineRunTimeSelectionTable(solarLoadAbsorptionEmissionModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarLoad::solarLoadAbsorptionEmissionModel::solarLoadAbsorptionEmissionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::solarLoad::solarLoadAbsorptionEmissionModel::~solarLoadAbsorptionEmissionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::solarLoad::solarLoadAbsorptionEmissionModel::a(const label bandI) const
{
    return aDisp(bandI) + aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::solarLoad::solarLoadAbsorptionEmissionModel::aCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "aCont",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::solarLoad::solarLoadAbsorptionEmissionModel::aDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "aDisp",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::solarLoad::solarLoadAbsorptionEmissionModel::e(const label bandI) const
{
    return eDisp(bandI) + eCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::solarLoad::solarLoadAbsorptionEmissionModel::eCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "eCont",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::solarLoad::solarLoadAbsorptionEmissionModel::eDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "eDisp",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::solarLoad::solarLoadAbsorptionEmissionModel::E(const label bandI) const
{
    return EDisp(bandI) + ECont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::solarLoad::solarLoadAbsorptionEmissionModel::ECont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ECont",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::solarLoad::solarLoadAbsorptionEmissionModel::EDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "EDisp",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}


Foam::label Foam::solarLoad::solarLoadAbsorptionEmissionModel::nBands() const
{
    return pTraits<label>::one;
}


const Foam::Vector2D<Foam::scalar>&
Foam::solarLoad::solarLoadAbsorptionEmissionModel::bands(const label n) const
{
    return Vector2D<scalar>::one;
}


bool Foam::solarLoad::solarLoadAbsorptionEmissionModel::isGrey() const
{
    return false;
}


Foam::tmp<Foam::volScalarField>
Foam::solarLoad::solarLoadAbsorptionEmissionModel::addIntensity
(
    const label rayI,
    const volScalarField& ILambda
) const
{
    return ILambda;
}


void Foam::solarLoad::solarLoadAbsorptionEmissionModel::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aj
) const
{
    a = this->a();
    aj[0] =  a;
}


// ************************************************************************* //
