/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "vegetationModel.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace vegetation
    {
        defineTypeNameAndDebug(vegetationModel, 0);
        defineRunTimeSelectionTable(vegetationModel, T);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::vegetation::vegetationModel::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "vegetationProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


void Foam::vegetation::vegetationModel::initialise()
{
    if (vegetation_)
    {
        solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vegetation::vegetationModel::vegetationModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "vegetationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    vegetation_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true)
{
    initialise();
}


Foam::vegetation::vegetationModel::vegetationModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary(createIOobject(T.mesh())),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    vegetation_(lookupOrDefault("vegetation", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true)
{
    if (readOpt() == IOobject::NO_READ)
    {
        vegetation_ = false;
    }

    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::vegetation::vegetationModel::~vegetationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vegetation::vegetationModel::read()
{
    if (regIOobject::read())
    {
        lookup("vegetation") >> vegetation_;
        coeffs_ = subOrEmptyDict(type() + "Coeffs");

        solverFreq_ = lookupOrDefault<label>("solverFreq", 1);
        solverFreq_ = max(1, solverFreq_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::vegetation::vegetationModel::correct
(
    volVectorField& U_,
    volScalarField& T_, 
    volScalarField& q_
)
{
    if (!vegetation_)
    {
        return;
    }

    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate(U_, T_, q_);
        firstIter_ = false;
    }
}

// ************************************************************************* //
