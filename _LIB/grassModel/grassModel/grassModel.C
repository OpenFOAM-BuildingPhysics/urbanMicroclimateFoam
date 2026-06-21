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

#include "grassModel.H"
#include "volFields.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace grass
    {
        defineTypeNameAndDebug(grassModel, 0);
        defineRunTimeSelectionTable(grassModel, T);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::grass::grassModel::createIOobject
(
    const fvMesh& mesh
) const
{
    // typeIOobject<IOdictionary>::headerOk() searches the global (case-level)
    // path in parallel; a plain IOobject probes processorN/constant and fails,
    // selecting NO_READ on every rank so the IOdictionary stays empty and the
    // <model>Coeffs subdict loses its keys, triggering "keyword ... undefined".
    typeIOobject<IOdictionary> io
    (
        IOobject
        (
            "grassProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (io.headerOk())
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


void Foam::grass::grassModel::initialise()
{
    if (grass_)
    {
        solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::grass::grassModel::grassModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "grassProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    grass_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true)
{
    initialise();
}


Foam::grass::grassModel::grassModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary(createIOobject(T.mesh())),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    grass_(lookupOrDefault("grass", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true)
{
    if (readOpt() == IOobject::NO_READ)
    {
        grass_ = false;
    }

    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::grass::grassModel::~grassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::grass::grassModel::read()
{
    if (regIOobject::read())
    {
        lookup("grass") >> grass_;
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


void Foam::grass::grassModel::correct
(
    const volScalarField& T_, 
    const volScalarField& w_,
    const volVectorField& U_
)
{
    if (!grass_)
    {
        return;
    }

    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate(T_, w_, U_);
        firstIter_ = false;
    }
}

// ************************************************************************* //
