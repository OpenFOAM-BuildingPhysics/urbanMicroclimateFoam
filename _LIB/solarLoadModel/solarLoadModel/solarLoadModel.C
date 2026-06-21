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

#include "solarLoadModel.H"
#include "solarLoadAbsorptionEmissionModel.H"
#include "solarLoadScatterModel.H"
#include "fvmSup.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace solarLoad
    {
        defineTypeNameAndDebug(solarLoadModel, 0);
        defineRunTimeSelectionTable(solarLoadModel, T);
        defineRunTimeSelectionTable(solarLoadModel, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::solarLoad::solarLoadModel::createIOobject
(
    const fvMesh& mesh
) const
{
    typeIOobject<IOdictionary> io
    (
        IOobject
        (
            "solarLoadProperties",
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


void Foam::solarLoad::solarLoadModel::initialise()
{
    if (solarLoad_)
    {
        solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

        solarLoadAbsorptionEmission_.reset
        (
            solarLoadAbsorptionEmissionModel::New(*this, mesh_).ptr()
        );

        solarLoadScatter_.reset(solarLoadScatterModel::New(*this, mesh_).ptr());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarLoad::solarLoadModel::solarLoadModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "solarLoadProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    solarLoad_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    solarLoadAbsorptionEmission_(nullptr),
    solarLoadScatter_(nullptr)
{}


Foam::solarLoad::solarLoadModel::solarLoadModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary(createIOobject(T.mesh())),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    solarLoad_(lookupOrDefault("solarLoad", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    solarLoadAbsorptionEmission_(nullptr),
    solarLoadScatter_(nullptr)
{
    if (readOpt() == IOobject::NO_READ)
    {
        solarLoad_ = false;
    }

    initialise();
}


Foam::solarLoad::solarLoadModel::solarLoadModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "solarLoadProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    solarLoad_(lookupOrDefault("solarLoad", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    solarLoadAbsorptionEmission_(nullptr),
    solarLoadScatter_(nullptr)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::solarLoad::solarLoadModel::~solarLoadModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solarLoad::solarLoadModel::read()
{
    if (regIOobject::read())
    {
        lookup("solarLoad") >> solarLoad_;
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


void Foam::solarLoad::solarLoadModel::correct()
{
    if (!solarLoad_)
    {
        return;
    }

    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate();
        firstIter_ = false;
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::solarLoad::solarLoadModel::Sh
(
    const basicThermo& thermo,
    const volScalarField& he
) const
{
    const volScalarField Cpv(thermo.Cpv());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cpv, he)
      - Rp()*T3*(T_ - 4.0*he/Cpv)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::solarLoad::solarLoadModel::ST
(
    const dimensionedScalar& rhoCp,
    volScalarField& T
) const
{
    return
    (
        Ru()/rhoCp
      - fvm::Sp(Rp()*pow3(T)/rhoCp, T)
    );
}


const Foam::solarLoad::solarLoadAbsorptionEmissionModel&
Foam::solarLoad::solarLoadModel::solarLoadAbsorptionEmission() const
{
    if (!solarLoadAbsorptionEmission_.valid())
    {
        FatalErrorInFunction
            << "Requested solarLoad solarLoadAbsorptionEmission model, but model is "
            << "not activate" << abort(FatalError);
    }

    return solarLoadAbsorptionEmission_();
}


// ************************************************************************* //
