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
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solarLoad::solarLoadModel>
Foam::solarLoad::solarLoadModel::New
(
    const volScalarField& T
)
{
    typeIOobject<IOdictionary> solIO
    (
        IOobject
        (
            "solarLoadProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    word modelType("none");
    if (solIO.headerOk())
    {
        IOdictionary(solIO).lookup("solarLoadModel") >> modelType;
    }
    else
    {
        Info<< "solarLoad model not active: solarLoadProperties not found"
            << endl;
    }

    Info<< "Selecting solarLoadModel " << modelType << endl;

    TConstructorTable::iterator cstrIter =
        TConstructorTablePtr_->find(modelType);

    if (cstrIter == TConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown solarLoadModel type "
            << modelType << nl << nl
            << "Valid solarLoadModel types are:" << nl
            << TConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<solarLoadModel>(cstrIter()(T));
}


Foam::autoPtr<Foam::solarLoad::solarLoadModel>
Foam::solarLoad::solarLoadModel::New
(
    const dictionary& dict,
    const volScalarField& T
)
{
    const word modelType(dict.lookup("solarLoadModel"));

    Info<< "Selecting solarLoadModel " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown solarLoadModel type "
            << modelType << nl << nl
            << "Valid solarLoadModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<solarLoadModel>(cstrIter()(dict, T));
}


// ************************************************************************* //
