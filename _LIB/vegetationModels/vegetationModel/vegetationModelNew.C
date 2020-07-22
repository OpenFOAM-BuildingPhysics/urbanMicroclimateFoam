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

#include "vegetationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::vegetation::vegetationModel>
Foam::vegetation::vegetationModel::New
(
    const volScalarField& T
)
{
    IOobject vegetationIO
    (
        "vegetationProperties",
        T.time().constant(),
        T.mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );

    word modelType("none");
    if (vegetationIO.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary(vegetationIO).lookup("vegetationModel") >> modelType;
    }
    else
    {
        Info<< "Vegetation model not active: vegetationProperties not found"
            << endl;
    }

    Info<< "Selecting vegetationModel " << modelType << endl;

    TConstructorTable::iterator cstrIter =
        TConstructorTablePtr_->find(modelType);

    if (cstrIter == TConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "vegetationModel::New(const volScalarField&)"
        )   << "Unknown vegetationModel type "
            << modelType << nl << nl
            << "Valid vegetationModel types are:" << nl
            << TConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<vegetationModel>(cstrIter()(T));
}


// ************************************************************************* //
