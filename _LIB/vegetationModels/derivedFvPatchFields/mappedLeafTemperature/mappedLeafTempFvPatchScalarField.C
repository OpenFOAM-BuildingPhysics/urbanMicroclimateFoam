/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "mappedLeafTempFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    mappedPatchFieldBase<scalar>(this->mapper(p, iF), *this)
{}


Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    mappedPatchFieldBase<scalar>(this->mapper(p, iF), *this, dict)
{}


Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const mappedLeafTempFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    mappedPatchFieldBase<scalar>(this->mapper(p, iF), *this, ptf)
{}


Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const mappedLeafTempFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    mappedPatchFieldBase<scalar>(ptf)
{}


Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const mappedLeafTempFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    mappedPatchFieldBase<scalar>(this->mapper(this->patch(), iF), *this, ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::mappedPatchBase& Foam::mappedLeafTempFvPatchScalarField::mapper
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
{
    if (!isA<mappedPatchBase>(p.patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.patch().name()
            << " of field " << iF.name()
            << " in file " << iF.objectPath()
            << exit(FatalError);
    }
    return refCast<const mappedPatchBase>(p.patch());
}


void Foam::mappedLeafTempFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //this->operator==(this->mappedField());

    //The default nearestCell mapping will sometimes map from a cell with Tl = 0
    //Replaced with the part below based on Foam::meshSearch::findNearestCellLinear function
    //but we are searching only within the cells where Tl is defined

    const fvMesh& airMesh = this->sampleField().mesh();
    const volScalarField& Tl = airMesh.lookupObject<volScalarField>("Tl");

    scalarField& Tp = *this;

    List<List<point>> vegCellCentres(Pstream::nProcs());
    List<List<scalar>> vegCellValues(Pstream::nProcs());
    forAll(Tl.internalField(), cellI)
    {
        if (Tl.internalField()[cellI] > 0)
        {
            vegCellCentres[Pstream::myProcNo()].append(airMesh.cellCentres()[cellI]);
            vegCellValues[Pstream::myProcNo()].append(Tl.internalField()[cellI]);
        }
    }
    Pstream::gatherList(vegCellCentres);
    Pstream::scatterList(vegCellCentres);

    Pstream::gatherList(vegCellValues);
    Pstream::scatterList(vegCellValues);

    List<Tuple2<scalar, scalar>> nearest(Tp.size());
    //Tuple2 comprising 0: sqr(distance), 1: value

    forAll(nearest, i)
    {
        //initialize
        nearest[i].first() = great;
        nearest[i].second() = great;

        //findNearer

        forAll (vegCellCentres, proci)
        {
            forAll(vegCellCentres[proci], pointi)
            {
                const point& location = this->patch().Cf()[i];
                scalar distSqr = magSqr(vegCellCentres[proci][pointi] - location);

                if (distSqr < nearest[i].first())
                {
                    nearest[i].first() = distSqr;
                    nearest[i].second() = vegCellValues[proci][pointi];
                }
            }
        }

        //update value
        Tp[i] = nearest[i].second();
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::mappedLeafTempFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    mappedPatchFieldBase<scalar>::write(os);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mappedLeafTempFvPatchScalarField
    );
}

// ************************************************************************* //
