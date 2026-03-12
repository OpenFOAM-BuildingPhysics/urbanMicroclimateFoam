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
//v8: #include "mappedFvPatchBaseBase.H"
#include "mappedInternalPatchBase.H" //v12: mappedInternal uses separate class hierarchy

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    fieldName_(iF.name())
    //v8: mapperPtr_(nullptr)
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
    fieldName_(dict.lookupOrDefault<word>("field", iF.name()))
    //v8: mapperPtr_(mappedPatchBase::specified(dict) ? new mappedPatchBase(...) : nullptr)
{
    //v12: mappedInternal uses mappedInternalPatchBase, not mappedPatchBase
    if (!isA<mappedInternalPatchBase>(p.patch()))
    {
        FatalIOErrorInFunction(dict)
            << "Field " << iF.name() << " on patch " << p.name()
            << " is not of mappedInternal type"
            << exit(FatalIOError);
    }
}


Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const mappedLeafTempFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_)
    //v8: mapperPtr_(ptf.mapperPtr_.valid() ? new mappedPatchBase(...) : nullptr)
{}


Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const mappedLeafTempFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    fieldName_(ptf.fieldName_)
    //v8: mapperPtr_(ptf.mapperPtr_.valid() ? new mappedPatchBase(...) : nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//v8: const Foam::mappedPatchBase& Foam::mappedLeafTempFvPatchScalarField::mapper() const
//v8: {
//v8:     return
//v8:         mapperPtr_.valid()
//v8:       ? mapperPtr_()
//v8:       : mappedPatchBase::getMap(this->patch().patch());
//v8: }

//v12: mappedInternal uses mappedInternalPatchBase (cell-based, not patch-based)
const Foam::mappedInternalPatchBase&
Foam::mappedLeafTempFvPatchScalarField::mapper() const
{
    return refCast<const mappedInternalPatchBase>(this->patch().patch());
}


void Foam::mappedLeafTempFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Get the mapper and the neighbouring mesh
    //v8: const mappedPatchBase& mpp = this->mapper();
    //v8: const polyMesh& nbrMesh = mpp.sampleMesh();
    //v8: const fvPatch& nbrPatch =
    //v8:     refCast<const fvMesh>(nbrMesh).boundary()[mpp.samplePolyPatch().index()];
    //v8: scalarField tnbrIntFld = nbrPatch.lookupPatchField<volScalarField, scalar>(fieldName_);
    //v8: mpp.distribute(tnbrIntFld);

    //v12: mappedInternal distributes internal cell values to this patch
    const mappedInternalPatchBase& mipb = this->mapper();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mipb.nbrMesh());
    const volScalarField& nbrField =
        nbrMesh.lookupObject<volScalarField>(fieldName_);

    // Distribute cell values from air region to vegetation patch faces
    tmp<scalarField> tnbrIntFld = mipb.distribute(nbrField);

    this->operator==(tnbrIntFld);

    const fvMesh& airMesh = db().time().lookupObject<fvMesh>("air");
    const volScalarField& Tl = airMesh.lookupObject<volScalarField>("Tl");
    const volScalarField& LAD = airMesh.lookupObject<volScalarField>("LAD");
    scalar Tl_avg = gSum(Tl.primitiveField()*LAD.primitiveField())/gSum(LAD.primitiveField());

    scalarField& Tp = *this;
    forAll(Tp, i)
    {
        if(Tp[i] < 1.0)
        {
            Tp[i] = Tl_avg;
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::mappedLeafTempFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "field", fieldName_);
    //v8: if (mapperPtr_.valid()) { mapperPtr_->write(os); }
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
