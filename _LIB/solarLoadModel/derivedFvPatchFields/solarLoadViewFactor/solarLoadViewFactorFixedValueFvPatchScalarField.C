/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "solarLoadViewFactorFixedValueFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
solarLoadViewFactorFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    solarRadiationCoupledBase(p, dict),
    qso_("qso", dimPower/dimArea, dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value",iF.dimensions(), dict, p.size())
        );

    }
    else
    {
         fvPatchScalarField::operator=(0.0);
    }
}


Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
solarLoadViewFactorFixedValueFvPatchScalarField
(
    const solarLoadViewFactorFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    solarRadiationCoupledBase
    (
        patch(),
        ptf.albedoMethod(),
        ptf.albedo_,
        mapper
    ),
    qso_(mapper(ptf.qso_))
{}


Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
solarLoadViewFactorFixedValueFvPatchScalarField
(
    const solarLoadViewFactorFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    solarRadiationCoupledBase
    (
        ptf.patch(),
        ptf.albedoMethod(),
        ptf.albedo_
    ),
    qso_(ptf.qso_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    fixedValueFvPatchScalarField::map(ptf, mapper);
    solarRadiationCoupledBase::map(ptf, mapper);
    const solarLoadViewFactorFixedValueFvPatchScalarField& mrptf =
        refCast<const solarLoadViewFactorFixedValueFvPatchScalarField>(ptf);

    mapper(qso_, mrptf.qso_);
}


void Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    fixedValueFvPatchScalarField::reset(ptf);
    solarRadiationCoupledBase::reset(ptf);
    const solarLoadViewFactorFixedValueFvPatchScalarField& mrptf =
        refCast<const solarLoadViewFactorFixedValueFvPatchScalarField>(ptf);

    qso_.reset(mrptf.qso_);
}


void Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Do nothing

    if (debug)
    {
        scalar Q = gSum((*this)*patch().magSf());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << " heat transfer rate:" << Q
            << " wall radiative heat flux "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::solarLoad::solarLoadViewFactorFixedValueFvPatchScalarField::
write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    solarRadiationCoupledBase::write(os);
    writeEntry(os, "qso", qso_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solarLoad
{
    makePatchTypeField
    (
        fvPatchScalarField,
        solarLoadViewFactorFixedValueFvPatchScalarField
    );
}
}


// ************************************************************************* //
