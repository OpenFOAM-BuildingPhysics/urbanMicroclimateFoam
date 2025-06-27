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

#include "CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fixedValueFvPatchFields.H"
#include "TableFile.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField::
CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    hcoeff_(),
    Tamb_()
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField::
CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
(
    const CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    hcoeff_(psf.hcoeff_),
    Tamb_(psf.Tamb_)
{}


CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField::
CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    hcoeff_(dict.lookupOrDefault<scalar>("hcoeff",0)),
    Tamb_(dict.lookupOrDefault<fileName>("Tamb", "none"))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField::
CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
(
    const CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    hcoeff_(psf.hcoeff_),
    Tamb_(psf.Tamb_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    scalarField& Tp = *this;
                          
    scalarField lambda_m(Tp.size(), 0.0);
        lambda_m = patch().lookupPatchField<volScalarField, scalar>("lambda_m");                               

    Time& time = const_cast<Time&>(this->patch().boundaryMesh().mesh().time());
    
    dictionary TambValueIO;
    TambValueIO.add(
        "file", 
        Tamb_
    );
    Function1s::TableFile<scalar> TambValue
    (
        "TambValue",
        TambValueIO
    );
    scalar TambValue_ = TambValue.value(time.value());

    refValue() = TambValue_;
    refGrad() = 0;

    const scalarField kappaDeltaCoeffs
    (
        lambda_m * patch().deltaCoeffs()
    );
    valueFraction() = hcoeff_ / (hcoeff_ + kappaDeltaCoeffs);

    mixedFvPatchScalarField::updateCoeffs(); 

    // Restore tag
    UPstream::msgType() = oldTag;

}


void CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("hcoeff")<< hcoeff_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tamb")<< Tamb_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
