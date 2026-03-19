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

#include "CFDHAMfluidMoistureCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
//#include "mappedPatchBase.H"
#include "mappedFvPatchBaseBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::
CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::
CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
(
    const CFDHAMfluidMoistureCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper)
{}


CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::
CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
    //v8: if (!isA<mappedPatchBase>(this->patch().patch()))
    if (!isA<mappedFvPatchBaseBase>(this->patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedFvPatchBaseBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

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


CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::
CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
(
    const CFDHAMfluidMoistureCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the mapper and the neighbouring patch
    //v8: const mappedPatchBase& mpp =
    //v8:     refCast<const mappedPatchBase>(patch().patch());
    //v8: const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    //v8: const fvPatch& nbrPatch =
    //v8:     refCast<const fvMesh>(nbrMesh).boundary()[mpp.samplePolyPatch().index()];
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvPatch& nbrPatch = mapper.nbrFvPatch();

    scalarField& Tp = *this;

    tmp<scalarField> TNbr = mapper.fromNeighbour(nbrPatch.lookupPatchField<volScalarField, scalar>("Ts"));
    tmp<scalarField> pcNbr = mapper.fromNeighbour(nbrPatch.lookupPatchField<volScalarField, scalar>("pc"));

    scalarField p(Tp.size(), 0.0);
        p = patch().lookupPatchField<volScalarField, scalar>("p");

    scalarField rhoair(Tp.size(), 0.0);
        rhoair = patch().lookupPatchField<volScalarField, scalar>("rho");            

    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);
    scalarField pvsat_s = exp(6.58094e1-7.06627e3/TNbr()-5.976*log(TNbr()));

    scalarField pv_s = pvsat_s*exp((pcNbr())/(rhol*Rv*TNbr()));

    valueFraction() = 1.0;
    refValue() = 0.62198*pv_s/p;
    refGrad() = 0.0;

    mixedFvPatchScalarField::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;

}


void CFDHAMfluidMoistureCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
