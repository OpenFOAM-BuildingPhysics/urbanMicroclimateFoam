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

#include "CFDHAMsolidMoistureCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "fixedValueFvPatchFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::
CFDHAMsolidMoistureCoupledMixedFvPatchScalarField
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


CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::
CFDHAMsolidMoistureCoupledMixedFvPatchScalarField
(
    const CFDHAMsolidMoistureCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper)
{}


CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::
CFDHAMsolidMoistureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
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


CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::
CFDHAMsolidMoistureCoupledMixedFvPatchScalarField
(
    const CFDHAMsolidMoistureCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchI = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];

    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);                        
    scalar Dm = 2.5e-5; scalar Sct = 0.7;

    scalarField& pcp = *this;

    const mixedFvPatchScalarField&
        fieldTs = refCast
            <const mixedFvPatchScalarField>
            (
                patch().lookupPatchField<volScalarField, scalar>("Ts")
            );
    const fvPatchScalarField&
        fieldpc = refCast
            <const fvPatchScalarField>
            (
                patch().lookupPatchField<volScalarField, scalar>("pc")
            );  
    const mixedFvPatchScalarField&
        nbrFieldw = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("w")
            );            
                        
    scalarField Ts(pcp.size(), 0.0);
        Ts = patch().lookupPatchField<volScalarField, scalar>("Ts"); 
    scalarField TNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("T");
        mpp.distribute(TNbr); 

    scalarField wcNbr(nbrFieldw.patchInternalField());
        mpp.distribute(wcNbr);
    scalarField wNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("w");
    scalarField rhoNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("rho");
    scalarField pv_o = wNbr*1e5/(0.621945*rhoNbr);
        mpp.distribute(wNbr);
        mpp.distribute(rhoNbr);
        mpp.distribute(pv_o);  
    scalarField pv_o_sat = exp(6.58094e1-7.06627e3/TNbr-5.976*log(TNbr));
    scalarField pc_o=log(pv_o/pv_o_sat)*rhol*Rv*TNbr; 

    scalarField gcrNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("gcr");
        mpp.distribute(gcrNbr);    

    scalarField Krel(pcp.size(), 0.0);
        Krel = patch().lookupPatchField<volScalarField, scalar>("Krel"); 

    scalarField K_v(pcp.size(), 0.0);
        K_v = patch().lookupPatchField<volScalarField, scalar>("K_v");             

    scalarField deltaCoeff_ = nbrPatch.deltaCoeffs(); 
        mpp.distribute(deltaCoeff_);
    scalarField nutNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("nut");
        mpp.distribute(nutNbr);
    
    scalarField pvsat_s = exp(6.58094e1-7.06627e3/Ts-5.976*log(Ts));
    scalarField pv_s = pvsat_s*exp((pcp)/(rhol*Rv*Ts));

    scalarField gl = ((gcrNbr*rhol)/(3600*1000));
    scalarField g_conv = rhoNbr*(Dm + nutNbr/Sct) * (wcNbr-(0.62198*pv_s/1e5)) *deltaCoeff_; 
//    scalarField g_cond = (Krel+K_v)*fieldpc.snGrad();  
    
    // term with temperature gradient:
    scalarField K_pt(pcp.size(), 0.0);
        K_pt = patch().lookupPatchField<volScalarField, scalar>("K_pt");                 
    scalarField X = K_pt*fieldTs.snGrad();
    //////////////////////////////////

    //-- Gravity flux --//
    //lookup gravity vector
    uniformDimensionedVectorField g = db().lookupObject<uniformDimensionedVectorField>("g");
    scalarField gn = g.value() & patch().nf();

    scalarField phiG = Krel*rhol*gn;
    //////////////////////////////////

    valueFraction() = 0.0;
    if(gMax(gl) > 0)
    {
        scalarField g_cond = (Krel+K_v)*(-10.0-fieldpc.patchInternalField())*patch().deltaCoeffs();        
        forAll(valueFraction(),faceI)
        {  
            //if(pcp[faceI] > -100.0 && (gl[faceI] > g_cond[faceI] - g_conv[faceI] - phiG[faceI] + X[faceI]) )    
            if( (gl[faceI] > g_cond[faceI] - g_conv[faceI] - phiG[faceI] + X[faceI]) )    
            {
                valueFraction()[faceI] = 1.0;
            }
        }
    }
   
    refGrad() = (g_conv + gl + phiG - X)/(Krel + K_v);
//    refValue() =  -100.0 + 1.0;
    refValue() =  -10.0;
        

    mixedFvPatchScalarField::updateCoeffs(); 

    // Restore tag
    UPstream::msgType() = oldTag;

}


void CFDHAMsolidMoistureCoupledMixedFvPatchScalarField::write
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
    CFDHAMsolidMoistureCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
