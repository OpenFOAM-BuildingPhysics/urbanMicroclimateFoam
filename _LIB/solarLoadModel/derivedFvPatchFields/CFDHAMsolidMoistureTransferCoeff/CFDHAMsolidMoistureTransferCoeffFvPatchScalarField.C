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

#include "CFDHAMsolidMoistureTransferCoeffFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
//#include "mappedPatchBase.H" //v8
#include "mappedFvPatchBaseBase.H" //v12
#include "fixedValueFvPatchFields.H"
#include "Function1.H"
#include "Table.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CFDHAMsolidMoistureTransferCoeffFvPatchScalarField::
CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    betacoeff_(),
    pv_o_()    
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


CFDHAMsolidMoistureTransferCoeffFvPatchScalarField::
CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
(
    const CFDHAMsolidMoistureTransferCoeffFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    betacoeff_(psf.betacoeff_),
    pv_o_(psf.pv_o_)    
{}


CFDHAMsolidMoistureTransferCoeffFvPatchScalarField::
CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    betacoeff_(dict.lookupOrDefault<scalar>("betacoeff",0)),
    pv_o_(dict.lookupOrDefault<fileName>("pv_o", "none"))    
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


CFDHAMsolidMoistureTransferCoeffFvPatchScalarField::
CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
(
    const CFDHAMsolidMoistureTransferCoeffFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    betacoeff_(psf.betacoeff_),
    pv_o_(psf.pv_o_)    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void CFDHAMsolidMoistureTransferCoeffFvPatchScalarField::updateCoeffs()
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
    const fvMesh& nbrMesh = mapper.nbrMesh();

    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);

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
    tmp<scalarField> TNbr = mapper.fromNeighbour(nbrPatch.lookupPatchField<volScalarField, scalar>("T"));

    tmp<scalarField> wcNbr = mapper.fromNeighbour(nbrFieldw.patchInternalField());
    scalarField wNbr_local = nbrPatch.lookupPatchField<volScalarField, scalar>("w");
    scalarField rhoNbr_local = nbrPatch.lookupPatchField<volScalarField, scalar>("rho");
    scalarField pv_o_local = wNbr_local*1e5/(0.621945*rhoNbr_local);
    tmp<scalarField> wNbr = mapper.fromNeighbour(wNbr_local);
    tmp<scalarField> rhoNbr = mapper.fromNeighbour(rhoNbr_local);
    tmp<scalarField> pv_o_nbr = mapper.fromNeighbour(pv_o_local);
    scalarField pv_o_sat = exp(6.58094e1-7.06627e3/TNbr()-5.976*log(TNbr()));
    scalarField pc_o=log(pv_o_nbr()/pv_o_sat)*rhol*Rv*TNbr();

    tmp<scalarField> gcrNbr = mapper.fromNeighbour(nbrPatch.lookupPatchField<volScalarField, scalar>("gcr"));

    scalarField Krel(pcp.size(), 0.0);
        Krel = patch().lookupPatchField<volScalarField, scalar>("Krel");

    scalarField K_v(pcp.size(), 0.0);
        K_v = patch().lookupPatchField<volScalarField, scalar>("K_v");

    tmp<scalarField> deltaCoeff_ = mapper.fromNeighbour(nbrPatch.deltaCoeffs());
    tmp<scalarField> nutNbr = mapper.fromNeighbour(nbrPatch.lookupPatchField<volScalarField, scalar>("nut"));

    scalarField pvsat_s = exp(6.58094e1-7.06627e3/Ts-5.976*log(Ts));
    scalarField pv_s = pvsat_s*exp((pcp)/(rhol*Rv*Ts));

    scalarField gl = ((gcrNbr()*rhol)/(3600*1000));

    Time& time = const_cast<Time&>(nbrMesh.time());

    dictionary pv_oValueIO;
    pv_oValueIO.add("type", "table");
    pv_oValueIO.add(
        "file",
        pv_o_
    );
    Function1s::Table<scalar> pv_oValue
    (
        "pv_oValue",
        dimTime,
        dimPressure,
        pv_oValueIO
    );
    scalar pv_oValue_ = pv_oValue.value(time.value());
    scalarField g_conv = betacoeff_*(pv_oValue_-pv_s); 
    
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


void CFDHAMsolidMoistureTransferCoeffFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("betacoeff")<< betacoeff_ << token::END_STATEMENT << nl;
    os.writeKeyword("pv_o")<< pv_o_ << token::END_STATEMENT << nl;    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
