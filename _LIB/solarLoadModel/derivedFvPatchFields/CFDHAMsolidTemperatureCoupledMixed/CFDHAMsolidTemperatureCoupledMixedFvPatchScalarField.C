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

#include "CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "fixedValueFvPatchFields.H"
#include "TableFile.H"
#include "uniformDimensionedFields.H"

#include "hashedWordList.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::
CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    qrNbrName_("undefined-qrNbr"),
    qsNbrName_("undefined-qsNbr"),
    qrNbr(0),
    qsNbr(0),
    timeOfLastRadUpdate(-1.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::
CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
(
    const CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    qrNbrName_(psf.qrNbrName_),
    qsNbrName_(psf.qsNbrName_),
    qrNbr(psf.qrNbr),
    qsNbr(psf.qsNbr),
    timeOfLastRadUpdate(psf.timeOfLastRadUpdate)
{}


CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::
CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    qrNbrName_(dict.lookupOrDefault<word>("qrNbr", "none")),
    qsNbrName_(dict.lookupOrDefault<word>("qsNbr", "none")),
    qrNbr(Zero),
    qsNbr(Zero),
    timeOfLastRadUpdate(-1.0)
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


CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::
CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
(
    const CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    qrNbrName_(psf.qrNbrName_),
    qsNbrName_(psf.qsNbrName_),
    qrNbr(psf.qrNbr),
    qsNbr(psf.qsNbr),
    timeOfLastRadUpdate(psf.timeOfLastRadUpdate)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::updateCoeffs()
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

    scalar cap_v = 1880; scalar Tref = 273.15; scalar L_v = 2.5e6; scalar cap_l = 4182;
    scalar cp = 1005; //specific heat of air [J/(kg K)]
    scalar muair = 1.8e-5; scalar Pr = 0.7;
    scalar Dm = 2.5e-5; scalar Sct = 0.7;
    scalar rhol=1.0e3; scalar Rv=8.31451*1000/(18.01534);

    scalarField& Tp = *this;

    const mixedFvPatchScalarField& //CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField&
        nbrField = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("T")
            );
    scalarField TcNbr(nbrField.patchInternalField()); 
        mpp.distribute(TcNbr);
    scalarField TNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("T");
        mpp.distribute(TNbr);

    const mixedFvPatchScalarField&
        nbrFieldw = refCast
            <const mixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("w")
            );
    scalarField wcNbr(nbrFieldw.patchInternalField());
        mpp.distribute(wcNbr);
    scalarField wNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("w");
    scalarField rhoNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("rho");
    scalarField pv_o = wNbr*1e5/(0.621945*rhoNbr);
        mpp.distribute(wNbr);
        mpp.distribute(rhoNbr);  
        mpp.distribute(pv_o); 

    const mixedFvPatchScalarField&
        fieldpc = refCast
            <const mixedFvPatchScalarField>
            (
                patch().lookupPatchField<volScalarField, scalar>("pc")
            );
    const fvPatchScalarField&
        fieldTs = refCast
            <const fvPatchScalarField>
            (
                patch().lookupPatchField<volScalarField, scalar>("Ts")
            );
                        
    scalarField pc(Tp.size(), 0.0);
        pc = patch().lookupPatchField<volScalarField, scalar>("pc");    
    scalarField K_pt(Tp.size(), 0.0);
        K_pt = patch().lookupPatchField<volScalarField, scalar>("K_pt"); 
    scalarField lambda_m(Tp.size(), 0.0);
        lambda_m = patch().lookupPatchField<volScalarField, scalar>("lambda_m");                               

    scalarField deltaCoeff_ = nbrPatch.deltaCoeffs();
        mpp.distribute(deltaCoeff_);
    scalarField alphatNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("alphat");
        mpp.distribute(alphatNbr);
    scalarField nutNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("nut");
        mpp.distribute(nutNbr); 
    
    scalarField q_conv = (muair/Pr + alphatNbr)*cp*(TcNbr-Tp)*deltaCoeff_; 
            
    scalarField pvsat_s = exp(6.58094e1-7.06627e3/Tp-5.976*log(Tp));
    scalarField pv_s = pvsat_s*exp((pc)/(rhol*Rv*Tp));
    
    scalarField g_conv = rhoNbr*(Dm + nutNbr/Sct) * (wcNbr-(0.62198*pv_s/1e5)) *deltaCoeff_; 
    scalarField LE = (cap_v*(Tp-Tref)+L_v)*g_conv;//Latent and sensible heat transfer due to vapor exchange   */

    scalarField K_v(Tp.size(), 0.0);
        K_v = patch().lookupPatchField<volScalarField, scalar>("K_v");  
    scalarField Krel(Tp.size(), 0.0);
        Krel = patch().lookupPatchField<volScalarField, scalar>("Krel");   

    scalarField gcrNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("gcr");
        mpp.distribute(gcrNbr); 

    scalarField gl = ((gcrNbr*rhol)/(3600*1000));

    // Set rain temperature //////////////////////////////////////////////////
    Time& time = const_cast<Time&>(nbrMesh.time());
    //label timestep = ceil( (time.value()/3600)-1E-6 ); timestep = timestep%24;

    fileName rainTempFile
    (
       nbrMesh.time().rootPath()
       /nbrMesh.time().globalCaseName()
       /"0/air/rainTemp"
    );
    scalar rainTemp = 293.15;
    if(isFile(rainTempFile))
    {
//        Info << "Found rainTemp file..." << endl;
        dictionary rainTempIO;
        rainTempIO.add(
            "file", 
            rainTempFile
        );
        Function1s::TableFile<scalar> rT
        (
            "rainTemp",
            rainTempIO
        );
        rainTemp = rT.value(time.value());
    }
    else
    {
//        Info << "Calculating rainTemp..." << endl;
        // Calculate rain temperature - approximation for wet-bulb temp///////////
        //obtain Tambient - can find a better way to import this value?
        dictionary TambientIO;
        TambientIO.add(
            "file", 
            fileName
            (
                nbrMesh.time().rootPath()
                /nbrMesh.time().globalCaseName()
                /"0/air/Tambient"
            )
        );
        Function1s::TableFile<scalar> Tambient
        (
            "Tambient",
            TambientIO
        );
        
        dictionary wambientIO;
        wambientIO.add(
            "file", 
            fileName
            (
                nbrMesh.time().rootPath()
                /nbrMesh.time().globalCaseName()
                /"0/air/wambient"
            )
        );
        Function1s::TableFile<scalar> wambient
        (
            "wambient",
            wambientIO
        );      
        ///////////
        scalar Tambient_ = Tambient.value(time.value());
        scalar wambient_ = wambient.value(time.value());
        scalar saturationPressure = 133.322*pow(10,(8.07131-(1730.63/(233.426+Tambient_-273.15))));
        scalar airVaporPressure = wambient_*1e5/0.621945;
        scalar relhum = airVaporPressure/saturationPressure*100;
        scalar dewPointTemp = Tambient_ - (100-relhum)/5;
        rainTemp = Tambient_ - (Tambient_-dewPointTemp)/3;
    }
    //////////////////////////////////////////////////////////////////////////

    //scalarField qrNbr(Tp.size(), 0.0);
    //scalarField qsNbr(Tp.size(), 0.0);
    dictionary controlDict_ = time.controlDict();
    const scalar deltaT_(readScalar(controlDict_.lookup("deltaT")));
    label moduloTest = int(time.value()/deltaT_);
    bool firstIter = false;
    if(time.value()/deltaT_ - moduloTest < SMALL)
    {
        if(timeOfLastRadUpdate != time.value())
        {
            firstIter = true; //check if first internal iteration
        }
    }
    bool radUpdateNow = false;
    if ((firstIter) or (time.value() - timeOfLastRadUpdate >= 600.0)) //update rad once at the beginning and every 600 s
    {
        radUpdateNow = true;
        timeOfLastRadUpdate = time.value();
    }

    //-- Access vegetation region and populate radiation if vegetation exists,
    //otherwise use radiation from air region --//
    regionProperties rp(time);
    const wordList vegNames(rp["vegetation"]);

    if (vegNames.size()>0)
    {
        if(radUpdateNow) //update qs and qr once at the beginning
        {
            const word& vegiRegion = "vegetation";
            const scalar mppVegDistance = 0;
     
            const polyMesh& vegiMesh =
                patch().boundaryMesh().mesh().time().lookupObject<polyMesh>(vegiRegion);
     
            const word& nbrPatchName = nbrPatch.name();
     
            const label patchi = vegiMesh.boundaryMesh().findPatchID(nbrPatchName);
        
            const fvPatch& vegiNbrPatch =
                refCast<const fvMesh>(vegiMesh).boundary()[patchi];
     
            const mappedPatchBase& mppVeg = mappedPatchBase(patch().patch(), vegiRegion, mpp.mode(), mpp.samplePatch(), mppVegDistance);
     
            if (qrNbrName_ != "none")
            {
                qrNbr = vegiNbrPatch.lookupPatchField<volScalarField, scalar>(qrNbrName_);
                mppVeg.distribute(qrNbr);
            }
            if (qsNbrName_ != "none")
            {
                qsNbr = vegiNbrPatch.lookupPatchField<volScalarField, scalar>(qsNbrName_);
                mppVeg.distribute(qsNbr);
            }
            timeOfLastRadUpdate = time.value();
        }
    }
    else
    {
        if(radUpdateNow) //update qs and qr once at the beginning
        {
            if (qrNbrName_ != "none")
            {
                qrNbr = nbrPatch.lookupPatchField<volScalarField, scalar>(qrNbrName_);
                mpp.distribute(qrNbr);
            }
            if (qsNbrName_ != "none")
            {
                qsNbr = nbrPatch.lookupPatchField<volScalarField, scalar>(qsNbrName_);
                mpp.distribute(qsNbr);
            }
            timeOfLastRadUpdate = time.value();
        }
    }
    //////////////////////////////


    //-- Grass adjustments --//
    IOdictionary grassProperties
    (
		IOobject
		(
		    "grassProperties",
		    nbrMesh.time().constant(),
		    nbrMesh,
		    IOobject::READ_IF_PRESENT,
		    IOobject::NO_WRITE
        )
    );

    if (grassProperties.typeHeaderOk<IOdictionary>(true))
    {
        word grassModel(grassProperties.lookup("grassModel"));
        if (grassModel != "none")
        {
            const dictionary& modelCoeffs = grassProperties.subDict(grassModel + "Coeffs");
            hashedWordList grassPatches = modelCoeffs.lookup("grassPatches");

            if (grassPatches.found(nbrPatch.name()))//if patch is covered with grass
            {
                if(radUpdateNow) //update qs and qr once at the beginning
                {                
                    scalarField TgNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("Tg");
                        mpp.distribute(TgNbr);
    
                    const dictionary& coeffs = grassProperties.subDict(grassModel + "Coeffs");
                    scalar LAI = coeffs.lookupOrDefault("LAI", 2.0);
                    scalar beta = coeffs.lookupOrDefault("beta", 0.78);
                    scalar albedoSoil = coeffs.lookupOrDefault("albedoSoil", 0.0);
                    qrNbr = 6*(TgNbr-Tp); //thermal radiation between grass and surface - Malys et al 2014
                                          //assuming external thermal radiation is fully absorbed with grass layer
                    qsNbr = qsNbr*exp(-beta*LAI)*(1-albedoSoil); //solar radiation transmitted through grass layer - solar radiation reflected from soil surface
                }
            }
        }
    }
    ///////////////////////////

    //-- Gravity-enthalpy flux --//
    //lookup gravity vector
    uniformDimensionedVectorField g = db().lookupObject<uniformDimensionedVectorField>("g");
    scalarField gn = g.value() & patch().nf();

    scalarField phiG = Krel*rhol*gn;
    scalarField phiGT = (cap_l*(Tp-Tref))*phiG;

    // term with capillary moisture gradient:                          
    scalarField X = ((cap_l*(Tp-Tref)*Krel)+(cap_v*(Tp-Tref)+L_v)*K_v)*fieldpc.snGrad();
    // moisture flux term with temperature gradient:               
    scalarField Xmoist = K_pt*fieldTs.snGrad();    
    //////////////////////////////////  

    scalarField CR(Tp.size(), 0.0);
    if(gMax(gl) > 0)
    {
        //scalarField g_cond = (Krel+K_v)*fieldpc.snGrad();
        scalarField g_cond = (Krel+K_v)*(-10.0-fieldpc.patchInternalField())*patch().deltaCoeffs();       
        forAll(CR,faceI)
        {
            scalar rainFlux = 0;
            //if(pc[faceI] > -100.0 && (gl[faceI] > g_cond[faceI] - g_conv[faceI] - phiG[faceI] + Xmoist[faceI]) )
            if( (gl[faceI] > g_cond[faceI] - g_conv[faceI] - phiG[faceI] + Xmoist[faceI]) )
            {
                rainFlux = g_cond[faceI] - g_conv[faceI] - phiG[faceI] + Xmoist[faceI];
            }
            else
            {
                rainFlux = gl[faceI];
            }
            CR[faceI] = rainFlux * cap_l*(rainTemp - Tref);
        }
    }         

    if(fieldpc.type() == "compressible::CFDHAMsolidMoistureCoupledImpermeable")
    {
        valueFraction() = 0;
        refValue() = 0;
        refGrad() = (q_conv + qrNbr + qsNbr)/(lambda_m);
    }
    else
    {
        valueFraction() = 0;
        refValue() = 0;
        refGrad() = (q_conv + LE + qrNbr + qsNbr + CR + phiGT -X)/(lambda_m+(cap_v*(Tp-Tref)+L_v)*K_pt);
    }

    mixedFvPatchScalarField::updateCoeffs(); 

    // Restore tag
    UPstream::msgType() = oldTag;

}


void CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("qrNbr")<< qrNbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("qsNbr")<< qsNbrName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
