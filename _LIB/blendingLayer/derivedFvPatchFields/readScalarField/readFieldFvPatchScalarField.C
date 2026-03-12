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

#include "readFieldFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::readFieldFvPatchScalarField::
readFieldFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    inputTimeStep(),
    Target_Field(p.size()),
    fieldName(iF.name())
{
}


Foam::readFieldFvPatchScalarField::
readFieldFvPatchScalarField
(
    const readFieldFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    inputTimeStep(ptf.inputTimeStep),
    Target_Field(ptf.Target_Field),
    fieldName(ptf.fieldName)
{}


Foam::readFieldFvPatchScalarField::
readFieldFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    inputTimeStep(dict.lookup<label>("inputTimeStep")),
    Target_Field(p.size()),
    fieldName(iF.name())
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::readFieldFvPatchScalarField::
readFieldFvPatchScalarField
(
    const readFieldFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    inputTimeStep(ptf.inputTimeStep),
    Target_Field(ptf.Target_Field),
    fieldName(ptf.fieldName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::readFieldFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar timeValue = this->db().time().value();
    scalar timeIndex = this->db().time().timeIndex();
    word boundaryName = this->patch().name();
      
    if (timeIndex == 1)
    {
        label moduloTest = int(timeValue/inputTimeStep); 
        word TargetFile = fieldName + "target_" + boundaryName + "/" + fieldName + "target_" + boundaryName +"_" + name(moduloTest*inputTimeStep);
        IOList<scalar> Target(
            IOobject
            (
                TargetFile,
                db().time().caseConstant(),
                db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )            
        );
        /////interpolate between two input files if necessary/////
        if (timeValue/inputTimeStep - moduloTest > 0)
        {
            word TargetFile_B = fieldName + "target_" + boundaryName + "/" + fieldName + "target_" + boundaryName +"_" + name(moduloTest*inputTimeStep+inputTimeStep);
            IOList<scalar> Target_B(
                IOobject
                (
                    TargetFile_B,
                    db().time().caseConstant(),
                    db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )            
            );
            scalar ratio = (timeValue-moduloTest*inputTimeStep)/inputTimeStep;
            Target = Target*(1-ratio) + Target_B*(ratio);       
        }
        //////////////////////////////////////////////////////////

        if (Pstream::parRun()) //in parallel runs, need to read correct values from global wdrFile
        {
            const labelIOList& localFaceProcAddr //global addresses for local faces, includes also internal faces
            (
                IOobject
                (
            	"faceProcAddressing",
            	this->patch().boundaryMesh().mesh().facesInstance(),
            	this->patch().boundaryMesh().mesh().meshSubDir,
            	this->patch().boundaryMesh().mesh(),
            	IOobject::MUST_READ,
            	IOobject::NO_WRITE
                )
            );

            label startFace = this->patch().start(); //local startFace for this patch
            label nFaces = this->patch().size(); //local total face number for this patch

            List<label> globalFaceAddr_;
            globalFaceAddr_.setSize(nFaces);
            forAll(globalFaceAddr_, i) //global address for local faces, only for this patch
            {
                globalFaceAddr_[i] = localFaceProcAddr[startFace + i] - 1; //subtracted 1 to get global index, as localFaceProcAddr starts from 1, not 0
            }
            label minGlobalFaceAddr_ = gMin(globalFaceAddr_); //get the minimum global address for this patch = startFace in global patch

            List<scalar> Target_;
            Target_.setSize(nFaces);
            forAll(Target_, i) //read correct values from global wdrFile
            {
                Target_[i] = Target[localFaceProcAddr[startFace + i] - 1 - minGlobalFaceAddr_]; //subtracted 1 to get global index, as localFaceProcAddr starts from 1, not 0
            }
            
            Target_Field = Target_;                 
        }
        else
        {
            Target_Field = Target;
        }
    }        

    operator==(Target_Field);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::readFieldFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "inputTimeStep", inputTimeStep);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        readFieldFvPatchScalarField
    );
}

// ************************************************************************* //
