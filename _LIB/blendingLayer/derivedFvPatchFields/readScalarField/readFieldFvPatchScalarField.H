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

Class
    Foam::
    readFieldFvPatchScalarField

Description
    Reads lists of boundary values with spatial and temporal variation based
    on mesoscale simulation

    Written by Aytac Kubilay    

    Example usage:
        myPatchName
        {
            type            readScalarField;
            inputTimeStep   3600;            
            value           uniform 300.0;
        }

SourceFiles
    readFieldFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef readFieldFvPatchScalarField_H
#define readFieldFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class readFieldFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class readFieldFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{        
    // Private data
        
        //mesoscale output timestep
        label inputTimeStep;

        //mesoscale target field
        scalarField Target_Field;
        
        //field name
        word fieldName;

public:

    //- Runtime type information
    TypeName("readScalarField");


    // Constructors

        //- Construct from patch and internal field
        readFieldFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        readFieldFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given readFieldFvPatchScalarField onto a new patch
        readFieldFvPatchScalarField
        (
            const readFieldFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        readFieldFvPatchScalarField
        (
            const readFieldFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        readFieldFvPatchScalarField
        (
            const readFieldFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new readFieldFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
