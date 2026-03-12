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
    readFieldFvPatchVectorField

Description
    Reads lists of boundary values with spatial and temporal variation based
    on mesoscale simulation

    Written by Aytac Kubilay    

    Example usage:
        myPatchName
        {
            type            readVectorField;
            inputTimeStep   3600;            
            value           uniform (1 1 0);
        }

SourceFiles
    readFieldFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef readFieldFvPatchVectorField_H
#define readFieldFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class readFieldFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class readFieldFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{        
    // Private data
        
        //mesoscale output timestep
        label inputTimeStep;

        //mesoscale target field
        vectorField Target_Field;
        
        //field name
        word fieldName;

public:

    //- Runtime type information
    TypeName("readVectorField");


    // Constructors

        //- Construct from patch and internal field
        readFieldFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        readFieldFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        readFieldFvPatchVectorField
        (
            const
            readFieldFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );
        
        readFieldFvPatchVectorField
        (
            const readFieldFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        readFieldFvPatchVectorField
        (
            const readFieldFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new readFieldFvPatchVectorField
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
