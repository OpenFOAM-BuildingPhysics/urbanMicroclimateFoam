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

Class
    Foam::mappedLeafTempFvPatchScalarField

Description
    This boundary condition maps the value at a set of cells or patch faces
    back to *this.

    The sample mode is set by the underlying mapping engine, provided by the
    mappedPatchBase class.

Usage
    \table
        Property     | Description             | Required    | Default value
        fieldName    | name of field to be mapped | no       | this field name
        setAverage   | flag to activate setting of average value | yes |
        average      | average value to apply if \c setAverage = yes | yes |
        interpolationScheme | type of interpolation scheme | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedLeafTemperature;
        fieldName       Tl;
        setAverage      no;
        average         0;
        interpolationScheme cell;
        value           uniform 285;
    }
    \endverbatim

    Boundary condition to map leaf temperature from air region to vegetation
    region, based on the default mapped boundary condition. 
    To be used with patch "air_to_vegetation" using "nearestCell" sampleMode.

    Written by Aytac Kubilay
    
SourceFiles
    mappedLeafTempFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedLeafTempFvPatchScalarField_H
#define mappedLeafTempFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
//v8: #include "mappedPatchBase.H"
#include "mappedInternalPatchBase.H" //v12: cell-based mapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mappedLeafTempFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mappedLeafTempFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected Data

        //- Name of field to map
        const word fieldName_;

        //v8: mutable autoPtr<mappedPatchBase> mapperPtr_;
        //v12: mapperPtr_ no longer needed - mapper comes from patch polyPatch


    // Protected Member Functions

        //v12: returns mappedInternalPatchBase (cell-based mapping)
        const mappedInternalPatchBase& mapper() const;


public:

    //- Runtime type information
    TypeName("mappedLeafTemperature");


    // Constructors

        //- Construct from patch and internal field
        mappedLeafTempFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedLeafTempFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mappedLeafTempFvPatchScalarField
        //  onto a new patch
        mappedLeafTempFvPatchScalarField
        (
            const mappedLeafTempFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        mappedLeafTempFvPatchScalarField
        (
            const mappedLeafTempFvPatchScalarField&
        ) = delete;

        //- Construct as copy setting internal field reference
        mappedLeafTempFvPatchScalarField
        (
            const mappedLeafTempFvPatchScalarField&,
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
                new mappedLeafTempFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

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
