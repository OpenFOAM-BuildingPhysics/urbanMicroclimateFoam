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
    compressible::
    CFDHAMfluidMoistureCoupledMixedFvPatchScalarField

Description
    Coupled boundary condition for moisture value at fluid boundaries
    Written by Aytac Kubilay    

    Example usage:
        myInterfacePatchName
        {
            type            compressible::CFDHAMfluidMoistureCoupledMixed;
            value           uniform 0.008;
        }

SourceFiles
    CFDHAMfluidMoistureCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef CFDHAMfluidMoistureCoupledMixedFvPatchScalarField_H
#define CFDHAMfluidMoistureCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class CFDHAMfluidMoistureCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{        
    // Private data                            

public:

    //- Runtime type information
    TypeName("compressible::CFDHAMfluidMoistureCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
        (
            const
            CFDHAMfluidMoistureCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
        (
            const CFDHAMfluidMoistureCoupledMixedFvPatchScalarField&,
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
                new CFDHAMfluidMoistureCoupledMixedFvPatchScalarField
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

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
