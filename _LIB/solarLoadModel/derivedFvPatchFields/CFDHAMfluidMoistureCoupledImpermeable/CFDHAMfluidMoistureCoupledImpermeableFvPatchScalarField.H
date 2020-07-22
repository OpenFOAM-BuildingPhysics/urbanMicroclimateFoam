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
    CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField

Description
    Coupled boundary condition for moisture value at fluid impermeable boundaries
    Written by Aytac Kubilay
    
    Value entry is only place-holder

    Example usage:
        myInterfacePatchName
        {
            type            compressible::CFDHAMfluidMoistureCoupledImpermeable;
            value           uniform 0.008;
        }

SourceFiles
    CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField_H
#define CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField
:
    public mixedFvPatchScalarField
{        
    // Private data                            

public:

    //- Runtime type information
    TypeName("compressible::CFDHAMfluidMoistureCoupledImpermeable");


    // Constructors

        //- Construct from patch and internal field
        CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleImpermeableFvPatchScalarField onto a
        //  new patch
        CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField
        (
            const
            CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField
        (
            const CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField&,
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
                new CFDHAMfluidMoistureCoupledImpermeableFvPatchScalarField
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
