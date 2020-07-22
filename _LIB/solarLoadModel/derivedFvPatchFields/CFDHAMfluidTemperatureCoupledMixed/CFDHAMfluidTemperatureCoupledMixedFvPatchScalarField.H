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
    CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField

Description
    Coupled boundary condition for temperature value at fluid boundaries
    Written by Aytac Kubilay    

    Example usage:
        myInterfacePatchName
        {
            type            compressible::CFDHAMfluidTemperatureCoupledMixed;
            value           uniform 300.0;
        }

SourceFiles
    CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField_H
#define CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{        
    // Private data       

        //- Name of field on the neighbour region
//        const word TnbrName_;                        

public:

    //- Runtime type information
    TypeName("compressible::CFDHAMfluidTemperatureCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField
        (
            const
            CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField
        (
            const CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField&,
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
                new CFDHAMfluidTemperatureCoupledMixedFvPatchScalarField
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
