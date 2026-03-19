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
    CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField

Description
    Boundary condition for heat flux at solid interior boundaries
    Written by Aytac Kubilay    

    Example usage:
        myInterfacePatchName
        {
            type            compressible::CFDHAMsolidTemperatureTransferCoeffInterior;
            hcoeff          0.405;
            Tamb            "$FOAM_CASE/0/regionName/T_room.dat";         
            value           uniform 293.15;
        }

SourceFiles
    CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField_H
#define CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
:
    public mixedFvPatchScalarField
{        
    // Private data
        
         //- Convective heat transfer coefficient [W/m2K]
        const scalar hcoeff_;
        
         //- Path to the file with the reference temperature values
        const fileName Tamb_;                               	                    

public:

    //- Runtime type information
    TypeName("compressible::CFDHAMsolidTemperatureTransferCoeffInterior");


    // Constructors

        //- Construct from patch and internal field
        CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
        (
            const
            CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
        (
            const CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField&
        ) = delete;

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
                (
                    *this,
                    this->internalField()
                )
            );
        }

        //- Construct as copy setting internal field reference
        CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
        (
            const CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField&,
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
                new CFDHAMsolidTemperatureTransferCoeffInteriorFvPatchScalarField
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
