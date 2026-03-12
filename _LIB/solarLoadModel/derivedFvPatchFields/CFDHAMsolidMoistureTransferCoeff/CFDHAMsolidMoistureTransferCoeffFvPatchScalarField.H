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
    CFDHAMsolidMoistureTransferCoeffFvPatchScalarField

Description
    Coupled boundary condition for moisture flux at solid boundaries
    Written by Aytac Kubilay    

    Example usage:
        myInterfacePatchName
        {
            type            compressible::CFDHAMsolidMoistureTransferCoeff;
            betacoeff       1e-8;
            pv_o            "$FOAM_CASE/0/pv_o.dat";
            value           uniform -1e7;
        }

SourceFiles
    CFDHAMsolidMoistureTransferCoeffFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef CFDHAMsolidMoistureTransferCoeffFvPatchScalarField_H
#define CFDHAMsolidMoistureTransferCoeffFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class CFDHAMsolidMoistureTransferCoeffFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
:
    public mixedFvPatchScalarField
{        
    // Private data
    
         //- Convective moisture transfer coefficient [s/m]
        const scalar betacoeff_;
        
         //- Path to the file with the reference exterior vapor pressure values
        const fileName pv_o_;      

public:

    //- Runtime type information
    TypeName("compressible::CFDHAMsolidMoistureTransferCoeff");


    // Constructors

        //- Construct from patch and internal field
        CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
        (
            const
            CFDHAMsolidMoistureTransferCoeffFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
        (
            const CFDHAMsolidMoistureTransferCoeffFvPatchScalarField&
        ) = delete;

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
                (
                    *this,
                    this->internalField()
                )
            );
        }

        //- Construct as copy setting internal field reference
        CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
        (
            const CFDHAMsolidMoistureTransferCoeffFvPatchScalarField&,
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
                new CFDHAMsolidMoistureTransferCoeffFvPatchScalarField
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
