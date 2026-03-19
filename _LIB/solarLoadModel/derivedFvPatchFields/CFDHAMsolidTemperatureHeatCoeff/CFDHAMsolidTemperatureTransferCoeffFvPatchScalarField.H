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
    CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField

Description
    Coupled boundary condition for heat flux at solid boundaries
    Written by Aytac Kubilay    

    Example usage:
        myInterfacePatchName
        {
            type            compressible::CFDHAMsolidTemperatureTransferCoeff;
            qrNbr           qr; //or none, name of thermal radiation from neighbour region
            qsNbr           qs; //or none, name of solar radiation from neighbour region
            hcoeff          8.0;
            Tamb            "$FOAM_CASE/0/Ta.dat";
            betacoeff       1e-8;
            pv_o            "$FOAM_CASE/0/pv_o.dat";            
            value           uniform 293.15;
        }

SourceFiles
    CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField_H
#define CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField
:
    public mixedFvPatchScalarField
{        
    // Private data

         //- Name of the thermal radiative heat flux in the neighbour region
        const word qrNbrName_;
        
         //- Name of the solar radiative heat flux in the neighbour region
        const word qsNbrName_;
        
         //- Convective heat transfer coefficient [W/m2K]
        const scalar hcoeff_;
        
         //- Path to the file with the reference temperature values
        const fileName Tamb_;   
        
         //- Convective moisture transfer coefficient [s/m]
        const scalar betacoeff_;
        
         //- Path to the file with the reference exterior vapor pressure values
        const fileName pv_o_;
        
        //veg qr
        scalarField qrNbr;
        //veg qs
        scalarField qsNbr;                               	                    

        scalar timeOfLastRadUpdate;

public:

    //- Runtime type information
    TypeName("compressible::CFDHAMsolidTemperatureTransferCoeff");


    // Constructors

        //- Construct from patch and internal field
        CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField
        (
            const
            CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField
        (
            const CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField&
        ) = delete;

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField
                (
                    *this,
                    this->internalField()
                )
            );
        }

        //- Construct as copy setting internal field reference
        CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField
        (
            const CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField&,
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
                new CFDHAMsolidTemperatureTransferCoeffFvPatchScalarField
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
