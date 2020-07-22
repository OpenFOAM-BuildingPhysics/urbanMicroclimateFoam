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
    CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField

Description
    Coupled boundary condition for heat flux at solid boundaries
    Written by Aytac Kubilay    

    Example usage:
        myInterfacePatchName
        {
            type            compressible::CFDHAMsolidTemperatureCoupledMixed;
            qrNbr           qr; //or none, name of thermal radiation from neighbour region
            qsNbr           qs; //or none, name of solar radiation from neighbour region
            value           uniform -1e7;
        }

SourceFiles
    CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField_H
#define CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{        
    // Private data

         //- Name of the thermal radiative heat flux in the neighbour region
        const word qrNbrName_;
        
         //- Name of the solar radiative heat flux in the neighbour region
        const word qsNbrName_;          	                    
                 

public:

    //- Runtime type information
    TypeName("compressible::CFDHAMsolidTemperatureCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
        (
            const
            CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
        (
            const CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField&,
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
                new CFDHAMsolidTemperatureCoupledMixedFvPatchScalarField
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
