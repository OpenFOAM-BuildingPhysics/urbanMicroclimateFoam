/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::blendingLayer

Description
    Blending layer for imposing boundary conditions from mesoscale model by
    imposing sink/source terms for the cells near lateral boundaries.
    
    - Aytac Kubilay, ETH Zurich, 2023

    Ref: Rodrigues, C.V., Palma, J.M.L.M., Rodrigues, Á.H., 2016. Atmospheric 
    Flow over a Mountainous Region by a One-Way Coupled Approach Based on 
    Reynolds-Averaged Turbulence Modelling. Boundary-Layer Meteorol 159, 407–437. 
    https://doi.org/10.1007/s10546-015-0116-7

SourceFiles
    blendingLayer.C

\*---------------------------------------------------------------------------*/

#ifndef blendingLayer_H
#define blendingLayer_H

#include "typeInfo.H"
#include "volFields.H"
#include "DimensionedField.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                         Class blendingLayer Declaration
\*---------------------------------------------------------------------------*/

class blendingLayer
{

protected:

    // Protected Data
    
        //- Reference to the mesh database
        const fvMesh& mesh_; 
        
        //- Reference to the time database
        const Time& time_;           

        // Field - blending layer boundary face Ids
        volScalarField bL_;
        
        // Field - blending layer U source
        volVectorField USource_;
        
        // Field - blending layer T source
        volScalarField TSource_;           

        // Dictionary - blending layer settings
        dictionary coeffs_;

        scalar dampingThickness;
        scalar alphaCoeffU;
        scalar alphaCoeffT;

    // Private Member Functions

        //- Disallow copy construct
        blendingLayer(const blendingLayer&);

        //- Disallow default bitwise assignment
        void operator=(const blendingLayer&);
        

public:

    // Static data members

        //- Run-time type information
        TypeName("blendingLayer");


    // Constructors

        //- Construct from components
        blendingLayer(const volVectorField& U, const volScalarField& T);


    //- Destructor
    virtual ~blendingLayer();


    // Member Functions

        //- Initialize blending layer
        virtual void initialize();
        
        //- Populate blending layer
        virtual void getValues(volVectorField& USource_, const volVectorField& U);

        //- Populate blending layer
        virtual void getValues(volScalarField& TSource_, const volScalarField& T);

        //- Return blending layer source for U
        virtual tmp<volVectorField> bL_USource(const volVectorField& U);

        //- Return blending layer source for T
        virtual tmp<volScalarField> bL_TSource(const volScalarField& T);
        
        //- Return blending layer field
        virtual tmp<volScalarField> bL() const
        {
            return bL_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
