/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Namespace
    Foam::buildingMaterialModels

Description
    A namespace for various buildingMaterialModel implementations.

Class
    Foam::buildingMaterialModel

Description
    An abstract base class for buildingMaterialModels

    Constitutive relations for buildingMaterials: theta(h) and K(h)

SourceFiles
    buildingMaterialModel.C
    newbuildingMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef buildingMaterialModel_H
#define buildingMaterialModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class buildingMaterialModel Declaration
\*---------------------------------------------------------------------------*/

class buildingMaterialModel
{

protected:

    // Protected data

        word name_;
        dictionary buildingMaterialDict_;
        word cellZoneModel_;

    // Private Member Functions

        //- Disallow copy construct
        buildingMaterialModel(const buildingMaterialModel&);

        //- Disallow default bitwise assignment
        void operator=(const buildingMaterialModel&);


public:

    //- Runtime type information
    TypeName("buildingMaterialModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            buildingMaterialModel,
            dictionary,
            (
                const word& name,
                const dictionary& buildingMaterialDict,
                const word& cellZoneModel
            ),
            (name, buildingMaterialDict, cellZoneModel)
        );


    // Selectors

        //- Return a reference to the selected buildingMaterial model
        static autoPtr<buildingMaterialModel> New
        (
            const word& name,
            const dictionary& buildingMaterialDict,
            const word& cellZoneModel
        );


    // Constructors

        //- Construct from components
        buildingMaterialModel
        (
            const word& name,
            const dictionary& buildingMaterialDict,
            const word& cellZoneModel
        );


    // Destructor

        virtual ~buildingMaterialModel()
        {}


    // Member Functions

        //- Return the buildingMaterial properties dictionary
        const dictionary& buildingMaterialDict() const
        {
            return buildingMaterialDict_;
        }

        //- Correct the buildingMaterial moisture content (cell)
        virtual void update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli) = 0;

        //- Correct the buildingMaterial liquid permeability (cell)
        virtual void update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli) = 0;

        //- Correct the buildingMaterial vapor permeability (cell)
        virtual void update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli) = 0;

        //- Correct the buildingMaterial vapor permeability (cell)
        virtual void update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
