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

\*---------------------------------------------------------------------------*/

#include "Impermeable.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(Impermeable, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        Impermeable,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::Impermeable::Impermeable
(
    const word& name,
    const dictionary& buildingMaterialDict,
    const word& cellZoneModel
)
:
    buildingMaterialModel(name, buildingMaterialDict, cellZoneModel)
{
    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Correct the buildingMaterial moisture content (cell)
void Foam::buildingMaterialModels::Impermeable::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    w[celli] = SMALL;
    Crel[celli] = GREAT;
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::Impermeable::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    Krel[celli] = SMALL;
}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::Impermeable::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{
    K_v[celli] = SMALL;
}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::Impermeable::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{
    K_pt[celli] = SMALL;
}


//*********************************************************** //
