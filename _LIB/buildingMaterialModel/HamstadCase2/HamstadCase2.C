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

#include "HamstadCase2.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(HamstadCase2, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        HamstadCase2,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::HamstadCase2::HamstadCase2
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
void Foam::buildingMaterialModels::HamstadCase2::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar T=293.15;
    
    scalar phi = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T));
    w[celli] = 116/(pow(1-(1/0.118*log(phi)),0.869));
    Crel[celli] = mag( (854.271)/ pow((rho_l*R_v*T)-8.47458*pc.internalField()[celli],1.869) );
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::HamstadCase2::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar T=293.15;
    
    scalar diffusivity = 6e-10;
    
    scalar Crel_tmp = mag( (854.271)/ pow((rho_l*R_v*T)-8.47458*pc.internalField()[celli],1.869) );
    
    Krel[celli] = diffusivity * Crel_tmp;
}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::HamstadCase2::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{
    /*
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar T=293.15;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T - 5.976*Foam::log(T)); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T)); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/1.57e2); 
    scalar delta = 2.61e-5 * tmp/(R_v*T*30*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for brick" [s]
    */
    K_v[celli] = 0;//(delta*p_vsat*relhum)/(rho_l*R_v*T);
}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::HamstadCase2::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    scalar dpsatdt = (7.06627e3/(T.internalField()[celli]*T.internalField()[celli]) - 5.976/T.internalField()[celli]) * p_vsat; // saturation vapour pressure [Pa]
        
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar delta = 1e-15;

    K_pt[celli] = delta*relhum*dpsatdt;
}

//*********************************************************** //
