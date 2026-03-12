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

#include "PorousAsphalt.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(PorousAsphalt, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        PorousAsphalt,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::PorousAsphalt::PorousAsphalt
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
void Foam::buildingMaterialModels::PorousAsphalt::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    List<scalar> reta; reta.setSize(3); reta[0]=-0.00283; reta[1]=-2.041e-3; reta[2]=-2.041e-8;
    List<scalar> retn; retn.setSize(3); retn[0]=8.2; retn[1]=1.4; retn[2]=1.4;
    List<scalar> retm; retm.setSize(3); retm[0]=0.8780; retm[1]=0.2857; retm[2]=0.2857;
    List<scalar> retw; retw.setSize(3); retw[0]=0.7; retw[1]=0.2; retw[2]=0.1;
    scalar w_tmp = 0; scalar tmp = 0; scalar C_tmp = 0; scalar tmp2 = 0;    
    for (int i=0; i<=2; i++)
    {
        tmp = pow( (reta[i]*pc.internalField()[celli]) , retn[i] );
        w_tmp = w_tmp + retw[i] / ( pow( (1 + tmp) , retm[i] ));
        tmp2 = pow( (1 + tmp) , retm[i] );
        C_tmp = C_tmp - retw[i]/tmp2 * retm[i]*retn[i]*tmp/((1 + tmp)*pc.internalField()[celli]); 
    } 
    w[celli] = w_tmp*48.8;   
    Crel[celli] = mag( C_tmp*48.8 );   
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::PorousAsphalt::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    scalar tmp=w.internalField()[celli]-73;
    tmp=-39.2619e0 +0.0704e0*tmp -1.742e-4*pow(tmp,2) -2.7953e-6*pow(tmp,3) -1.1566e-7*pow(tmp,4) +2.5969e-9*pow(tmp,5);
    Krel[celli] = pow(10,tmp*0.4342944819e0);

}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::PorousAsphalt::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/48.8); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*2*(0.89*tmp*tmp + 0.11)); 
    
    K_v[celli] = (delta*p_vsat*relhum)/(rho_l*R_v*T.internalField()[celli]);
}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::PorousAsphalt::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar L_v = 2.5e6;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    //scalar dpsatdt = (7.06627e3/(T.internalField()[celli]*T.internalField()[celli]) - 5.976/T.internalField()[celli]) * p_vsat; // saturation vapour pressure [Pa]
        
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/48.8); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*2*(0.89*tmp*tmp + 0.11)); // Water vapour diffusion coefficient "for concrete" [s]

    K_pt[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) ) * (rho_l*L_v - pc.internalField()[celli]);
}

//*********************************************************** //
