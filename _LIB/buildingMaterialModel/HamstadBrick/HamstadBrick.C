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

#include "HamstadBrick.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(HamstadBrick, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        HamstadBrick,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::HamstadBrick::HamstadBrick
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
void Foam::buildingMaterialModels::HamstadBrick::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    List<scalar> reta; reta.setSize(2); reta[0]=-1.25e-5; reta[1]=-1.80e-5;
    List<scalar> retn; retn.setSize(2); retn[0]=1.65e0; retn[1]=6.00e0;
    List<scalar> retm; retm.setSize(2); retm[0]=0.39394e0; retm[1]=0.83333e0;
    List<scalar> retw; retw.setSize(2); retw[0]=0.300e0; retw[1]=0.700e0;
    scalar w_tmp = 0; scalar tmp = 0; scalar C_tmp = 0; scalar tmp2 = 0;    
    for (int i=0; i<=1; i++)
    {
        tmp = pow( (reta[i]*pc.internalField()[celli]) , retn[i] );
        w_tmp = w_tmp + retw[i] / ( pow( (1 + tmp) , retm[i] ));
        tmp2 = pow( (1 + tmp) , retm[i] );
        C_tmp = C_tmp - retw[i]/tmp2 * retm[i]*retn[i]*tmp/((1 + tmp)*pc.internalField()[celli]);   
    }
    w[celli] = w_tmp*157;
    Crel[celli] = mag( C_tmp*157 );   
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::HamstadBrick::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    scalar Ks=1.907E-9; scalar tau=-1.631;
    List<scalar> reta; reta.setSize(3); reta[0]=2.96E-5; reta[1]=4.17E-7; reta[2]=1.09E-6;
    List<scalar> retn; retn.setSize(3); retn[0]=6.62; retn[1]=1.17; retn[2]=2.04;
    List<scalar> retm; retm.setSize(3); retm[0]=0.84894; retm[1]=0.14530; retm[2]=0.50980;
    List<scalar> retw; retw.setSize(3); retw[0]=0.891; retw[1]=0.500E-3; retw[2]=0.1085;
    scalar dum1=0; scalar dum2=0; scalar dum3=0; scalar dum4=0;  
    for (int i=0; i<=2; i++)
    {
        dum1=pow( (-reta[i]*pc.internalField()[celli]) , retn[i]);
        dum2=dum2 + retw[i]*(pow( 1+dum1 , -retm[i]));
        dum3=dum3 + retw[i]*reta[i]*(1-pow( (dum1/(1+dum1)) , retm[i]));
        dum4=dum4 + retw[i]*reta[i];
    }    
    
    Krel[celli] = Ks*(pow( dum2 , tau))*(pow( (dum3/dum4) , 2));
}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::HamstadBrick::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/1.57e2); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*30*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for brick" [s]
    
    K_v[celli] = (delta*p_vsat*relhum)/(rho_l*R_v*T.internalField()[celli]);
}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::HamstadBrick::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar L_v = 2.5e6;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    //scalar dpsatdt = (7.06627e3/(T.internalField()[celli]*T.internalField()[celli]) - 5.976/T.internalField()[celli]) * p_vsat; // saturation vapour pressure [Pa]
        
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/1.57e2); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*30*(0.503*tmp*tmp + 0.497)); // Water vapour diffusion coefficient "for brick" [s]

    K_pt[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) ) * (rho_l*L_v - pc.internalField()[celli]);
}

//*********************************************************** //
