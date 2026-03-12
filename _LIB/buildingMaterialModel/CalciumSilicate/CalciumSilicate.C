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

#include "CalciumSilicate.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(CalciumSilicate, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        CalciumSilicate,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::CalciumSilicate::CalciumSilicate
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
void Foam::buildingMaterialModels::CalciumSilicate::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    scalar A = 0.004342; scalar n = 0.741839763;
    scalar rhol = 1000; scalar Rv = 8.31451*1000/(18.01534); scalar T = 293.15;
    
    scalar rh = Foam::exp(pc.internalField()[celli]/(rhol*Rv*T));
    scalar wcap = 793;
    
    w[celli] = wcap*pow( 1-log(rh)/A , (-1/n));
    
    scalar rh2 = Foam::exp((pc.internalField()[celli]+100)/(rhol*Rv*T));
    scalar w2 = wcap*pow( 1-log(rh2)/A , (-1/n));
    Crel[celli] = (w2-w.internalField()[celli])/100;   
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::CalciumSilicate::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    scalar logpc = log10(-pc.internalField()[celli]);
    scalar logKl = 0;
    int i;
    double logpc_M[]={2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
         3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
         4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9,
         5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9,
         6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9,
         7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9,
         8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9,
         9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9,
         10}; 
    double logKl_M[]={-7.861221861,-7.861450776,-7.861738933,-7.862101654,-7.862558218,-7.863132879,-7.863856146,-7.864766387,-7.865911841,-7.867353134,
         -7.869166430,-7.871447360,-7.874315914,-7.877922517,-7.882455540,-7.888150548,-7.895301619,-7.904275063,-7.915525911,-7.929617429,
         -7.947243788,-7.969255703,-7.996688363,-8.030790085,-8.073048923,-8.125212606,-8.189294878,-8.267558515,-8.362462600,-8.476560026,
         -8.612332437,-8.771956239,-8.957007275,-9.168133833,-9.404754592,-9.664860432,-9.945002580,-10.24052138,-10.54600899,-10.85592415,
         -11.16522125,-11.46984897,-11.76702093,-12.05523832,-12.33411469,-12.60409074,-12.86612561,-13.12142686,-13.37124856,-13.61676173,
         -13.85898588,-14.09876470,-14.33676965,-14.57351869,-14.80940090,-15.04470187,-15.27962674,-15.51431980,-15.74888026,-15.98337462,
         -16.21784611,-16.45232164,-16.68681694,-16.92134027,-17.15589498,-17.39048136,-17.62509786,-17.85974196,-18.09441068,-18.32910092,
         -18.56380973,-18.79853435,-19.03327231,-19.26802145,-19.50277989,-19.73754602,-19.97231850,-20.20709618,-20.44187811,-20.67666351,
         -20.92450516};

    if (logpc < scalar(2))
    {
        i = 0;
        logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
    }
    else if (logpc >= scalar(10))
    {
        i = 79;
        logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
    }
    else
    {
        for (i=0; i<=79; ++i)
        {
            if ( (logpc_M[i] <= logpc) && (logpc < logpc_M[i+1]) )
            {
                logKl = logKl_M[i] + (((logKl_M[i+1] - logKl_M[i])/(logpc_M[i+1] - logpc_M[i]))*(logpc - logpc_M[i])) ;
                break;
            }
        }
    }
    Krel[celli] = pow(10,logKl);
}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::CalciumSilicate::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/793); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*3.8*(0.5*tmp*tmp + 0.5)); // Water vapour diffusion coefficient "for brick" [s]
    
    K_v[celli] = (delta*p_vsat*relhum)/(rho_l*R_v*T.internalField()[celli]);
}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::CalciumSilicate::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{
    scalar rho_l = 1.0e3; 
    scalar R_v = 8.31451*1000/(18.01534); 
    scalar L_v = 2.5e6;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    //scalar dpsatdt = (7.06627e3/(T.internalField()[celli]*T.internalField()[celli]) - 5.976/T.internalField()[celli]) * p_vsat; // saturation vapour pressure [Pa]
        
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]
    
    scalar tmp = 1 - (w.internalField()[celli]/793); 
    scalar delta = 2.61e-5 * tmp/(R_v*T.internalField()[celli]*3.8*(0.5*tmp*tmp + 0.5)); // Water vapour diffusion coefficient "for brick" [s]

    K_pt[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) ) * (rho_l*L_v - pc.internalField()[celli]);
}

//*********************************************************** //
