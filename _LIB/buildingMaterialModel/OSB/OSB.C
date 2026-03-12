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

#include "OSB.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(OSB, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        OSB,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::OSB::OSB
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
void Foam::buildingMaterialModels::OSB::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
    List<scalar> reta; reta.setSize(2); reta[0] = -6.78e-6; reta[1] = -1.66e-7;
    List<scalar> retn; retn.setSize(2); retn[0] = 3.2225; retn[1] = 1.478;
    List<scalar> retm; retm.setSize(2); retm[0] = 0.6897; retm[1] = 0.3236;
    List<scalar> retw; retw.setSize(2); retw[0] = 0.4676; retw[1] = 0.5324;
    scalar w_tmp = 0; scalar tmp = 0; scalar C_tmp = 0; scalar tmp2 = 0;    
    for (int i=0; i<=1; i++)
    {
        tmp = pow( (reta[i]*(pc.internalField()[celli])) , retn[i] );
        w_tmp = w_tmp + retw[i] / ( pow( (1 + tmp) , retm[i] ));
        tmp2 = pow( (1 + tmp) , retm[i] );
        C_tmp = C_tmp - retw[i]/tmp2 * retm[i]*retn[i]*tmp/((1 + tmp)*(pc.internalField()[celli])); 
    }
    w[celli] = w_tmp*263;
    Crel[celli] = mag( C_tmp*263);   
}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::OSB::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    scalar logpc = log10(-pc.internalField()[celli]);
    scalar logKl = 0;
    int i;
    List<scalar> logpc_M = { 2.13615, 2.43718, 3.43719, 3.73822, 3.91432, 4.03926, 
                        4.13618, 4.43723, 4.61334, 4.7383, 4.7797, 4.80821, 4.83523, 
                        4.91444, 4.9814, 5.03942, 5.06576, 5.09059, 5.11408, 5.13637, 
                        5.17779, 5.2156, 5.21564, 5.4221, 5.56146, 5.6668, 5.75151, 5.82236, 
                        6.22123, 6.42558, 6.44159, 6.61987, 7.05723, 7.08472, 7.15884, 7.2428, 
                        7.31462, 7.37759, 7.43381, 7.48474, 7.53141, 7.5746, 7.61487, 7.6527, 
                        7.68843, 7.72236, 7.75474, 7.78577, 7.81562, 7.84443, 7.87233, 7.89945, 7.92587, 7.95168, 7.97698, 8.00183 };
    
    List<scalar> logKl_M = { -9.44224, -9.44224, -9.44224, -9.44224, -9.44224, -9.44224,
                        -9.44224, -9.44224, -9.44224, -9.44224, -9.44224, -9.44224, -9.44224,
                        -9.44224, -9.44224, -9.44224, -9.44224, -9.44224, -9.44523, -9.45335,
                        -9.48228, -9.52417, -9.52422, -9.9593, -10.61, -11.467, -11.984, -12.3339,
                        -13.409, -13.7882, -13.8231, -14.2848, -15.7508, -15.8401, -16.0724, -16.3215,
                        -16.5245, -16.6975, -16.8504, -16.9898, -17.1202, -17.2447, -17.3659, -17.4859, 
                        -17.6064, -17.7291, -17.8557, -17.9879, -18.1274, -18.2759, -18.4346, -18.6026, -18.7735, -18.9283, -19.0296, -19.0414 };

    if (logpc < logpc_M.first())
    {
        i = 0;
        logKl = logKl_M[i] + (((logKl_M[i + 1] - logKl_M[i]) / (logpc_M[i + 1] - logpc_M[i]))*(logpc - logpc_M[i]));
    }
    else if (logpc >= logpc_M.last())
    {
        i = logpc_M.size()-2;
        logKl = logKl_M[i] + (((logKl_M[i + 1] - logKl_M[i]) / (logpc_M[i + 1] - logpc_M[i]))*(logpc - logpc_M[i]));
    }
    else
    {
        for (i = 0; i < logpc_M.size()-1; ++i)
        {
            if ((logpc_M[i] <= logpc) && (logpc < logpc_M[i + 1]))
            {
                logKl = logKl_M[i] + (((logKl_M[i + 1] - logKl_M[i]) / (logpc_M[i + 1] - logpc_M[i]))*(logpc - logpc_M[i]));
                break;
            }
        }
    }
    Krel[celli] = pow(10, logKl);
}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::OSB::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{
    scalar rho_l = 1.0e3;
    scalar R_v = 8.31451 * 1000 / (18.01534);
    scalar a = 2.61e-5;
    scalar b = 264;
    scalar c = 0.503;
    scalar d = 0.497;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]

    scalar tmp = 1 - (w.internalField()[celli] / 263);
    scalar delta = a * tmp / (R_v*T.internalField()[celli] * b * (c*tmp*tmp + d)); // Water vapour diffusion coefficient "for OSB" [s]

    K_v[celli] = (delta*p_vsat*relhum) / (rho_l*R_v*T.internalField()[celli]);
}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::OSB::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{
    scalar rho_l = 1.0e3;
    scalar R_v = 8.31451 * 1000 / (18.01534);
    scalar L_v = 2.5e6;
    scalar a = 2.61e-5;
    scalar b = 264;
    scalar c = 0.503;
    scalar d = 0.497;

    scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]
    //scalar dpsatdt = (7.06627e3/(T.internalField()[celli]*T.internalField()[celli]) - 5.976/T.internalField()[celli]) * p_vsat; // saturation vapour pressure [Pa]
        
    scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]

    scalar tmp = 1 - (w.internalField()[celli] / 263);
    scalar delta = a * tmp / (R_v*T.internalField()[celli] * b * (c*tmp*tmp + d)); // Water vapour diffusion coefficient "for OSB" [s]

    K_pt[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) ) * (rho_l*L_v - pc.internalField()[celli]);
}

//*********************************************************** //
