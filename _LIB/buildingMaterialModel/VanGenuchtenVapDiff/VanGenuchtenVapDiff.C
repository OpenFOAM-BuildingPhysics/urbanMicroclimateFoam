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

#include "VanGenuchtenVapDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(VanGenuchtenVapDiff, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        VanGenuchtenVapDiff,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::VanGenuchtenVapDiff::VanGenuchtenVapDiff
(
    const word& name,
    const dictionary& buildingMaterialDict,
    const word& cellZoneModel
)
:
    buildingMaterialModel(name, buildingMaterialDict, cellZoneModel),
    VanGenuchtenVapDiffCoeffs_(buildingMaterialDict),
    wcap_
    (
       readScalar
       (
          VanGenuchtenVapDiffCoeffs_.lookup("wcap")
       )
    ),

    n_
    (
       readScalar
       (
          VanGenuchtenVapDiffCoeffs_.lookup("n")
       )
    ),

    alpha_
    (
       readScalar
       (
          VanGenuchtenVapDiffCoeffs_.lookup("alpha")
       )
    ),

   Ks_
   (
       readScalar
       (
          VanGenuchtenVapDiffCoeffs_.lookup("Ks")
       )
    ),

    muDry_
   (
       readScalar
       (
          VanGenuchtenVapDiffCoeffs_.lookup("muDry")
       )
   ),
   
   A_
   (
       readScalar
       (
          VanGenuchtenVapDiffCoeffs_.lookup("A")
       )
   ),
   
   wr_
   (
       readScalar
       (
          VanGenuchtenVapDiffCoeffs_.lookup("wr")
       )
   ),
   minCrel_(VanGenuchtenVapDiffCoeffs_.lookupOrDefault<scalar>("minCrel",VSMALL))
   



{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Correct the buildingMaterial moisture content (cell)
void Foam::buildingMaterialModels::VanGenuchtenVapDiff::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{
  
   scalar pci = pc.internalField()[celli];

   scalar m_ = 1.0 - 1.0/n_;
   scalar tmp = pow(-alpha_*pci, n_);
   w[celli] = (wcap_ - wr_)*pow(1+tmp,-m_) + wr_;
   scalar tmp2 = 1+tmp;

   Crel[celli] = max(mag(-(wcap_-wr_)*m_*n_*alpha_*pow(tmp2,-1-m_)*pow(-alpha_*pci,n_-1)),minCrel_);

}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::VanGenuchtenVapDiff::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    scalar m_ = 1.0 - 1.0/n_;

    scalar tmp = (w.internalField()[celli]-wr_)/(wcap_-wr_);
    scalar tmp1 = pow(tmp,1/m_);
    scalar tmp2 = pow(1-tmp1, m_);
    Krel[celli] = Ks_*(Foam::sqrt(tmp))*pow(1-tmp2,2);

}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::VanGenuchtenVapDiff::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{

    if(muDry_ == 0.0)
        {
           Info << "Specify mudry != 0.0 or use VanGenuchten" << endl;
           Foam::FatalError.exit();
        }
    else
        {
           scalar rho_l = 1.0e3;
           scalar R_v = 8.31451*1000/(18.01534);
           A_ = max(1.0,A_);
           scalar B_ = 1.0 - A_;

           scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]

           scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]

           scalar tmp3 = 1 - ((w.internalField()[celli]-wr_)/(wcap_-wr_));
           scalar delta = 2.61e-5 * tmp3/(R_v*T.internalField()[celli]*muDry_*(A_*tmp3*tmp3 + B_)); // Water vapour diffusion coefficient
           K_v[celli] = (delta*p_vsat*relhum)/(rho_l*R_v*T.internalField()[celli]);
    }

}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::VanGenuchtenVapDiff::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{

        if(muDry_ == 0.0)
        {
           Info << "Specify mudry != 0.0 or use VanGenuchten" << endl;
           Foam::FatalError.exit();
        }
        else
        {
           scalar rho_l = 1.0e3;
           scalar R_v = 8.31451*1000/(18.01534);
           scalar L_v = 2.5e6;
           A_ = max(1.0,A_);
           scalar B_ = 1.0 - A_;

           scalar p_vsat = Foam::exp(6.58094e1 - 7.06627e3/T.internalField()[celli] - 5.976*Foam::log(T.internalField()[celli])); // saturation vapour pressure [Pa]

           scalar relhum = Foam::exp(pc.internalField()[celli]/(rho_l*R_v*T.internalField()[celli])); // relative humidity [-]

           scalar tmp3 = 1 - ((w.internalField()[celli]-wr_)/(wcap_-wr_));
           scalar delta = 2.61e-5 * tmp3/(R_v*T.internalField()[celli]*muDry_*(A_*tmp3*tmp3 + B_)); // Water vapour diffusion coefficient
           K_pt[celli] = ( (delta*p_vsat*relhum)/(rho_l*R_v*pow(T.internalField()[celli],2)) )*(rho_l*L_v - pc.internalField()[celli]);
        }

}

//*********************************************************** //
