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

#include "VanGenuchten.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{
    defineTypeNameAndDebug(VanGenuchten, 0);

    addToRunTimeSelectionTable
    (
        buildingMaterialModel,
        VanGenuchten,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingMaterialModels::VanGenuchten::VanGenuchten
(
    const word& name,
    const dictionary& buildingMaterialDict,
    const word& cellZoneModel
)
:
    buildingMaterialModel(name, buildingMaterialDict, cellZoneModel),
    VanGenuchtenCoeffs_(buildingMaterialDict),
    wcap_
    (
       readScalar
       (
          VanGenuchtenCoeffs_.lookup("wcap")
       )
   ),

    n_
    (
       readScalar
       (
          VanGenuchtenCoeffs_.lookup("n")
       )
   ),

    alpha_
    (
       readScalar
       (
          VanGenuchtenCoeffs_.lookup("alpha")
       )
   ),

   Ks_
   (
       readScalar
       (
          VanGenuchtenCoeffs_.lookup("Ks")
       )
   )
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Correct the buildingMaterial moisture content (cell)
void Foam::buildingMaterialModels::VanGenuchten::update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli)
{

   scalar m_ = 1.0 - 1.0/n_;
   scalar tmp = pow(-alpha_*pc.internalField()[celli], n_);
   w[celli] = wcap_*pow(1+tmp,-m_);
   scalar tmp2 = 1+tmp;
   Crel[celli] = mag(-wcap_*m_*n_*alpha_*pow(tmp2,-1-m_)*pow(-alpha_*pc.internalField()[celli],n_-1));

}

//- Correct the buildingMaterial liquid permeability (cell)
void Foam::buildingMaterialModels::VanGenuchten::update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli)
{
    scalar m_ = 1.0 - 1.0/n_;

    scalar tmp = w.internalField()[celli]/wcap_;
    scalar tmp2 = pow(1-pow(tmp,1/m_), m_);
    Krel[celli] = Ks_*(Foam::sqrt(tmp))*pow(1-tmp2,2);

}

//- Correct the buildingMaterial vapor permeability (cell)
void Foam::buildingMaterialModels::VanGenuchten::update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli)
{
    scalar m_ = 1.0 - 1.0/n_;

    scalar tmp = w.internalField()[celli]/wcap_;
    scalar tmp2 = pow(1-pow(tmp,1/m_),2*m_);
    K_v[celli] = Ks_*(Foam::sqrt(1-tmp))*tmp2;
}

//- Correct the buildingMaterial K_pt (cell)
void Foam::buildingMaterialModels::VanGenuchten::update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli)
{

        scalar m_ = 1.0 - 1.0/n_;
        scalar rho_l = 1.0e3;
        scalar L_v = 2.5e6;

        scalar tmp = w.internalField()[celli]/wcap_;
        scalar tmp2 = pow(1-pow(tmp,1/m_),2*m_);
        scalar Kv = Ks_*(Foam::sqrt(1-tmp))*tmp2;
        K_pt[celli] = (Kv/T.internalField()[celli]) * (rho_l*L_v - pc.internalField()[celli]);
}

//*********************************************************** //
