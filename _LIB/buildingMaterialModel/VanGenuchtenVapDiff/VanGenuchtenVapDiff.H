
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

Class
    Foam::buildingMaterialModels::VanGenuchtenVapDiff

Description
    buildingMaterial constitutive relationship given by Van Genuchten relationship
    for water retention curve and liquid permeability. Vapor permeability is given 
    in terms of vapor diffusion coefficient
    (for single pore system material)

    m = 1 - 1/n

    water retention
    w = (wcap-wr) * (1+(alpha*pc)^n)^-m + wr

    liquid permeability
    We = (w-wcap)/(wcap-wr)
    Krel = Ks * (We)^0.5 * (1-(1-(We)^1/m)^m)^2

    vapor permeability
    Kv = (2.61e-5/(muDry*Rv*T))*(1-We)/(A*(1-We)^2 + B)

    B = 1 - A;
    muDry = vapor resistance factor of dry material [-]

    *requires definition of the following additional entries:
         -wcap:  saturated moisture content [kg/m^3]
         -wr:    residual moisture content [kg/m^3]
         -n:     VG exponent [-]
         -alpha: VG coefficient [Pa^-1]
         -Ks:    saturated hydraulic conductivity divided gravity [s]
         -muDry: vapor resistance factor of dry material [-]
         -A:     coefficient [-]

         -minCrel: minimum  Crel value (if not specified it takes VSMALL)


    Author: Andrea Ferrari

SourceFiles
    VanGenuchtenVapDiff.C

\*---------------------------------------------------------------------------*/

#ifndef VanGenuchtenVapDiff_H
#define VanGenuchtenVapDiff_H

#include "buildingMaterialModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{

/*---------------------------------------------------------------------------*\
                           Class VanGenuchtenVapDiff Declaration
\*---------------------------------------------------------------------------*/

class VanGenuchtenVapDiff
:
    public buildingMaterialModel
{
    // Private data

        dictionary VanGenuchtenVapDiffCoeffs_;
        scalar wcap_;
        scalar n_;
        scalar alpha_;
        scalar Ks_;
        scalar muDry_;
        scalar A_;
        scalar wr_;
        scalar minCrel_;

        bool vaporDiffusion_;
    // Private Member Functions

public:

    //- Runtime type information
    TypeName("VanGenuchtenVapDiff");


    // Constructors

        //- Construct from components
        VanGenuchtenVapDiff
        (
            const word& name, 
            const dictionary& buildingMaterialDict,
            const word& cellZoneModel
        );


    // Destructor

        ~VanGenuchtenVapDiff()
        {}


    // Member Functions

        //- Correct the buildingMaterial moisture content (cell)
        void update_w_C_cell(const volScalarField& pc, volScalarField& w, volScalarField& Crel, label& celli);

        //- Correct the buildingMaterial liquid permeability (cell)
        void update_Krel_cell(const volScalarField& pc, const volScalarField& w, volScalarField& Krel, label& celli);

        //- Correct the buildingMaterial vapor permeability (cell)
        void update_Kv_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_v, label& celli);

        //- Correct the buildingMaterial K_pt (cell)
        void update_Kpt_cell(const volScalarField& pc, const volScalarField& w, const volScalarField& T, volScalarField& K_pt, label& celli);      

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace buildingMaterialModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
