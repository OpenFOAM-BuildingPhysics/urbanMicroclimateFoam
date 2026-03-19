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
    Foam::buildingMaterialModels::HamstadBrick

Description
    buildingMaterial constitutive relationship as in HamstadBrick et al. (1977).
   
    Reference:
    HamstadBrick R, Vauclin M, Touma J, Wierenga P, Vachaud G. Comparision of numerical
    simulation models for one dimensional infiltration. buildingMaterial Sci. Soc. Am J. 1977,
    41:285-294

SourceFiles
    HamstadBrick.C

\*---------------------------------------------------------------------------*/

#ifndef HamstadBrick_H
#define HamstadBrick_H

#include "buildingMaterialModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace buildingMaterialModels
{

/*---------------------------------------------------------------------------*\
                           Class HamstadBrick Declaration
\*---------------------------------------------------------------------------*/

class HamstadBrick
:
    public buildingMaterialModel
{
    // Private data

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("HamstadBrick");


    // Constructors

        //- Construct from components
        HamstadBrick
        (
            const word& name, 
            const dictionary& buildingMaterialDict,
            const word& cellZoneModel
        );


    // Destructor

        ~HamstadBrick()
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
