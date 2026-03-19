/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::grass::simpleGrass

Description
    Grass model intended for grass surfaces and green roof/facades
    Written by Aytac Kubilay

SourceFiles
    simpleGrass.C

\*---------------------------------------------------------------------------*/

#ifndef grassModelsimpleGrass_H
#define grassModelsimpleGrass_H

#include "grassModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace grass
{

/*---------------------------------------------------------------------------*\
                           Class simpleGrass Declaration
\*---------------------------------------------------------------------------*/

class simpleGrass
:
    public grassModel
{
    // Private data

        //- grass leaf temperature [K]
        volScalarField Tg_;

        volScalarField Sw_;
        volScalarField Sh_;
        volScalarField Cf_;     

        //- Selected patches
        labelList selectedPatches_;

        scalar nEvapSides_; // number of sides, leaf evaporates from
        scalar Cd_; // leaf drag coefficient
        scalar beta_; // extinction coefficient for short-wave radiation
        scalar betaLW_; // extinction coefficient for long-wave radiation
        scalar LAI_; // leaf area index
        scalar l_; // characteristic leaf size
        scalar albedoSoil_; // short-wave reflectivity for soil
        scalar emissivitySoil_; // long-wave reflectivity for soil

        scalar p_; // air pressure - Pa
        scalar rhoa; // density of air - kg/m3
        scalar cpa; // specific heat of air at constant pressure - J/(kgK)
        scalar Ra; //gas constant for air
        scalar Rv; //gas constant for vapour

        scalar rs; //stomatal resistance - s/m
        scalar ra; //aerodynamic resistance - s/m
        
        //- debug switch
        bool debug_;

        //- Grass patch ID
        label grassPatchID;

        scalar Tg_relax;
        scalar Tg_residualControl;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- calculate saturation vapour poressure
        scalarField calc_pvsat(const scalarField& T_);

        //- Disallow default bitwise copy construct
        simpleGrass(const simpleGrass&);

        //- Disallow default bitwise assignment
        void operator=(const simpleGrass&);


public:

    //- Runtime type information
    TypeName("simpleGrass");


    // Constructors

        //- Construct from components
        simpleGrass(const volScalarField& T);


    //- Destructor
    virtual ~simpleGrass();


    // Member functions

        // Edit

            //- Solve grass equation(s)
            void calculate
            (
                const volScalarField& T_, 
                const volScalarField& w_,
                const volVectorField& U_
            );
            
            // -----------------------------------------------
            // source terms

            // return vegetation energy source
            virtual tmp<volScalarField> Sh() const;

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Cf() const;

            // return vegetation specific humidity source
            virtual tmp<volScalarField> Sw() const;

            // -----------------------------------------------            

            //- Read grass properties dictionary
            bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace grass
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
