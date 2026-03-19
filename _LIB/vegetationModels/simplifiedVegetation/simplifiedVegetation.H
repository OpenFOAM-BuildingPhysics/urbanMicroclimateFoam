/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::simplifiedVegetation

Description
    Simplified Vegetation model

Author
    Lento Manickathan. manickathan@arch.ethz.ch

SourceFiles
    vegetationModel.C

\*---------------------------------------------------------------------------*/

#ifndef vegetationModelsimplifiedVegetation_H
#define vegetationModelsimplifiedVegetation_H

#include "vegetationModel.H"

#include "fvm.H"
#include "IOdictionary.H"
#include "bound.H"
#include "scalarListIOList.H"
#include "scalarIOList.H"
#include "vectorIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam 
{
namespace vegetation
{
/*---------------------------------------------------------------------------*\
                         Class simplifiedVegetation Declaration
\*---------------------------------------------------------------------------*/

class simplifiedVegetation
:
    public vegetationModel
{

    // Private data

        //- debug switch
        bool debug_;

		// Model coefficients
		dimensionedScalar a1_;     // stomatal resistance specific model parameters
		dimensionedScalar a2_;     // stomatal resistance specific model parameters
		dimensionedScalar a3_;     // stomatal resistance specific model parameters
		dimensionedScalar C_;       // proportionality factor of ra
		dimensionedScalar D0_;      // vapour pressure corresponding to minimal resistance
		dimensionedScalar nEvapSides_; // number of sides, leaf evaporates from
		dimensionedScalar H_;       // Height of vegetation
		dimensionedScalar kc_;      // extinction coefficient of radiation
		dimensionedScalar l_;       // characteristic length of leaf
		dimensionedScalar Rg0_;     // global radiation at top of canopy
		dimensionedScalar Rl0_;     // long-wave radiation at top of canopy
		dimensionedScalar rsMin_;   // minimum stomatal resistance
		dimensionedScalar TlMin_;   // lower limit of T
		dimensionedScalar UMin_;    // lower limit of U
                scalar Cd_;
		dimensionedScalar rhoa_;    // denisty of air
		dimensionedScalar cpa_;     // specific heat of air at constant pressure
		dimensionedScalar lambda_;  // latent heat of vaporization


		scalarListIOList divqrsw;
		// -----------------------------------------------
		// Model fields
		volScalarField E_;      // transpiration rate
		volScalarField ev_;     // water vapor pressure
		volScalarField evsat_;  // saturated water vapor pressure
		volScalarField LAD_;    // leaf area density
                volScalarField Tl_;     // leaf temperature
		volScalarField qsat_;   // saturated specific humidity
		volScalarField Qlat_;     // latent heat flux
		volScalarField Qsen_;     // sensible heat flux
		volScalarField ra_;     // aerodynamic resistance
		volScalarField rs_;     // stomatal resistance
		volScalarField rhosat_; // saturated density field
		volScalarField Rg_;     // global radiation
		volScalarField Rn_;     // global radiation density in volume
		volScalarField VPD_;    // vapor pressure deficit
		
		scalar Tl_relax;
		scalar Tl_residualControl;

    // Private Member Functions
        //- Initialise
        void initialise();

        //- calculate saturation vapour pressure
        scalar calc_evsat(double& T);

        //- calculate saturated density of water vapour
        scalar calc_rhosat(double& T);

        //- Disallow default bitwise copy construct
        simplifiedVegetation(const simplifiedVegetation&);

        //- Disallow default bitwise assignment
        void operator=(const simplifiedVegetation&);

public:
    // runtime type info
    TypeName("simplifiedVegetation");

    // constructors
        simplifiedVegetation
        (
            const volScalarField& T
        );

    // destructor
        virtual ~simplifiedVegetation();

    // Member Functions

        // -----------------------------------------------

        // calc radiation
        void radiation();

        // calc aerodynamic, stomatal resistances
        void resistance(volScalarField& magU, volScalarField& T, volScalarField& q, volScalarField& Tl);

        // solve all
        void calculate(volVectorField&U, volScalarField& T, volScalarField& q);

        // -----------------------------------------------
        // source terms

        // return vegetation energy source
        virtual tmp<volScalarField> Sh() const;

        //- return leaf drag coefficient
        virtual tmp<volScalarField> Cf() const;

        // return vegetation specific humidity source
        virtual tmp<volScalarField> Sq() const;

        // -----------------------------------------------

        // read vegetationProperties dictionary
        bool read();
};


} // end namepsace Foam
} // End namespace Foam

#endif


// ************************************************************************* //
