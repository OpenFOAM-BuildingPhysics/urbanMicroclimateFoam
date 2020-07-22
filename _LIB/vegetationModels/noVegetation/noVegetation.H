/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::vegetation::noVegetation

Description
    No vegetation - does nothing to energy equation source terms
    (returns zeros)

SourceFiles
    noVegetation.C

\*---------------------------------------------------------------------------*/

#ifndef noVegetation_H
#define noVegetation_H

#include "vegetationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vegetation
{

/*---------------------------------------------------------------------------*\
                         Class noVegetation Declaration
\*---------------------------------------------------------------------------*/

class noVegetation
:
    public vegetationModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        noVegetation(const noVegetation&);

        //- Disallow default bitwise assignment
        void operator=(const noVegetation&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noVegetation(const volScalarField& T);


    //- Destructor
    virtual ~noVegetation();


    // Member functions

        // Edit

            //- Solve vegetation equation(s)
            void calculate
            (
                volVectorField& U_,
                volScalarField& T_, 
                volScalarField& q_
            );

            //- Read vegetationProperties dictionary
            bool read();

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Sh() const;

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Cf() const;

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Sq() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vegetation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
