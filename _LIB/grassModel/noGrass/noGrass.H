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
    Foam::grass::noGrass

Description
    No grass - does nothing to energy equation source terms
    (returns zeros)

SourceFiles
    noGrass.C

\*---------------------------------------------------------------------------*/

#ifndef noGrass_H
#define noGrass_H

#include "grassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace grass
{

/*---------------------------------------------------------------------------*\
                         Class noGrass Declaration
\*---------------------------------------------------------------------------*/

class noGrass
:
    public grassModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        noGrass(const noGrass&);

        //- Disallow default bitwise assignment
        void operator=(const noGrass&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noGrass(const volScalarField& T);


    //- Destructor
    virtual ~noGrass();


    // Member functions

        // Edit

            //- Solve grass equation(s)
            void calculate
            (
                const volScalarField& T_, 
                const volScalarField& w_,
                const volVectorField& U_
            );

            //- Read grassProperties dictionary
            bool read();
            
            //- return leaf drag coefficient
            virtual tmp<volScalarField> Sh() const;

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Cf() const;

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Sw() const;            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace grass
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
