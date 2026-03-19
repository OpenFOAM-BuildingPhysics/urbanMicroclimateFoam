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
    Foam::radiation::noRadiation

Description
    No radiation - does nothing to energy equation source terms
    (returns zeros)

SourceFiles
    noRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef noSolarLoad_H
#define noSolarLoad_H

#include "solarLoadModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solarLoad
{

/*---------------------------------------------------------------------------*\
                         Class noSolarLoad Declaration
\*---------------------------------------------------------------------------*/

class noSolarLoad
:
    public solarLoadModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        noSolarLoad(const noSolarLoad&);

        //- Disallow default bitwise assignment
        void operator=(const noSolarLoad&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noSolarLoad(const volScalarField& T);

        //- Construct from components
        noSolarLoad(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~noSolarLoad();


    // Member functions

        // Edit

            //- Solve solarLoad equation(s)
            void calculate();

            //- Read solarLoadProperties dictionary
            bool read();

            //- Source term component (for power of T^4)
            tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            tmp<volScalarField::Internal> Ru() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solarLoad
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
