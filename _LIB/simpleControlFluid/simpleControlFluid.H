/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::simpleControlFluid

Description
    Simple control class. Provides time-loop control methods which exit the
    simulation once convergence criteria have been reached. Example usage:

    \verbatim
    simpleControlFluid simple(mesh);

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // solve ...
    }
    \endverbatim

SourceFiles
    simpleControlFluid.C

\*---------------------------------------------------------------------------*/

#ifndef simpleControlFluid_H
#define simpleControlFluid_H

#include "fluidSolutionControl.H"
#include "singleRegionConvergenceControl.H"

#define SIMPLE_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class simpleControlFluid Declaration
\*---------------------------------------------------------------------------*/

class simpleControlFluid
:
    public fluidSolutionControl,
    public singleRegionConvergenceControl
{
protected:

    // Protected Data

        //- Initialised flag
        bool initialised_;

        virtual bool read();

public:

    // Static data members

        //- Run-time type information
        TypeName("simpleControlFluid");


    // Constructors

        //- Construct from a mesh and the name of the algorithm
        simpleControlFluid(fvMesh& mesh, const word& algorithmName="SIMPLE");


    //- Destructor
    virtual ~simpleControlFluid();


    // Member Functions

        // Evolution

            //- Time run loop
            bool run(Time& time);

            bool loop(Time& time);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
