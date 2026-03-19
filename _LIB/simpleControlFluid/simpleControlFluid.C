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

\*---------------------------------------------------------------------------*/

#include "simpleControlFluid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleControlFluid, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::simpleControlFluid::read()
{
    return
        fluidSolutionControl::read()
     && singleRegionConvergenceControl::read();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleControlFluid::simpleControlFluid(fvMesh& mesh, const word& algorithmName)
:
    fluidSolutionControl(mesh, algorithmName),
    singleRegionConvergenceControl
    (
        static_cast<singleRegionSolutionControl&>(*this)
    ),
    initialised_(false)
{
    read();
    printResidualControls();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleControlFluid::~simpleControlFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::simpleControlFluid::run(Time& time)
{
    read();

    Time& runTime = const_cast<Time&>(mesh().time());

    if (initialised_)
    {
        scalar timeValue = time.value();
        label timeIndex = time.timeIndex();
        const convergenceData cd(criteriaSatisfied());

        if (cd.checked && cd.satisfied)
        {
            time.setTime(timeValue,timeIndex+1); //iteration counter needs to continue in case minFluidIteration is not yet reached
            return false;
        }
        else
        {
            storePrevIterFields();
            time.setDeltaT(1); //this is related to the calculation of timestep continuity error
            runTime.loop(); //this is needed to use functionObjects with runTime writeControl, e.g. probes in controlDict
            time.setTime(timeValue,timeIndex+1); //necessary for iter().stream() to get the correct iteration for convergence control
            //time value must stay the same, otherwise wrong ambient value is read
            return true;
        }
    }
    else
    {
        initialised_ = true;
        storePrevIterFields();
    }

    return true;
}

bool Foam::simpleControlFluid::loop(Time& time)
{
    if (!endIfConverged(time))
    {
        storePrevIterFields();
    }

    return time.loop();
}

// ************************************************************************* //
