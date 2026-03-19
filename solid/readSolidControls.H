/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Global
    readSolidControls

Description
    Read the control parameters used in the solid

\*---------------------------------------------------------------------------*/

scalar increase_factor =
    runTime.controlDict().lookupOrDefault<scalar>("increase_factor", 1.2);

scalar decrease_factor =
    runTime.controlDict().lookupOrDefault<scalar>("decrease_factor", 0.5);    

scalar initialSolidTimestepFactor =
    runTime.controlDict().lookupOrDefault<scalar>("initialSolidTimestepFactor", 0.1);    

int nInternalIterMax =
    runTime.controlDict().lookupOrDefault<int>("nInternalIterMax", 10);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", 600.0); 

scalar minDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", 1.0);  

scalar PicardTolerancews =
    runTime.controlDict().lookupOrDefault<scalar>("PicardTolerancews", 0.01);

scalar PicardToleranceTs =
    runTime.controlDict().lookupOrDefault<scalar>("PicardToleranceTs", 0.01); 

word pcEqnForm =
    runTime.controlDict().lookupOrDefault<word>("pcEqnForm", "pc-based"); 

scalar minCrel =
    runTime.controlDict().lookupOrDefault<scalar>("minCrel", VSMALL); 

// ************************************************************************* //
