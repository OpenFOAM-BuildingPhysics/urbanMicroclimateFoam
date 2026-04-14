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
    setMultiRegionDeltaT

Description
    Reset the timestep to maintain a constant maximum courant and
    diffusion Numbers. Reduction of time-step is immediate, but
    increase is damped to avoid unstable oscillations.

\*---------------------------------------------------------------------------*/


if(timeStepDecrease == false && runTime.deltaT().value() <= maxDeltaT) //time step is increased (slowly)
{
  Info << "Increasing time step slowly...deltaT = " 
       <<  min
           (
             min
             (
              runTime.deltaT().value()*increase_factor,
              maxDeltaT
             ),
             timeToOutput
           )
       << endl;

  runTime.setDeltaT
  (
        min(
            min
            (
              runTime.deltaT().value()*increase_factor,
              maxDeltaT
            ),
            timeToOutput
        )
  );
}
else if(timeStepDecrease == true) //time step is decreased
{
  if (runTime.deltaT().value() - minDeltaT < VSMALL)
  {
    volScalarField pc_dump
    (
        IOobject
        (
            "pc_dump",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pc
    );
    pc_dump.write();

    volScalarField Ts_dump
    (
        IOobject
        (
            "Ts_dump",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ts
    );
    Ts_dump.write();

    volScalarField ws_dump
    (
        IOobject
        (
            "ws_dump",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ws
    );
    ws_dump.write();

    volScalarField Crel_dump
    (
        IOobject
        (
            "Crel_dump",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Crel
    );
    Crel_dump.write();

    FatalError
      << "Cannot decrease time step further. Dumping fields and exiting..." 
      << abort(FatalError);
  }
  Info << "Decreasing time step...deltaT = " 
       << min
            (
              max
              (
                  runTime.deltaT().value()*decrease_factor,
                  minDeltaT
              ),
              timeToOutput
            )
       << endl;
     

  runTime.setDeltaT
  (
    min
    (
        max
        (
          runTime.deltaT().value()*decrease_factor,
          minDeltaT
        ),
        timeToOutput
    )
  );
}

// ************************************************************************* //
