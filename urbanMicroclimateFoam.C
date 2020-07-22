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

Application
    urbanMicroclimateFoam

Description
    Solves for air flow (CFD) and transport in porous building materials (HAM)
    Written by Aytac Kubilay, December 2015, ETH Zurich/Empa
    
    Contributions:
    Aytac Kubilay, akubilay@ethz.ch
    Andrea Ferrari, andferra@ethz.ch
    Lento Manickathan, lento.manickathan@empa.ch

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "buildingMaterialModel.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "solarLoadModel.H"
#include "grassModel.H"
#include "simpleControlFluid.H"
#include "fvOptions.H"

#include "vegetationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    regionProperties rp(runTime);

    #include "createFluidMeshes.H"
    #include "createSolidMeshes.H"
    #include "createVegMeshes.H"

    #include "createFluidFields.H"
    #include "createSolidFields.H"
    #include "createVegFields.H"

    #include "initContinuityErrs.H"
    #include "initSolidContinuityErrs.H"
    #include "readFluidControls.H"
    #include "readSolidControls.H"

    while (runTime.loop())
    {
        Info<< nl << "Time = " << runTime.timeName() << endl;

        forAll(vegRegions, i)
        {
			Info<< "\nVegetation region found..." << endl;
			#include "setRegionVegFields.H"
			#include "solveVeg.H"
        }

        forAll(fluidRegions, i)
        {
            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
            #include "setRegionFluidFields.H"
            #include "readFluidMultiRegionSIMPLEControls.H"
            #include "solveFluid.H"
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        forAll(solidRegions, i)
        {
            Info<< "\nSolving for solid region "
                << solidRegions[i].name() << endl;
            #include "setRegionSolidFields.H"
            #include "solveSolid.H"
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
