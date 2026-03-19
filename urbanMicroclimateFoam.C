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

    Updated for OpenFOAM v12 - regions now defined in controlDict

    Contributions:
    Clément Nevers, clement.nevers@usherbrooke.ca
    Aytac Kubilay, akubilay@ethz.ch
    Andrea Ferrari, andferra@ethz.ch
    Lento Manickathan, lento.manickathan@empa.ch

\*---------------------------------------------------------------------------*/

// Core OpenFOAM includes (replacing fvCFD.H which no longer exists in v12)
#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "findRefCell.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "OSspecific.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "IOMRFZoneList.H"
#include "linear.H"

// Thermo and transport
#include "rhoThermo.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModel.H"
#include "fluidThermophysicalTransportModel.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidThermoThermophysicalTransportModel.H"

// Custom libraries
#include "buildingMaterialModel.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "noRadiation.H"
#include "solarLoadModel.H"
#include "grassModel.H"
#include "simpleControlFluid.H"
#include "blendingLayer.H"
#include "vegetationModel.H"

// Pressure and constraints
#include "pressureReference.H"
#include "fvConstraints.H"
#include "fvModels.H"

// FVC operations
#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcFlux.H"
#include "fvcVolumeIntegrate.H"

// FVM operations
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

// Patch fields
#include "mixedFvPatchFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // Read regions from controlDict
    // Format in controlDict:
    //   regions
    //   {
    //       fluid  (air);           // list of fluid region names
    //       solid  (building);      // list of solid region names
    //       vegetation (vegetation); // list of vegetation region names
    //   }

    const dictionary& regionsDict =
        runTime.controlDict().subDict("regions");

    wordList fluidNames;
    wordList solidNames;
    wordList vegNames;

    if (regionsDict.found("fluid"))
    {
        fluidNames = wordList(regionsDict.lookup("fluid"));
    }

    if (regionsDict.found("solid"))
    {
        solidNames = wordList(regionsDict.lookup("solid"));
    }
    if (regionsDict.found("vegetation"))
    {
        vegNames = wordList(regionsDict.lookup("vegetation"));
    }

    Info<< "Fluid regions: " << fluidNames << endl;
    Info<< "Solid regions: " << solidNames << endl;
    Info<< "Vegetation regions: " << vegNames << endl;

    #include "createFluidMeshes.H"
    #include "createSolidMeshes.H"
    #include "createVegMeshes.H"
    
    Info<< "all create meshes" <<endl;
    
    #include "createFluidFields.H"
    #include "createSolidFields.H"
    #include "createVegFields.H"

    Info<< "all create fiels" <<endl;

    #include "initContinuityErrs.H"
    #include "initSolidContinuityErrs.H"
    #include "readFluidControls.H"
    #include "readSolidControls.H"

    Info<< "before loop" <<endl;

    while (runTime.loop())
    {
        Info<< nl << "Time = " << runTime.name() << endl;

        forAll(fluidRegions, i)
        {
            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
            #include "setRegionFluidFields.H"
            #include "readFluidMultiRegionSIMPLEControls.H"
            #include "solveFluid.H"
        }

        forAll(vegRegions, i)
        {
			Info<< "\nVegetation region found..." << endl;
			#include "setRegionVegFields.H"
			#include "solveVeg.H"
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
