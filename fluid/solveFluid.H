simpleControlFluid simpleControlFluid(fluidRegions[i]);

runTime.setTime(runTime.value(),1);
TimeState pts(runTime); //store time state

while ((simpleControlFluid.run(runTime) || runTime.timeIndex() <= minFluidIteration) && runTime.timeIndex() <= maxFluidIteration)
{
//  Pressure-velocity SIMPLE corrector
    Info << "Internal iteration number: " << runTime.timeIndex() << endl;

    if(!simple.lookupOrDefault<bool>("frozenFlow", 0))
    {
        p_rgh.storePrevIter();
        rho.storePrevIter();

        //v12
        thermo.T().correctBoundaryConditions();

        //v12: predict before equations (updates alphat = rho*nut/Prt)
        turb.predict();
        fluidThermophys.predict();

        {
            #include "UEqn.H"
            #include "EEqn.H"
            #include "pEqn.H"
            #include "wEqn.H"
        }

        gra.correct(thermo.T(), w, U);

        // Solve vegetation energy balance
        veg.correct(U, thermo.T(), w);

        turb.correct();
        fluidThermophys.correct();
    }

    Info << nl;
}

Info<< nl << fluidRegions[i].name() << " region solution ended in "
    << runTime.timeIndex()-1 << " iterations" << nl << endl;

runTime.TimeState::operator=(pts); //restore time state

//v12
{
    thermo.T().correctBoundaryConditions();
}
Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
rad.correct();
Info << "Updating short-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
sol.correct();

gcr.correctBoundaryConditions();
