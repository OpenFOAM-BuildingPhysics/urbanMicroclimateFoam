    // fluidNames is defined in urbanMicroclimateFoam.C from controlDict regions dictionary

    PtrList<fvMesh> fluidRegions(fluidNames.size());

    forAll(fluidNames, i)
    {
        Info<< "Create fluid mesh for region " << fluidNames[i]
            << " for time = " << runTime.name() << nl << endl;

        fluidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    fluidNames[i],
                    runTime.name(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
