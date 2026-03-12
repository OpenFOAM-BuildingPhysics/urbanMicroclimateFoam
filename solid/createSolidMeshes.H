    // solidNames is defined in urbanMicroclimateFoam.C from controlDict regions dictionary

    PtrList<fvMesh> solidRegions(solidNames.size());

    forAll(solidNames, i)
    {
        Info<< "Create solid mesh for region " << solidNames[i]
            << " for time = " << runTime.name() << nl << endl;

        solidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    solidNames[i],
                    runTime.name(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        // Force calculation of geometric properties to prevent it being done
        // later in e.g. some boundary evaluation
        //(void)solidRegions[i].weights();
        //(void)solidRegions[i].deltaCoeffs();
    }
