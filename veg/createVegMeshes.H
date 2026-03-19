    // vegNames is defined in urbanMicroclimateFoam.C from controlDict regions dictionary

    PtrList<fvMesh> vegRegions(vegNames.size());

    forAll(vegNames, i)
    {
        Info<< "Create vegetation mesh for region " << vegNames[i]
            << " for time = " << runTime.name() << nl << endl;

        vegRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    vegNames[i],
                    runTime.name(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

    }
