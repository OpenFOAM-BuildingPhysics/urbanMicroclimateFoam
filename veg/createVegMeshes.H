    const wordList vegNames(rp["vegetation"]);

    PtrList<fvMesh> vegRegions(vegNames.size());

    forAll(vegNames, i)
    {
        Info<< "Create vegetation mesh for region " << vegNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        vegRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    vegNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

    }
