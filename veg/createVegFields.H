// Initialise vegetation field pointer lists
PtrList<radiationModel> radiation2(vegRegions.size());
PtrList<solarLoad::solarLoadModel> solarLoad2(vegRegions.size());
PtrList<volScalarField> TVeg(vegRegions.size());

// Populate vegetation field pointer lists
forAll(vegRegions, i)
{
    Info<< "*** Reading vegetation mesh properties for region "
        << vegRegions[i].name() << nl << endl;

    Info<< "    Adding to TVeg\n" << endl;
    TVeg.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.name(),
                vegRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            vegRegions[i]
        )
    );

    // Check if radiationProperties exists; if not, use noRadiation
    {
        typeIOobject<IOdictionary> radIO
        (
            IOobject
            (
                "radiationProperties",
                runTime.constant(),
                vegRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        if (radIO.headerOk())
        {
            radiation2.set
            (
                i,
                radiationModel::New(TVeg[i])
            );
        }
        else
        {
            Info<< "    Radiation model not active for veg region: "
                << "radiationProperties not found" << endl;
            radiation2.set
            (
                i,
                new radiationModels::noRadiation(TVeg[i])
            );
        }
    }

    solarLoad2.set
    (
        i,
        solarLoad::solarLoadModel::New(TVeg[i])
    );

}

