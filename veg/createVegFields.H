// Initialise fluid field pointer lists
PtrList<rhoThermo> thermoVeg(vegRegions.size());
PtrList<radiation::radiationModel> radiation2(vegRegions.size());
PtrList<solarLoad::solarLoadModel> solarLoad2(vegRegions.size());
PtrList<volScalarField> TVeg(vegRegions.size());

// Populate fluid field pointer lists
forAll(vegRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
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
                runTime.timeName(),
                vegRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),  
            vegRegions[i]
        )   
    );  

    radiation2.set
    (
        i,
        radiation::radiationModel::New(TVeg[i])
    );
    
    solarLoad2.set
    (
        i,
        solarLoad::solarLoadModel::New(TVeg[i])
    );

}

