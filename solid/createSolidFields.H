// Initialise solid field pointer lists
PtrList<IOdictionary> solidTransportProperties(solidRegions.size());
PtrList<volScalarField> rho_mSolid(solidRegions.size());
PtrList<volScalarField> cap_mSolid(solidRegions.size());
PtrList<volScalarField> lambda_mSolid(solidRegions.size());
PtrList<volScalarField> TsSolid(solidRegions.size());
PtrList<volScalarField> wsSolid(solidRegions.size());
PtrList<volScalarField> pcSolid(solidRegions.size());
PtrList<volScalarField> KrelSolid(solidRegions.size());
PtrList<volScalarField> K_vSolid(solidRegions.size());
PtrList<volScalarField> K_ptSolid(solidRegions.size());
PtrList<volScalarField> CrelSolid(solidRegions.size());
PtrList<uniformDimensionedVectorField> gSolid(solidRegions.size());


// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "    Adding solidTransportProperties\n" << endl;
    solidTransportProperties.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to rho_mSolid\n" << endl;
    rho_mSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho_m",
                runTime.name(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("rho_m", dimensionSet(1,-3,0,0,0,0,0), scalar(1E-3)),
            ("zeroGradient")
        )
    );      

    Info<< "    Adding to cap_mSolid\n" << endl;
    cap_mSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "cap_m",
                runTime.name(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("cap_m", dimensionSet(0,2,-2,-1,0,0,0), scalar(1E-3)),
            ("zeroGradient")
        )
    );    

    Info<< "    Adding to lambda_mSolid\n" << endl;
    lambda_mSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "lambda_m",
                runTime.name(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("lambda_m", dimensionSet(1,1,-3,-1,0,0,0), scalar(1E-3)),
            ("zeroGradient")
        )
    );

    Info<< "    Adding to TsSolid\n" << endl;
    TsSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Ts",
                runTime.name(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    ); 

    Info<< "    Adding to wsSolid\n" << endl;
    wsSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "ws",
                runTime.name(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );                               

    Info<< "    Adding to pcSolid\n" << endl;
    pcSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "pc",
                runTime.name(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );  

    Info<< "    Adding to KrelSolid\n" << endl;
    KrelSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Krel",
                runTime.name(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("Krel", dimensionSet(0, 0, 1, 0, 0, 0, 0), scalar(1E-3)),
            ("zeroGradient")
        )
    );  

    Info<< "    Adding to K_vSolid\n" << endl;
    K_vSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "K_v",
                runTime.name(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("K_v", dimensionSet(0, 0, 1, 0, 0, 0, 0), scalar(1E-3)),
            ("zeroGradient")
        )
    ); 

    Info<< "    Adding to K_ptSolid\n" << endl;
    K_ptSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "K_pt",
                runTime.name(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("K_pt", dimensionSet(1, -1, -1, -1, 0, 0, 0), scalar(1E-3)),
            ("zeroGradient")
        )
    ); 

    Info<< "    Adding to CrelSolid\n" << endl;
    CrelSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Crel",
                runTime.name(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("Crel", dimensionSet(0, -2, 2, 0, 0, 0, 0), scalar(1E-3)),
            ("zeroGradient")
        )
    );

    Info<< "    Adding to gSolid\n" << endl;
    gSolid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ) 
        )
    );

}
