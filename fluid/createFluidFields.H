// Initialise fluid field pointer lists
PtrList<fluidThermo> thermoFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<compressible::momentumTransportModel> turbulence(fluidRegions.size());
PtrList<fluidThermoThermophysicalTransportModel> fluidThermophysicalTransport(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<radiationModel> radiation(fluidRegions.size());
PtrList<solarLoad::solarLoadModel> solarLoad(fluidRegions.size());
PtrList<grass::grassModel> grass(fluidRegions.size());
PtrList<vegetation::vegetationModel> vegetation(fluidRegions.size());
PtrList<volScalarField> wFluid(fluidRegions.size());
PtrList<volScalarField> gcrFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());

PtrList<pressureReference> pressureReferenceFluid(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());

PtrList<volVectorField> netForceFluid(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{

    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoFluid\n" << endl;

    thermoFluid.set
    (
        i,
        fluidThermo::New(fluidRegions[i]).ptr()
    );

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.name(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()
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
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        if (radIO.headerOk())
        {
            radiation.set
            (
                i,
                radiationModel::New(thermoFluid[i].T())
            );
        }
        else
        {
            Info<< "    Radiation model not active: radiationProperties not found"
                << endl;
            radiation.set
            (
                i,
                new radiationModels::noRadiation(thermoFluid[i].T())
            );
        }
    }

    solarLoad.set
    (
        i,
        solarLoad::solarLoadModel::New(thermoFluid[i].T())
    );

    vegetation.set
    (
        i,
        vegetation::vegetationModel::New(thermoFluid[i].T())
    );

    grass.set
    (
        i,
        grass::grassModel::New(thermoFluid[i].T())
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.name(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to wFluid\n" << endl;
    wFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "w",
                runTime.name(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to gcrFluid\n" << endl;
    gcrFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "gcr",
                runTime.name(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.name(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to gFluid\n" << endl;
    gFluid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to turbulence\n" << endl;
    turbulence.set
    (
        i,
        compressible::momentumTransportModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to thermophysical transport\n" << endl;
    fluidThermophysicalTransport.set
    (
        i,
        fluidThermoThermophysicalTransportModel::New
        (
            turbulence[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField("gh", gFluid[i] & fluidRegions[i].C())
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField("ghf", gFluid[i] & fluidRegions[i].Cf())
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.name(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p
    p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

    fluidRegions[i].schemes().setFluxRequired(p_rghFluid[i].name());

    Info<< "    Adding to netForceFluid\n" << endl;
    netForceFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "netForce",
                runTime.name(),
                fluidRegions[i]
            ),
            fvc::reconstruct
            (
                (
                  - ghfFluid[i]*fvc::snGrad(rhoFluid[i])
                  - fvc::snGrad(p_rghFluid[i])
                )*fluidRegions[i].magSf()
            )
        )
    );

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    pressureReferenceFluid.set
    (
        i,
        new pressureReference
        (
            thermoFluid[i].p(),
            p_rghFluid[i],
            fluidRegions[i].solution().dict().subDict("SIMPLE")
        )
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

}


