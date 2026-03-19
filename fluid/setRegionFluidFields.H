    fvMesh& mesh = fluidRegions[i];

    fluidThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    volScalarField& w = wFluid[i];
    volScalarField& gcr = gcrFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::momentumTransportModel& turb = turbulence[i];
    fluidThermoThermophysicalTransportModel& fluidThermophys = fluidThermophysicalTransport[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    IOMRFZoneList& MRF = MRFfluid[i];
    Foam::fvModels& fvModels = Foam::fvModels::New(mesh);
    Foam::fvConstraints& fvConstraints = Foam::fvConstraints::New(mesh);

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    pressureReference& pRef = pressureReferenceFluid[i];

    radiationModel& rad = radiation[i];
    solarLoad::solarLoadModel& sol = solarLoad[i];
    grass::grassModel& gra = grass[i];
    vegetation::vegetationModel& veg = vegetation[i];

    volScalarField& p_rgh = p_rghFluid[i];

    volVectorField& netForce = netForceFluid[i];

    const dimensionedVector& g = gFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];


    bool blending =
        runTime.controlDict().lookupOrDefault<bool>("blending", false);
    blendingLayer bL(U, thermo.T());
    if (blending)
    {
        bL.initialize();
    }


