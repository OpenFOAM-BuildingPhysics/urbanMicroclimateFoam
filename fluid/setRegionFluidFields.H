    const fvMesh& mesh = fluidRegions[i];

    rhoThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    volScalarField& w = wFluid[i];    
    volScalarField& gcr = gcrFluid[i];  
    surfaceScalarField& phi = phiFluid[i]; 

    volScalarField& Sh = ShFluid[i];
    volScalarField& Sw = SwFluid[i];   
    volScalarField& Cf = CfFluid[i]; 

    compressible::turbulenceModel& turb = turbulence[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    radiation::radiationModel& rad = radiation[i];
    solarLoad::solarLoadModel& sol = solarLoad[i];
    grass::grassModel& gra = grass[i];
    vegetation::vegetationModel& veg = vegetation[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

    volScalarField& p_rgh = p_rghFluid[i];

    const dimensionedVector& g = gFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];





