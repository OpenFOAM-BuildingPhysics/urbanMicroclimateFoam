fvMesh& mesh = solidRegions[i];

Info<< "Reading buildingMaterial Information\n" << endl;

PtrList<dictionary> Materials(solidTransportProperties[i].lookup("buildingMaterials"));

volScalarField& rho_m = rho_mSolid[i];
volScalarField& cap_m = cap_mSolid[i];
volScalarField& lambda_m = lambda_mSolid[i];
volScalarField& Ts = TsSolid[i];
volScalarField& ws = wsSolid[i];
volScalarField& pc = pcSolid[i];
volScalarField& Krel = KrelSolid[i];
volScalarField& K_v = K_vSolid[i];
volScalarField& K_pt = K_ptSolid[i];
volScalarField& Crel = CrelSolid[i];
const uniformDimensionedVectorField& g = gSolid[i];

