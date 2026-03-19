//Solid-region continuity errors

dimensionedScalar rhol
(
   "rhol",
   dimMass/dimVolume,
   scalar(1000)
);

//calc continuity errors
dimensionedScalar totalWs = fvc::domainIntegrate(ws); // [Kg]
surfaceScalarField phiK = fvc::interpolate(Krel+K_v,"Krel")*fvc::snGrad(pc)*mesh.magSf();
surfaceScalarField phiT = fvc::interpolate(K_pt,"Krel")*fvc::snGrad(Ts)*mesh.magSf();
surfaceScalarField phiG = (fvc::interpolate(Krel, "Krel")*rhol*g) & mesh.Sf();
surfaceScalarField phi = phiK + phiT - phiG; //[Kg/s]
volScalarField contErr = (ws-ws_old) - fvc::div(phi)*runTime.deltaT(); //[Kg/m³]

scalar sumLocalContErr = (fvc::domainIntegrate(mag(contErr))/totalWs).value(); //[-]
scalar globalContErr = (fvc::domainIntegrate(contErr)/totalWs).value();
cumulativeSolidContErr += globalContErr;

Info << "Time step continuity errors: sum local = " << sumLocalContErr
     << ", global = " << globalContErr
     << ", cumulative = " << cumulativeSolidContErr
     << endl;


