//Construct gravity flux
dimensionedScalar rhol
(
   "rhol",
   dimMass/dimVolume,
   scalar(1000)
);

surfaceScalarField phiG =  (fvc::interpolate(Krel,"Krel")*rhol*g) & mesh.Sf();


if(pcEqnForm == "mixed")
{
    //using mixed form for moisture equation - see suGWFoam
    dimensionedScalar rDeltaT = 1.0/runTime.deltaT();
    volScalarField pc_ss = rDeltaT*(ws_old - ws_n); //additional source term appearing due to dw/dt

    pc.oldTimeRef() = pc_n; //Tricking OpenFOAM to use pc value from previous Picard iteration

    fvScalarMatrix pcEqn
    (
        Crel*fvm::ddt(pc)
        ==
        fvm::laplacian(Krel+K_v,pc,"laplacian(Krel,pc)")
        +fvc::laplacian(K_pt,Ts,"laplacian(Krel,pc)")
        -fvc::div(phiG)
        +pc_ss
    );
    pcEqn.solve();

    pc.oldTimeRef() = pc_old; //restoring pc.oldTime()
}
else
{
    fvScalarMatrix pcEqn
    (
        Crel*fvm::ddt(pc)
        ==
        fvm::laplacian(Krel+K_v,pc,"laplacian(Krel,pc)")
        +fvc::laplacian(K_pt,Ts,"laplacian(Krel,pc)")
        -fvc::div(phiG)
    );
    pcEqn.solve();
}
