dimensionedScalar cap_l("cap_l",dimensionSet(0,2,-2,-1,0,0,0),scalar(4182));
dimensionedScalar cap_v("cap_v",dimensionSet(0,2,-2,-1,0,0,0),scalar(1880));
dimensionedScalar Tref("Tref",dimensionSet(0,0,0,1,0,0,0),scalar(273.15));
dimensionedScalar L_v("L_v",dimensionSet(0,2,-2,0,0,0,0), 2.5e6);//   +(cap_l.value()-cap_v.value())*Tref.value());

volScalarField C_t = rho_m*cap_m + ws*cap_l;

// enthalpy-gravity flux
surfaceScalarField phiGT = (cap_l*fvc::interpolate(Ts-Tref)*fvc::interpolate(Krel,"Krel")*rhol*g) & mesh.Sf();

if(pcEqnForm == "mixed")
{
   dimensionedScalar rDeltaT = 1.0/runTime.deltaT();
   volScalarField Ts_sn = Ts_n*(rho_m*cap_m + ws*cap_l);
   volScalarField Ts_st = Ts_old*(rho_m*cap_m + ws*cap_l);
   volScalarField Ts_ss = rDeltaT*(Ts_st - Ts_sn);

   Ts.oldTimeRef() = Ts_n;  //using Ts from previous Picard iteration

   fvScalarMatrix TsEqn
   (
     C_t*fvm::ddt(Ts)
     ==
     fvm::laplacian(lambda_m,Ts,"laplacian(Krel,pc)")
     + fvm::laplacian(((Ts-Tref)*cap_v+L_v)*K_pt,Ts,"laplacian(Krel,pc)")
     + fvc::laplacian((Ts-Tref)*cap_l*Krel,pc,"laplacian(Krel,pc)")
     + fvc::laplacian(((Ts-Tref)*cap_v+L_v)*K_v,pc,"laplacian(Krel,pc)")
     - fvc::div(phiGT)
     + Ts_ss
   );

   TsEqn.solve();
   Ts.oldTimeRef() = Ts_old;

}
else
{
   fvScalarMatrix TsEqn
   (
     C_t*fvm::ddt(Ts)
     ==
     fvm::laplacian(lambda_m,Ts,"laplacian(Krel,pc)")
     + fvm::laplacian(((Ts-Tref)*cap_v+L_v)*K_pt,Ts,"laplacian(Krel,pc)")
     + fvc::laplacian((Ts-Tref)*cap_l*Krel,pc,"laplacian(Krel,pc)")
     + fvc::laplacian(((Ts-Tref)*cap_v+L_v)*K_v,pc,"laplacian(Krel,pc)")
     - fvc::div(phiGT)
   );

   TsEqn.solve();


}
