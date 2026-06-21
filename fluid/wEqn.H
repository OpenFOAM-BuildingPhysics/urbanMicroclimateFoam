{
    dimensionedScalar Dm("Dm",dimensionSet(0,2,-1,0,0,0,0),scalar(2.5e-5)); 
    scalar Sct = 0.7;

    fvScalarMatrix wEqn
    (
        fvm::ddt(rho,w)
      + fvm::div(phi, w)
      - fvm::laplacian(rho*Dm + rho*turb.nut()/Sct, w)  // rho*nut() = mut() in v12
      ==
        gra.Sw()
      + veg.Sq() // vegetation energy source term
      + fvModels.source(rho, w)
    );
    wEqn.relax();

    fvConstraints.constrain(wEqn);

    wEqn.solve();

    fvConstraints.constrain(w);

}
