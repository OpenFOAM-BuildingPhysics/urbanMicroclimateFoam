{
    dimensionedScalar Dm("Dm",dimensionSet(0,2,-1,0,0,0,0),scalar(2.5e-5)); 
    scalar Sct = 0.7;

    fvScalarMatrix wEqn
    (
        fvm::ddt(rho,w)
      + fvm::div(phi, w)
      - fvm::laplacian(rho*Dm + turb.mut()/Sct, w)
      ==
        Sw
      + veg.Sq() // vegetation energy source term
    );
    wEqn.relax();

    wEqn.solve();

}
