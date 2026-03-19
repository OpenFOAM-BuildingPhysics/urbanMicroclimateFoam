{

    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      + fluidThermophys.divq(he)
     ==
        rho*(U&g)
      + rad.Sh(thermo, he)
      + sol.Sh(thermo, he)
      + gra.Sh()
      + veg.Sh() // vegetation energy source term
      + fvModels.source(rho, he)
    );

    if (blending)
    {
        EEqn -= rho*thermo.Cp()*bL.bL_TSource(thermo.T());
    }

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
