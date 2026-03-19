    // Solve the Momentum equation

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turb.divDevTau(U)
     ==
        fvModels.source(rho, U)
      - fvm::SuSp(gra.Cf()*rho*mag(U), U)
      - fvm::SuSp(veg.Cf()*rho*mag(U), U) // vegetation momentum source
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    if (blending)
    {
        UEqn -= rho*bL.bL_USource(U);
    }

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    solve
    (
        UEqn
     ==
        netForce
    );

    fvConstraints.constrain(U);
