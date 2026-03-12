{
    rho = thermo.rho();
    rho.relax();

    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    tUEqn.clear();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);

    // Buoyancy flux (kept separate for reuse in velocity correction)
    const surfaceScalarField ghGradRhof
    (
        -ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    phiHbyA += rhorAUf*ghGradRhof;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

    fvScalarMatrix p_rghEqn(p_rgh, dimMass/dimTime);

    // Solve pressure
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        p_rghEqn =
            fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p_rgh)
          ==
            fvModels.source(psi, p_rgh);

        p_rghEqn.setReference
        (
            pRef.refCell(),
            pRef.refValue()
        );

        fvConstraints.constrain(p_rghEqn);

        p_rghEqn.solve();
    }

    phi = phiHbyA + p_rghEqn.flux();

    // Calculate and relax the net pressure-buoyancy force
    netForce.relax
    (
        fvc::reconstruct(ghGradRhof + p_rghEqn.flux()/rhorAUf),
        p_rgh.relaxationFactor()
    );

    // Correct the momentum source with the relaxed net force
    U = HbyA + rAU*netForce;
    U.correctBoundaryConditions();
    fvConstraints.constrain(U);

    p = p_rgh + rho*gh;

    #include "continuityErrs.H"

    // For closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (closedVolume && !thermo.incompressible())
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /fvc::domainIntegrate(psi);
        p_rgh = p - rho*gh;
        p.correctBoundaryConditions();
    }

    // Optionally relax pressure for the thermophysics
    p.relax(); //v12

    rho = thermo.rho();
    rho.relax();
}
