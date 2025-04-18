{
    rho = thermo.rho();

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

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    phiHbyA += phig;

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
            fvOptions(psi, p_rgh, rho.name());
            
        p_rghEqn.setReference
        (
            pControl.refCell(),
            pControl.refValue()
        );

        p_rghEqn.solve();
    }
    
    phi = phiHbyA + p_rghEqn.flux();

    p = p_rgh + rho*gh;

    #include "incompressible/continuityErrs.H"
    
    // Explicitly relax pressure for momentum corrector
    p_rgh.relax();
    
    // Correct the momentum source with the pressure gradient flux
    // calculated from the relaxed pressure
    U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    
    pControl.limit(p);

    // For closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (closedVolume && !thermo.incompressible())
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /fvc::domainIntegrate(psi);
        p_rgh = p - rho*gh;
        p.correctBoundaryConditions();
    }

    rho = thermo.rho();

    rho.relax();

}
