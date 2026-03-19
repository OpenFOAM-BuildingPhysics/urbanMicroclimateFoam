Info << "Updating T boundary fields..." << endl;
vegT.correctBoundaryConditions();

Info << "Updating long-wave radiation heat transfer for region: " << vegRegions[i].name() << endl;
rad.correct();
Info << "Updating short-wave radiation heat transfer for region: " << vegRegions[i].name() << endl;
sol.correct();

