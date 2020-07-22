
    // Read vegetation+air mesh
    fvMesh vegMesh
    (
      IOobject
      (
          "vegetation",//"radiation",
          runTime.timeName(),
          runTime,
          IOobject::MUST_READ
      )
    );

    volScalarField vegQr
    (
      IOobject
      (
          "qr",
          runTime.timeName(),
          vegMesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      ),
      vegMesh
    );

    scalarListIOList kcLAIboundaryList
    (
        IOobject
        (
            "kcLAIboundary",
            vegMesh.facesInstance(),
            vegMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        sunPosVector.size()
    );

    labelListIOList vegFinalAgglom
    (
        IOobject
        (
            "finalAgglom",
            vegMesh.facesInstance(),
            vegMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    singleCellFvMesh vegCoarseMesh
    (
        IOobject
        (
            vegMesh.name(),
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vegMesh,
        vegFinalAgglom
    );

    Info << "done !!! " << endl;

    const polyBoundaryMesh& vegPatches = vegMesh.boundaryMesh();
    const polyBoundaryMesh& vegCoarsePatches = vegCoarseMesh.boundaryMesh();

    const volScalarField::Boundary& vegQrb = vegQr.boundaryField();

    labelList vegViewFactorsPatches(vegPatches.size());
    labelList vegHowManyCoarseFacesPerPatch(vegPatches.size());
    label count = 0;
    label countAll = 0;
    label vegNCoarseFaces = 0;
    label vegNCoarseFacesAll = 0;   //Also includes non-wall faces with greyDiffusive boundary

    //forAll(LADb, patchI)
    forAll(vegQrb, patchI)
    {
        const polyPatch& pp = vegPatches[patchI];
        //const fvPatchScalarField& LADpI = LADb[patchI];
        const fvPatchScalarField& vegQrpI = vegQrb[patchI];

    		if ((isA<wallFvPatch>(vegMesh.boundary()[patchI])) && (pp.size() > 0))
        {
            //viewFactorsPatches[count] = LADpI.patch().index();
            vegViewFactorsPatches[count] = vegQrpI.patch().index();
            count ++;
            vegHowManyCoarseFacesPerPatch[countAll] = vegCoarsePatches[patchI].size();
            vegNCoarseFaces += vegCoarsePatches[patchI].size();
            vegNCoarseFacesAll += vegCoarsePatches[patchI].size();
        }
        else if ((isA<fixedValueFvPatchScalarField>(vegQrpI)) && (pp.size() > 0))
		    {
            vegHowManyCoarseFacesPerPatch[countAll] = vegCoarsePatches[patchI].size();
            vegNCoarseFacesAll += vegCoarsePatches[patchI].size();
        }
        else
        {
            vegHowManyCoarseFacesPerPatch[countAll] = 0;
        }
        countAll ++;
    }

    vegViewFactorsPatches.resize(count--);

    DynamicList<point> vegLocalCoarseCf(vegNCoarseFaces);

    forAll (vegViewFactorsPatches, i)
    {
        const label patchID = vegViewFactorsPatches[i];

        const polyPatch& pp = vegPatches[patchID];
        const labelList& agglom = vegFinalAgglom[patchID];

        label nAgglom = max(agglom)+1;
        labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
        const labelList& vegCoarsePatchFace = vegCoarseMesh.patchFaceMap()[patchID];

        const pointField& vegCoarseCf = vegCoarseMesh.Cf().boundaryField()[patchID];
        //const pointField& coarseSf = coarseMesh.Sf().boundaryField()[patchID];

        labelHashSet includePatches;
        includePatches.insert(patchID);

        forAll(vegCoarseCf, faceI)
        {
            point cf = vegCoarseCf[faceI];

            const label vegCoarseFaceI = vegCoarsePatchFace[faceI];
            const labelList& fineFaces = coarseToFine[vegCoarseFaceI];
            // const label agglomI =
            //     agglom[fineFaces[0]] + coarsePatches[patchID].start();

            // Construct single face
            uindirectPrimitivePatch upp
            (
                UIndirectList<face>(pp, fineFaces),
                pp.points()
            );

            List<point> availablePoints
            (
                upp.faceCentres().size()
              + upp.localPoints().size()
            );

            SubList<point>
            (
                availablePoints,
                upp.faceCentres().size()
            ) = upp.faceCentres();

            SubList<point>
            (
                availablePoints,
                upp.localPoints().size(),
                upp.faceCentres().size()
            ) = upp.localPoints();

            point cfo = cf;
            scalar dist = GREAT;
            forAll(availablePoints, iPoint)
            {
                point cfFine = availablePoints[iPoint];
                if(mag(cfFine-cfo) < dist)
                {
                    dist = mag(cfFine-cfo);
                    cf = cfFine;
                }
            }
            vegLocalCoarseCf.append(cf);
        }
    }

    //Info << "localCoarseCf " <<  localCoarseCf << endl;
