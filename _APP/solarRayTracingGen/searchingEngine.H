Random rndGen(653213);

// Determine mesh bounding boxes:
List<treeBoundBox> meshBb
(
    1,
    treeBoundBox
    (
        boundBox(coarseMesh.points(), false)
    ).extend(1e-3)
);

// Dummy bounds dictionary
dictionary dict;
dict.add("bounds", meshBb);
dict.add
(
    "distributionType",
    distributedTriSurfaceMeshBugFix::distributionTypeNames_
    [
        distributedTriSurfaceMeshBugFix::FROZEN
    ]
);
dict.add("mergeDistance", SMALL);

labelHashSet includePatches;
forAll(patches, patchi)
{
    const polyPatch& pp = patches[patchi];
    if 
    (
        !pp.coupled() 
     && !isA<cyclicAMIPolyPatch>(pp)
     && !isA<emptyPolyPatch>(pp)
    )
    {
        includePatches.insert(patchi);
    }
}

labelList triSurfaceToAgglom(5*nFineFacesTotal);

const triSurface localSurface = triangulate
(
    patches,
    includePatches,
    finalAgglom,
    triSurfaceToAgglom,
    globalNumbering,
    coarsePatches
);

distributedTriSurfaceMeshBugFix surfacesMesh
(
    IOobject
    (
        "wallSurface.stl",
        runTime.constant(),     // directory
        "triSurface",           // instance
        runTime,                // registry
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    localSurface,
    dict
);

triSurfaceToAgglom.resize(surfacesMesh.size());

//surfacesMesh.searchableSurface::write();

surfacesMesh.setField(triSurfaceToAgglom);
