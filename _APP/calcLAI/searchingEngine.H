
const polyBoundaryMesh& patches = mesh.boundaryMesh();

// Determine mesh bounding boxes:
List<treeBoundBox> meshBb
(
    1,
    treeBoundBox
    (
        boundBox(mesh.points(), false)
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
     && !isA<cyclicTransform>(pp)
     && !isA<emptyPolyPatch>(pp)
    )
    {
        includePatches.insert(patchi);
    }
}

const triSurface localSurface = Foam::triSurfaceTools::triangulate
(
    patches,
    includePatches
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
//surfacesMesh.searchableSurface::write();

