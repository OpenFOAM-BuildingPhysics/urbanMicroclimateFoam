
const polyBoundaryMesh& patches = mesh.boundaryMesh();
const polyBoundaryMesh& coarsePatches = coarseMesh.boundaryMesh();

label nCoarseFaces = 0;      //total number of coarse faces
label nCoarseFacesAll = 0;   //Also includes non-wall faces with greyDiffusive boundary
label nFineFacesTotal = 0;        //total number of fine faces including non-fixedValueFvPatchScalarField patches following advise from bug report

//const volScalarField::Boundary& qrb = qr.boundaryField();
const volScalarField::Boundary& LADb = LAD.boundaryField();

forAll(LADb, patchI)
//forAll(qrb, patchI)
{
    const polyPatch& pp = patches[patchI];
    //const fvPatchScalarField& LADI = LADb[patchI];

    if ((isA<wallFvPatch>(mesh.boundary()[patchI])) && (pp.size() > 0))
    {
      nCoarseFaces += coarsePatches[patchI].size();
      nCoarseFacesAll += coarsePatches[patchI].size();
      nFineFacesTotal += patches[patchI].size();
    }
    else if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp) && !isA<emptyPolyPatch>(pp))
    //else if ((isA<fixedValueFvPatchScalarField>(qrpI)) && (pp.size() > 0))
    {
        nCoarseFacesAll += coarsePatches[patchI].size();
        nFineFacesTotal += patches[patchI].size();
    }
/*    else
    {
        nFineFacesTotal += patches[patchI].size(); Info << "AAA" << endl;
    }*/

}

globalIndex globalNumbering(nCoarseFaces);

Random rndGen(653213);

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
    distributedTriSurfaceMesh::distributionTypeNames_
    [
        distributedTriSurfaceMesh::FROZEN
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

distributedTriSurfaceMesh surfacesMesh
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

surfacesMesh.setField(triSurfaceToAgglom);

