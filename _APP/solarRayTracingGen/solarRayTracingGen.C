/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    solarRayTracingGen

Description
    Aytac Kubilay, 2015, Empa
    Based on viewFactorsGen

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "distributedTriSurfaceMeshBugFix.H"
#include "cyclicTransform.H"
#include "triSurfaceTools.H"
#include "distributionMap.H"

#include "OFstream.H"
#include "meshTools.H"
#include "plane.H"
#include "uindirectPrimitivePatch.H"
#include "DynamicField.H"
#include "IFstream.H"
#include "unitConversion.H"

#include "mathematicalConstants.H"
#include "scalarMatrices.H"
#include "CompactListList.H"
#include "labelIOList.H"
#include "labelListIOList.H"
#include "scalarListIOList.H"
#include "scalarIOList.H"
#include "vectorIOList.H"

#include "singleCellFvMesh.H"
#include "IOdictionary.H"
#include "fixedValueFvPatchFields.H"
#include "wallFvPatch.H"
#include "mappedInternalFvPatch.H"

#include "unitConversion.H"
#include "dimensionSets.H"
#include "polygonTriangulate.H"
#include "randomGenerator.H"

#include "Table.H"

using namespace Foam;

triSurface triangulate
(
    const polyBoundaryMesh& bMesh,
    const labelHashSet& includePatches,
    const labelListIOList& finalAgglom,
    labelList& triSurfaceToAgglom,
    const globalIndex& globalNumbering,
    const polyBoundaryMesh& coarsePatches
)
{
    const polyMesh& mesh = bMesh.mesh();

    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles
    (
        mesh.nFaces() - mesh.nInternalFaces()
    );

    label newPatchi = 0;
    label localTriFacei = 0;

    polygonTriangulate triEngine;

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const label patchi = iter.key();
        const polyPatch& patch = bMesh[patchi];
        const pointField& points = patch.points();

        forAll(patch, patchFacei)
        {
            const face& f = patch[patchFacei];

            triEngine.triangulate(UIndirectList<point>(points, f));

            forAll(triEngine.triPoints(), triFacei)
            {
                triangles.append
                (
                    labelledTri
                    (
                        triEngine.triPoints(triFacei, f),
                        newPatchi
                    )
                );

                triSurfaceToAgglom[localTriFacei++] = globalNumbering.toGlobal
                (
                    Pstream::myProcNo(),
                    finalAgglom[patchi][patchFacei]
                  + coarsePatches[patchi].start()
                );
            }
        }

        newPatchi++;
    }

    triSurfaceToAgglom.resize(localTriFacei);

    triangles.shrink();

    // Create globally numbered tri surface
    triSurface rawSurface(triangles, mesh.points());

    // Create locally numbered tri surface
    triSurface surface
    (
        rawSurface.localFaces(),
        rawSurface.localPoints()
    );

    // Add patch names to surface
    surface.patches().setSize(newPatchi);

    newPatchi = 0;

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const label patchi = iter.key();
        const polyPatch& patch = bMesh[patchi];

        surface.patches()[newPatchi].index() = patchi;
        surface.patches()[newPatchi].name() = patch.name();
        surface.patches()[newPatchi].geometricType() = patch.type();

        newPatchi++;
    }

    return surface;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createRegionMeshNoChangers.H"
    
    // Read view factor dictionary
    IOdictionary viewFactorDict
    (
       IOobject
       (
            "viewFactorsDict",
            runTime.system(),  // read from system/<region>/ for consistency with viewFactorsGen
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
       )
    );

    vector skyPos = viewFactorDict.lookup("skyPosVector");

    // Read sunPosVector list
    dictionary sunPosVectorIO;
    sunPosVectorIO.add(
        "file",
        fileName
        (
            mesh.time().constant()
            /"sunPosVector"
        )
    );
    Function1s::Table<vector> sunPosVector
    (
        "sunPosVector",
        dimTime,
        dimless,
        sunPosVectorIO
    );
    // Read solar radiation intensity flux
    dictionary IDNIO;
    IDNIO.add(
        "file",
        fileName
        (
            mesh.time().constant()
            /"IDN"
        )
    );
    Function1s::Table<scalar> IDN
    (
        "IDN",
        dimTime,
        dimless,
        IDNIO
    );
    // Read diffuse solar radiation intensity flux
    dictionary IdifIO;
    IdifIO.add(
        "file",
        fileName
        (
            mesh.time().constant()
            /"Idif"
        )
    );
    Function1s::Table<scalar> Idif
    (
        "Idif",
        dimTime,
        dimless,
        IdifIO
    );
    
    vectorField sunPosVector_y = sunPosVector.y();
    scalarField IDN_y = IDN.y();
    scalarField Idif_y = Idif.y();

    const label debug = viewFactorDict.lookupOrDefault<label>("debug", 0);

    volScalarField qr
    (
        IOobject
        (
            "qr",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Read agglomeration map
    labelListIOList finalAgglom
    (
        IOobject
        (
            "finalAgglom",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // Create the coarse mesh  using agglomeration
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Info << "\nCreating single cell mesh..." << endl;
    }

    singleCellFvMesh coarseMesh
    (
        IOobject
        (
            "coarse:" + mesh.name(),
            runTime.name(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        finalAgglom
    );

    if (debug)
    {
        Pout << "\nCreated single cell mesh..." << endl;
    }

    // Calculate total number of fine and coarse faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nCoarseFaces = 0;      //total number of coarse faces
    label nCoarseFacesAll = 0;   //Also includes non-wall faces with greyDiffusive boundary
    label nFineFaces = 0;        //total number of fine faces
    label nFineFacesTotal = 0;        //total number of fine faces including non-fixedValueFvPatchScalarField patches following advise from bug report

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const polyBoundaryMesh& coarsePatches = coarseMesh.boundaryMesh();

    labelList viewFactorsPatches(patches.size());
    labelList howManyCoarseFacesPerPatch(patches.size());
    labelList howManyFineFacesPerPatch(patches.size());
    DynamicList<label> sunskyMap_(nCoarseFaces);
    const volScalarField::Boundary& qrb = qr.boundaryField();

    label count = 0;
    label countAll = 0;
    label countForMapping = 0;
    forAll(qrb, patchi)
    {
        const polyPatch& pp = patches[patchi];
        const fvPatchScalarField& qrpI = qrb[patchi];

        //if ((isA<fixedValueFvPatchScalarField>(QrpI)) && (pp.size() > 0))
        // v12: also include mappedInternal patches (e.g. air_to_vegetation)
        // which were mappedWall in v8 but changed to mappedInternal in v12
        if
        (
            (
                isA<wallFvPatch>(mesh.boundary()[patchi])
             || isA<mappedInternalFvPatch>(mesh.boundary()[patchi])
            )
         && (pp.size() > 0)
        )
        {
            viewFactorsPatches[count] = qrpI.patch().index();
            nCoarseFaces += coarsePatches[patchi].size();
            nCoarseFacesAll += coarsePatches[patchi].size();
            nFineFaces += patches[patchi].size();
            count ++;
            
            howManyCoarseFacesPerPatch[countAll] = coarsePatches[patchi].size();
            howManyFineFacesPerPatch[countAll] = patches[patchi].size();

            label i = 0;
            for (; i < howManyCoarseFacesPerPatch[countAll]; i++)
            {
                sunskyMap_.append(countForMapping);
                countForMapping ++;
            }
            nFineFacesTotal += patches[patchi].size();             
        }
        else if ((isA<fixedValueFvPatchScalarField>(qrpI)) && (pp.size() > 0))
        {
            nCoarseFacesAll += coarsePatches[patchi].size();
            
            howManyCoarseFacesPerPatch[countAll] = coarsePatches[patchi].size();  
            howManyFineFacesPerPatch[countAll] = patches[patchi].size();

            label i = 0;
            for (; i < howManyCoarseFacesPerPatch[countAll]; i++)
            {
                sunskyMap_.append(countForMapping);
                countForMapping ++;
            }        

            nFineFacesTotal += patches[patchi].size();        
        }
        else 
        {
            howManyCoarseFacesPerPatch[countAll] = 0;  
            howManyFineFacesPerPatch[countAll] = 0;

            nFineFacesTotal += patches[patchi].size();        
        }
        countAll ++;
    }
    viewFactorsPatches.resize(count);
    Info << "howManyCoarseFacesPerPatch: " << howManyCoarseFacesPerPatch << endl;
    Info << "howManyFineFacesPerPatch: " << howManyFineFacesPerPatch << endl;

    List<labelField> sunskyMap__(Pstream::nProcs());
    sunskyMap__[Pstream::myProcNo()] = sunskyMap_;  
    Pstream::gatherList(sunskyMap__);
    Pstream::scatterList(sunskyMap__);

    label sunskyMapCounter = 0;
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    { 
        if (procI > 0)
        {
            sunskyMap__[procI] += sunskyMapCounter;
        }   
        sunskyMapCounter += sunskyMap__[procI].size();           
    }

    sunskyMap_ = sunskyMap__[Pstream::myProcNo()];

    labelIOList sunskyMap
    (
        IOobject
        (
            "sunskyMap",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        sunskyMap_.size()
    );
    sunskyMap = sunskyMap_;
    sunskyMap.write();
    
    // total number of coarse faces
    label totalNCoarseFaces = nCoarseFaces;

    reduce(totalNCoarseFaces, sumOp<label>());

    if (Pstream::master())
    {
        Info << "\nTotal number of coarse faces: "<< totalNCoarseFaces << endl;
    }

    if (Pstream::master() && debug)
    {
        Pout << "\nView factor patches included in the calculation : "
             << viewFactorsPatches << endl;
    }

    // Collect local Cf and Sf on coarse mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<point> localCoarseCf(nCoarseFaces); DynamicList<point> localFINECf(nFineFaces);
    DynamicList<point> localCoarseSf(nCoarseFaces); DynamicList<point> localFINESf(nFineFaces);

    DynamicList<label> localAgg(nCoarseFaces);

    forAll (viewFactorsPatches, i)
    {
        const label patchID = viewFactorsPatches[i];
 
        localFINECf.append(mesh.Cf().boundaryField()[patchID]);
        localFINESf.append(mesh.Sf().boundaryField()[patchID]);

        const polyPatch& pp = patches[patchID];
        const labelList& agglom = finalAgglom[patchID];
        
        label nAgglom = max(agglom)+1;
        labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
        const labelList& coarsePatchFace = coarseMesh.patchFaceMap()[patchID];

        const pointField& coarseCf = coarseMesh.Cf().boundaryField()[patchID];
        const pointField& coarseSf = coarseMesh.Sf().boundaryField()[patchID];

        labelHashSet includePatches;
        includePatches.insert(patchID);

        forAll(coarseCf, faceI)
        {
            point cf = coarseCf[faceI];

            const label coarseFaceI = coarsePatchFace[faceI];
            const labelList& fineFaces = coarseToFine[coarseFaceI];
            const label agglomI =
                agglom[fineFaces[0]] + coarsePatches[patchID].start();

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

            point sf = coarseSf[faceI];
            localCoarseCf.append(cf);
            localCoarseSf.append(sf);
            localAgg.append(agglomI);            
        }
    }

    globalIndex globalNumbering(nCoarseFaces);       
    
    // Set up searching engine for obstacles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #include "searchingEngine.H"    


    // Determine rays between coarse face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DynamicList<label> rayStartFace(nCoarseFaces + 0.001*nCoarseFaces); DynamicList<label> rayStartFaceFINE(nFineFaces + 0.001*nFineFaces);

    DynamicList<label> rayEndFace(rayStartFace.size()); DynamicList<label> rayEndFaceFINE(rayStartFaceFINE.size()); 


    // Find the bounding box of the domain
    // ######################################
    List<point> minList(Pstream::nProcs());
    List<point> maxList(Pstream::nProcs());
    minList[Pstream::myProcNo()] = Foam::min(mesh.points());
    maxList[Pstream::myProcNo()] = Foam::max(mesh.points());
    Pstream::gatherList(minList);
    Pstream::gatherList(maxList);
    Pstream::scatterList(minList);
    Pstream::scatterList(maxList);

    point min_(point::max);
    point max_(point::min);
    for (label i = 0; i < minList.size(); i++)
    {
        min_ = ::Foam::min(min_, minList[i]);
        max_ = ::Foam::max(max_, maxList[i]);
    }    

    // Find the Solar Ray Start Points within domain
    // ######################################   
    List<point> solarStart(localCoarseCf); List<point> solarStartFINE(localFINECf);
    List<point> solarEnd(solarStart.size()); List<point> solarEndFINE(solarStartFINE.size()); 

    // Number of visible faces from local index
    labelListList nVisibleFaceFacesList(sunPosVector_y.size()); labelListList nVisibleFaceFacesListFINE(sunPosVector_y.size()); 
    labelListList visibleFaceFaces(nCoarseFaces); labelListList visibleFaceFacesFINE(nFineFaces); 

    forAll(sunPosVector_y, vectorId)
    {   
        labelList nVisibleFaceFaces(nCoarseFaces, 0); labelList nVisibleFaceFacesFINE(nFineFaces, 0);

        vector sunPos = sunPosVector_y[vectorId];

        //List<pointIndexHit> hitInfo(1);
        forAll(solarStart, pointI)
        {
            scalar i1 = 0; scalar i2 = 0; scalar i3 = 0;

            if (sunPos.x() > 0.0)
            {
                i1 = (max_.x() - solarStart[pointI].x())/sunPos.x();
            } 
            else if (sunPos.x() < 0.0)
            {
                i1 = (min_.x() - solarStart[pointI].x())/sunPos.x();
            } 
            else {i1 = VGREAT;}

            if (sunPos.y() > 0.0)
            {
                i2 = (max_.y() - solarStart[pointI].y())/sunPos.y();
            } 
            else if (sunPos.y() < 0.0)
            {
                i2 = (min_.y() - solarStart[pointI].y())/sunPos.y();
            }
            else{i2 = VGREAT;}

            if (sunPos.z() > 0.0)
            {
                i3 = (max_.z() - solarStart[pointI].z())/sunPos.z();
            } 
            else if (sunPos.z() < 0.0)
            {
                i3 = (min_.z() - solarStart[pointI].z())/sunPos.z();
            }
            else{i3 = VGREAT;}

            scalar i = min(i1, min(i2, i3));
            point solarEndPoint = i*point(sunPos.x(),sunPos.y(),sunPos.z())+solarStart[pointI];
            solarEnd[pointI] = solarEndPoint;
        }  

        forAll(solarStartFINE, pointI)
        {
            scalar i1 = 0; scalar i2 = 0; scalar i3 = 0;

            if (sunPos.x() > 0.0)
            {
                i1 = (max_.x() - solarStartFINE[pointI].x())/sunPos.x();
            } 
            else if (sunPos.x() < 0.0)
            {
                i1 = (min_.x() - solarStartFINE[pointI].x())/sunPos.x();
            } 
            else {i1 = VGREAT;}

            if (sunPos.y() > 0.0)
            {
                i2 = (max_.y() - solarStartFINE[pointI].y())/sunPos.y();
            } 
            else if (sunPos.y() < 0.0)
            {
                i2 = (min_.y() - solarStartFINE[pointI].y())/sunPos.y();
            }
            else{i2 = VGREAT;}

            if (sunPos.z() > 0.0)
            {
                i3 = (max_.z() - solarStartFINE[pointI].z())/sunPos.z();
            } 
            else if (sunPos.z() < 0.0)
            {
                i3 = (min_.z() - solarStartFINE[pointI].z())/sunPos.z();
            }
            else{i3 = VGREAT;}

            scalar i = min(i1, min(i2, i3));
            point solarEndPointFINE = i*point(sunPos.x(),sunPos.y(),sunPos.z())+solarStartFINE[pointI];
            solarEndFINE[pointI] = solarEndPointFINE;
        }

        //Info << "solarStart: " << solarStart << endl;
        //Info << "solarEnd: " << solarEnd << endl;

        // Collect Cf and Sf on coarse mesh
        // #############################################

        List<pointField> remoteCoarseCf_(Pstream::nProcs());
        remoteCoarseCf_[Pstream::myProcNo()] = solarEnd;
        
        List<pointField> localCoarseCf_(Pstream::nProcs());
        localCoarseCf_[Pstream::myProcNo()] = solarStart;

        List<pointField> localCoarseSf_(Pstream::nProcs());
        localCoarseSf_[Pstream::myProcNo()] = localCoarseSf;         

        List<pointField> remoteFINECf_(Pstream::nProcs());
        remoteFINECf_[Pstream::myProcNo()] = solarEndFINE;
        
        List<pointField> localFINECf_(Pstream::nProcs());
        localFINECf_[Pstream::myProcNo()] = solarStartFINE;

        List<pointField> localFINESf_(Pstream::nProcs());
        localFINESf_[Pstream::myProcNo()] = localFINESf;

        // Collect remote Cf and Sf on fine mesh
        // #############################################

        /*List<pointField> remoteFineCf(Pstream::nProcs());
        List<pointField> remoteFineSf(Pstream::nProcs());

        remoteCoarseCf[Pstream::myProcNo()] = solarEnd;
        remoteCoarseSf[Pstream::myProcNo()] = localCoarseSf;*/
        
        // Distribute local coarse Cf and Sf for shooting rays
        // #############################################

        Pstream::gatherList(remoteCoarseCf_);
        Pstream::scatterList(remoteCoarseCf_);
        Pstream::gatherList(localCoarseCf_);
        Pstream::scatterList(localCoarseCf_);
        Pstream::gatherList(localCoarseSf_);
        Pstream::scatterList(localCoarseSf_);   

        Pstream::gatherList(remoteFINECf_);
        Pstream::scatterList(remoteFINECf_);
        Pstream::gatherList(localFINECf_);
        Pstream::scatterList(localFINECf_);
        Pstream::gatherList(localFINESf_);
        Pstream::scatterList(localFINESf_);

        // Return rayStartFace in local index and rayEndFace in global index
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //#include "shootRays.H"
        #include "shootRaysFINE.H"

        forAll(rayStartFaceFINE, i)
        {
            nVisibleFaceFacesFINE[rayStartFaceFINE[i]]++;
        }

        /*label nViewFactors = 0;
        forAll(nVisibleFaceFaces, faceI)
        {
            visibleFaceFaces[faceI].setSize(nVisibleFaceFaces[faceI]);
            nViewFactors += nVisibleFaceFaces[faceI];
        }*/

        //Info << "nVisibleFaceFaces: " << nVisibleFaceFaces << endl;
        nVisibleFaceFacesList[vectorId] = nVisibleFaceFaces;
        nVisibleFaceFacesListFINE[vectorId] = nVisibleFaceFacesFINE;

        rayStartFace.clear(); rayStartFaceFINE.clear();
        rayEndFace.clear(); rayEndFaceFINE.clear();

    }   

    Info << "nVisibleFaceFacesList: " << nVisibleFaceFacesList << endl;
    Info << "nVisibleFaceFacesListFINE: " << nVisibleFaceFacesListFINE << endl; 

    // Fill local view factor matrix
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarListIOList solarLoadFineFaces
    (
        IOobject
        (
            "solarLoadFineFaces",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        sunPosVector_y.size()
    );   
    scalarListIOList sunViewCoeff
    (
        IOobject
        (
            "sunViewCoeff",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        sunPosVector_y.size()
    );
    scalarListIOList skyViewCoeff
    (
        IOobject
        (
            "skyViewCoeff",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        sunPosVector_y.size()
    );    
    
    scalarList init(nCoarseFacesAll, 0.0);
    scalarList initFine(nFineFaces, 0.0);
    forAll(sunViewCoeff, vectorId)
    {
        sunViewCoeff[vectorId] = init;
        skyViewCoeff[vectorId] = init;
        solarLoadFineFaces[vectorId] = initFine;
    }

    scalar cosPhi = 0;
    scalar radAngleBetween = 0;
    scalar degAngleBetween = 0;
    
    label faceNo = 0;
    label fineFaceNo = 0;     
    label patchIDall = 0;
    label j = 0;
    label faceNoAll = 0;

    // Read vegetation region names from controlDict (v12 style)
    wordList vegNames;
    if (runTime.controlDict().found("regions"))
    {
        const dictionary& regionsDict = runTime.controlDict().subDict("regions");
        if (regionsDict.found("vegetation"))
        {
            vegNames = wordList(regionsDict.lookup("vegetation"));
        }
    } 
    scalarListIOList kcLAIboundaryList
    (
        IOobject
        (
            "kcLAIboundary",
            mesh.facesInstance(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );  
    if (vegNames.size()>0 && !kcLAIboundaryList.headerOk())
    {
        FatalErrorInFunction
            << "File kcLAIboundary not found! Did you not run calcLAI before?"
            << exit(FatalError);
    }      

    if (vegNames.size()>0)
    {
        if (kcLAIboundaryList.size() != sunPosVector_y.size())
        {
            FatalErrorInFunction
                << "kcLAIboundary holds " << kcLAIboundaryList.size()
                << " sun directions, expected " << sunPosVector_y.size()
                << ". Re-run calcLAI with the same sunPosVector."
                << exit(FatalError);
        }
        forAll(kcLAIboundaryList, vectorId)
        {
            if (kcLAIboundaryList[vectorId].size() != nCoarseFacesAll)
            {
                FatalErrorInFunction
                    << "kcLAIboundary[" << vectorId << "] holds "
                    << kcLAIboundaryList[vectorId].size()
                    << " coarse faces, expected " << nCoarseFacesAll
                    << ". calcLAI patch selection or decomposition does not "
                    << "match this utility."
                    << exit(FatalError);
            }
        }
    }

    forAll(sunPosVector_y, vectorId)
    {    
        vector sunPos = sunPosVector_y[vectorId];

        forAll(viewFactorsPatches, patchID)
        {
            while (patchIDall < viewFactorsPatches[patchID])
            {
                while (j < howManyCoarseFacesPerPatch[patchIDall])
                {
                    faceNoAll++;
                    j++;
                }
                j = 0;
                patchIDall++;
            } 
            
            while (j < howManyCoarseFacesPerPatch[patchIDall])
            {

                /////////////////////////////////////////////////////////////////////////
                const labelList& agglom = finalAgglom[patchIDall];
                label nAgglom = max(agglom)+1;
                labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
                const labelList& coarsePatchFace = coarseMesh.patchFaceMap()[patchIDall];
                const label coarseFaceI = coarsePatchFace[j];
                const labelList& fineFaces = coarseToFine[coarseFaceI];
                
                forAll(fineFaces, k)
                {
                     label faceI = fineFaces[k];
                     cosPhi = (localFINESf[fineFaceNo+faceI] & sunPos)/(mag(localFINESf[fineFaceNo+faceI])*mag(sunPos) + SMALL);
                     solarLoadFineFaces[vectorId][fineFaceNo+faceI] = nVisibleFaceFacesListFINE[vectorId][fineFaceNo+faceI]*mag(cosPhi) * IDN_y[vectorId];                                                  
                }
                                
                // Denominator as  the SUM of fine areas (Sigma|Sf_fine|)
                scalar coarseFineArea = 0;
                forAll(fineFaces,fineFaceI)
                {
                    coarseFineArea +=
                        mesh.magSf().boundaryField()[patchIDall][fineFaces[fineFaceI]];
                }
                scalar nVisibleFaceFacesListFINE_avg = 0;
                forAll(fineFaces,fineFaceI)
                {
                    nVisibleFaceFacesListFINE_avg += (nVisibleFaceFacesListFINE[vectorId][fineFaceNo+fineFaces[fineFaceI]])
                                                    * (mesh.magSf().boundaryField()[patchIDall][fineFaces[fineFaceI]])
                                                    / (coarseFineArea + SMALL);
                }                                  
                
                /////////////////////////////////////////////////////////////////////////

                cosPhi = (localCoarseSf[faceNo] & sunPos)/(mag(localCoarseSf[faceNo])*mag(sunPos) + SMALL);                
                //sunViewCoeff[vectorId][faceNoAll] = nVisibleFaceFacesList[vectorId][faceNo]*mag(cosPhi) * IDN[vectorId].second();
                sunViewCoeff[vectorId][faceNoAll] = nVisibleFaceFacesListFINE_avg*mag(cosPhi) * IDN_y[vectorId];

                if (vegNames.size()>0)
                {
                    if (kcLAIboundaryList[vectorId][faceNoAll]-0>SMALL && cosPhi < 0) 
                    //if LAIboundary value is positive and if the surface is looking towards the sun, update sunViewCoeff
                    //nVisibleFaceFacesListFINE_avg indicates the ratio of coarseFace that see the sun
                    {
                        sunViewCoeff[vectorId][faceNoAll] += (1 - nVisibleFaceFacesListFINE_avg) * mag(cosPhi) * IDN_y[vectorId] * Foam::exp(-kcLAIboundaryList[vectorId][faceNoAll]); // beer-lambert law
                        forAll(fineFaces, k)
                        {
                             label faceI = fineFaces[k];
                             cosPhi = (localFINESf[fineFaceNo+faceI] & sunPos)/(mag(localFINESf[fineFaceNo+faceI])*mag(sunPos) + SMALL);
                             if(!nVisibleFaceFacesListFINE[vectorId][fineFaceNo+faceI])
                             {
                                 solarLoadFineFaces[vectorId][fineFaceNo+faceI] = mag(cosPhi) * IDN_y[vectorId] * Foam::exp(-kcLAIboundaryList[vectorId][faceNoAll]); // beer-lambert law;
                             }
                        }
                    }
                }
                
                cosPhi = (localCoarseSf[faceNo] & skyPos)/(mag(localCoarseSf[faceNo])*mag(skyPos) + SMALL);
                radAngleBetween = Foam::acos( min(max(cosPhi, -1), 1) );
                degAngleBetween = radToDeg(radAngleBetween);
                if (degAngleBetween > 90 && degAngleBetween <= 180){degAngleBetween=90 - (degAngleBetween-90);}
                skyViewCoeff[vectorId][faceNoAll] = (1-0.5*(degAngleBetween/90)) * Idif_y[vectorId];               
 
                faceNoAll++;
                j++;
                faceNo++;
            }
            fineFaceNo += howManyFineFacesPerPatch[patchIDall];
        }
        patchIDall = 0;
        j = 0;
        faceNoAll = 0;
        faceNo = 0;
        fineFaceNo = 0;
    }
    
    Info << "sunViewCoeff: " << sunViewCoeff << endl;    
    Info << "skyViewCoeff: " << skyViewCoeff << endl;    

    sunViewCoeff.write();    
    skyViewCoeff.write();
    solarLoadFineFaces.write();  

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
