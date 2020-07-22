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
#include "cyclicAMIPolyPatch.H"
#include "triSurfaceTools.H"
#include "mapDistribute.H"

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

#include "unitConversion.H"

#include "regionProperties.H"

#include "interpolationTable.H"

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

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const label patchi = iter.key();
        const polyPatch& patch = bMesh[patchi];
        const pointField& points = patch.points();

        label nTriTotal = 0;

        forAll(patch, patchFacei)
        {
            const face& f = patch[patchFacei];

            faceList triFaces(f.nTriangles(points));

            label nTri = 0;

            f.triangles(points, nTri, triFaces);

            forAll(triFaces, triFacei)
            {
                const face& f = triFaces[triFacei];

                triangles.append(labelledTri(f[0], f[1], f[2], newPatchi));

                nTriTotal++;

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
    #include "createNamedMesh.H"
    
    // Read view factor dictionary
    IOdictionary viewFactorDict
    (
       IOobject
       (
            "viewFactorsDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
       )
    );

    vector skyPos = viewFactorDict.lookup("skyPosVector");

    // Read sunPosVector list
    interpolationTable<vector> sunPosVector
    (
        runTime.time().rootPath()
        /runTime.time().globalCaseName()
        /runTime.time().constant()
        /"sunPosVector"
    );  
    interpolationTable<scalar> IDN // direct solar radiation intensity flux
    (
        runTime.time().rootPath()
        /runTime.time().globalCaseName()
        /runTime.time().constant()
        /"IDN"
    );    
    interpolationTable<scalar> Idif // diffuse solar radiation intensity flux
    (
        runTime.time().rootPath()
        /runTime.time().globalCaseName()
        /runTime.time().constant()
        /"Idif"
    );

    const label debug = viewFactorDict.lookupOrDefault<label>("debug", 0);

    volScalarField qr
    (
        IOobject
        (
            "qr",
            runTime.timeName(),
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
            runTime.timeName(),
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
        if ((isA<wallFvPatch>(mesh.boundary()[patchi])) && (pp.size() > 0))
        {
            viewFactorsPatches[count] = qrpI.patch().index();
            nCoarseFaces += coarsePatches[patchi].size();
            nCoarseFacesAll += coarsePatches[patchi].size();
            nFineFaces += patches[patchi].size();
            count ++;
            
            howManyCoarseFacesPerPatch[countAll] = coarsePatches[patchi].size();

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

            nFineFacesTotal += patches[patchi].size();        
        }
        countAll ++;
    }
    viewFactorsPatches.resize(count);
    Info << "howManyCoarseFacesPerPatch: " << howManyCoarseFacesPerPatch << endl;

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

    DynamicList<point> localCoarseCf(nCoarseFaces);
    DynamicList<point> localCoarseSf(nCoarseFaces);

    DynamicList<label> localAgg(nCoarseFaces);

    forAll (viewFactorsPatches, i)
    {
        const label patchID = viewFactorsPatches[i];

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
    DynamicList<label> rayStartFace(nCoarseFaces + 0.01*nCoarseFaces);

    DynamicList<label> rayEndFace(rayStartFace.size());  


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

    point min_(point::zero);
    point max_(point::zero);
    for (label i = 0; i < minList.size(); i++)
    {
        min_ = ::Foam::min(min_, minList[i]);
        max_ = ::Foam::max(max_, maxList[i]);
    }    

    // Find the Solar Ray Start Points within domain
    // ######################################   
    List<point> solarStart(localCoarseCf);
    List<point> solarEnd(solarStart.size());   

    // Number of visible faces from local index
    labelListList nVisibleFaceFacesList(sunPosVector.size());    
    labelListList visibleFaceFaces(nCoarseFaces); 

    forAll(sunPosVector, vectorId)
    {   
        labelList nVisibleFaceFaces(nCoarseFaces, 0);

        vector sunPos = sunPosVector[vectorId].second();

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

        // Return rayStartFace in local index and rayEndFace in global index
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #include "shootRays.H"


        forAll(rayStartFace, i)
        {
            nVisibleFaceFaces[rayStartFace[i]]++;
        }

        label nViewFactors = 0;
        forAll(nVisibleFaceFaces, faceI)
        {
            visibleFaceFaces[faceI].setSize(nVisibleFaceFaces[faceI]);
            nViewFactors += nVisibleFaceFaces[faceI];
        }

        //Info << "nVisibleFaceFaces: " << nVisibleFaceFaces << endl;
        nVisibleFaceFacesList[vectorId] = nVisibleFaceFaces;

        rayStartFace.clear();
        rayEndFace.clear();

    }   

    Info << "nVisibleFaceFacesList: " << nVisibleFaceFacesList << endl;
    

    // Fill local view factor matrix
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListIOList sunVisibleOrNot
    (
        IOobject
        (
            "sunVisibleOrNot",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        sunPosVector.size()
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
        sunPosVector.size()
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
        sunPosVector.size()
    );    
    
    labelList dummy(nCoarseFacesAll, -1);
    scalarList dummy2(nCoarseFacesAll, 0.0);
    forAll(sunVisibleOrNot, vectorId)
    {
        sunVisibleOrNot[vectorId] = dummy;
        sunViewCoeff[vectorId] = dummy2;
        skyViewCoeff[vectorId] = dummy2;
    }

    scalar cosPhi = 0;
    scalar radAngleBetween = 0;
    scalar degAngleBetween = 0;
    
    label faceNo = 0;
    label i = 0;
    label j = 0;
    label k = 0;

    regionProperties rp(runTime); 
    const wordList vegNames(rp["vegetation"]); 
    scalarListIOList kcLAIboundaryList
    (
        IOobject
        (
            "kcLAIboundary",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );   
    if (vegNames.size()>0 && !kcLAIboundaryList.typeHeaderOk<IOList<scalarList>>())
    {
        FatalErrorInFunction
            << "File kcLAIboundary not found! Did you not run calcLAI before?"
            << exit(FatalError);
    }        

    forAll(sunPosVector, vectorId)
    {    
        vector sunPos = sunPosVector[vectorId].second();

        forAll(viewFactorsPatches, patchID)
        {
            while (i < viewFactorsPatches[patchID])
            {
                while (j < howManyCoarseFacesPerPatch[i])
                {
                    sunVisibleOrNot[vectorId][k] = 0;
                    k++;
                    j++;
                }
                j = 0;
                i++;
            } 
            
            while (j < howManyCoarseFacesPerPatch[i])
            {
                sunVisibleOrNot[vectorId][k] = nVisibleFaceFacesList[vectorId][faceNo];

                cosPhi = (localCoarseSf[faceNo] & sunPos)/(mag(localCoarseSf[faceNo])*mag(sunPos) + SMALL);                
                sunViewCoeff[vectorId][k] = nVisibleFaceFacesList[vectorId][faceNo]*mag(cosPhi) * IDN[vectorId].second();
                if (vegNames.size()>0)
                {
                    if (kcLAIboundaryList[vectorId][k]-0>SMALL && cosPhi < 0) //if LAIboundary value is nonzero and if the surface is looking towards the sun, update sunViewCoeff
                    {
                        sunViewCoeff[vectorId][k] = mag(cosPhi) * IDN[vectorId].second() * Foam::exp(-kcLAIboundaryList[vectorId][k]); // beer-lambert law
                    }
                }
                
                cosPhi = (localCoarseSf[faceNo] & skyPos)/(mag(localCoarseSf[faceNo])*mag(skyPos) + SMALL);
                radAngleBetween = Foam::acos( min(max(cosPhi, -1), 1) );
                degAngleBetween = radToDeg(radAngleBetween);
                if (degAngleBetween > 90 && degAngleBetween <= 180){degAngleBetween=90 - (degAngleBetween-90);}
                skyViewCoeff[vectorId][k] = (1-0.5*(degAngleBetween/90)) * Idif[vectorId].second();            
                
                k++;
                j++;
                faceNo++;
            }
        }
        i = 0;
        j = 0;
        k = 0;
        faceNo = 0;
    }
    Info << "sunVisibleOrNot: " << sunVisibleOrNot << endl;    

    Info << "localCoarseCf: " << localCoarseCf << endl;    
    Info << "localCoarseSf: " << localCoarseSf << endl;
    
    Info << "sunViewCoeff: " << sunViewCoeff << endl;    
    Info << "skyViewCoeff: " << skyViewCoeff << endl;    

    sunVisibleOrNot.write();
    sunViewCoeff.write();    
    skyViewCoeff.write();    

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
