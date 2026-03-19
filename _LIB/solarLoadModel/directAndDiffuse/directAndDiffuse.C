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

\*---------------------------------------------------------------------------*/

#include "directAndDiffuse.H"
#include "surfaceFields.H"
#include "constants.H"
#include "solarLoadViewFactorFixedValueFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "typeInfo.H"
#include "Time.H"

#include "vectorIOList.H"

#include "Function1.H"
#include "Table.H"
#include "IFstream.H"
#include "OFstream.H"
#include "OSspecific.H"

#include "mappedPatchBase.H"
#include "mappedInternalFvPatch.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace solarLoad
    {
        defineTypeNameAndDebug(directAndDiffuse, 0);
        addToSolarLoadRunTimeSelectionTables(directAndDiffuse);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solarLoad::directAndDiffuse::initialise()
{
    const polyBoundaryMesh& coarsePatches = coarseMesh_.boundaryMesh();
    const volScalarField::Boundary& qsp = qs_.boundaryField();

    label count = 0;
    forAll(qsp, patchI)
    {
        //const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        const fvPatchScalarField& qsPatchI = qsp[patchI];

        if ((isA<fixedValueFvPatchScalarField>(qsPatchI)))
        {
            selectedPatches_[count] = qsPatchI.patch().index();
            nLocalCoarseFaces_ += coarsePatches[patchI].size();
            
            // v12: mappedInternal patches (e.g. air_to_vegetation) need
            // wall-like treatment for solar radiation
            if
            (
                isA<wallFvPatch>(mesh_.boundary()[patchI])
             || isA<mappedInternalFvPatch>(mesh_.boundary()[patchI])
            )
            {
                wallPatchOrNot_[count] = 1;
                nLocalWallCoarseFaces_ += coarsePatches[patchI].size();
                nLocalFineFaces_ += qsPatchI.patch().size();
            }

            count++;
        }
    }
    Info<< "Selected patches:" << selectedPatches_ << endl;
    Info<< "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    Info << "wallPatchOrNot_: " << wallPatchOrNot_ << endl;
    Info << "nLocalWallCoarseFaces_: " << nLocalWallCoarseFaces_ << endl;
    
    selectedPatches_.resize(count);
    wallPatchOrNot_.resize(count);

    Info<< "Selected patches:" << selectedPatches_ << endl;
    Info<< "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    Info << "wallPatchOrNot_: " << wallPatchOrNot_ << endl;
    Info << "nLocalWallCoarseFaces_: " << nLocalWallCoarseFaces_ << endl;

    if (debug)
    {
        Pout<< "Selected patches:" << selectedPatches_ << endl;
        Pout<< "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    }

    totalNCoarseFaces_ = nLocalCoarseFaces_;
    reduce(totalNCoarseFaces_, sumOp<label>());
    totalNFineFaces_ = nLocalFineFaces_;
    reduce(totalNFineFaces_, sumOp<label>());    

    if (Pstream::master())
    {
        Info<< "Total number of clusters : " << totalNCoarseFaces_ << endl;
        Info<< "Total number of fine faces : " << totalNFineFaces_ << endl;
    }

    labelListIOList subMap
    (
        IOobject
        (
            "subMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelListIOList constructMap
    (
        IOobject
        (
            "constructMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    IOList<label> consMapDim
    (
        IOobject
        (
            "constructMapDim",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );    

    map_.reset
    (
        new distributionMap
        (
            consMapDim[0],
            move(subMap),
            move(constructMap)
        )
    );

    scalarListIOList FmyProc
    (
        IOobject
        (
            "F",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
    scalarListIOList solarLoadFineFacesmyProc
    (
        IOobject
        (
            "solarLoadFineFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    solarLoadFineFacesSize = solarLoadFineFacesmyProc.size();       
    
    scalarListIOList skyViewCoeffmyProc
    (
        IOobject
        (
            "skyViewCoeff",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    skyViewCoeffSize = skyViewCoeffmyProc.size();    
    
    scalarListIOList sunViewCoeffmyProc
    (
        IOobject
        (
            "sunViewCoeff",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    sunViewCoeffSize = sunViewCoeffmyProc.size();

    labelIOList sunskyMapmyProc
    (
        IOobject
        (
            "sunskyMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );    

    labelListIOList globalFaceFaces
    (
        IOobject
        (
            "globalFaceFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ); 

    List<labelList> sunskyMap(Pstream::nProcs());
    sunskyMap[Pstream::myProcNo()] = sunskyMapmyProc;
    Pstream::gatherList(sunskyMap);

    List<labelListList> globalFaceFacesProc(Pstream::nProcs());
    globalFaceFacesProc[Pstream::myProcNo()] = globalFaceFaces;
    Pstream::gatherList(globalFaceFacesProc);

    List<scalarListList> F(Pstream::nProcs());
    F[Pstream::myProcNo()] = FmyProc;
    Pstream::gatherList(F);
    
    List<scalarListList> solarLoadFineFaces(Pstream::nProcs());
    solarLoadFineFaces[Pstream::myProcNo()] = solarLoadFineFacesmyProc;
    Pstream::gatherList(solarLoadFineFaces);         
    
    List<scalarListList> skyViewCoeff(Pstream::nProcs());
    skyViewCoeff[Pstream::myProcNo()] = skyViewCoeffmyProc;
    Pstream::gatherList(skyViewCoeff);    
    
    List<scalarListList> sunViewCoeff(Pstream::nProcs());
    sunViewCoeff[Pstream::myProcNo()] = sunViewCoeffmyProc;
    Pstream::gatherList(sunViewCoeff);        

    globalIndex globalNumbering(nLocalCoarseFaces_);
    globalIndex globalNumberingFine(nLocalFineFaces_);

    solarLoadFineFacesGlobal_.reset
    (
        new scalarListList(solarLoadFineFacesSize)
    );    
    forAll(solarLoadFineFacesGlobal_(), vectorId)
    {
        scalarList globalCoeffs(totalNFineFaces_, 0.0);   
        solarLoadFineFacesGlobal_()[vectorId] = globalCoeffs;
    }
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        insertScalarListListElements
        (
            globalNumberingFine,
            procI,
            sunskyMap,
            globalFaceFacesProc[procI],
            solarLoadFineFaces[procI],
            solarLoadFineFacesGlobal_(),
            "fine"
        );
    }
 
    skyViewCoeffGlobal_.reset
    (
        new scalarListList(sunViewCoeffSize)
    );    
    forAll(skyViewCoeffGlobal_(), vectorId)
    {
        scalarList globalCoeffs(totalNCoarseFaces_, 0.0);   
        skyViewCoeffGlobal_()[vectorId] = globalCoeffs;
    } 
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        insertScalarListListElements
        (
            globalNumbering,
            procI,
            sunskyMap,
            globalFaceFacesProc[procI],
            skyViewCoeff[procI],
            skyViewCoeffGlobal_(),
            "coarse"
        );
    }

    sunViewCoeffGlobal_.reset
    (
        new scalarListList(sunViewCoeffSize)
    );
    forAll(sunViewCoeffGlobal_(), vectorId)
    {
        scalarList globalCoeffs(totalNCoarseFaces_, 0.0);   
        sunViewCoeffGlobal_()[vectorId] = globalCoeffs;
    }        
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        insertScalarListListElements
        (
            globalNumbering,
            procI,
            sunskyMap,
            globalFaceFacesProc[procI],
            sunViewCoeff[procI],
            sunViewCoeffGlobal_(),
            "coarse"
        );          
    }

    if (Pstream::master())
    {
        Fmatrix_.reset
        (
            new scalarSquareMatrix(totalNCoarseFaces_, 0.0)
        );    

        Info<< "Insert elements in the matrix..." << endl;

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            insertMatrixElements
            (
                globalNumbering,
                procI,
                globalFaceFacesProc[procI],
                F[procI],
                Fmatrix_()
            );
        }

        bool smoothing = readBool(coeffs_.lookup("smoothing"));
        if (smoothing)
        {
            Info<< "Smoothing the matrix..." << endl;

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                scalar sumF = 0.0;
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    sumF += Fmatrix_()(i, j);
                }
                scalar delta = sumF - 1.0;
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    Fmatrix_()(i, j) *= (1.0 - delta/sumF);
                }
            }
        }

        constAlbedo_ = readBool(coeffs_.lookup("constantAlbedo"));
        if (constAlbedo_)
        {
            CLU_.reset
            (
                new scalarSquareMatrix
                (
                    totalNCoarseFaces_,
                    0.0
                )
            );

            pivotIndices_.setSize(CLU_().m());
        }
    
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarLoad::directAndDiffuse::directAndDiffuse(const volScalarField& T)
:
    solarLoadModel(typeName, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    map_(),
    coarseMesh_
    (
        IOobject
        (
            mesh_.name(),
            mesh_.polyMesh::instance(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        finalAgglom_
    ),
    qs_
    (
        IOobject
        (
            "qs",
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    solarLoadFineFacesGlobal_(),
    skyViewCoeffGlobal_(),   
    sunViewCoeffGlobal_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    wallPatchOrNot_(mesh_.boundary().size(), 0),    
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    nLocalWallCoarseFaces_(0),
    nLocalFineFaces_(0),        
    totalNFineFaces_(0),
    constAlbedo_(false),
    timestepsInADay_(24),
    iterCounter_(0),
    pivotIndices_(0)
{
    initialise();
}


Foam::solarLoad::directAndDiffuse::directAndDiffuse
(
    const dictionary& dict,
    const volScalarField& T
)
:
    solarLoadModel(typeName, dict, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    map_(),
    coarseMesh_
    (
        IOobject
        (
            mesh_.name(),
            mesh_.polyMesh::instance(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        finalAgglom_
    ),
    qs_
    (
        IOobject
        (
            "Qs",
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    solarLoadFineFacesGlobal_(),    
    skyViewCoeffGlobal_(),
    sunViewCoeffGlobal_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    wallPatchOrNot_(mesh_.boundary().size(), 0),    
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    nLocalWallCoarseFaces_(0),
    nLocalFineFaces_(0),        
    totalNFineFaces_(0),
    constAlbedo_(false),
    timestepsInADay_(24),
    iterCounter_(0),
    pivotIndices_(0)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solarLoad::directAndDiffuse::~directAndDiffuse()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solarLoad::directAndDiffuse::read()
{
    if (solarLoadModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::solarLoad::directAndDiffuse::insertMatrixElements
(
    const globalIndex& globalNumbering,
    const label procI,
    const labelListList& globalFaceFaces,
    const scalarListList& viewFactors,
    scalarSquareMatrix& Fmatrix
)
{
    forAll(viewFactors, faceI)
    {
        const scalarList& vf = viewFactors[faceI];
        const labelList& globalFaces = globalFaceFaces[faceI];

        label globalI = globalNumbering.toGlobal(procI, faceI);
        forAll(globalFaces, i)
        {
            Fmatrix[globalI][globalFaces[i]] = vf[i];
        }
    }

    //Info << "Fmatrix: " << Fmatrix << endl;
}

void Foam::solarLoad::directAndDiffuse::insertScalarListListElements
(
    const globalIndex& globalNumbering,
    const label procI,
    const labelListList& sunskyMap,
    const labelListList& globalFaceFaces,
    const scalarListList& localCoeffs,
    scalarListList& globalCoeffs,
    const word& coarseOrFine
)
{
    forAll(localCoeffs, vectorId)
    {
        const scalarList& vf = localCoeffs[vectorId];
        
        forAll(vf, faceI)
        {
            if (coarseOrFine == "coarse")
            {               
                globalCoeffs[vectorId][sunskyMap[procI][faceI]] = vf[faceI];
            }
            else
            {
                label globalI = globalNumbering.toGlobal(procI, faceI);
                globalCoeffs[vectorId][globalI] = vf[faceI];
            }
        }        
    }
}

void Foam::solarLoad::directAndDiffuse::calculate()
{
    // Store previous iteration
    qs_.storePrevIter();

    scalarField compactCoarseA(map_->constructSize(), 0.0);
    scalarField compactCoarseHo(map_->constructSize(), 0.0);

    globalIndex globalNumbering(nLocalCoarseFaces_);
    globalIndex globalNumberingFine(nLocalFineFaces_);    

    // Fill local averaged Albedo(A) and external heatFlux(Ho)
    DynamicList<scalar> localCoarseAave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseHoave(nLocalCoarseFaces_);

    volScalarField::Boundary& qsBf = qs_.boundaryFieldRef();

    forAll(selectedPatches_, i)
    {
        label patchID = selectedPatches_[i];

        const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

        fvPatchScalarField& qsPatch = qsBf[patchID];

        solarLoadViewFactorFixedValueFvPatchScalarField& qsp =
            refCast
            <
                solarLoadViewFactorFixedValueFvPatchScalarField
            >(qsPatch);

        const scalarList ab = qsp.albedo();

        const scalarList& Hoi = qsp.qso();

        const polyPatch& pp = coarseMesh_.boundaryMesh()[patchID]; 
        const labelList& coarsePatchFace = coarseMesh_.patchFaceMap()[patchID]; 

        scalarList Aave(pp.size(), 0.0);
        scalarList Hoiave(Aave.size(), 0.0);

        if (pp.size() > 0)
        {
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom) + 1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            forAll(coarseToFine, coarseI)
            {
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                UIndirectList<scalar> fineSf
                (
                    sf,
                    fineFaces
                );
                const scalar area = sum(fineSf());
                // albedo and external flux area weighting
                forAll(fineFaces, j)
                {
                    label faceI = fineFaces[j];
                    Aave[coarseI] += (ab[faceI]*sf[faceI])/area;
                    Hoiave[coarseI] += (Hoi[faceI]*sf[faceI])/area;
                }
            }
        }

        //localCoarseTave.append(Tave);
        localCoarseAave.append(Aave);
        localCoarseHoave.append(Hoiave);
    }

    // Fill the local values to distribute
    SubList<scalar>(compactCoarseA,nLocalCoarseFaces_) = localCoarseAave;
    SubList<scalar>(compactCoarseHo,nLocalCoarseFaces_) = localCoarseHoave;

    // Distribute data
    map_->distribute(compactCoarseA);
    map_->distribute(compactCoarseHo);

    // Distribute local global ID
    labelList compactGlobalIds(map_->constructSize(), 0.0);

    labelList localGlobalIds(nLocalCoarseFaces_);

    for(label k = 0; k < nLocalCoarseFaces_; k++)
    {
        localGlobalIds[k] = globalNumbering.toGlobal(Pstream::myProcNo(), k);
    }

    SubList<label>
    (
        compactGlobalIds,
        nLocalCoarseFaces_
    ) = localGlobalIds;

    map_->distribute(compactGlobalIds);

    // Create global size vectors
    scalarField A(totalNCoarseFaces_, 0.0);
    scalarField qsExt(totalNCoarseFaces_, 0.0);

    // Fill lists from compact to global indexes.
    forAll(compactCoarseA, i)
    {
        A[compactGlobalIds[i]] = compactCoarseA[i];
        qsExt[compactGlobalIds[i]] = compactCoarseHo[i];
    }

    Pstream::listCombineGather(A, maxEqOp<scalar>());
    Pstream::listCombineGather(qsExt, maxEqOp<scalar>());

    Pstream::listCombineScatter(A);
    Pstream::listCombineScatter(qsExt);

    // Net solarLoad
    scalarField q(totalNCoarseFaces_, 0.0);
    
    Time& time = const_cast<Time&>(mesh_.time());   
    // Read sunPosVector list
    dictionary sunPosVectorIO;
    sunPosVectorIO.add("type", "table");
    sunPosVectorIO.add(
        "file",
        fileName
        (
            mesh_.time().constant()
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
    // look for the correct range
    label lo = 0;
    label hi = 0;
    scalarField sunPosVector_x = sunPosVector.x();
    forAll(sunPosVector_x, i)
    {
        if (time.value() >= sunPosVector_x[i])
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }   
    }
    scalar hi_fraction = 0; 
    if (lo != hi) //if timestep is between two time values in sunPosVector
    {
        hi_fraction = (time.value() - sunPosVector_x[lo]) / (sunPosVector_x[hi] - sunPosVector_x[lo]);
    }  

    if (Pstream::master())
    {
        // Variable Albedo
        if (!constAlbedo_) //this is not tested - aytac
        {
            scalarSquareMatrix C(totalNCoarseFaces_, 0.0);

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    //scalar invEj = 1.0/E[j];
                    scalar Isol = (skyViewCoeffGlobal_()[lo][j] + sunViewCoeffGlobal_()[lo][j]);
                    if (i==j)
                    {
                        C(i, j) = (1/(1-A[j]))-(A[j]/(1-A[j]))*Fmatrix_()(i, j);
                        q[i] += (- 1.0)*(-Isol) - qsExt[j];
                    }
                    else
                    {
                        C(i, j) = -(A[j]/(1-A[j]))*Fmatrix_()(i, j);
                        q[i] += 0 - qsExt[j];
                    }

                }
            }

            Info<< "\nSolving view factor equations..." << endl;
            // Negative coming into the fluid
            LUsolve(C, q);
        }
        else //Constant albedo
        {
            // Initial iter calculates CLU and chaches it
            if (iterCounter_ == 0)
            {
                for (label i=0; i<totalNCoarseFaces_; i++)
                {
                    for (label j=0; j<totalNCoarseFaces_; j++)
                    { 
                        //scalar invEj = 1/E[j];
                        if (i==j)
                        {
                            CLU_()(i, j) = (1/(1-A[j]))-(A[j]/(1-A[j]))*Fmatrix_()(i, j);
                        }
                        else
                        {
                            CLU_()(i, j) = -(A[j]/(1-A[j]))*Fmatrix_()(i, j);
                        }
                    }
                }
                
                fileName fileCLU
                (
                   mesh_.time().rootPath()
                   /mesh_.time().globalCaseName()
                   /"processor0/CLU_qs"
                ); //under processor0 to avoid keeping CLU file for future uses by mistake
                // Check if file already exists
                IFstream is(fileCLU);
                label testCLU = -1;
                if (is.good())
                {
                    is >> testCLU;
                    if (testCLU == totalNCoarseFaces_)
                    {
                        is >> CLU_() >> pivotIndices_;
                        Info << "Read decomposed C matrix from existing file!" << endl;
                    }
                    else
                    {
                        testCLU = -1;
                        Info << "Warning: File for decomposed C matrix does not match totalNCoarseFaces! Will decompose C matrix again..." << endl;
                    }
                }                
                if (testCLU == -1)
                {                                                                
                    Info<< "\nDecomposing C matrix..." << endl;
                    LUDecompose(CLU_(), pivotIndices_);
                    
                    if (Pstream::nProcs() > 1)
                    {
                        // Write file - only in parallel cases
                        OFstream os(fileCLU);
                        os << totalNCoarseFaces_ << endl;
                        os << CLU_() << endl;
                        os << pivotIndices_ << endl;
                    }
                }                    
            }
            
            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    scalar Isol = (skyViewCoeffGlobal_()[lo][j]*(1-hi_fraction) + skyViewCoeffGlobal_()[hi][j]*(hi_fraction)
                                 + sunViewCoeffGlobal_()[lo][j]*(1-hi_fraction) + sunViewCoeffGlobal_()[hi][j]*(hi_fraction));
                    if (i==j)
                    {
                        q[i] += (- 1.0)*(-Isol) - qsExt[j];
                    }
                    else
                    {
                        q[i] += 0 - qsExt[j];
                    }
                }
            }

            Info<< "\nLU Back substitute C matrix.." << endl;
            LUBacksubstitute(CLU_(), pivotIndices_, q);
            iterCounter_ ++;
        }
    }

    // Scatter q and fill qs
    Pstream::listCombineScatter(q);
    Pstream::listCombineGather(q, maxEqOp<scalar>());
    
    Pstream::listCombineScatter(A);
    Pstream::listCombineGather(A, maxEqOp<scalar>());    

    label globCoarseId = 0;
    //label globFineId = 0;    
    label fineFaceNo = 0;
    forAll(selectedPatches_, i)
    {
        const label patchID = selectedPatches_[i];
        const polyPatch& pp = mesh_.boundaryMesh()[patchID];
        if (pp.size() > 0)
        {
            scalarField& qsp = qsBf[patchID];
            const scalarField& sf = mesh_.magSf().boundaryField()[patchID];
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom)+1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            const labelList& coarsePatchFace =
                coarseMesh_.patchFaceMap()[patchID];

            scalar heatFlux = 0.0;
            forAll(coarseToFine, coarseI)
            {
                label globalCoarse =
                    globalNumbering.toGlobal(Pstream::myProcNo(), globCoarseId);
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                forAll(fineFaces, k)
                {
                    label faceI = fineFaces[k];

                    qsp[faceI] = q[globalCoarse];
                    if
                    (
                        isA<wallFvPatch>(mesh_.boundary()[patchID])
                     || isA<mappedInternalFvPatch>(mesh_.boundary()[patchID])
                    )
                    {
                        label globalFine =
                            globalNumberingFine.toGlobal(Pstream::myProcNo(), fineFaceNo+faceI);   
                        qsp[faceI] -= (sunViewCoeffGlobal_()[lo][globalCoarse]*(1-hi_fraction) + sunViewCoeffGlobal_()[hi][globalCoarse]*(hi_fraction)) * (1-A[globalCoarse]);
                        qsp[faceI] += (solarLoadFineFacesGlobal_()[lo][globalFine]*(1-hi_fraction) + solarLoadFineFacesGlobal_()[hi][globalFine]*(hi_fraction)) * (1-A[globalCoarse]);
                    }
                    heatFlux += qsp[faceI]*sf[faceI];
                }
                globCoarseId ++;
            }
        }
        if
        (
            isA<wallFvPatch>(mesh_.boundary()[patchID])
         || isA<mappedInternalFvPatch>(mesh_.boundary()[patchID])
        )
        {
            fineFaceNo += pp.size();
        }
    }

    if (debug)
    {
        forAll(qsBf, patchID)
        {
            const scalarField& qsp = qs_.boundaryField()[patchID];
            const scalarField& magSf = mesh_.magSf().boundaryField()[patchID];
            const scalar heatFlux = gSum(qsp*magSf);
            Info<< "Total heat transfer rate at patch: "
                << patchID << " "
                << heatFlux << endl;
        }
    }

    // Relax qs if necessary
    qs_.relax();
}


Foam::tmp<Foam::volScalarField> Foam::solarLoad::directAndDiffuse::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimMass/pow3(dimTime)/dimLength/pow4(dimTemperature),
                0.0
            )
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::solarLoad::directAndDiffuse::Ru() const
{
    return tmp<volScalarField::Internal>
    (
        new volScalarField::Internal 
        (
            IOobject
            (
                "Ru",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}

// ************************************************************************* //
