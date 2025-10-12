/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "viewFactorSky.H"
#include "surfaceFields.H"
#include "constants.H"
#include "greyDiffusiveViewFactorFixedValueFvPatchScalarField.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"

#include "wallFvPatch.H"
#include "TableFile.H"

#include "mappedPatchBase.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
    defineTypeNameAndDebug(viewFactorSky, 0);
    addToRadiationRunTimeSelectionTables(viewFactorSky);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiationModels::viewFactorSky::initialise()
{
    const polyBoundaryMesh& coarsePatches = coarseMesh_.boundaryMesh();
    const volScalarField::Boundary& qrp = qr_.boundaryField();

    label count = 0;
    forAll(qrp, patchi)
    {
        // const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const fvPatchScalarField& qrPatchi = qrp[patchi];

        if ((isA<fixedValueFvPatchScalarField>(qrPatchi)))
        {
            selectedPatches_[count] = qrPatchi.patch().index();
            nLocalCoarseFaces_ += coarsePatches[patchi].size();
            count++;
        }
    }

    selectedPatches_.resize(count--);

    if (debug)
    {
        Pout<< "radiationModels::viewFactorSky::initialise() "
            << "Selected patches:" << selectedPatches_ << endl;
        Pout<< "radiationModels::viewFactorSky::initialise() "
            << "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    }

    totalNCoarseFaces_ = nLocalCoarseFaces_;
    reduce(totalNCoarseFaces_, sumOp<label>());

    if (debug && Pstream::master())
    {
        InfoInFunction
            << "Total number of clusters : " << totalNCoarseFaces_ << endl;
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
        new mapDistribute
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

    List<labelListList> globalFaceFacesProc(Pstream::nProcs());
    globalFaceFacesProc[Pstream::myProcNo()] = globalFaceFaces;
    Pstream::gatherList(globalFaceFacesProc);

    List<scalarListList> F(Pstream::nProcs());
    F[Pstream::myProcNo()] = FmyProc;
    Pstream::gatherList(F);

    globalIndex globalNumbering(nLocalCoarseFaces_);

    if (Pstream::master())
    {
        Fmatrix_.reset
        (
            new scalarSquareMatrix(totalNCoarseFaces_, 0.0)
        );

        if (debug)
        {
            InfoInFunction
                << "Insert elements in the matrix..." << endl;
        }

        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            insertMatrixElements
            (
                globalNumbering,
                proci,
                globalFaceFacesProc[proci],
                F[proci],
                Fmatrix_()
            );
        }


        bool smoothing = readBool(coeffs_.lookup("smoothing"));
        if (smoothing)
        {
            if (debug)
            {
                InfoInFunction
                    << "Smoothing the matrix..." << endl;
            }

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                scalar sumF = 0.0;
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    sumF += Fmatrix_()(i, j);
                }

                const scalar delta = sumF - 1.0;
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    Fmatrix_()(i, j) *= (1.0 - delta/(sumF + 0.001));
                }
            }
        }

        constEmissivity_ = readBool(coeffs_.lookup("constantEmissivity"));
        if (constEmissivity_)
        {
            CLU_.reset
            (
                new scalarSquareMatrix(totalNCoarseFaces_, 0.0)
            );

            pivotIndices_.setSize(CLU_().m());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::viewFactorSky::viewFactorSky(const volScalarField& T)
:
    radiationModel(typeName, T),
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
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    constEmissivity_(false),
    iterCounter_(0),
    pivotIndices_(0)
{
    initialise();
}


Foam::radiationModels::viewFactorSky::viewFactorSky
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
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
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    constEmissivity_(false),
    iterCounter_(0),
    pivotIndices_(0)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::viewFactorSky::~viewFactorSky()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiationModels::viewFactorSky::read()
{
    if (radiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiationModels::viewFactorSky::insertMatrixElements
(
    const globalIndex& globalNumbering,
    const label proci,
    const labelListList& globalFaceFaces,
    const scalarListList& viewFactors,
    scalarSquareMatrix& Fmatrix
)
{
    forAll(viewFactors, facei)
    {
        const scalarList& vf = viewFactors[facei];
        const labelList& globalFaces = globalFaceFaces[facei];

        label globalI = globalNumbering.toGlobal(proci, facei);
        forAll(globalFaces, i)
        {
            Fmatrix[globalI][globalFaces[i]] = vf[i];
        }
    }
}


void Foam::radiationModels::viewFactorSky::calculate()
{
    // Store previous iteration
    qr_.storePrevIter();

    scalarField compactCoarseT4(map_->constructSize(), 0.0);
    scalarField compactCoarseE(map_->constructSize(), 0.0);
    scalarField compactCoarseHo(map_->constructSize(), 0.0);

    globalIndex globalNumbering(nLocalCoarseFaces_);

    // Fill local averaged(T), emissivity(E) and external heatFlux(Ho)
    DynamicList<scalar> localCoarseT4ave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseEave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseHoave(nLocalCoarseFaces_);

    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();
    
    //////////////////////////////////////////////////////////////////////////
    //obtain Tambient to calculate Tsky - can find a better way to import Tambient?
    Time& time = const_cast<Time&>(mesh_.time());
    //label timestep = ceil( (time.value()/3600)-1E-6 ); timestep = timestep%24;
    
    dictionary TambientIO;
    TambientIO.add(
        "file", 
        fileName
        (
            mesh_.time().rootPath()
            /mesh_.time().globalCaseName()
            /"0/air/Tambient"
        )
    );
    Function1s::TableFile<scalar> Tambient
    (
        "Tambient",
        TambientIO
    );
    //////////////////////////////////////////////////////////////////////////
    fileName cloudCoverFile
    (
       mesh_.time().rootPath()
       /mesh_.time().globalCaseName()
       /"0/air/cloudCover"
    );
    
    scalar cc = 0; //cloud cover
    if(isFile(cloudCoverFile))
    {
        Info << "Reading cloud cover values..." << endl;
        dictionary cloudCoverIO;
        cloudCoverIO.add(
            "file", 
            cloudCoverFile
        );
        Function1s::TableFile<scalar> cloudCover
        (
            "cloudCover",
            cloudCoverIO
        );        
        cc = cloudCover.value(time.value());
    }
    else
    {
        Info << "Constant cloud cover of 0 is being used..." << endl;
    }
    
    //////////////////////////////////////////////////////////////////////////
    //is grass model activated?
    const polyMesh& airMesh =
        mesh_.time().lookupObject<polyMesh>("air");
    IOdictionary grassProperties
    (
        IOobject
        (
            "grassProperties",
            airMesh.time().constant(),
            airMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    if (grassProperties.typeHeaderOk<IOdictionary>(true))
    {
        word grassModel(grassProperties.lookup("grassModel"));
        if (grassModel != "none")
        {
            const dictionary& modelCoeffs = grassProperties.subDict(grassModel + "Coeffs");
            const List<word>& grassPatches_ = modelCoeffs.lookup("grassPatches");
            grassPatches = grassPatches_;
        }
    }
    //////////////////////////////////////////////////////////////////////////

    forAll(selectedPatches_, i)
    {
        label patchID = selectedPatches_[i];

        const scalarField& Tp = T_.boundaryField()[patchID];
        const scalarField& sf = mesh_.magSf().boundaryField()[patchID];
        scalarField Tg_(T_.size(), -1.0);
        if (grassPatches.found(mesh_.boundary()[patchID].name()))
        {
            if(mesh_.name() == "vegetation")
            {
                // Get the coupling information from the mappedPatchBase
                const mappedPatchBase& mpp =
                    refCast<const mappedPatchBase>(mesh_.boundary()[patchID].patch());
                const polyMesh& nbrMesh = mpp.sampleMesh();
                const label samplePatchI = mpp.samplePolyPatch().index();
                const fvPatch& nbrPatch =
                    refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];
                scalarField TgNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("Tg");
                    mpp.distribute(TgNbr);
                Tg_ = TgNbr;
            }
            else
            {
                Tg_ = mesh_.boundary()[patchID].lookupPatchField<volScalarField, scalar>("Tg");
            }
        }

        fvPatchScalarField& qrPatch = qrBf[patchID];

        greyDiffusiveViewFactorFixedValueFvPatchScalarField& qrp =
            refCast
            <
                greyDiffusiveViewFactorFixedValueFvPatchScalarField
            >(qrPatch);

        const scalarList eb = qrp.emissivity();

        const scalarList& Hoi = qrp.qro();

        const polyPatch& pp = coarseMesh_.boundaryMesh()[patchID];
        const labelList& coarsePatchFace = coarseMesh_.patchFaceMap()[patchID];

        scalarList T4ave(pp.size(), 0.0);
        scalarList Eave(pp.size(), 0.0);
        scalarList Hoiave(pp.size(), 0.0);

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

                // Temperature, emissivity and external flux area weighting
                forAll(fineFaces, j)
                {
                    label facei = fineFaces[j];
                    if (!isA<wallFvPatch>(mesh_.boundary()[patchID])) // added to take into account sky temperature
                    {
                        scalar Tambient_ = Tambient.value(time.value());
                        scalar ec = (1-0.84*cc)*(0.527 + 0.161*Foam::exp(8.45*(1-273/Tambient_))) +0.84*cc; //cloud emissivity
                        scalar Tsky = pow(9.365574E-6*(1-cc)*pow(Tambient_,6) + pow(Tambient_,4)*cc*ec ,0.25); // Swinbank model (1963, Cole 1976)
                        T4ave[coarseI] += (pow4(Tsky)*sf[facei])/area;
                    }
                    else
                    {
                        if (grassPatches.found(mesh_.boundary()[patchID].name()))//use Tg if patch is covered with grass
                        {
                            T4ave[coarseI] += (pow4(Tg_[facei])*sf[facei])/area;
                        }
                        else//otherwise use T wall temperature
                        {
                            T4ave[coarseI] += (pow4(Tp[facei])*sf[facei])/area;
                        }
                    }
                    Eave[coarseI] += (eb[facei]*sf[facei])/area;
                    Hoiave[coarseI] += (Hoi[facei]*sf[facei])/area;
                }
            }
        }

        localCoarseT4ave.append(T4ave);
        localCoarseEave.append(Eave);
        localCoarseHoave.append(Hoiave);
    }

    // Fill the local values to distribute
    SubList<scalar>(compactCoarseT4, nLocalCoarseFaces_) = localCoarseT4ave;
    SubList<scalar>(compactCoarseE, nLocalCoarseFaces_) = localCoarseEave;
    SubList<scalar>(compactCoarseHo, nLocalCoarseFaces_) = localCoarseHoave;

    // Distribute data
    map_->distribute(compactCoarseT4);
    map_->distribute(compactCoarseE);
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
    scalarField T4(totalNCoarseFaces_, 0.0);
    scalarField E(totalNCoarseFaces_, 0.0);
    scalarField qrExt(totalNCoarseFaces_, 0.0);

    // Fill lists from compact to global indexes.
    forAll(compactCoarseT4, i)
    {
        T4[compactGlobalIds[i]] = compactCoarseT4[i];
        E[compactGlobalIds[i]] = compactCoarseE[i];
        qrExt[compactGlobalIds[i]] = compactCoarseHo[i];
    }

    Pstream::listCombineGather(T4, maxEqOp<scalar>());
    Pstream::listCombineGather(E, maxEqOp<scalar>());
    Pstream::listCombineGather(qrExt, maxEqOp<scalar>());

    Pstream::listCombineScatter(T4);
    Pstream::listCombineScatter(E);
    Pstream::listCombineScatter(qrExt);

    // Net radiation
    scalarField q(totalNCoarseFaces_, 0.0);

    if (Pstream::master())
    {
        // Variable emissivity
        if (!constEmissivity_)
        {
            scalarSquareMatrix C(totalNCoarseFaces_, 0.0);

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    const scalar invEj = 1.0/E[j];
                    const scalar sigmaT4 = physicoChemical::sigma.value()*T4[j];

                    if (i==j)
                    {
                        C(i, j) = invEj - (invEj - 1.0)*Fmatrix_()(i, j);
                        q[i] += (Fmatrix_()(i, j) - 1.0)*sigmaT4 - qrExt[j];
                    }
                    else
                    {
                        C(i, j) = (1.0 - invEj)*Fmatrix_()(i, j);
                        q[i] += Fmatrix_()(i, j)*sigmaT4;
                    }

                }
            }

            Info<< "\nSolving view factor equations..." << endl;

            // Negative coming into the fluid
            LUsolve(C, q);
        }
        else // Constant emissivity
        {
            // Initial iter calculates CLU and chaches it
            if (iterCounter_ == 0)
            {
                for (label i=0; i<totalNCoarseFaces_; i++)
                {
                    for (label j=0; j<totalNCoarseFaces_; j++)
                    {
                        const scalar invEj = 1.0/E[j];
                        if (i==j)
                        {
                            CLU_()(i, j) = invEj-(invEj-1.0)*Fmatrix_()(i, j);
                        }
                        else
                        {
                            CLU_()(i, j) = (1.0 - invEj)*Fmatrix_()(i, j);
                        }
                    }
                }

                fileName fileCLU
                (
                    mesh_.time().rootPath()
                    /mesh_.time().globalCaseName()
                    /"processor0/CLU_qr"
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
                    const scalar sigmaT4 =
                        constant::physicoChemical::sigma.value()*T4[j];

                    if (i==j)
                    {
                        q[i] += (Fmatrix_()(i, j) - 1.0)*sigmaT4  - qrExt[j];
                    }
                    else
                    {
                        q[i] += Fmatrix_()(i, j)*sigmaT4;
                    }
                }
            }

            Info<< "\nLU Back substitute C matrix.." << endl;
            LUBacksubstitute(CLU_(), pivotIndices_, q);
            iterCounter_ ++;
        }
    }

    // Scatter q and fill qr
    Pstream::listCombineScatter(q);
    Pstream::listCombineGather(q, maxEqOp<scalar>());

    label globCoarseId = 0;
    forAll(selectedPatches_, i)
    {
        const label patchID = selectedPatches_[i];
        const polyPatch& pp = mesh_.boundaryMesh()[patchID];
        if (pp.size() > 0)
        {
            scalarField& qrp = qrBf[patchID];
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
                    label facei = fineFaces[k];

                    qrp[facei] = q[globalCoarse];
                    heatFlux += qrp[facei]*sf[facei];
                }
                globCoarseId ++;
            }
        }
    }

    if (debug)
    {
        forAll(qrBf, patchID)
        {
            const scalarField& qrp = qrBf[patchID];
            const scalarField& magSf = mesh_.magSf().boundaryField()[patchID];
            const scalar heatFlux = gSum(qrp*magSf);

            InfoInFunction
                << "Total heat transfer rate at patch: "
                << patchID << " "
                << heatFlux << endl;
        }
    }

    // Relax qr if necessary
    qr_.relax();
}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::viewFactorSky::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        mesh_,
        dimensionedScalar
        (
            dimMass/pow3(dimTime)/dimLength/pow4(dimTemperature),
            0
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::viewFactorSky::Ru() const
{
    return volScalarField::Internal::New
    (
        "Ru",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}

// ************************************************************************* //
