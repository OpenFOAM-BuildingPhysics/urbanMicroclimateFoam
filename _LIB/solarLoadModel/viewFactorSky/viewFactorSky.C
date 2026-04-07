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
#include "Function1.H"
#include "Table.H"
#include "IFstream.H"
#include "OFstream.H"
#include "OSspecific.H"

//#include "mappedPatchBase.H" //v8
#include "mappedFvPatchBaseBase.H" //v12
#include "mappedInternalFvPatch.H"

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

    globalFaceFaces_.setSize(globalFaceFaces.size());
    viewFactors_.setSize(FmyProc.size());

    forAll(globalFaceFaces_, facei)
    {
        globalFaceFaces_[facei] = globalFaceFaces[facei];
        viewFactors_[facei] = FmyProc[facei];
    }

    bool smoothing = readBool(coeffs_.lookup("smoothing"));
    if (smoothing)
    {
        if (debug)
        {
            Pout<< "Smoothing local sparse view-factor rows..." << endl;
        }

        forAll(viewFactors_, facei)
        {
            scalarList& vf = viewFactors_[facei];

            scalar sumF = 0.0;
            forAll(vf, i)
            {
                sumF += vf[i];
            }

            const scalar delta = sumF - 1.0;
            const scalar scale = 1.0 - delta/(sumF + 0.001);
            forAll(vf, i)
            {
                vf[i] *= scale;
            }
        }
    }

    constEmissivity_ = readBool(coeffs_.lookup("constantEmissivity"));

    globalIndex globalNumbering(nLocalCoarseFaces_);
    localGlobalIds_.setSize(nLocalCoarseFaces_);
    forAll(localGlobalIds_, i)
    {
        localGlobalIds_[i] = globalNumbering.toGlobal(Pstream::myProcNo(), i);
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
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    globalFaceFaces_(),
    viewFactors_(),
    localGlobalIds_(),
    qPrev_(),
    qPrevValid_(false),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    constEmissivity_(false),
    grassPatches()
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
            mesh_.time().name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    globalFaceFaces_(),
    viewFactors_(),
    localGlobalIds_(),
    qPrev_(),
    qPrevValid_(false),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    constEmissivity_(false),
    grassPatches()
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


void Foam::radiationModels::viewFactorSky::assembleGlobal
(
    scalarField& x
)
const
{
    Pstream::listCombineGather(x, plusEqOp<scalar>());
    Pstream::listCombineScatter(x);
}


Foam::scalar Foam::radiationModels::viewFactorSky::localDot
(
    const scalarField& a,
    const scalarField& b
)
const
{
    scalar result = 0.0;
    scalar correction = 0.0;
    forAll(localGlobalIds_, i)
    {
        const label globalI = localGlobalIds_[i];
        const scalar y = a[globalI]*b[globalI] - correction;
        const scalar t = result + y;
        correction = (t - result) - y;
        result = t;
    }

    reduce(result, sumOp<scalar>());
    return result;
}


void Foam::radiationModels::viewFactorSky::multiply
(
    const scalarField& x,
    const scalarField& E,
    scalarField& Ax
)
const
{
    Ax = scalarField(totalNCoarseFaces_, 0.0);

    forAll(viewFactors_, facei)
    {
        const label globalI = localGlobalIds_[facei];
        const scalar invEi = 1.0/E[globalI];
        scalar Axi = invEi*x[globalI];
        scalar correction = 0.0;

        const scalarList& vf = viewFactors_[facei];
        const labelList& globalFaces = globalFaceFaces_[facei];

        forAll(globalFaces, i)
        {
            const label globalJ = globalFaces[i];
            const scalar invEj = 1.0/E[globalJ];

            if (globalJ == globalI)
            {
                const scalar y =
                  - (invEj - 1.0)*vf[i]*x[globalJ] - correction;
                const scalar t = Axi + y;
                correction = (t - Axi) - y;
                Axi = t;
            }
            else
            {
                const scalar y =
                    (1.0 - invEj)*vf[i]*x[globalJ] - correction;
                const scalar t = Axi + y;
                correction = (t - Axi) - y;
                Axi = t;
            }
        }

        Ax[globalI] = Axi;
    }

    assembleGlobal(Ax);
}


Foam::scalarField Foam::radiationModels::viewFactorSky::solveViewFactorSystem
(
    const scalarField& b,
    const scalarField& E
)
const
{
    const label maxIter =
        coeffs_.lookupOrDefault<label>("viewFactorMaxIter", 10000);
    const scalar tolerance =
        coeffs_.lookupOrDefault<scalar>("viewFactorTolerance", 1e-14);
    const scalar relTol =
        coeffs_.lookupOrDefault<scalar>("viewFactorRelTol", 1e-14);
    const bool usePreviousSolution =
        coeffs_.lookupOrDefault<bool>("viewFactorUsePreviousSolution", true);
    const label residualReplacementInterval =
        coeffs_.lookupOrDefault<label>
        (
            "viewFactorResidualReplacementInterval",
            1
        );

    scalarField diag(totalNCoarseFaces_, 0.0);
    forAll(viewFactors_, facei)
    {
        const label globalI = localGlobalIds_[facei];
        const scalar invEi = 1.0/E[globalI];
        diag[globalI] = invEi;

        const scalarList& vf = viewFactors_[facei];
        const labelList& globalFaces = globalFaceFaces_[facei];
        forAll(globalFaces, i)
        {
            if (globalFaces[i] == globalI)
            {
                diag[globalI] -= (invEi - 1.0)*vf[i];
                break;
            }
        }
    }
    assembleGlobal(diag);

    scalarField x(totalNCoarseFaces_, 0.0);
    if
    (
        usePreviousSolution
     && qPrevValid_
     && qPrev_.size() == totalNCoarseFaces_
    )
    {
        x = qPrev_;
    }

    scalarField Ax(totalNCoarseFaces_, 0.0);
    multiply(x, E, Ax);

    scalarField r(b);
    forAll(r, i)
    {
        r[i] -= Ax[i];
    }
    scalarField r0(r);
    scalarField p(totalNCoarseFaces_, 0.0);
    scalarField v(totalNCoarseFaces_, 0.0);
    scalarField s(totalNCoarseFaces_, 0.0);
    scalarField t(totalNCoarseFaces_, 0.0);
    scalarField pHat(totalNCoarseFaces_, 0.0);
    scalarField sHat(totalNCoarseFaces_, 0.0);

    const scalar normB = Foam::sqrt(localDot(b, b));
    const scalar stop =
        max(tolerance, relTol*max(normB, VSMALL));

    scalar residual = Foam::sqrt(localDot(r, r));
    if (Pstream::master())
    {
        Info<< "\nSolving sparse view factor equations, initial residual = "
            << residual << ", stop = " << stop << endl;
    }

    scalar rho = 1.0;
    scalar alpha = 1.0;
    scalar omega = 1.0;

    label iter = 0;
    for (iter = 0; iter < maxIter && residual > stop; iter++)
    {
        const scalar rhoNew = localDot(r0, r);
        if (mag(rhoNew) < VSMALL || mag(omega) < VSMALL)
        {
            WarningInFunction
                << "BiCGStab breakdown while solving sparse view factors"
                << endl;
            break;
        }

        const scalar beta = (rhoNew/rho)*(alpha/omega);

        forAll(p, i)
        {
            p[i] = r[i] + beta*(p[i] - omega*v[i]);
        }

        pHat = scalarField(totalNCoarseFaces_, 0.0);
        forAll(localGlobalIds_, i)
        {
            const label globalI = localGlobalIds_[i];
            pHat[globalI] = p[globalI]/stabilise(diag[globalI], VSMALL);
        }
        assembleGlobal(pHat);

        multiply(pHat, E, v);

        const scalar r0v = localDot(r0, v);
        if (mag(r0v) < VSMALL)
        {
            WarningInFunction
                << "BiCGStab alpha breakdown while solving sparse view factors"
                << endl;
            break;
        }

        alpha = rhoNew/r0v;

        forAll(s, i)
        {
            s[i] = r[i] - alpha*v[i];
        }

        residual = Foam::sqrt(localDot(s, s));
        if (residual <= stop)
        {
            forAll(x, i)
            {
                x[i] += alpha*pHat[i];
            }
            break;
        }

        sHat = scalarField(totalNCoarseFaces_, 0.0);
        forAll(localGlobalIds_, i)
        {
            const label globalI = localGlobalIds_[i];
            sHat[globalI] = s[globalI]/stabilise(diag[globalI], VSMALL);
        }
        assembleGlobal(sHat);

        multiply(sHat, E, t);

        const scalar tt = localDot(t, t);
        if (mag(tt) < VSMALL)
        {
            WarningInFunction
                << "BiCGStab omega breakdown while solving sparse view factors"
                << endl;
            break;
        }

        omega = localDot(t, s)/tt;

        forAll(x, i)
        {
            x[i] += alpha*pHat[i] + omega*sHat[i];
            r[i] = s[i] - omega*t[i];
        }

        rho = rhoNew;
        if
        (
            residualReplacementInterval > 0
         && ((iter + 1) % residualReplacementInterval) == 0
        )
        {
            Ax = scalarField(totalNCoarseFaces_, 0.0);
            multiply(x, E, Ax);
            forAll(r, i)
            {
                r[i] = b[i] - Ax[i];
            }
        }

        residual = Foam::sqrt(localDot(r, r));
    }

    Ax = scalarField(totalNCoarseFaces_, 0.0);
    multiply(x, E, Ax);

    scalarField trueR(b);
    forAll(trueR, i)
    {
        trueR[i] -= Ax[i];
    }

    const scalar trueResidual = Foam::sqrt(localDot(trueR, trueR));
    const scalar trueRelResidual = trueResidual/max(normB, VSMALL);

    qPrev_ = x;
    qPrevValid_ = true;

    if (Pstream::master())
    {
        Info<< "Sparse view factor solve completed in " << iter
            << " iterations, final residual = " << residual
            << ", true residual = " << trueResidual
            << ", true relative residual = " << trueRelResidual << endl;
    }

    return x;
}


void Foam::radiationModels::viewFactorSky::referenceCheck
(
    const scalarField& q,
    const scalarField& b,
    const scalarField& E
)
const
{
    if (!coeffs_.lookupOrDefault<bool>("viewFactorReferenceCheck", false))
    {
        return;
    }

    const label maxFaces =
        coeffs_.lookupOrDefault<label>("viewFactorReferenceCheckMaxFaces", 25000);

    if (totalNCoarseFaces_ > maxFaces)
    {
        if (Pstream::master())
        {
            Info<< "Skipping dense thermal view-factor reference check: "
                << totalNCoarseFaces_ << " faces exceeds "
                << maxFaces << endl;
        }
        return;
    }

    List<labelListList> globalFacesByProc(Pstream::nProcs());
    globalFacesByProc[Pstream::myProcNo()] = globalFaceFaces_;
    Pstream::gatherList(globalFacesByProc);

    List<scalarListList> viewFactorsByProc(Pstream::nProcs());
    viewFactorsByProc[Pstream::myProcNo()] = viewFactors_;
    Pstream::gatherList(viewFactorsByProc);

    List<labelList> localGlobalIdsByProc(Pstream::nProcs());
    localGlobalIdsByProc[Pstream::myProcNo()] = localGlobalIds_;
    Pstream::gatherList(localGlobalIdsByProc);

    if (!Pstream::master())
    {
        return;
    }

    scalarSquareMatrix C(totalNCoarseFaces_, 0.0);

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        const labelListList& globalFaces = globalFacesByProc[proci];
        const scalarListList& viewFactors = viewFactorsByProc[proci];
        const labelList& globalIds = localGlobalIdsByProc[proci];

        forAll(viewFactors, facei)
        {
            const label globalI = globalIds[facei];
            const scalarList& vf = viewFactors[facei];
            const labelList& faces = globalFaces[facei];

            C(globalI, globalI) = 1.0/E[globalI];

            forAll(faces, i)
            {
                const label globalJ = faces[i];
                const scalar invEj = 1.0/E[globalJ];

                if (globalI == globalJ)
                {
                    C(globalI, globalJ) -= (invEj - 1.0)*vf[i];
                }
                else
                {
                    C(globalI, globalJ) = (1.0 - invEj)*vf[i];
                }
            }
        }
    }

    scalarField qRef(b);
    LUsolve(C, qRef);

    scalar maxDiff = 0.0;
    scalar sumDiff = 0.0;
    label maxI = -1;

    forAll(qRef, i)
    {
        const scalar diff = mag(q[i] - qRef[i]);
        sumDiff += diff;
        if (diff > maxDiff)
        {
            maxDiff = diff;
            maxI = i;
        }
    }

    Info<< "Dense thermal view-factor reference check: max|dq| = "
        << maxDiff << " at coarse face " << maxI
        << ", mean|dq| = " << sumDiff/max(label(1), qRef.size())
        << endl;

    if (coeffs_.lookupOrDefault<bool>("viewFactorReferenceDump", false))
    {
        scalarField rowSum(totalNCoarseFaces_, 0.0);

        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            const scalarListList& viewFactors = viewFactorsByProc[proci];
            const labelList& globalIds = localGlobalIdsByProc[proci];

            forAll(viewFactors, facei)
            {
                const label globalI = globalIds[facei];
                rowSum[globalI] = sum(viewFactors[facei]);
            }
        }

        static label dumpCounter = 0;
        const fileName dumpFile
        (
            mesh_.time().rootPath()
           /mesh_.time().globalCaseName()
           /"processor0"
           /("viewFactorSky_sparseDump_"
           + Foam::name(dumpCounter++)
           + "_"
           + Time::timeName(mesh_.time().value(), 12))
        );

        OFstream os(dumpFile);
        os.precision(17);
        os << "n " << totalNCoarseFaces_ << nl;
        os << "time " << mesh_.time().value() << nl;
        os << "b " << b << nl;
        os << "q " << q << nl;
        os << "E " << E << nl;
        os << "rowSum " << rowSum << nl;

        Info<< "Wrote sparse thermal view-factor reference dump: "
            << dumpFile << endl;
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
    TambientIO.add("type", "table");
    TambientIO.add(
        "file",
        fileName
        (
            "$FOAM_CASE/0/air/Tambient"
        )
    );
    Function1s::Table<scalar> Tambient
    (
        "Tambient",
        dimTime,
        dimTemperature,
        TambientIO
    );
    //////////////////////////////////////////////////////////////////////////
    fileName cloudCoverFile
    (
       "$FOAM_CASE/0/air/cloudCover"
    );
    
    scalar cc = 0; //cloud cover
    if(isFile(cloudCoverFile))
    {
        Info << "Reading cloud cover values..." << endl;
        dictionary cloudCoverIO;
        cloudCoverIO.add("type", "table");
        cloudCoverIO.add(
            "file",
            cloudCoverFile
        );
        Function1s::Table<scalar> cloudCover
        (
            "cloudCover",
            dimTime,
            dimless,
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
    if (grassProperties.headerOk())
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
                // Get the mapper and the neighbouring patch
                //v8: const mappedPatchBase& mpp =
                //v8:     refCast<const mappedPatchBase>(mesh_.boundary()[patchID].patch());
                //v8: const polyMesh& nbrMesh = mpp.sampleMesh();
                //v8: const fvPatch& nbrPatch =
                //v8:     refCast<const fvMesh>(nbrMesh).boundary()[mpp.samplePolyPatch().index()];
                //v8: scalarField TgNbr = nbrPatch.lookupPatchField<volScalarField, scalar>("Tg");
                //v8: mpp.distribute(TgNbr);
                const mappedFvPatchBaseBase& mapperVeg =
                    mappedFvPatchBaseBase::getMap(mesh_.boundary()[patchID]);
                const fvPatch& nbrPatch = mapperVeg.nbrFvPatch();
                tmp<scalarField> TgNbr = mapperVeg.fromNeighbour(
                    nbrPatch.lookupPatchField<volScalarField, scalar>("Tg"));
                Tg_ = TgNbr();
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
                    // v12: mappedInternal patches are physical surfaces, not sky
                    if
                    (
                        !isA<wallFvPatch>(mesh_.boundary()[patchID])
                     && !isA<mappedInternalFvPatch>(mesh_.boundary()[patchID])
                    )
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
    scalarField b(totalNCoarseFaces_, 0.0);
    forAll(viewFactors_, facei)
    {
        const label globalI = localGlobalIds_[facei];
        scalar bi =
            - constant::physicoChemical::sigma.value()*T4[globalI]
          - qrExt[globalI];
        scalar correction = 0.0;

        const scalarList& vf = viewFactors_[facei];
        const labelList& globalFaces = globalFaceFaces_[facei];
        forAll(globalFaces, i)
        {
            const label globalJ = globalFaces[i];
            const scalar y =
                vf[i]*constant::physicoChemical::sigma.value()*T4[globalJ]
              - correction;
            const scalar t = bi + y;
            correction = (t - bi) - y;
            bi = t;
        }

        b[globalI] = bi;
    }
    assembleGlobal(b);

    scalarField q(solveViewFactorSystem(b, E));
    referenceCheck(q, b, E);

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
