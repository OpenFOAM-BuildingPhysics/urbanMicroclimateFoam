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

    globalIndex globalNumbering(nLocalCoarseFaces_);

    globalFaceFaces_.setSize(globalFaceFaces.size());
    viewFactors_.setSize(FmyProc.size());
    forAll(globalFaceFaces_, faceI)
    {
        globalFaceFaces_[faceI] = globalFaceFaces[faceI];
        viewFactors_[faceI] = FmyProc[faceI];
    }

    bool smoothing = readBool(coeffs_.lookup("smoothing"));
    if (smoothing)
    {
        Info<< "Smoothing local sparse view-factor rows..." << endl;
        forAll(viewFactors_, faceI)
        {
            scalarList& vf = viewFactors_[faceI];

            scalar sumF = 0.0;
            forAll(vf, i)
            {
                sumF += vf[i];
            }

            const scalar scale = sumF > VSMALL ? 1.0/sumF : 1.0;
            forAll(vf, i)
            {
                vf[i] *= scale;
            }
        }
    }

    constAlbedo_ = readBool(coeffs_.lookup("constantAlbedo"));

    localGlobalIds_.setSize(nLocalCoarseFaces_);
    forAll(localGlobalIds_, i)
    {
        localGlobalIds_[i] = globalNumbering.toGlobal(Pstream::myProcNo(), i);
    }

    sunskyMap_.setSize(sunskyMapmyProc.size());
    forAll(sunskyMap_, i)
    {
        sunskyMap_[i] = sunskyMapmyProc[i];
    }

    solarLoadFineFaces_.setSize(solarLoadFineFacesmyProc.size());
    skyViewCoeff_.setSize(skyViewCoeffmyProc.size());
    sunViewCoeff_.setSize(sunViewCoeffmyProc.size());
    forAll(solarLoadFineFaces_, i)
    {
        solarLoadFineFaces_[i] = solarLoadFineFacesmyProc[i];
    }
    forAll(skyViewCoeff_, i)
    {
        skyViewCoeff_[i] = skyViewCoeffmyProc[i];
    }
    forAll(sunViewCoeff_, i)
    {
        sunViewCoeff_[i] = sunViewCoeffmyProc[i];
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
    globalFaceFaces_(),
    viewFactors_(),
    localGlobalIds_(),
    qPrev_(),
    qPrevValid_(false),
    solarLoadFineFaces_(),
    skyViewCoeff_(),
    sunViewCoeff_(),
    sunskyMap_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    wallPatchOrNot_(mesh_.boundary().size(), 0),    
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    nLocalWallCoarseFaces_(0),
    nLocalFineFaces_(0),        
    totalNFineFaces_(0),
    constAlbedo_(false),
    timestepsInADay_(24)
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
    globalFaceFaces_(),
    viewFactors_(),
    localGlobalIds_(),
    qPrev_(),
    qPrevValid_(false),
    solarLoadFineFaces_(),
    skyViewCoeff_(),
    sunViewCoeff_(),
    sunskyMap_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    wallPatchOrNot_(mesh_.boundary().size(), 0),    
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    nLocalWallCoarseFaces_(0),
    nLocalFineFaces_(0),        
    totalNFineFaces_(0),
    constAlbedo_(false),
    timestepsInADay_(24)
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


void Foam::solarLoad::directAndDiffuse::assembleGlobal
(
    scalarField& x
)
const
{
    Pstream::listCombineGather(x, plusEqOp<scalar>());
    Pstream::listCombineScatter(x);
}


Foam::scalar Foam::solarLoad::directAndDiffuse::localDot
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


void Foam::solarLoad::directAndDiffuse::multiply
(
    const scalarField& x,
    const scalarField& A,
    scalarField& Ax
)
const
{
    Ax = scalarField(totalNCoarseFaces_, 0.0);

    forAll(viewFactors_, faceI)
    {
        const label globalI = localGlobalIds_[faceI];
        scalar Axi = x[globalI]/(1.0 - A[globalI]);
        scalar correction = 0.0;

        const scalarList& vf = viewFactors_[faceI];
        const labelList& globalFaces = globalFaceFaces_[faceI];

        forAll(globalFaces, i)
        {
            const label globalJ = globalFaces[i];
            const scalar y =
              - (A[globalJ]/(1.0 - A[globalJ]))*vf[i]*x[globalJ]
              - correction;
            const scalar t = Axi + y;
            correction = (t - Axi) - y;
            Axi = t;
        }

        Ax[globalI] = Axi;
    }

    assembleGlobal(Ax);
}


Foam::scalarField Foam::solarLoad::directAndDiffuse::solveViewFactorSystem
(
    const scalarField& b,
    const scalarField& A
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
    forAll(viewFactors_, faceI)
    {
        const label globalI = localGlobalIds_[faceI];
        diag[globalI] = 1.0/(1.0 - A[globalI]);

        const scalarList& vf = viewFactors_[faceI];
        const labelList& globalFaces = globalFaceFaces_[faceI];
        forAll(globalFaces, i)
        {
            if (globalFaces[i] == globalI)
            {
                diag[globalI] -=
                    (A[globalI]/(1.0 - A[globalI]))*vf[i];
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
    multiply(x, A, Ax);

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
        Info<< "\nSolving sparse solar view factor equations, "
            << "initial residual = " << residual
            << ", stop = " << stop << endl;
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
                << "BiCGStab breakdown while solving sparse solar view factors"
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

        multiply(pHat, A, v);

        const scalar r0v = localDot(r0, v);
        if (mag(r0v) < VSMALL)
        {
            WarningInFunction
                << "BiCGStab alpha breakdown while solving sparse solar view factors"
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

        multiply(sHat, A, t);

        const scalar tt = localDot(t, t);
        if (mag(tt) < VSMALL)
        {
            WarningInFunction
                << "BiCGStab omega breakdown while solving sparse solar view factors"
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
            multiply(x, A, Ax);
            forAll(r, i)
            {
                r[i] = b[i] - Ax[i];
            }
        }

        residual = Foam::sqrt(localDot(r, r));
    }

    Ax = scalarField(totalNCoarseFaces_, 0.0);
    multiply(x, A, Ax);

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
        Info<< "Sparse solar view factor solve completed in " << iter
            << " iterations, final residual = " << residual
            << ", true residual = " << trueResidual
            << ", true relative residual = " << trueRelResidual << endl;
    }

    return x;
}


void Foam::solarLoad::directAndDiffuse::referenceCheck
(
    const scalarField& q,
    const scalarField& b,
    const scalarField& A
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
            Info<< "Skipping dense solar view-factor reference check: "
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

        forAll(viewFactors, faceI)
        {
            const label globalI = globalIds[faceI];
            const scalarList& vf = viewFactors[faceI];
            const labelList& faces = globalFaces[faceI];

            C(globalI, globalI) = 1.0/(1.0 - A[globalI]);

            forAll(faces, i)
            {
                const label globalJ = faces[i];
                const scalar coeff =
                    -(A[globalJ]/(1.0 - A[globalJ]))*vf[i];

                if (globalI == globalJ)
                {
                    C(globalI, globalJ) += coeff;
                }
                else
                {
                    C(globalI, globalJ) = coeff;
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

    Info<< "Dense solar view-factor reference check: max|dq| = "
        << maxDiff << " at coarse face " << maxI
        << ", mean|dq| = " << sumDiff/max(label(1), qRef.size())
        << endl;
}


Foam::scalarField Foam::solarLoad::directAndDiffuse::globalCoarseCoeff
(
    const scalarListList& localCoeffs,
    const label vectorId
)
const
{
    scalarField coeff(totalNCoarseFaces_, 0.0);

    if (vectorId >= 0 && vectorId < localCoeffs.size())
    {
        const scalarList& localCoeff = localCoeffs[vectorId];
        forAll(localCoeff, faceI)
        {
            coeff[sunskyMap_[faceI]] = localCoeff[faceI];
        }
    }

    assembleGlobal(coeff);
    return coeff;
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

    const scalarField skyLo(globalCoarseCoeff(skyViewCoeff_, lo));
    const scalarField skyHi(globalCoarseCoeff(skyViewCoeff_, hi));
    const scalarField sunLo(globalCoarseCoeff(sunViewCoeff_, lo));
    const scalarField sunHi(globalCoarseCoeff(sunViewCoeff_, hi));

    scalarField Isol(totalNCoarseFaces_, 0.0);
    forAll(Isol, i)
    {
        Isol[i] =
            skyLo[i]*(1 - hi_fraction) + skyHi[i]*hi_fraction
          + sunLo[i]*(1 - hi_fraction) + sunHi[i]*hi_fraction;
    }

    scalarField b(totalNCoarseFaces_, 0.0);
    const scalar qsExtSum = sum(qsExt);
    forAll(localGlobalIds_, i)
    {
        const label globalI = localGlobalIds_[i];
        b[globalI] = Isol[globalI] - qsExtSum;
    }
    assembleGlobal(b);

    q = solveViewFactorSystem(b, A);
    referenceCheck(q, b, A);

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
                        const label localFine = fineFaceNo + faceI;
                        const scalar fineSolar =
                            solarLoadFineFaces_[lo][localFine]*(1 - hi_fraction)
                          + solarLoadFineFaces_[hi][localFine]*hi_fraction;

                        qsp[faceI] -=
                            (sunLo[globalCoarse]*(1 - hi_fraction)
                           + sunHi[globalCoarse]*hi_fraction)
                           *(1 - A[globalCoarse]);

                        qsp[faceI] += fineSolar*(1 - A[globalCoarse]);
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
