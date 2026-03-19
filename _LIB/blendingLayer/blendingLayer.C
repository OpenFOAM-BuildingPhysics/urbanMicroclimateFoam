/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "blendingLayer.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blendingLayer, 0);
    
    //- Private class for finding nearest
    //  Comprising:
    //  - sqr(distance)
    //  - face id (local)
    typedef Tuple2<scalar, label> nearInfo;
    
    class nearestEqOp
    {
    
    public:
        void operator()(nearInfo& x, const nearInfo& y) const
        {
            if (y.first() < x.first())
            {
                x = y;
            }
        }    
    };
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingLayer::blendingLayer
(
    const volVectorField& U,
    const volScalarField& T
)
:
    mesh_(U.mesh()),
    time_(U.time()),
    bL_
    (
        IOobject
        (
            "bL",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("bL",dimensionSet(0,0,0,0,0,0,0),-1)
    ), 
    USource_
    (
        IOobject
        (
            "USource_",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("USource_",dimensionSet(0,1,-2,0,0,0,0),vector::zero)
    ),
    TSource_
    (
        IOobject
        (
            "TSource_",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("TSource_",dimensionSet(0,0,-1,1,0,0,0),scalar(0))
    ),
    coeffs_(dictionary::null),
    dampingThickness(60),
    alphaCoeffU(0.3),
    alphaCoeffT(0.1)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingLayer::~blendingLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::blendingLayer::initialize()
{
    coeffs_ = time_.controlDict().subDict("blendingCoeffs");
    dampingThickness = readScalar(coeffs_.lookup("dampingThickness"));
    alphaCoeffU = readScalar(coeffs_.lookup("alphaCoeffU"));
    alphaCoeffT = readScalar(coeffs_.lookup("alphaCoeffT"));

    Info << "blendingLayer dampingThickness: " << dampingThickness << endl; 
    Info << "blendingLayer alphaCoeffU: " << alphaCoeffU << endl; 
    Info << "blendingLayer alphaCoeffT: " << alphaCoeffT << endl; 

    meshSearch ms(mesh_);
    
    const vectorField &centres = mesh_.C();
    scalar minX = mesh_.bounds().min().x();
    scalar maxX = mesh_.bounds().max().x();
    scalar minY = mesh_.bounds().min().y();
    scalar maxY = mesh_.bounds().max().y();
    
    word patches[] = {"west", "east", "north", "south"};
    List<label> patchID_WENS(4);
    
    List<List<label>> startFace(Pstream::nProcs());
    List<List<label>> nFaces(Pstream::nProcs());  

    forAll(patchID_WENS,i)
    {
        patchID_WENS[i] = mesh_.boundaryMesh().findIndex(patches[i]);
        
        startFace[Pstream::myProcNo()].append(mesh_.boundary()[patchID_WENS[i]].start()); //local startFace for this patch
        nFaces[Pstream::myProcNo()].append(mesh_.boundary()[patchID_WENS[i]].size()); //local total face number for this patch        
    }
    
    Pstream::gatherList(startFace);
    Pstream::scatterList(startFace); 
    
    Pstream::gatherList(nFaces);
    Pstream::scatterList(nFaces); 


    DynamicList<Tuple2<label, point>> pointsOnBoundaries;
    DynamicList<label> blendingCells;

    forAll (centres, cellI)
    {
        const vector& cell = centres[cellI];
        if ((cell.x() <=  minX + dampingThickness) && (cell.x()-minX < cell.y()-minY) && (cell.x()-minX < maxY-cell.y())) //WEST
        {
            point boundaryLoc(minX, cell.y(), cell.z());
            Tuple2<label, point> pointOnBoundary(0, boundaryLoc); //patchID_WENS and position on boundary
            pointsOnBoundaries.append(pointOnBoundary);
            
            blendingCells.append(cellI);
        }
        else if ((cell.x() >=  maxX - dampingThickness) && (maxX-cell.x() < cell.y()-minY) && (maxX-cell.x() < maxY-cell.y())) //EAST
        {
            point boundaryLoc(maxX, cell.y(), cell.z());
            Tuple2<label, point> pointOnBoundary(1, boundaryLoc); //patchID_WENS and position on boundary
            pointsOnBoundaries.append(pointOnBoundary);
            
            blendingCells.append(cellI);
        }
        else if ((cell.y() >=  maxY - dampingThickness) && (maxY-cell.y() < cell.x()-minX) && (maxY-cell.y() < maxX-cell.x())) //NORTH
        {
            point boundaryLoc(cell.x(), maxY, cell.z());
            Tuple2<label, point> pointOnBoundary(2, boundaryLoc); //patchID_WENS and position on boundary
            pointsOnBoundaries.append(pointOnBoundary);
            
            blendingCells.append(cellI);
        }        
        else if ((cell.y() <=  minY + dampingThickness) && (cell.y()-minY < cell.x()-minX) && (cell.y()-minY < maxX-cell.x())) //SOUTH
        {
            point boundaryLoc(cell.x(), minY, cell.z());
            Tuple2<label, point> pointOnBoundary(3, boundaryLoc); //patchID_WENS and position on boundary
            pointsOnBoundaries.append(pointOnBoundary);
            
            blendingCells.append(cellI);
        }           
    }

    blendingCells.shrink();
    pointsOnBoundaries.shrink();

    List<List<Tuple2<label, point>>> pointsOnBoundaries_All(Pstream::nProcs());    
    pointsOnBoundaries_All[Pstream::myProcNo()] = pointsOnBoundaries;
    Pstream::gatherList(pointsOnBoundaries_All);
    Pstream::scatterList(pointsOnBoundaries_All);
    
    List<List<label>> blendingCells_All(Pstream::nProcs());
    blendingCells_All[Pstream::myProcNo()] = blendingCells;
    Pstream::gatherList(blendingCells_All);
    Pstream::scatterList(blendingCells_All);    

    List<List<nearInfo>> nearest(Pstream::nProcs());

    forAll(nearest, proci)
    {
        List<nearInfo>& nearest_ = nearest[proci];
              
        nearest_.setSize(blendingCells_All[proci].size());
       
        forAll(nearest_, i)
        {
            nearest_[i].first() = great;
            nearest_[i].second() = labelMax;
        }
        forAll(nearest_, i)
        {
            label found = ms.findNearestBoundaryFace(pointsOnBoundaries_All[proci][i].second());
            
            label patchId = pointsOnBoundaries_All[proci][i].first();
            label startFace_ = startFace[Pstream::myProcNo()][patchId];
            nearest_[i].first() = magSqr(pointsOnBoundaries_All[proci][i].second() - 
                                            mesh_.boundary()[patchID_WENS[patchId]].Cf()[found-startFace_]);
            nearest_[i].second() = found;
        }

        Pstream::listCombineGather(nearest_, nearestEqOp());
        Pstream::listCombineScatter(nearest_);
    }
    
    bool blendingWarning(0);
    forAll(blendingCells, i)
    {
        label cellId = blendingCells[i];

        label found = nearest[Pstream::myProcNo()][i].second(); //local face id
        label patchId = pointsOnBoundaries[i].first();

        List<label> startFace_n = startFace[Pstream::myProcNo()];
        List<label> nFaces_n = nFaces[Pstream::myProcNo()];
        if ((found >= startFace_n[patchId]) && (found < startFace_n[patchId] + nFaces_n[patchId]))
        {
            label faceId = found - startFace_n[patchId];
            bL_[cellId] = faceId;
        }
        else
        {
            bL_[cellId] = -1000;
            blendingWarning = 1;
        }
    }

    reduce(blendingWarning, orOp<bool>());
    if(blendingWarning)
    {
        Info << "Warning: blendingLayer: boundary face could not be found for some blending cells. Maybe the terrain is not flat everywhere within the blending layer?" << endl;
    }
    Info << "Blending layer initialized" << endl;
    
    return;
}

void Foam::blendingLayer::getValues(volVectorField& USource_, const volVectorField& U)
{
    const vectorField &centres = mesh_.C();
    scalar minX = mesh_.bounds().min().x();
    scalar maxX = mesh_.bounds().max().x();
    scalar minY = mesh_.bounds().min().y();
    scalar maxY = mesh_.bounds().max().y();
    
    label patchId = -1;
    
    patchId = mesh_.boundaryMesh().findIndex("west");
    List<vector> UTarget_W = U.boundaryField()[patchId];
    
    patchId = mesh_.boundaryMesh().findIndex("east");
    List<vector> UTarget_E = U.boundaryField()[patchId];
    
    patchId = mesh_.boundaryMesh().findIndex("north");
    List<vector> UTarget_N = U.boundaryField()[patchId];
        
    patchId = mesh_.boundaryMesh().findIndex("south");
    List<vector> UTarget_S = U.boundaryField()[patchId];

    forAll (centres, cellI)
    {
        const vector& cell = centres[cellI];
        if ((cell.x() <=  minX + dampingThickness) && (cell.x()-minX < cell.y()-minY) && (cell.x()-minX < maxY-cell.y())) //WEST
        {
            label faceId = bL_.internalField()[cellI];
            if(faceId >= 0)
            {
                vector UTarget = UTarget_W[faceId];
                scalar distance = cell.x() - minX;
                scalar sinusInput = (dampingThickness - distance)/dampingThickness;
                USource_[cellI] = (UTarget - U.internalField()[cellI]) * 
                            alphaCoeffU * pow(sin( (Foam::constant::mathematical::pi/2)*sinusInput ), 2.0);
            }
        }
        else if ((cell.x() >=  maxX - dampingThickness) && (maxX-cell.x() < cell.y()-minY) && (maxX-cell.x() < maxY-cell.y())) //EAST
        {
            label faceId = bL_.internalField()[cellI];
            if(faceId >= 0)
            {
                vector UTarget = UTarget_E[faceId];
                scalar distance = maxX - cell.x();
                scalar sinusInput = (dampingThickness - distance)/dampingThickness;
                USource_[cellI] = (UTarget - U.internalField()[cellI]) * 
                            alphaCoeffU * pow(sin( (Foam::constant::mathematical::pi/2)*sinusInput ), 2.0);
            }
        }
        else if ((cell.y() >=  maxY - dampingThickness) && (maxY-cell.y() < cell.x()-minX) && (maxY-cell.y() < maxX-cell.x())) //NORTH
        {
            label faceId = bL_.internalField()[cellI];
            if(faceId >= 0)
            {
                vector UTarget = UTarget_N[faceId];
                scalar distance = maxY - cell.y();
                scalar sinusInput = (dampingThickness - distance)/dampingThickness;
                USource_[cellI] = (UTarget - U.internalField()[cellI]) * 
                            alphaCoeffU * pow(sin( (Foam::constant::mathematical::pi/2)*sinusInput ), 2.0);
            }
        }  
        else if ((cell.y() <=  minY + dampingThickness) && (cell.y()-minY < cell.x()-minX) && (cell.y()-minY < maxX-cell.x())) //SOUTH
        {
            label faceId = bL_.internalField()[cellI];
            if(faceId >= 0)
            {
                vector UTarget = UTarget_S[faceId];
                scalar distance = cell.y() - minY;
                scalar sinusInput = (dampingThickness - distance)/dampingThickness;
                USource_[cellI] = (UTarget - U.internalField()[cellI]) * 
                            alphaCoeffU * pow(sin( (Foam::constant::mathematical::pi/2)*sinusInput ), 2.0);
            }
        }
    }
}

void Foam::blendingLayer::getValues(volScalarField& TSource_, const volScalarField& T)
{
    const vectorField &centres = mesh_.C();
    scalar minX = mesh_.bounds().min().x();
    scalar maxX = mesh_.bounds().max().x();
    scalar minY = mesh_.bounds().min().y();
    scalar maxY = mesh_.bounds().max().y();
    
    label patchId = -1;
    
    patchId = mesh_.boundaryMesh().findIndex("west");
    List<scalar> TTarget_W = T.boundaryField()[patchId];
    
    patchId = mesh_.boundaryMesh().findIndex("east");
    List<scalar> TTarget_E = T.boundaryField()[patchId];
    
    patchId = mesh_.boundaryMesh().findIndex("north");
    List<scalar> TTarget_N = T.boundaryField()[patchId];
        
    patchId = mesh_.boundaryMesh().findIndex("south");
    List<scalar> TTarget_S = T.boundaryField()[patchId];

    forAll (centres, cellI)
    {
        const vector& cell = centres[cellI];
        if ((cell.x() <=  minX + dampingThickness) && (cell.x()-minX < cell.y()-minY) && (cell.x()-minX < maxY-cell.y())) //WEST
        {
            label faceId = bL_.internalField()[cellI];
            if(faceId >= 0)
            {
                scalar TTarget = TTarget_W[faceId];
                scalar distance = cell.x() - minX;
                scalar sinusInput = (dampingThickness - distance)/dampingThickness;
                TSource_[cellI] = (TTarget - T.internalField()[cellI]) * 
                            alphaCoeffT * pow(sin( (Foam::constant::mathematical::pi/2)*sinusInput ), 2.0);
            }
        }
        else if ((cell.x() >=  maxX - dampingThickness) && (maxX-cell.x() < cell.y()-minY) && (maxX-cell.x() < maxY-cell.y())) //EAST
        {
            label faceId = bL_.internalField()[cellI];
            if(faceId >= 0)
            {
                scalar TTarget = TTarget_E[faceId];
                scalar distance = maxX - cell.x();
                scalar sinusInput = (dampingThickness - distance)/dampingThickness;
                TSource_[cellI] = (TTarget - T.internalField()[cellI]) * 
                            alphaCoeffT * pow(sin( (Foam::constant::mathematical::pi/2)*sinusInput ), 2.0);
            }
        }
        else if ((cell.y() >=  maxY - dampingThickness) && (maxY-cell.y() < cell.x()-minX) && (maxY-cell.y() < maxX-cell.x())) //NORTH
        {
            label faceId = bL_.internalField()[cellI];
            if(faceId >= 0)
            {
                scalar TTarget = TTarget_N[faceId];
                scalar distance = maxY - cell.y();
                scalar sinusInput = (dampingThickness - distance)/dampingThickness;
                TSource_[cellI] = (TTarget - T.internalField()[cellI]) * 
                            alphaCoeffT * pow(sin( (Foam::constant::mathematical::pi/2)*sinusInput ), 2.0);
            }
        }  
        else if ((cell.y() <=  minY + dampingThickness) && (cell.y()-minY < cell.x()-minX) && (cell.y()-minY < maxX-cell.x())) //SOUTH
        {
            label faceId = bL_.internalField()[cellI];
            if(faceId >= 0)
            {
                scalar TTarget = TTarget_S[faceId];
                scalar distance = cell.y() - minY;
                scalar sinusInput = (dampingThickness - distance)/dampingThickness;
                TSource_[cellI] = (TTarget - T.internalField()[cellI]) * 
                            alphaCoeffT * pow(sin( (Foam::constant::mathematical::pi/2)*sinusInput ), 2.0);
            }
        }
    }
}

Foam::tmp<Foam::volVectorField> Foam::blendingLayer::bL_USource(const volVectorField& U)
{
    getValues(USource_, U);

    return tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "USource",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            USource_
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::blendingLayer::bL_TSource(const volScalarField& T)
{
    getValues(TSource_, T);    

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "TSource",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            TSource_
        )
    );
}

// ************************************************************************* //
