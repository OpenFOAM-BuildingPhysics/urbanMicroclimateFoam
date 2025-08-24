/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    calcLAI

Description
    Original version by Lento Manickathan

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "fvc.H"
#include "fvCFD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "distributedTriSurfaceMeshBugFix.H"
#include "cyclicAMIPolyPatch.H"
#include "triSurfaceTools.H"
#include "mapDistribute.H"
//#include "regionProperties.H"
#include "OFstream.H"
#include "meshTools.H"
#include "meshSearch.H"
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
#include "vectorListIOList.H"
#include "scalarIOList.H"
#include "vectorIOList.H"

#include "singleCellFvMesh.H"
#include "interpolation.H"
#include "IOdictionary.H"
#include "fixedValueFvPatchFields.H"
#include "wallFvPatch.H"
//#include "treeDataFace.H"
#include "unitConversion.H"
//#include "fvIOoptionList.H"

#include "TableFile.H"

using namespace Foam;

// calculate the end point for a ray hit check
point calcEndPoint
(
    const point &start,
    const point &n2,
    const point &pminO,
    const point &pmaxO
)
{
  scalar ix = 0; scalar iy = 0; scalar iz = 0;

  if (n2.x() > 0.0)
    ix = (pmaxO.x() - start.x())/n2.x();
  else if (n2.x() < 0.0)
    ix = (pminO.x() - start.x())/n2.x();
  else
    ix = VGREAT;

  if (n2.y() > 0.0)
    iy = (pmaxO.y() - start.y())/n2.y();
  else if (n2.y() < 0.0)
    iy = (pminO.y() - start.y())/n2.y();
  else iy = VGREAT;

  if (n2.z() > 0.0)
    iz = (pmaxO.z() - start.z())/n2.z();
  else if (n2.z() < 0.0)
    iz = (pminO.z() - start.z())/n2.z();
  else
    iz = VGREAT;

  // closest edg direction
  scalar i = min(ix, min(iy, iz));

  return 0.999*i*n2 + start;
}


// trilinear interpolation
scalar interp3D
(
    const point &ptemp,
    const scalarField &valInterp,
    const int &nx,
    const int &ny,
    const point &pmin,
    const point &dp
)
{

  // point &dp, int &i0, int &j0, int &k0, int &nx, int &ny, int &nz)
  //point pmin = pInterp[0];
  //point dp   = pInterp[(nx*ny) + nx + 1] - pInterp[0];

  // Offset from p min.
  double xp = ptemp.x()-pmin.x();
  double yp = ptemp.y()-pmin.y();
  double zp = ptemp.z()-pmin.z();

  // Determine index of lower bound
  int i0 = floor(xp/dp.x());
  int j0 = floor(yp/dp.y());
  int k0 = max(0, floor(zp/dp.z())); // if (k0 < 0) k0 = 0;


  // indices
  int i000 = (nx*ny)*k0 + j0*nx + i0;
  int i100 = (nx*ny)*k0 + j0*nx + i0+1;
  int i010 = (nx*ny)*k0 + (j0+1)*nx + i0;
  int i110 = (nx*ny)*k0 + (j0+1)*nx + i0+1;
  int i001 = (nx*ny)*(k0+1) + j0*nx + i0;
  int i101 = (nx*ny)*(k0+1) + j0*nx + i0+1;
  int i011 = (nx*ny)*(k0+1) + (j0+1)*nx + i0;
  int i111 = (nx*ny)*(k0+1) + (j0+1)*nx + i0+1;

  point pInterp000;// = pInterp[i000];
  pInterp000.x() = pmin.x() + i0*dp.x();
  pInterp000.y() = pmin.y() + j0*dp.y();
  pInterp000.z() = pmin.z() + k0*dp.z();

  // Bilinear interpolation
  double xd = (ptemp.x()-pInterp000.x())/dp.x();
  double yd = (ptemp.y()-pInterp000.y())/dp.y();
  double zd = (ptemp.z()-pInterp000.z())/dp.z();

  // Interpolation in x-dir
  double c00 = valInterp[i000]*(1-xd) + valInterp[i100]*xd;
  double c01 = valInterp[i001]*(1-xd) + valInterp[i101]*xd;
  double c10 = valInterp[i010]*(1-xd) + valInterp[i110]*xd;
  double c11 = valInterp[i011]*(1-xd) + valInterp[i111]*xd;

  // Interpolation in y-dir
  double c0 = c00*(1.0-yd) + c10*yd;
  double c1 = c01*(1.0-yd) + c11*yd;

  // Interpolate in z-dir
  return c0*(1.0-zd) + c1*zd;

}

void calcVegBBOX
(
    const pointField &pmeshC,
    const volScalarField &LAD,
    point &pmin,
    point &pmax
)
{
    point ptemp;

    forAll(LAD, cellI)
    {
        // where vegetation is present
        if (LAD[cellI] > 10*SMALL)
        {
            ptemp = pmeshC[cellI];
            pmin = min(pmin,ptemp);
            pmax = max(pmax,ptemp);
        }
    }

    List<point> pmin_(Pstream::nProcs());
    List<point> pmax_(Pstream::nProcs());

    pmin_[Pstream::myProcNo()] = pmin;
    pmax_[Pstream::myProcNo()] = pmax;
    Pstream::gatherList(pmin_);
    Pstream::scatterList(pmin_);
    Pstream::gatherList(pmax_);
    Pstream::scatterList(pmax_);

    pmin = gMin(pmin_);
    pmax = gMax(pmax_);

}

void interpfvMeshToCartesian
(
    const fvMesh& mesh,
    const volScalarField &LAD,
    point &pmin,
    point &pmax,
    scalarField &LADInterp,
    int &nx,
    int &ny,
    int &nz,
    point &dp,
    const scalar &minCellSizeFactor
)
{

    // Define search mesh
    meshSearch ms(mesh);

    // Define cartesian interpolation grid
    scalar minCellV = gMin(mesh.V()); // Cartesian mesh resolution (determine from minimum cell size)
    scalar minCellL = Foam::pow(minCellV, 1.0/3.0)*minCellSizeFactor;
    Info << "minCellSizeFactor = " << minCellSizeFactor << ", minCellL = " << minCellL << endl;

    dp = vector(minCellL,minCellL,minCellL); // grid spacing

    // Extend the cartesian grid to include vegetation
    pmin -= 5*dp;
    pmax += 5*dp;

    // Define cartesian grid size
    nx = ceil( (pmax.x()-pmin.x()) / dp.x()) + 1;
    ny = ceil( (pmax.y()-pmin.y()) / dp.y()) + 1;
    nz = ceil( (pmax.z()-pmin.z()) / dp.z()) + 1;

    // Generate cartesian interpolation grid
    // coordinates

    DynamicList<scalar> LADInterpDyn;
    DynamicList<label> pIndexDyn;
    ////////////////////////////////////////////////////////////////////////////

    /////////////// Interpolate LAD onto cartesian interpolation mesh

    int cellIndex;
    int pIndex;
    point ptemp;

    for (int k=0; k < nz; k++)
    {
      for (int j=0; j < ny; j++)
      {
        for (int i=0; i < nx; i++)
        {
          // index of node p
          pIndex = (nx*ny)*k + j*nx + i;

          // coordinate of point in original coordinate system
          ptemp.x() = pmin.x() + i*dp.x();
          ptemp.y() = pmin.y() + j*dp.y();
          ptemp.z() = pmin.z() + k*dp.z();

          pmax = max(pmax, ptemp);

          // Find intersecting cell
          cellIndex = ms.findCell(ptemp,-1,true); // faster

          // if point is inside domain
          if (cellIndex != -1)
          {
              if(LAD[cellIndex] > 0)
              {
                  pIndexDyn.append(pIndex);
                  LADInterpDyn.append(LAD[cellIndex]);
              }
          }
        }
      }
    }

    List<DynamicList<scalar>> LADInterpDyn_(Pstream::nProcs());
    LADInterpDyn_[Pstream::myProcNo()] = LADInterpDyn;
    Pstream::gatherList(LADInterpDyn_);

    List<DynamicList<label>> pIndexDyn_(Pstream::nProcs());
    pIndexDyn_[Pstream::myProcNo()] = pIndexDyn;
    Pstream::gatherList(pIndexDyn_);
    
    if (Pstream::master())
    {
        LADInterp.setSize(nx*ny*nz, pTraits<scalar>::zero);
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            List<scalar> LADInterpLocal = LADInterpDyn_[procI];
            List<label> pIndexLocal = pIndexDyn_[procI];
            forAll (LADInterpLocal, i)
            {
                LADInterp[pIndexLocal[i]] = LADInterpLocal[i];
            }
        }
    }    

    //Pstream::scatter(LADInterp);

    // update maximum point
    reduce(pmax, maxOp<point>());

}

void interpcartesianToRotCartesian
(
    point &pminRot,
    point &pmaxRot,
    scalarField &LADInterpRot,
    int &nxRot,
    int &nyRot,
    int &nzRot,
    const tensor &Tinv,
    const scalarField &LADInterp,
    const point& pmin,
    const point& pmax,
    const int &nx,
    const int &ny,
    const point &dp
)
{
    pminRot -= 5*dp;
    pmaxRot += 5*dp;

    // Define rotated cartesian grid size
    nxRot = ceil( (pmaxRot.x()-pminRot.x()) / dp.x()) + 1;
    nyRot = ceil( (pmaxRot.y()-pminRot.y()) / dp.y()) + 1;
    nzRot = ceil( (pmaxRot.z()-pminRot.z()) / dp.z()) + 1;

    ////////////////////////////////////////////////////////////////////
    // Interpolate onto rotated cartesian grid
    int pIndex;
    point ptemp;
    point ptempRot;
    // keeping only on proc0
    if (Pstream::master())
    {
        // Generate rotated cartesian interpolation grid coordinates
        LADInterpRot.setSize(nxRot*nyRot*nzRot, pTraits<scalar>::zero);
        
        for (int k=0; k < nzRot; k++)
        {
          for (int j=0; j < nyRot; j++)
          {
            for (int i=0; i < nxRot; i++)
            {
              // index of node p
              pIndex = (nxRot*nyRot)*k + j*nxRot + i;

              // x,y,z coordinates in rotated coordinate system
              ptempRot.x() = pminRot.x() + i*dp.x();
              ptempRot.y() = pminRot.y() + j*dp.y();
              ptempRot.z() = pminRot.z() + k*dp.z();

              pmaxRot = max(pmaxRot, ptempRot);
    
              // coordinate of point in original coordinate system
              ptemp = transform(Tinv, ptempRot);

              // If point is within the bbox of original cartesian grid
              if ( (ptemp.x() >= pmin.x()) && (ptemp.x() <= pmax.x()) &&
                   (ptemp.y() >= pmin.y()) && (ptemp.y() <= pmax.y()) &&
                   (ptemp.z() >= pmin.z()) && (ptemp.z() <= pmax.z()) )
              LADInterpRot[pIndex] = interp3D(ptemp, LADInterp, nx, ny, pmin, dp);

            }
          }
        }

    }

    // update maximum point
    reduce(pmaxRot, maxOp<point>());

}

void integrateLAD
(
    const scalarField &LADInterpRot,
    const int &nxRot,
    const int &nyRot,
    const int &nzRot,
    const point &dp,
    scalarField &LAIInterpRot
)
{
    int pIndex, pIndexkp1;

    // keeping only on proc0
    if (Pstream::master())
    {
        LAIInterpRot.setSize(nxRot*nyRot*nzRot, pTraits<scalar>::zero);

        for (int i=0; i < nxRot; i++)
        {
          for (int j=0; j < nyRot; j++)
          {
            for (int k=(nzRot-2); k>=0; k--)
            {
              // lower and upper row index
              pIndex = (nxRot*nyRot)*k + j*nxRot + i;
              pIndexkp1 = (nxRot*nyRot)*(k+1) + j*nxRot + i;
              // trapezoidal integration
              LAIInterpRot[pIndex] = LAIInterpRot[pIndexkp1] + 0.5*(LADInterpRot[pIndex]+LADInterpRot[pIndexkp1])*dp.z();

            }
          }
        }

    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    // start timer
    clock_t tstart = std::clock();
    //clock_t tstartlocal = std::clock();
    clock_t tstartStep;

    timeSelector::addOptions();
    #include "addRegionOption.H"

    Foam::argList::addOption
    (
         "writeFields",
         "",
         "write LAI volScalarFields of all time steps"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    Info << "timeDirs: " << timeDirs << endl;
    runTime.setTime(timeDirs[0], 0);

    #include "createNamedMesh.H"

    volScalarField LAD
    (
      IOobject
      (
          "LAD",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      ),
      mesh
    );

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
    Function1s::TableFile<vector> sunPosVector
    (
        "sunPosVector",
        sunPosVectorIO
    );
    scalarField sunPosVector_x = sunPosVector.x();
    vectorField sunPosVector_y = sunPosVector.y();

    scalarListIOList divqrswList
    (
        IOobject
        (
            "divqrsw",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        sunPosVector_y.size()
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
    Function1s::TableFile<scalar> IDN
    (
        "IDN",
        IDNIO
    ); 
    scalarField IDN_y = IDN.y();

    IOdictionary vegetationProperties
    (
        IOobject
        (
            "vegetationProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word vegModel = vegetationProperties.lookup("vegetationModel");
    dictionary coeffs = vegetationProperties.subDict(vegModel + "Coeffs");

    #include "readGravitationalAcceleration.H"
    Info << "Gravity is = " << g << endl;
    const vector ez = - g.value()/mag(g.value());
    Info << "Vertical vector : " << ez << endl;

    scalar kc = dimensioned<scalar>::lookupOrDefault("kc", coeffs, 0.5).value();
    scalar minCellSizeFactor = dimensioned<scalar>::lookupOrDefault("minCellSizeFactor", coeffs, 10).value(); 

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // Global setup
    Info << "Creating interpolation setup...\n" << endl;

    // Mesh setup
    int nMeshCells = mesh.cells().size();

    vector n1(0,0,1); // original vector
    n1 /= mag(n1);

    // Mesh cell centers
    const pointField& pmeshC = mesh.C();

    // mesh bounding box
    point pminO = gMin(pmeshC);
    point pmaxO = gMax(pmeshC);


    // Set up searching engine for obstacles
    #include "searchingEngine.H"

    ////// Determine BBOX of vegetation (original coordinate system)
    point pmin = gMax(pmeshC);
    point pmax = gMin(pmeshC);

    // calculate vegetation bbox - decompose
    calcVegBBOX(pmeshC, LAD, pmin, pmax);

    ////// Define LAD interpolator to arbitrary locations
    tstartStep = std::clock();

    Info << "Interpolation from fvMesh to Cartesian mesh...";

    scalarField LADInterp;
    int nx, ny, nz;
    point dp;

    interpfvMeshToCartesian(mesh, LAD, pmin, pmax, LADInterp, nx, ny, nz, dp, minCellSizeFactor);

    Info << " took " << (std::clock()-tstartStep) / (double)CLOCKS_PER_SEC
         << " second(s)."<< endl;


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    // Coarse mesh faces
    #include "findCoarseMeshFaces.H"


    ////////////////////////////////////////////////////////////////////////////
    // Status Info
    // Info << "Info: Time, Initialization took: "
    //      << (std::clock()-tstart) / (double)CLOCKS_PER_SEC
    //      <<" second(s)."<< endl;

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    Info << "Interpolation from cartesian to rotated cartesian mesh...\n" << endl;

    // Iterate for each sun ray
    int iter = 0;
    label faceNo = 0;
    label insideFaceI = 0;
  	label i = 0;
  	label j = 0;
  	label k = 0;

    forAll(sunPosVector_y, vectorID)
    {
        // start clock
        tstartStep = std::clock();

        ////////////////////////////////////////////////////////////////////////
        // Setup for each sun ray

        // Reference
        scalarList &kcLAIboundary = kcLAIboundaryList[vectorID];
        scalarList &divqrsw = divqrswList[vectorID];

        // Initialize LAI
        scalarList LAI(nMeshCells, 0.0);
        kcLAIboundary.setSize(vegNCoarseFacesAll, 0.0);
        divqrsw.setSize(nMeshCells, 0.0);

        // sunPosVector i
        vector n2 = sunPosVector_y[vectorID];
        n2 /= mag(n2);

        // only if sun is above the horizon
        if ( (n2 & ez) > 0)
        {
            ////////////////////////////////////////////////////////////////////
            // Setup

            // Define rotation matrix
            tensor T(rotationTensor(n2,n1));       // from n1 to n2
            tensor Tinv(rotationTensor(n1,n2));    // from n2 back to n1

            /////////////// (Step 1.d) Determine bbox of vegetation (rotated coordinate system)
            // Calculated in the rotated coordinate system

            // Mesh cell centers (rotated coordinate system)
            pointField pmeshCRot = transform(T,pmeshC);
            // Minimum point
            point pmeshMinRot = gMin(pmeshCRot);
            ////// Calculate BBOX of rotated vegetation
            point pminRot = gMax(pmeshCRot);
            point pmaxRot = gMin(pmeshCRot);

            calcVegBBOX(pmeshCRot, LAD, pminRot, pmaxRot);

            // Generate rotated cartesian interpolation grid
            scalarField LADInterpRot; // interpolated LAD
            int nxRot, nyRot, nzRot;

            interpcartesianToRotCartesian(pminRot, pmaxRot, LADInterpRot, nxRot, nyRot, nzRot, Tinv, LADInterp, pmin, pmax, nx, ny, dp);

            // Info << "gMin(LADInterpRot): " << gMin(LADInterpRot) << endl;
            // Info << "gMax(LADInterpRot): " << gMax(LADInterpRot) << endl;

            ////////////////////////////////////////////////////////////////////
            // Integrate LAD on the rotated cartesian grid

            scalarField LAIInterpRot;
            integrateLAD(LADInterpRot, nxRot, nyRot, nzRot, dp, LAIInterpRot);
            
            // Info << "gMin(LAIInterpRot): " << gMin(LAIInterpRot) << endl;
            // Info << "gMax(LAIInterpRot): " << gMax(LAIInterpRot) << endl;

            ////////////////////////////////////////////////////////////////////
            // Calculated short-wave radiation intensity
            // keeping only on proc0
            tmp<scalarField> qrswInterpRot;
            if (Pstream::master())
            {
                qrswInterpRot = IDN_y[vectorID]*Foam::exp(-kc*LAIInterpRot);
            }
            // Info << "qrswInterpRot: gMin: " << gMin(qrswInterpRot) << endl;
            // Info << "qrswInterpRot: gMax: " << gMax(qrswInterpRot) << endl;

            ////////////////////////////////////////////////////////////////////
            // Calculated divergence of short-wave radiation intensity - forward differencing
            scalarField divqrswInterpRot;
            int pIndex, pIndexkp1;

            // calculating only on proc0
            if (Pstream::master())
            {
                divqrswInterpRot.setSize(nxRot*nyRot*nzRot, pTraits<scalar>::zero);
                const scalarField& qrswInterpRot_ = qrswInterpRot;
                for (int kiter = 0; kiter < (nzRot-1); kiter++)
                {
                    for (int jiter = 0; jiter < nyRot; jiter++)
                    {
                        for (int iiter = 0; iiter < nxRot; iiter++)
                        {
                            pIndex = (nxRot*nyRot)*kiter + jiter*nxRot + iiter; // k
                            pIndexkp1 = (nxRot*nyRot)*(kiter+1) + jiter*nxRot + iiter; // p+1 index (forward)
                            divqrswInterpRot[pIndex] = -(qrswInterpRot_[pIndexkp1] - qrswInterpRot_[pIndex])/dp.z();
                        }
                    }
                }
            }         

            //calcDiv(qrswInterpRot, divqrswInterpRot, nx, ny, nz, dp);
            // Info << "divqrswInterpRot: " << gMax(divqrswInterpRot) << endl;
            // Info << "divqrswInterpRot: " << gMin(divqrswInterpRot) << endl;


            ////////////////////////////////////////////////////////////////////
            // Interpolate LAI from rotated cartesian grid onto original grid

            int nCellsinVegetationBBOX = 0;
            point ptemp;
            forAll(LAD, cellI)
            {
                // Cell center point (rotated coordinate system)
                ptemp = pmeshCRot[cellI];

                if ( (ptemp.x() >= pminRot.x()) && (ptemp.x() <= pmaxRot.x()) &&
                     (ptemp.y() >= pminRot.y()) && (ptemp.y() <= pmaxRot.y()) &&
                     (ptemp.z() >= pmeshMinRot.z()) && (ptemp.z() <= pmaxRot.z()) )
                {
                    nCellsinVegetationBBOX++;
                }
            }

            DynamicField<point> startList(nCellsinVegetationBBOX);
            DynamicField<point> endList(nCellsinVegetationBBOX);
            List<pointIndexHit> pHitList(nCellsinVegetationBBOX);
            DynamicField<label> insideCellIList(nCellsinVegetationBBOX);

            forAll(LAD, cellI)
            {
                // Cell center point (rotated coordinate system)
                ptemp = pmeshCRot[cellI];

                if ( (ptemp.x() >= pminRot.x()) && (ptemp.x() <= pmaxRot.x()) &&
                     (ptemp.y() >= pminRot.y()) && (ptemp.y() <= pmaxRot.y()) &&
                     (ptemp.z() >= pmeshMinRot.z()) && (ptemp.z() <= pmaxRot.z()) )
                {

                    // Check if cell in building shadow shadow
                    point starti = transform(Tinv, ptemp);
                    point endi = calcEndPoint(starti, n2, pminO, pmaxO);
                    const vector& d = endi - starti;
                    startList.append(starti + 0.001*d);
                    endList.append(endi);
                    insideCellIList.append(cellI);
                }
            }

            surfacesMesh.findLine(startList, endList, pHitList);

            DynamicList<point> ptempDyn;
            DynamicList<label> ptempIndexDyn;

            // Updated LAI fields
            forAll(pHitList, rayI)
            {
                label cellI = insideCellIList[rayI];

                if (!pHitList[rayI].hit())
                {
                    ptemp = pmeshCRot[cellI];
                    ptempDyn.append(ptemp);
                    ptempIndexDyn.append(cellI);
                }
                else if (LAD[cellI] > 10*SMALL)
                {
                    LAI[cellI] = 1000; // large enough to ensure qr = exp(-LAI*k)  \approx 0
                }
            }
            
            List<DynamicList<point>> ptempDyn_(Pstream::nProcs());
            ptempDyn_[Pstream::myProcNo()] = ptempDyn;
            Pstream::gatherList(ptempDyn_);
            
            List<List<scalar>> LAI_(Pstream::nProcs());
            LAI_[Pstream::myProcNo()].setSize(ptempDyn.size(), 0.0);
            Pstream::gatherList(LAI_);
            
            List<List<scalar>> divqrsw_(Pstream::nProcs());
            divqrsw_[Pstream::myProcNo()].setSize(ptempDyn.size(), 0.0);
            Pstream::gatherList(divqrsw_);            
            
            if (Pstream::master())
            {
                for (label procI = 0; procI < Pstream::nProcs(); procI++)
                {
                    List<scalar>& LAILocal = LAI_[procI];
                    List<scalar>& divqrswLocal = divqrsw_[procI];
                    forAll (LAILocal, i)
                    {
                        LAILocal[i] = interp3D(ptempDyn_[procI][i], LAIInterpRot, nxRot, nyRot, pminRot, dp);
                        divqrswLocal[i] = interp3D(ptempDyn_[procI][i], divqrswInterpRot, nxRot, nyRot, pminRot, dp);
                    }
                }
            }

            Pstream::listCombineScatter(LAI_);
            Pstream::listCombineScatter(divqrsw_);
            forAll (LAI_[Pstream::myProcNo()], i)
            {
                label cellI = ptempIndexDyn[i];
                LAI[cellI] = LAI_[Pstream::myProcNo()][i];
                if (LAD[cellI] > 10*SMALL)
                {
                    divqrsw[cellI] = divqrsw_[Pstream::myProcNo()][i];
                }
            }
            
            ptempDyn.clear();
            ptempIndexDyn.clear();
            ptempDyn_.clear();
            LAI_.clear();
            divqrsw_.clear();

            //Info << "gMin(LAI): " << gMin(LAI) << endl;
            //Info << "gMax(LAI): " << gMax(LAI) << endl;

            //clock_t tstartlocal = std::clock();

            ////////////////////////////////////////////////////////////////////
            // Interpolate LAI onto coarse mesh faces

            Info << "vegLocalCoarseCf " << vegLocalCoarseCf.size() << endl;

            int nFacesinVegetationBBOX = 0;
            forAll(vegLocalCoarseCf, faceI)
            {
                // Cell center point (rotated coordinate system)
                ptemp = transform(T, vegLocalCoarseCf[faceI]);

                if ( (ptemp.x() >= pminRot.x()) && (ptemp.x() <= pmaxRot.x()) &&
                     (ptemp.y() >= pminRot.y()) && (ptemp.y() <= pmaxRot.y()) &&
                     (ptemp.z() <= pmaxRot.z()) )
                {
                    nFacesinVegetationBBOX++;
                }
            }


            DynamicField<point> vegCoarseFaceStartList(nFacesinVegetationBBOX);
            DynamicField<point> vegCoarseFaceEndList(nFacesinVegetationBBOX);
            List<pointIndexHit> vegCoarseFacePHitList(nFacesinVegetationBBOX);
            DynamicField<label> vegCoarseFaceInsideFaceIList(nFacesinVegetationBBOX);
            List<bool> vegCoarseFaceInsideList(vegLocalCoarseCf.size(), false);

            forAll(vegLocalCoarseCf, faceI)
            {
                // Cell center point (rotated coordinate system)
                ptemp = transform(T, vegLocalCoarseCf[faceI]);

                if ( (ptemp.x() >= pminRot.x()) && (ptemp.x() <= pmaxRot.x()) &&
                     (ptemp.y() >= pminRot.y()) && (ptemp.y() <= pmaxRot.y()) &&
                     (ptemp.z() <= pmaxRot.z()) )
                {

                    // Check if cell in building shadow shadow
                    point starti = vegLocalCoarseCf[faceI];
                    point endi = calcEndPoint(starti, n2, pminO, pmaxO);
                    const vector& d = endi - starti;
                    vegCoarseFaceStartList.append(starti + 0.001*d);
                    vegCoarseFaceEndList.append(endi);
                    vegCoarseFaceInsideFaceIList.append(faceI);
                    vegCoarseFaceInsideList[faceI] = true;
                }
            }

            surfacesMesh.findLine(vegCoarseFaceStartList, vegCoarseFaceEndList, vegCoarseFacePHitList);

            // Updated LAI boundary fields
            // forAll(coarseFacePHitList, rayI)
            // {
            //     label faceI = coarseFaceInsideFaceIList[rayI];
            //
            //     if (!coarseFacePHitList[rayI].hit())
            //     {
            //         ptemp = transform(T, localCoarseCf[faceI]);
            //
            //         LAIboundary[faceI] = interp3D(ptemp, pInterpRot, LAIInterpRot, nxRot, nyRot);
            //     }
            // }

            forAll(vegViewFactorsPatches, patchID)
          	{
                while (i < vegViewFactorsPatches[patchID])
          		  {
          			     while (j < vegHowManyCoarseFacesPerPatch[i])
          			     {
          				         k++;
          				         j++;
          			     }
          			     j = 0;
          			     i++;
          		  }
          		  while (j < vegHowManyCoarseFacesPerPatch[i])
          		  {
                    // if face inside bbox of vegetation
                    if (vegCoarseFaceInsideList[faceNo])
                    {
                        // if face did not hit any wall
                        if (!vegCoarseFacePHitList[insideFaceI].hit())
                        {
                            ptemp = transform(T, vegLocalCoarseCf[faceNo]);
                            ptempDyn.append(ptemp);
                            ptempIndexDyn.append(k);
                            //kcLAIboundary[k] = kc*interp3D(ptemp, LAIInterpRot, nxRot, nyRot, pminRot, dp);
                            //Info << "k = " << k << ", faceNo = " << faceNo << ", insideFaceNo = " << insideFaceNo << endl;
                        }
                        else
                        {
                            kcLAIboundary[k] = 1000;
                        }
                        insideFaceI++;
                    }
          			    k++;
          			    j++;
          			    faceNo++;
          		  }
            }
            i = 0;
            j = 0;
            k = 0;
            faceNo = 0;
            insideFaceI = 0;

            ptempDyn_.setSize(Pstream::nProcs());
            ptempDyn_[Pstream::myProcNo()] = ptempDyn;
            Pstream::gatherList(ptempDyn_);

            LAI_.setSize(Pstream::nProcs());
            LAI_[Pstream::myProcNo()].setSize(ptempDyn.size(), 0.0);
            Pstream::gatherList(LAI_);

            if (Pstream::master())
            {
                for (label procI = 0; procI < Pstream::nProcs(); procI++)
                {
                    List<scalar>& LAILocal = LAI_[procI];
                    forAll (LAILocal, i)
                    {
                        LAILocal[i] = kc*interp3D(ptempDyn_[procI][i], LAIInterpRot, nxRot, nyRot, pminRot, dp);
                    }
                }
            }
                       
            Pstream::listCombineScatter(LAI_);
            forAll (LAI_[Pstream::myProcNo()], i)
            {
                label faceI = ptempIndexDyn[i];
                kcLAIboundary[faceI] = LAI_[Pstream::myProcNo()][i];
            }
            
            ////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////
            iter +=1;
        }

        // Export --- temp.
        if (args.optionFound("writeFields")==true)
        {
            Info << "Info: Exporting step " << vectorID << endl;

            //word nameLAIi("LAI" + name(vectorID));

             volScalarField LAIi
             (
               IOobject
               (
                  "LAI",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
               ),
               mesh,
               dimensionedScalar("0", dimensionSet(0,0,0,0,0,0,0), 0.0),
               zeroGradientFvPatchScalarField::typeName
             );

            volScalarField divqrswi
            (
                IOobject
                (
                   "divqrsw",
                   runTime.timeName(),
                   mesh,
                   IOobject::NO_READ
                ),
                mesh,
                dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0),
                zeroGradientFvPatchScalarField::typeName
            );

            forAll(LAIi, cellI)
            {
                LAIi[cellI] = LAI[cellI];
                divqrswi[cellI] = divqrsw[cellI];
            }
            LAIi.correctBoundaryConditions();
            divqrswi.correctBoundaryConditions();
            LAIi.write();
            divqrswi.write();

            //runTime++;
            if (vectorID + 1 < sunPosVector_y.size())
            {
                runTime.setTime(sunPosVector_x[vectorID+1],runTime.timeIndex()+1);
            }
        }

        ////////////////////////////////////////////////////////////////////////
        // Status Info
        Info << "Solar ray direction " << vectorID
             << ", It took " << (std::clock()-tstartStep) / (double)CLOCKS_PER_SEC
             << " second(s)."<< endl;

        ////////////////////////////////////////////////////////////////////////
        // if (iter >= 2)
        //     break;

    }

    Info << "\nWriting fields: kcLAI boundary" << endl;
    kcLAIboundaryList.write();

    Info << "\nWriting fields: div qrsw" << endl;
    divqrswList.write();

    // // Status Info
    Info << "\nTotal time took: "
         << (std::clock()-tstart) / (double)CLOCKS_PER_SEC
         <<" second(s).\n"<< endl;

    Info<< "End\n" << endl;
    return 0;

}


// ************************************************************************* //
