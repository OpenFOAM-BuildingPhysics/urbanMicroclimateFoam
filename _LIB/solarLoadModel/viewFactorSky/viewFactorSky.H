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

Class
    Foam::radiationModels::viewFactorSky

Description
    View factor radiation model.
    Sky temperature is used for non-wall patches, e.g. side and top boundaries.
    Written by Aytac Kubilay based on the original viewFactor model.
    
    The system solved is: C q = b
    where:
            Cij  = deltaij/Ej - (1/Ej - 1)Fij
            q    = heat flux
            b    = A eb - Ho
    and:
            eb   = sigma*T^4
            Ej   = emissivity
            Aij  = deltaij - Fij
            Fij  = view factor matrix


SourceFiles
    viewFactorSky.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelviewFactorSky_H
#define radiationModelviewFactorSky_H

#include "radiationModel.H"
#include "singleCellFvMesh.H"
#include "scalarMatrices.H"
#include "globalIndex.H"
#include "scalarListIOList.H"
#include "polyDistributionMap.H"
#include "volFields.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                           Class viewFactorSky Declaration
\*---------------------------------------------------------------------------*/

class viewFactorSky
:
    public radiationModel
{
    // Private Data

        //- Agglomeration List
        labelListIOList finalAgglom_;

        //- Map distributed
        autoPtr<distributionMap> map_;

        //- Coarse mesh
        singleCellFvMesh coarseMesh_;

        //- Net radiative heat flux [W/m^2]
        volScalarField qr_;

        //- Local sparse view-factor column ids
        labelListList globalFaceFaces_;

        //- Local sparse view-factor coefficients
        scalarListList viewFactors_;

        //- Local row to global coarse-face ids
        labelList localGlobalIds_;

        //- Previous sparse solve result, used as Krylov initial guess
        mutable scalarField qPrev_;

        //- Whether qPrev_ contains a valid previous solution
        mutable bool qPrevValid_;

        //- Selected patches
        labelList selectedPatches_;

        //- Total global coarse faces
        label totalNCoarseFaces_;

        //- Total local coarse faces
        label nLocalCoarseFaces_;

        //- Constant emissivity
        bool constEmissivity_;

        //- List of grass patches
        hashedWordList grassPatches;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Assemble locally-owned global vector entries across processors
        void assembleGlobal(scalarField& x) const;

        //- Dot product over locally-owned rows
        scalar localDot
        (
            const scalarField& a,
            const scalarField& b
        ) const;

        //- Sparse thermal view-factor matrix-vector multiply
        void multiply
        (
            const scalarField& x,
            const scalarField& E,
            scalarField& Ax
        ) const;

        //- Solve the sparse thermal view-factor system
        scalarField solveViewFactorSystem
        (
            const scalarField& b,
            const scalarField& E
        ) const;

        //- Compare sparse solution with a rank-0 dense solve for validation
        void referenceCheck
        (
            const scalarField& q,
            const scalarField& b,
            const scalarField& E
        ) const;


public:

    //- Runtime type information
    TypeName("viewFactorSky");


    // Constructors

        //- Construct from components
        viewFactorSky(const volScalarField& T);

        //- Construct from components
        viewFactorSky(const dictionary& dict, const volScalarField& T);

        //- Disallow default bitwise copy construction
        viewFactorSky(const viewFactorSky&) = delete;


    //- Destructor
    virtual ~viewFactorSky();


    // Member Functions

        // Edit

            //- Solve system of equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const;


    // Access

        //- Const access to total radiative heat flux field
        inline const volScalarField& qr() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const viewFactorSky&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "viewFactorSkyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
