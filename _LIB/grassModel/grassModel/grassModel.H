/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

Namespace
    Foam::grass

Description
    Namespace for grass modelling

Class
    Foam::grass::grassModel

Description
    Top level model for grass modelling

SourceFiles
    grassModel.C
    grassModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef grassModel_H
#define grassModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "fvMatricesFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidThermo;
class fvMesh;

namespace grass
{

/*---------------------------------------------------------------------------*\
                       Class grassModel Declaration
\*---------------------------------------------------------------------------*/

class grassModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the temperature field
        const volScalarField& T_;

        //- Radiation model on/off flag
        Switch grass_;

        //- Radiation model dictionary
        dictionary coeffs_;

        //- Radiation solver frequency - number flow solver iterations per
        //  grass solver iteration
        label solverFreq_;

        //- Flag to enable grass model to be evaluated on first iteration
        bool firstIter_;


private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        grassModel(const grassModel&);

        //- Disallow default bitwise assignment
        void operator=(const grassModel&);


public:

    //- Runtime type information
    TypeName("grassModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            grassModel,
            T,
            (
                const volScalarField& T
            ),
            (T)
        );


    // Constructors

        //- Null constructor
        grassModel(const volScalarField& T);

        //- Construct from components
        grassModel(const word& type, const volScalarField& T);


    // Selectors

        //- Return a reference to the selected grass model
        static autoPtr<grassModel> New(const volScalarField& T);


    //- Destructor
    virtual ~grassModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct
            (
                const volScalarField& T_, 
                const volScalarField& w_,
                const volVectorField& U_
            );

            //- Solve grass equation(s)
            virtual void calculate
            (
                const volScalarField& T_, 
                const volScalarField& w_,
                const volVectorField& U_
            ) = 0;

            //- Read grassProperties dictionary
            virtual bool read() = 0;
            
            //- return leaf drag coefficient
            virtual tmp<volScalarField> Sh() const = 0;

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Cf() const = 0;

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Sw() const = 0;            

        // Access

            //- Radiation model on/off flag
            const Switch grass() const
            {
                return grass_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToGrassRunTimeSelectionTables(model)                           \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        grassModel,                                                       \
        model,                                                                \
        T                                                                     \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace grass
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
