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
    Foam::vegetation

Description
    Namespace for vegetation modelling

Class
    Foam::vegetation::vegetationModel

Description
    Top level model for vegetation modelling

SourceFiles
    vegetationModel.C
    vegetationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef vegetationModel_H
#define vegetationModel_H

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

namespace vegetation
{

/*---------------------------------------------------------------------------*\
                       Class vegetationModel Declaration
\*---------------------------------------------------------------------------*/

class vegetationModel
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
        Switch vegetation_;

        //- Radiation model dictionary
        dictionary coeffs_;

        //- Radiation solver frequency - number flow solver iterations per
        //  vegetation solver iteration
        label solverFreq_;

        //- Flag to enable vegetation model to be evaluated on first iteration
        bool firstIter_;


private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        vegetationModel(const vegetationModel&);

        //- Disallow default bitwise assignment
        void operator=(const vegetationModel&);


public:

    //- Runtime type information
    TypeName("vegetationModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            vegetationModel,
            T,
            (
                const volScalarField& T
            ),
            (T)
        );


    // Constructors

        //- Null constructor
        vegetationModel(const volScalarField& T);

        //- Construct from components
        vegetationModel(const word& type, const volScalarField& T);


    // Selectors

        //- Return a reference to the selected vegetation model
        static autoPtr<vegetationModel> New(const volScalarField& T);


    //- Destructor
    virtual ~vegetationModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct
            (
                volVectorField& U_,
                volScalarField& T_, 
                volScalarField& q_
            );

            //- Solve vegetation equation(s)
            virtual void calculate
            (
                volVectorField& U_,
                volScalarField& T_, 
                volScalarField& q_
            ) = 0;

            //- Read vegetationProperties dictionary
            virtual bool read() = 0;

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Sh() const = 0;

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Cf() const = 0;

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Sq() const = 0;

        // Access

            //- Radiation model on/off flag
            const Switch vegetation() const
            {
                return vegetation_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToVegetationRunTimeSelectionTables(model)                           \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        vegetationModel,                                                       \
        model,                                                                \
        T                                                                     \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vegetation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
