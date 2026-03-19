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

Namespace
    Foam::radiation

Description
    Namespace for radiation modelling

Class
    Foam::radiation::radiationModel

Description
    Top level model for radiation modelling

SourceFiles
    radiationModel.C
    radiationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef solarLoadModel_H
#define solarLoadModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "DimensionedField.H"
#include "fvMatricesFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;
class fvMesh;

namespace solarLoad
{

// Forward declaration of classes
class solarLoadAbsorptionEmissionModel;
class solarLoadScatterModel;

/*---------------------------------------------------------------------------*\
                       Class solarLoadModel Declaration
\*---------------------------------------------------------------------------*/

class solarLoadModel
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

        //- solarLoad model on/off flag
        Switch solarLoad_;

        //- solarLoad model dictionary
        dictionary coeffs_;

        //- solarLoad solver frequency - number flow solver iterations per
        //  solarLoad solver iteration
        label solverFreq_;

        //- Flag to enable solarLoad model to be evaluated on first iteration
        bool firstIter_;


        // References to the solarLoad sub-models

            //- Absorption/emission model
            autoPtr<solarLoadAbsorptionEmissionModel> solarLoadAbsorptionEmission_;

            //- Scatter model
            autoPtr<solarLoadScatterModel> solarLoadScatter_;


private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        solarLoadModel(const solarLoadModel&);

        //- Disallow default bitwise assignment
        void operator=(const solarLoadModel&);


public:

    //- Runtime type information
    TypeName("solarLoadModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solarLoadModel,
            T,
            (
                const volScalarField& T
            ),
            (T)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            solarLoadModel,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& T
            ),
            (dict, T)
        );


    // Constructors

        //- Null constructor
        solarLoadModel(const volScalarField& T);

        //- Construct from components
        solarLoadModel(const word& type, const volScalarField& T);

        //- Construct from components
        solarLoadModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& T
        );


    // Selectors

        //- Return a reference to the selected solarLoad model
        static autoPtr<solarLoadModel> New(const volScalarField& T);

        //- Return a reference to the selected solarLoad model
        static autoPtr<solarLoadModel> New
        (
            const dictionary& dict,
            const volScalarField& T
        );


    //- Destructor
    virtual ~solarLoadModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve solarLoad equation(s)
            virtual void calculate() = 0;

            //- Read solarLoadProperties dictionary
            virtual bool read() = 0;


        // Access

            //- solarLoad model on/off flag
            const Switch solarLoad() const
            {
                return solarLoad_;
            }

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const = 0;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const = 0;

            //- Energy source term
            virtual tmp<fvScalarMatrix> Sh
            (
                const basicThermo& thermo,
                const volScalarField& he
            ) const;

            //- Temperature source term
            virtual tmp<fvScalarMatrix> ST
            (
                const dimensionedScalar& rhoCp,
                volScalarField& T
            ) const;

            //- Access to absorptionEmissionModel
            const solarLoadAbsorptionEmissionModel& solarLoadAbsorptionEmission() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToSolarLoadRunTimeSelectionTables(model)                           \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        solarLoadModel,                                                       \
        model,                                                                \
        dictionary                                                            \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        solarLoadModel,                                                       \
        model,                                                                \
        T                                                                     \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solarLoad
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
