/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::RASModels::porousrealizableKE

Description
    Realizable k-epsilon turbulence model for incompressible and compressible
    flows.

    References:
    \verbatim
        Shih, T. H., Liou, W. W., Shabbir, A., Yang, Z., & Zhu, J. (1994).
        A new k-epsilon eddy viscosity model for high Reynolds number
        turbulent flows: Model development and validation.
        NASA STI/Recon Technical Report N, 95, 11442.

        Shih, T. H., Liou, W. W., Shabbir, A., Yang, Z., & Zhu, J. (1995).
        A New k-epsilon Eddy Viscosity Model for High Reynolds Number
        Turbulent Flows.
        Computers and Fluids, 24(3), 227-238.
    \endverbatim

    The default model coefficients are
    \verbatim
        porousrealizableKECoeffs
        {
            A0          4.0;
            C2          1.9;
            sigmak      1.0;
            sigmaEps    1.2;
        }
    \endverbatim

SourceFiles
    porousrealizableKE.C

\*---------------------------------------------------------------------------*/

#ifndef porousrealizableKE_H
#define porousrealizableKE_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class porousrealizableKE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class porousrealizableKE
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar A0_;
            dimensionedScalar C2_;
            dimensionedScalar C4_;
            dimensionedScalar C5_;  
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar betaP_;
            dimensionedScalar betaD_;  

            scalar Cd_;

        // Fields

            volScalarField k_;
            volScalarField epsilon_;
            volScalarField LAD_;


   // Protected Member Functions

        tmp<volScalarField> rCmu
        (
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS
        );

        virtual void correctNut
        (
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS
        );

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("porousrealizableKE");

    // Constructors

        //- Construct from components
        porousrealizableKE
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );


    //- Destructor
    virtual ~porousrealizableKE()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    (this->nut_/sigmak_ + this->nu())
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    (this->nut_/sigmaEps_ + this->nu())
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "porousrealizableKE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
