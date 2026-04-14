/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Realizable k-epsilon turbulence model for compressible flows with added
    vegetation source terms.

SourceFiles
    porousrealizableKE.C

\*---------------------------------------------------------------------------*/

#ifndef porousrealizableKE_H
#define porousrealizableKE_H

#include "RASModel.H"
#include "eddyViscosity.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
class porousrealizableKE
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
protected:

    dimensionedScalar A0_;
    dimensionedScalar C2_;
    dimensionedScalar C4_;
    dimensionedScalar C5_;
    dimensionedScalar sigmak_;
    dimensionedScalar sigmaEps_;
    dimensionedScalar betaP_;
    dimensionedScalar betaD_;

    scalar Cd_;

    volScalarField k_;
    volScalarField epsilon_;
    volScalarField LAD_;

    void boundEpsilon();

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

    TypeName("porousrealizableKE");

    porousrealizableKE
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity,
        const word& type = typeName
    );

    virtual ~porousrealizableKE()
    {}

    virtual bool read();

    tmp<volScalarField> DkEff() const
    {
        return volScalarField::New
        (
            "DkEff",
            (this->nut_/sigmak_ + this->nu())
        );
    }

    tmp<volScalarField> DepsilonEff() const
    {
        return volScalarField::New
        (
            "DepsilonEff",
            (this->nut_/sigmaEps_ + this->nu())
        );
    }

    virtual tmp<volScalarField> k() const
    {
        return k_;
    }

    virtual tmp<volScalarField> epsilon() const
    {
        return epsilon_;
    }

    virtual tmp<volScalarField> omega() const
    {
        return volScalarField::New
        (
            "omega",
            epsilon_/(0.09*k_)
        );
    }

    virtual void correct();
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "porousrealizableKE.C"
#endif

#endif

