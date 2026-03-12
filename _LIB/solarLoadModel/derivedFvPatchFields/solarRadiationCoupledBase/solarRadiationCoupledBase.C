/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "solarRadiationCoupledBase.H"
#include "volFields.H"
#include "fieldMapper.H"
#include "mappedFvPatchBaseBase.H"
#include "radiationModel.H"
#include "opaqueSolid.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solarRadiationCoupledBase, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::solarRadiationCoupledBase::albedoMethodType,
        2
    >::names[] =
    {
        "solidRadiation",
        "lookup"
    };
}


const Foam::NamedEnum<Foam::solarRadiationCoupledBase::albedoMethodType, 2>
    Foam::solarRadiationCoupledBase::albedoMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarRadiationCoupledBase::solarRadiationCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const scalarField& albedo
)
:
    patch_(patch),
    method_(albedoMethodTypeNames_[calculationType]),
    albedo_(albedo)
{}

Foam::solarRadiationCoupledBase::solarRadiationCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const scalarField& albedo,
    const fieldMapper& mapper
)
:
    patch_(patch),
    method_(albedoMethodTypeNames_[calculationType]),
    albedo_(mapper(albedo))
{}

Foam::solarRadiationCoupledBase::solarRadiationCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(albedoMethodTypeNames_.read(dict.lookup("albedoMode")))
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            // if (!isA<mappedPatchBase>(patch_.patch()))
            // {
            //     FatalIOErrorInFunction
            //     (
            //         dict
            //     )   << "\n    patch type '" << patch_.type()
            //         << "' not type '" << mappedPatchBase::typeName << "'"
            //         << "\n    for patch " << patch_.name()
            //         << exit(FatalIOError);
            // }

            albedo_ = scalarField(patch_.size(), 0.0);
        }
        break;

        case LOOKUP:
        {
            if (!dict.found("albedo"))
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "\n    albedo key does not exist for patch "
                    << patch_.name()
                    << exit(FatalIOError);
            }
            else
            {
                albedo_ = scalarField("albedo", unitFraction,dict, patch_.size());
            }
        }
        break;
    }
}

// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::solarRadiationCoupledBase::~solarRadiationCoupledBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::solarRadiationCoupledBase::albedo() const
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            // Get the mapper and the neighbouring mesh and patch
            const mappedFvPatchBaseBase& mapper =
                mappedFvPatchBaseBase::getMap(patch_);
            const fvMesh& nbrMesh = mapper.nbrMesh();
            const fvPatch& nbrPatch = mapper.nbrFvPatch();

            const radiationModels::opaqueSolid& radiation =
                nbrMesh.lookupObject<radiationModels::opaqueSolid>
                (
                    "radiationProperties"
                );

            // NOTE: for an opaqueSolid the absorptionEmission model returns the
            // emissivity of the surface rather than the emission coefficient
            // and the input specification MUST correspond to this.
            return
                mapper.fromNeighbour
            (
                radiation.absorptionEmission().e()().boundaryField()
                [
                    nbrPatch.index()
                ]
            );
        }
        break;

        case LOOKUP:
        {
            // return local value
            return albedo_;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << method_ << endl
                << "Please set 'albedo' to one of "
                << albedoMethodTypeNames_.toc()
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


void Foam::solarRadiationCoupledBase::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    const solarRadiationCoupledBase& mrptf =
        refCast<const solarRadiationCoupledBase>(ptf);

    mapper(albedo_, mrptf.albedo_);
}


void Foam::solarRadiationCoupledBase::reset
(
    const fvPatchScalarField& ptf
)
{
    const solarRadiationCoupledBase& mrptf =
        refCast<const solarRadiationCoupledBase>(ptf);

    albedo_.reset(mrptf.albedo_);
}


void Foam::solarRadiationCoupledBase::write(Ostream& os) const
{
    writeEntry(os, "albedoMode", albedoMethodTypeNames_[method_]);
    writeEntry(os, "albedo", albedo_);
}


// ************************************************************************* //
