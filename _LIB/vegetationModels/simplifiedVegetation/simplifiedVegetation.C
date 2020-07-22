/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

################
Vegetation model implemented by L. Manickathan, Empa, February 2017
################

\*---------------------------------------------------------------------------*/

#include "simplifiedVegetation.H"
#include "addToRunTimeSelectionTable.H"

#include "interpolationTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace vegetation
    {
        defineTypeNameAndDebug(simplifiedVegetation, 0);
        addToVegetationRunTimeSelectionTables(simplifiedVegetation);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vegetation::simplifiedVegetation::initialise()
{  
// do nothing
}

Foam::scalar Foam::vegetation::simplifiedVegetation::calc_evsat(double& T)
{
    // saturated vapor pressure pws - ASHRAE 1.2
    return exp( - 5.8002206e3/T
                + 1.3914993
                - 4.8640239e-2*T
                + 4.1764768e-5*pow(T,2)
                - 1.4452093e-8*pow(T,3)
                + 6.5459673*log(T) );
}

// calc saturated density of water vapour
Foam::scalar Foam::vegetation::simplifiedVegetation::calc_rhosat(double& T)
{
    return calc_evsat(T)/(461.5*T);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vegetation::simplifiedVegetation::simplifiedVegetation
(
    const volScalarField& T
):
    vegetationModel(typeName, T),
    a1_(coeffs_.lookup("a1")),
    a2_(coeffs_.lookup("a2")),
    a3_(coeffs_.lookup("a3")),
    C_(coeffs_.lookup("C")),
    D0_(coeffs_.lookup("D0")),
    nEvapSides_(coeffs_.lookup("nEvapSides")),
    H_(coeffs_.lookup("H")),
    kc_(coeffs_.lookup("kc")),
    l_(coeffs_.lookup("l")),
    Rg0_(coeffs_.lookup("Rg0")),
    Rl0_(coeffs_.lookup("Rl0")),
    rsMin_(coeffs_.lookup("rsMin")),
    TlMin_("TlMin", dimTemperature, SMALL),
    UMin_("UMin", dimVelocity, SMALL),
    Cd_(coeffs_.lookupOrDefault("Cd", 0.2)),
    rhoa_(dimensionedScalar(
        "rhoa",
        dimensionSet(1, -3, 0, 0, 0, 0, 0),
        1.225)),
    cpa_(dimensionedScalar(
        "cpa",
        dimensionSet(0,2,-2,-1,0,0,0),
        1003.5)),
    lambda_(dimensionedScalar(
        "lambda",
        dimensionSet(0,2,-2,0,0,0,0),
        2500000)),
    divqrsw
    (
        IOobject
        (
            "divqrsw",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    ),
    ev_
    (
        IOobject
        (
            "ev",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),
    evsat_
    (
        IOobject
        (
            "evsat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),
    LAD_
    (
        IOobject
        (
            "LAD",
            "0",//runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    Tl_
    (
        IOobject
        (
            "Tl",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        T*pos(LAD_)
    ),
    qsat_
    (
        IOobject
        (
            "qsat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    Qlat_
    (
        IOobject
        (
            "Qlat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    Qsen_
    (
        IOobject
        (
            "Qsen",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    ra_
    (
        IOobject
        (
            "ra",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,-1,1,0,0,0,0), 0.0)
    ),
    rs_
    (
        IOobject
        (
            "rs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,-1,1,0,0,0,0), 0.0)
    ),
    rhosat_
    (
        IOobject
        (
            "rhosat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-3,0,0,0,0,0), 0.0)
    ),
    Rg_
    (
        IOobject
        (
            "Rg",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
    ),
    Rn_
    (
        IOobject
        (
            "Rn",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),
    VPD_
    (
        IOobject
        (
            "VPD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    )
    {
        // Bounding parameters
		//bound(Tl_, TlMin_);
        initialise();
        Info << " Defined simplifiedVegetation model" << endl;
    }

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vegetation::simplifiedVegetation::~simplifiedVegetation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// solve radiation
void Foam::vegetation::simplifiedVegetation::radiation()
{
    const fvMesh& vegiMesh = mesh_.time().lookupObject<fvMesh>("vegetation");
    const label patchi = vegiMesh.boundaryMesh().findPatchID("air_to_vegetation");
    const fvPatch& vegiPatch = vegiMesh.boundary()[patchi];

    scalarField vegiPatchQr = vegiPatch.lookupPatchField<volScalarField, scalar>("qr");
    scalar integrateQr = gSum(vegiPatch.magSf() * vegiPatchQr);

    scalar vegiVolume = gSum(pos(LAD_.primitiveField() - 10*SMALL)*mesh_.V().field());

//    label timestepsInADay_ = divqrsw.size(); //readLabel(coeffs_.lookup("timestepsInADay"));
    Time& time = const_cast<Time&>(mesh_.time());
    // Read sunPosVector list
    interpolationTable<vector> sunPosVector
    (
        mesh_.time().rootPath()
        /mesh_.time().globalCaseName()
        /mesh_.time().constant()
        /"sunPosVector"
    );
    // look for the correct range    
    label lo = 0;
    label hi = 0;        
    for (label i = 0; i < sunPosVector.size(); ++i)
    {
        if (time.value() >= sunPosVector[i].first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }   
    }
    label timestep = lo;            
    /*
    label timestep = ceil( (time.value()/(86400/timestepsInADay_))-0.5 );
    //Info << ", 1 timestep: " << timestep;
    timestep = timestep % timestepsInADay_;
    //Info << ", 2 timestep: " << timestep << endl;
    */
    scalarList divqrswi =  divqrsw[timestep]; // [W/m3]

    // radiation density inside vegetation
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > 10*SMALL)
        {
            Rn_[cellI] = -divqrswi[cellI] + (integrateQr)/(vegiVolume); // [W/m3]
            Rg_[cellI] = -divqrswi[cellI]/LAD_[cellI]; // [W/m2]
        }
    }

    Rn_.correctBoundaryConditions();
    Rg_.correctBoundaryConditions();
    //Rn_.write();

}

// solve aerodynamic resistance
void Foam::vegetation::simplifiedVegetation::resistance(volScalarField& magU, volScalarField& T, volScalarField& q, volScalarField& Tl_)
{
    const double p_ = 101325;

    // Calculate magnitude of velocity and bounding above Umin
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > 10*SMALL)
        {
            //Aerodynamic resistance
            ra_[cellI] = C_.value()*pow(l_.value()/magU[cellI], 0.5);

            // Calculate vapor pressure of air
            ev_[cellI] = p_*q[cellI]/(0.621945+q[cellI]);

            // Calculate sat. vapor pressure of air
            evsat_[cellI] = calc_evsat(T[cellI]);

            // Vapor pressure deficit - kPa
            // VPD_[cellI] = (calc_evsat(T[cellI]) - (q[cellI]*rhoa_.value()*T[cellI]*461.5))/1000.0; // kPa
            //VPD_[cellI] = ev_[cellI] - evsat_[cellI];
            VPD_[cellI] = evsat_[cellI] - ev_[cellI];


            // Stomatal resistance - type 1
            // rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]); // type 1
            //rs_[cellI] = rsMin_.value()*(31.0 + Rn_[cellI])*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rn_[cellI]);
            // rs_[cellI] = rsMin_.value()*(31.0 + Rg_[cellI].component(2))*(1.0+0.016*pow((T[cellI]-16.4-273.15),2))/(6.7+Rg_[cellI].component(2));


            // Stomatal resistance - type 2
            // rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2)); // type 2
            //if ((VPD_[cellI]/1000.0) < D0_.value())
            //    rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()));
            //else
            //    rs_[cellI] = rsMin_.value()*((a1_.value() + Rg0_.value())/(a2_.value() + Rg0_.value()))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2));
/*            if ((VPD_[cellI]/1000.0) < D0_.value())
                rs_[cellI] = rsMin_.value();//((a1_.value() + mag(Rg_[cellI]))/(a2_.value() + mag(Rg_[cellI])));
            else
                rs_[cellI] = rsMin_.value();//((a1_.value() + mag(Rg_[cellI]))/(a2_.value() + mag(Rg_[cellI])))*(1.0 + a3_.value()*pow(VPD_[cellI]/1000.0-D0_.value(),2));
*/
            scalar f1 = 7.119*exp(-0.05004*Rg_[cellI]) + 0.6174*exp(0.0006336*Rg_[cellI]);   
            scalar f2 = 1;
            if(VPD_[cellI] < 0)
            {                   
                f2 = 0.4372;
            }
            else
            {
                f2 = 0.4372*pow((VPD_[cellI]+1),0.204);
            }                     
            rs_[cellI] = rsMin_.value()*f1*f2;
        }
    }
    ev_.correctBoundaryConditions();
    evsat_.correctBoundaryConditions();
    VPD_.correctBoundaryConditions();
    ra_.correctBoundaryConditions();
    rs_.correctBoundaryConditions();
}

// solve vegetation model
void Foam::vegetation::simplifiedVegetation::calculate(volVectorField& U, volScalarField& T, volScalarField& q)
{
    // solve radiation within vegetation
    radiation();

    const double p_ = 101325;

    // Magnitude of velocity
    volScalarField magU("magU", mag(U));
    // Bounding velocity
    bound(magU, UMin_);

    volScalarField new_Tl("new_Tl", Tl_);

    // info
    Info << "    max leaf temp tl=" << max(T.internalField())
         << "k, iteration i=0" << endl;

//    const fvMesh& vegiMesh = mesh_.time().lookupObject<fvMesh>("vegetation");
//    const label patchi = vegiMesh.boundaryMesh().findPatchID("air_to_vegetation");
//    const fvPatch& vegiPatch = vegiMesh.boundary()[patchi];
//    scalarField vegiPatchQs = vegiPatch.lookupPatchField<volScalarField, scalar>("qs");
//    scalar integrateQs = gSum(vegiPatch.magSf() * vegiPatchQs);

    scalar maxError, maxRelError;
    int i;

    // solve leaf temperature, iteratively.
    int maxIter = 500;
    for (i=1; i<=maxIter; i++)
    {
        // Solve aerodynamc, stomatal resistance
        resistance(magU, T, q, new_Tl);

        forAll(LAD_, cellI)
        {
            if (LAD_[cellI] > 10*SMALL)
            {
                // Initial leaf temperature
                if (i==1)
                    Tl_[cellI];// = T[cellI];//*0. + 300.;//T[cellI];

                // Calculate saturated density, specific humidity
                rhosat_[cellI] = calc_rhosat(Tl_[cellI]);
                evsat_[cellI] = calc_evsat(Tl_[cellI]);
                qsat_[cellI] = 0.621945*(evsat_[cellI]/(p_-evsat_[cellI])); // ASHRAE 1, eq.23

                // Calculate transpiration rate
                //no transpiration at night when solar radiation is not >0
             	E_[cellI] = pos(Rg_[cellI]-SMALL)*nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]);

                // Calculate latent heat flux
                Qlat_[cellI] = lambda_.value()*E_[cellI];

                // Calculate new leaf temperature
                new_Tl[cellI] = T[cellI] + (Rn_[cellI] - Qlat_[cellI])*(ra_[cellI]/(2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]));

            }
        }

        // info
        Info << "    max leaf temp tl=" << gMax(new_Tl.internalField())
             << " K, iteration i="   << i << endl;

        // Check rel. L-infinity error
        maxError = gMax(mag(new_Tl.primitiveField()-Tl_.primitiveField()));
        maxRelError = maxError/gMax(mag(new_Tl.primitiveField()));

        // read relaxation factor for Tl - aytac
        dictionary relaxationDict = mesh_.solutionDict().subDict("relaxationFactors");
        scalar Tl_relax = relaxationDict.lookupOrDefault<scalar>("Tl", 0.5);

        // update leaf temp.
        forAll(Tl_, cellI)
            Tl_[cellI] = (1-Tl_relax)*Tl_[cellI]+(Tl_relax)*new_Tl[cellI];

         // convergence check
         if (maxRelError < 1e-8)
             break;
    }
    Tl_.correctBoundaryConditions();

    // Iteration info
    Info << "Vegetation model:  Solving for Tl, Final residual = " << maxError
         << ", Final relative residual = " << maxRelError
         << ", No Iterations " << i << endl;

    Info << "temperature parameters: max Tl = " << gMax(Tl_)
         << ", min T = " << gMin(T) << ", max T = " << gMax(T) << endl;

    Info << "resistances: max rs = " << gMax(rs_)
         << ", max ra = " << gMax(ra_) << endl;

    // Final: Solve aerodynamc, stomatal resistance
    resistance(magU, T, q, Tl_);

    // Final: Update sensible and latent heat flux
    forAll(LAD_, cellI)
    {
        if (LAD_[cellI] > 10*SMALL)
        {
            // Calculate saturated density, specific humidity
            rhosat_[cellI] = calc_rhosat(Tl_[cellI]);
            evsat_[cellI] = calc_evsat(Tl_[cellI]);
            qsat_[cellI] = 0.621945*(evsat_[cellI]/(p_-evsat_[cellI])); // ASHRAE 1, eq.23

            // Calculate transpiration rate
            E_[cellI] = pos(Rg_[cellI]-SMALL)*nEvapSides_.value()*LAD_[cellI]*rhoa_.value()*(qsat_[cellI]-q[cellI])/(ra_[cellI]+rs_[cellI]); // todo: implement switch for double or single side

            // Calculate latent heat flux
            Qlat_[cellI] = lambda_.value()*E_[cellI];

            // Calculate sensible heat flux
            Qsen_[cellI] = 2.0*rhoa_.value()*cpa_.value()*LAD_[cellI]*(Tl_[cellI]-T[cellI])/ra_[cellI];
        }
    }
    rhosat_.correctBoundaryConditions();
    qsat_.correctBoundaryConditions();
    E_.correctBoundaryConditions();
    Qlat_.correctBoundaryConditions();
    Qsen_.correctBoundaryConditions();

    // Iteration info
    // Info << "              Vegetation model:  max. Rn = " << max(mag(Rn_))
    //      << "; max. Qlat = " << max(mag(Qlat_))
    //      << "; max. Qsen = " << max(mag(Qsen_))
    //      << "; error: max. Esum = " << max(mag(Rn_.internalField() - Qsen_.internalField()- Qlat_.internalField())) << endl;

}

// -----------------------------------------------------------------------------

// return energy source term
Foam::tmp<Foam::volScalarField> Foam::vegetation::simplifiedVegetation::Sh() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            Qsen_
        )
    );
}

// solve & return momentum source term (explicit)
Foam::tmp<Foam::volScalarField> Foam::vegetation::simplifiedVegetation::Cf() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Cf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            Cd_*LAD_
        )
    );
}


// return specific humidity source term
Foam::tmp<Foam::volScalarField> Foam::vegetation::simplifiedVegetation::Sq() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Sq",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            E_
        )
    );
}


// -----------------------------------------------------------------------------

bool Foam::vegetation::simplifiedVegetation::read()
{
    if (vegetationModel::read())
    {
        // nothing to read

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

