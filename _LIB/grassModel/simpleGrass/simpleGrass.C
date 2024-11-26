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

\*---------------------------------------------------------------------------*/

#include "simpleGrass.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

#include "mixedFvPatchFields.H"
#include "mappedPatchBase.H"

#include "regionProperties.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace grass
    {
        defineTypeNameAndDebug(simpleGrass, 0);
        addToGrassRunTimeSelectionTables(simpleGrass);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::grass::simpleGrass::initialise()
{  
    nEvapSides_ = coeffs_.lookupOrDefault("nEvapSides", 1);
    Cd_ = coeffs_.lookupOrDefault("Cd", 0.2);
    beta_ = coeffs_.lookupOrDefault("beta", 0.78);
    betaLW_ = coeffs_.lookupOrDefault("betaLW", 0.83);
    LAI_ = coeffs_.lookupOrDefault("LAI", 2.0);
    l_ = coeffs_.lookupOrDefault("l", 0.1);
    albedoSoil_ = coeffs_.lookupOrDefault("albedoSoil", 0.2366);
    emissivitySoil_ = coeffs_.lookupOrDefault("emissivitySoil", 0.95);

    p_ = 101325;
    rhoa = 1.225;
    cpa = 1003.5;
    Ra = 287.042;
    Rv = 461.524;

    rs = coeffs_.lookupOrDefault("rs", 200);
    ra = coeffs_.lookupOrDefault("ra", -1); //if set to -1, it will be calculated below
    debug_ = coeffs_.lookupOrDefault("debug", 0);

    List<word> grassPatches(coeffs_.lookup("grassPatches"));  

    label count = 0;
    forAll(grassPatches, i)
    {
        grassPatchID = mesh_.boundaryMesh().findPatchID(grassPatches[i]);
        if (grassPatchID < 0)
        {
            FatalErrorInFunction
                << "Grass patch named " << grassPatches[i] << " not found." << nl
                << abort(FatalError);
        }
        else
        {
            selectedPatches_[count] = grassPatchID;
            count++;
        }
    }
    selectedPatches_.resize(count--);
}

Foam::scalarField Foam::grass::simpleGrass::calc_pvsat(const scalarField& T_)
{
     scalarField pvsat_ = exp( - 5.8002206e3/T_ // saturated vapor pressure pws - ASHRAE 1.2
            + 1.3914993
            - 4.8640239e-2*T_
            + 4.1764768e-5*pow(T_,2)
            - 1.4452093e-8*pow(T_,3)
            + 6.5459673*log(T_) );
     return pvsat_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::grass::simpleGrass::simpleGrass(const volScalarField& T)
:
    grassModel(typeName, T),
    Tg_
    (
        IOobject
        (
            "Tg",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        T*0
    ),
    Sw_
    (
        IOobject
        (
            "Sw",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    ),
    Sh_
    (
        IOobject
        (
            "Sh",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0)
    ),   
    Cf_
    (
        IOobject
        (
            "Cf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0,-1,0,0,0,0,0), 0.0)
    ),         
    selectedPatches_(mesh_.boundary().size(), -1)
{
    initialise();

    // read relaxation factor for Tg - aytac
    dictionary relaxationDict = mesh_.solutionDict().subDict("relaxationFactors");
    Tg_relax = relaxationDict.lookupOrDefault<scalar>("Tg", 0.5);

    dictionary residualControlDict = mesh_.solutionDict().subDict("SIMPLE").subDict("residualControl");
    Tg_residualControl = residualControlDict.lookupOrDefault<scalar>("Tg", 1e-8);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::grass::simpleGrass::~simpleGrass()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::grass::simpleGrass::read()
{
    if (grassModel::read())
    {
        // nothing to read

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::grass::simpleGrass::calculate
(
    const volScalarField& T_, 
    const volScalarField& w_,
    const volVectorField& U_
)
{    

    forAll(selectedPatches_, patchi)
    {
        label patchID = selectedPatches_[patchi];
        const fvPatch& thisPatch = mesh_.boundary()[patchID];

        const labelUList& patchInternalLabels = thisPatch.faceCells();
        scalarField Tg = thisPatch.patchInternalField(Tg_);
        scalarField Tc = thisPatch.patchInternalField(T_);
        scalarField wc = thisPatch.patchInternalField(w_);        
        scalarField deltaCoeffs = thisPatch.deltaCoeffs();

        scalarField Ts = thisPatch.lookupPatchField<volScalarField, scalar>("T");

        scalarField qr(thisPatch.size(), 0.0);
        scalarField qs(thisPatch.size(), 0.0);
        //-- Access vegetation region and populate radiation if vegetation exists,
        //otherwise use radiation from air region --//
        regionProperties rp(mesh_.time());
        const wordList vegNames(rp["vegetation"]);
        if (vegNames.size()>0)
        {
            const word& vegiRegion = "vegetation";
            const scalar mppVegDistance = 0; 
               
            const polyMesh& vegiMesh =
            	thisPatch.boundaryMesh().mesh().time().lookupObject<polyMesh>(vegiRegion);
     
            const word& nbrPatchName = thisPatch.name();

            const label patchi = vegiMesh.boundaryMesh().findPatchID(nbrPatchName);

            const fvPatch& vegiNbrPatch =
                refCast<const fvMesh>(vegiMesh).boundary()[patchi];

            // Get the coupling information from the mappedPatchBase
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(thisPatch.patch());

            const mappedPatchBase& mppVeg = mappedPatchBase(thisPatch.patch(), vegiRegion, mpp.mode(), thisPatch.name(), mppVegDistance);
            //const mappedPatchBase& mppVeg =
            //    refCast<const mappedPatchBase>(patch().patch(), vegiRegion, mpp.mode(), mpp.samplePatch(), mppVegDistance);
     
            qs = vegiNbrPatch.lookupPatchField<volScalarField, scalar>("qs");
            mppVeg.distribute(qs);

            qr = vegiNbrPatch.lookupPatchField<volScalarField, scalar>("qr");
            mppVeg.distribute(qr);           
        }
        else
        {
            qs = thisPatch.lookupPatchField<volScalarField, scalar>("qs");
            qr = thisPatch.lookupPatchField<volScalarField, scalar>("qr");
        }
        
        scalarField ra_(thisPatch.size(),ra);

        if(ra < 0) //calculate ra if it is not given in grassProperties
        {
            vectorField Uc = thisPatch.patchInternalField(U_);
            scalar C_ = 131.035; // proportionality factor
            scalarField magU = mag(Uc);
            magU = max(magU, SMALL);
            ra_ = C_*pow(l_/magU, 0.5); //aerodynamic resistance
        }    

        scalarField h_ch = (2*rhoa*cpa)/ra_; //convective heat transfer coefficient
        scalarField h_cm = (rhoa*Ra)/(p_*Rv*(ra_+rs)); //convective mass transfer coefficient

        //scalarField wsat = 0.621945*(pvsat/(p_-pvsat)); // saturated specific humidity, ASHRAE 1, eq.23
        scalarField pv = p_ * wc / (Ra/Rv + wc); // vapour pressure

        scalarField Qs_abs = qs*(1-exp(-beta_*LAI_)+albedoSoil_*exp(-beta_*LAI_));

        scalarField E(scalarField(Qs_abs.size(),0.0));

        ////calculate grass leaf temperature///////////
        label maxIter = 100;
        for (label i=1; i<=maxIter; i++)
        {
            if (i==1 && gMin(Tg) < SMALL)
            {
                Tg = Tc; //initialize if necessary
            }
            scalarField pvsat = calc_pvsat(Tg); //saturation vapour pressure
            E = pos(Qs_abs-SMALL)*nEvapSides_*h_cm*(pvsat-pv); //initialize transpiration rate [kg/(m2s)]
            //scalarField E = pos(Qs-SMALL)*nEvapSides_*rhoa*(wsat-wc)/(rs+ra);
            //no transpiration at night when Qs_abs is not >0

            scalar lambda = 2500000; // latent heat of vaporization of water J/kg
            scalarField Qlat = lambda*E*LAI_; //latent heat flux

            scalarField Qr2surrounding = qr;
            scalarField Qr2substrate = 6*(Ts-Tg); //thermal radiation between grass and surface - Malys et al 2014

            scalarField Tg_new = Tc + (Qr2surrounding + Qr2substrate + Qs_abs - Qlat)/ (h_ch*LAI_);

            scalar Tg_min = 250.0;
            scalar Tg_max = 400.0;
            bool boundTg = false;
            if((min(Tg_new) < Tg_min) or (max(Tg_new) > Tg_max))
            {
                boundTg = true;
                Tg_new = min
                (
                    Tg_new,
                    Tg_max
                );
                Tg_new = max
                (
                    Tg_new,
                    Tg_min
                );
            }
            reduce(boundTg, orOp<bool>());
            if(boundTg)
            {
                Info << "Warning, bounding Tg..." << endl;
                boundTg = false;
            }

            // info
            Info << " max grass leaf temp Tg=" << gMax(Tg_new)
                 << " K, iteration i="   << i << endl;

            // Check rel. L-infinity error
            scalar maxError = gMax(mag(Tg_new-Tg));
            scalar maxRelError = maxError/gMax(mag(Tg_new));

            // convergence check
            if (maxRelError < Tg_residualControl)
            {
                if(debug_)
                {
                    scalarField Qsen = h_ch*(Tc-Tg)*LAI_;
                    Info << " Qs_abs: " << gSum(thisPatch.magSf()*Qs_abs)/gSum(thisPatch.magSf()) << endl;
                    Info << " Qlat: " << gSum(thisPatch.magSf()*-Qlat)/gSum(thisPatch.magSf()) << endl;
                    Info << " Qsen: " << gSum(thisPatch.magSf()*Qsen)/gSum(thisPatch.magSf()) << endl;
                    Info << " Qr2surrounding: " << gSum(thisPatch.magSf()*Qr2surrounding)/gSum(thisPatch.magSf()) << endl;
                    Info << " Qr2substrate: " << gSum(thisPatch.magSf()*Qr2substrate)/gSum(thisPatch.magSf()) << endl;
                }
                break; 
            }
            else
            {
                Tg = (1-Tg_relax)*Tg+(Tg_relax)*Tg_new;// update leaf temp. 
            }                    
        }
        /////////////////////////////////////////

        ////update fields////////////////////////
        scalarField& Tg_Internal = Tg_.ref();
        volScalarField::Boundary& Tg_Bf = Tg_.boundaryFieldRef();
        scalarField& Tg_p = Tg_Bf[patchID];

        forAll(patchInternalLabels, i)
        {
            Tg_Internal[patchInternalLabels[i]] = Tg[i];
            Tg_p[i] = Tg[i];
            Sh_[patchInternalLabels[i]] = h_ch[i]*(Tg[i]-Tc[i])*(LAI_*(deltaCoeffs[i]/2));
            Sw_[patchInternalLabels[i]] = E[i]*(LAI_*(deltaCoeffs[i]/2));
            Cf_[patchInternalLabels[i]] = Cd_*(LAI_*(deltaCoeffs[i]/2));
        }
        /////////////////////////////////////////
    }
}

// return energy source term
Foam::tmp<Foam::volScalarField> Foam::grass::simpleGrass::Sh() const
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
            Sh_
        )
    );
}

// solve & return momentum source term (explicit)
Foam::tmp<Foam::volScalarField> Foam::grass::simpleGrass::Cf() const
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
            Cf_
        )
    );
}


// return specific humidity source term
Foam::tmp<Foam::volScalarField> Foam::grass::simpleGrass::Sw() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Sw",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            Sw_
        )
    );
}

// ************************************************************************* //

