# urbanMicroclimateFoam

An open-source solver for coupled physical processes modeling urban microclimate based on OpenFOAM.

urbanMicroclimateFoam is a multi-region solver consisting of an air subdomain together with subdomains for porous urban building materials. A computational fluid dynamics (CFD) model solves the turbulent, convective air flow and heat and moisture transport in the air subdomain. A coupled heat and moisture (HAM) transport model solves the absorption, transport and storage of heat and moisture in the porous building materials. A radiation model determines the net longwave and shortwave radiative heat fluxes for each surface using a view factor approach.

| Daily variation of surface temperature in a street canyon | Daily variation of air temperature and wind speed at pedestrian height in Münsterhof, Zürich. |
|:---:|:---:|
| <img src="https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/uploads/9fb2efac6b6827fa7604c3f58960093f/img1_out.gif"  width="75%"> | <img src="https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/uploads/95d6a8f84991e20b1223c307ca814fc9/img2b_out.gif"  width="45%"> &nbsp; <img src="https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/uploads/876fcd8fb7d6b18077f9ddeab44db013/img2c_out.gif"  width="45%"> |

The solver is tested for the following OpenFOAM versions:

* OpenFOAM-org (OpenFOAM Foundation) v6, v7, v8

### Usage

You can compile the solver for a specific OpenFOAM version by checking out the commit with corresponding tag. For example, for OpenFOAM v8:

	git clone https://github.com/OpenFOAM-BuildingPhysics/urbanMicroclimateFoam.git
	cd urbanmicroclimatefoam
	git checkout tags/of-org_v8.0
	./Allwmake

See the list of tags for different versions [here](https://github.com/OpenFOAM-BuildingPhysics/urbanMicroclimateFoam/tags)
	
### Documentation

Read [urbanMicroclimateFoam wiki](https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/home) for documentation.

### Tutorial cases

Tutorial case for the solver can be found [here](https://github.com/OpenFOAM-BuildingPhysics/urbanMicroclimateFoam-tutorials).

More information at the Chair of Building Physics: https://carmeliet.ethz.ch
