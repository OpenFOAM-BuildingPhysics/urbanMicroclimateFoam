# urbanMicroclimateFoam

An open-source solver for coupled physical processes modeling urban microclimate based on OpenFOAM.

urbanMicroclimateFoam is a multi-region solver consisting of an air subdomain together with subdomains for porous urban building materials. A computational fluid dynamics (CFD) model solves the turbulent, convective air flow and heat and moisture transport in the air subdomain. A coupled heat and moisture (HAM) transport model solves the absorption, transport and storage of heat and moisture in the porous building materials. A radiation model determines the net longwave and shortwave radiative heat fluxes for each surface using a view factor approach.

* Daily variation of surface temperature in a street canyon:

<img src="https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/uploads/9fb2efac6b6827fa7604c3f58960093f/img1_out.gif"  width="60%">

* Daily variation of air temperature and wind speed at pedestrian height in Münsterhof, Zürich.

<img src="https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/uploads/95d6a8f84991e20b1223c307ca814fc9/img2b_out.gif"  width="60%"><img src="https://gitlab.ethz.ch/openfoam-cbp/solvers/urbanmicroclimatefoam/-/wikis/uploads/876fcd8fb7d6b18077f9ddeab44db013/img2c_out.gif"  width="60%">

The solver is tested for OpenFOAM v6.

More information: https://carmeliet.ethz.ch/research/downloads
