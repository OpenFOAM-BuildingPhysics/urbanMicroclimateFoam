pc=pc_old;//pc.prevIter(); 
volScalarField::Boundary& pcBf = pc.boundaryFieldRef();
forAll(pcBf, patchi)
{ 
    forAll(pcBf[patchi], facei)
    {
        pcBf[patchi][facei] = pc_old.boundaryField()[patchi][facei];
    }
}
Ts=Ts_old;//.prevIter();
volScalarField::Boundary& TsBf = Ts.boundaryFieldRef();
forAll(TsBf, patchi)
{ 
    forAll(TsBf[patchi], facei)
    {
        TsBf[patchi][facei] = Ts_old.boundaryField()[patchi][facei];
    }
}
pcBf.updateCoeffs();
TsBf.updateCoeffs();
