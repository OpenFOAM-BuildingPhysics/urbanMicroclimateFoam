scalar maxChangews = gMax(mag(ws.primitiveField()-ws_n.primitiveField())); 
scalar maxChangeTs = gMax(mag(Ts.primitiveField()-Ts_n.primitiveField()));

if(maxChangews < PicardTolerancews && maxChangeTs < PicardToleranceTs && nInternalIter>=2) //force at least 1 internal iteration
{
    pc_n = pc; //pc_n.boundaryFieldRef().updateCoeffs();
    ws_n = ws; //ws_n.boundaryFieldRef().updateCoeffs();
    Ts_n = Ts; //Ts_n.boundaryFieldRef().updateCoeffs();
    Info << "Total internal iterations: " << nInternalIter << endl;
    timeStepDecrease = false;
    break;
}
else if (nInternalIter == nInternalIterMax) //nonlinear iteration reached maximum
{
    Info
       << "Nonlinear iteration didn't converge !"
       << endl
       << "maxChangews: " << maxChangews << endl
       << "maxChangeT: " << maxChangeTs << endl
       << "Total internal iterations: " << nInternalIter << endl;
    timeStepDecrease = true;
    #include "setSolidRegionDeltaT.H"                    
    #include "revertValues.H"
}   
else //not converged nor reached the maximum iteration yet, continue
{
    pc_n = pc; //pc_n.boundaryFieldRef().updateCoeffs();
    ws_n = ws; //ws_n.boundaryFieldRef().updateCoeffs();
    Ts_n = Ts; //Ts_n.boundaryFieldRef().updateCoeffs();
    Info << "maxChangews: " << maxChangews << ", maxChangeTs: " << maxChangeTs << endl;
}       
