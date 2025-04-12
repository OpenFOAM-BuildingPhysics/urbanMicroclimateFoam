#include "buildingMaterialModel.H"

labelList cellType(mesh.nCells(), -1);

forAll(Materials, MaterialsI)
{
    const dictionary& dict = Materials[MaterialsI];

    const word cellZoneName(dict.lookup("name"));
    const word cellZoneModel(dict.lookup("buildingMaterialModel"));
    const scalar rho_(readScalar(dict.lookup("rho")));
    const scalar cap_(readScalar(dict.lookup("cap")));
    const scalar lambda1_(readScalar(dict.lookup("lambda1")));
    const scalar lambda2_(readScalar(dict.lookup("lambda2")));

    //create the buildingMaterial model
    autoPtr<buildingMaterialModel> buildingMaterial
    (
       buildingMaterialModel::New("buildingMaterial", dict, cellZoneModel)
    );
    
    label cellZoneID = mesh.cellZones().findZoneID(cellZoneName);
    if (cellZoneID == -1)
    {
        FatalErrorInFunction
        << "Something is wrong, cannot find at least one of the necessary material cellZones!"
        << exit(FatalError);
    }
    const labelList& cells = mesh.cellZones()[cellZoneID];
//    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(cells, cellsI)
    {
        label celli = cells[cellsI];
        if (cellType[celli] > -1)
        {
            FatalErrorInFunction
            << "In solid region " << solidRegions[i].name() << ", a cell is assigned to more than one material cellZone!"
            << exit(FatalError);
        }
        else
        {
            cellType[celli] = MaterialsI;
        }

        buildingMaterial->update_w_C_cell(pc,ws,Crel,celli);
        buildingMaterial->update_Krel_cell(pc,ws,Krel,celli);
        buildingMaterial->update_Kv_cell(pc,ws,Ts,K_v,celli);
        rho_m.ref()[celli] = rho_;
        cap_m.ref()[celli] = cap_;
        lambda_m.ref()[celli] = lambda1_ + lambda2_*ws[celli];
        buildingMaterial->update_Kpt_cell(pc,ws,Ts,K_pt,celli);
    }
}
if (min(cellType) == -1)
{
    FatalErrorInFunction
    << "In solid region " << solidRegions[i].name() << ", not all cells are assigned to a material cellZone!"
    << exit(FatalError);
}

dimensionedScalar minCrel_("minCrel_", dimensionSet(0, -2, 2, 0, 0, 0, 0), minCrel);
Crel = max(Crel, minCrel_);

ws.correctBoundaryConditions();
Crel.correctBoundaryConditions();
Krel.correctBoundaryConditions();
K_v.correctBoundaryConditions();
rho_m.correctBoundaryConditions();
cap_m.correctBoundaryConditions();
lambda_m.correctBoundaryConditions();
K_pt.correctBoundaryConditions();
