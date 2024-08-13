#include "CompressibleMomentumTransportModel.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "makeMomentumTransportModel.H"

#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define createBaseTurbulenceModel(                                             \
    Alpha, Rho, baseModel, BaseModel, Transport)                               \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef BaseModel<Transport>                                           \
            Transport##BaseModel;                                              \
        typedef RASModel<Transport##BaseModel>                                 \
            RAS##Transport##BaseModel;                                         \
        typedef LESModel<Transport##BaseModel>                                 \
            LES##Transport##BaseModel;                                         \
    }

createBaseTurbulenceModel
(
    geometricOneField,
    volScalarField,
    compressibleMomentumTransportModel,
    CompressibleMomentumTransportModel,
    fluidThermo
);
#define makeRASModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                        \
    (fluidThermoCompressibleMomentumTransportModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                        \
    (fluidThermoCompressibleMomentumTransportModel, LES, Type)

#include "porousrealizableKE.H"
makeRASModel(porousrealizableKE);
