#include "cantera/thermo/ThermoPhase.h"

#define ARRAY_FUNC(PREFIX, CLASS_NAME, FUNC_NAME) \
    void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, double* data) \
    { object->FUNC_NAME(data); }

#define THERMO_1D(FUNC_NAME) ARRAY_FUNC(thermo, ThermoPhase, FUNC_NAME)

THERMO_1D(getMassFractions)
THERMO_1D(setMassFractions)
THERMO_1D(getMoleFractions)
THERMO_1D(setMoleFractions)
THERMO_1D(getConcentrations)
THERMO_1D(setConcentrations)

THERMO_1D(getMolecularWeights)
THERMO_1D(getChemPotentials)
THERMO_1D(getElectrochemPotentials)
THERMO_1D(getPartialMolarEnthalpies)
THERMO_1D(getPartialMolarEntropies)
THERMO_1D(getPartialMolarIntEnergies)
THERMO_1D(getPartialMolarCp)
THERMO_1D(getPartialMolarVolumes)
THERMO_1D(getEnthalpy_RT)
THERMO_1D(getEntropy_R)
THERMO_1D(getIntEnergy_RT)
THERMO_1D(getGibbs_RT)
THERMO_1D(getCp_R)
