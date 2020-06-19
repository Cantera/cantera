// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/logger.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/kinetics/Kinetics.h"

#include "Python.h"

// Wrappers for preprocessor defines
std::string get_cantera_version()
{
    return std::string(CANTERA_VERSION);
}

int get_sundials_version()
{
    return CT_SUNDIALS_VERSION;
}

class PythonLogger : public Cantera::Logger
{
public:
    virtual void write(const std::string& s) {
        // 1000 bytes is the maximum size permitted by PySys_WriteStdout
        static const size_t N = 999;
        for (size_t i = 0; i < s.size(); i+=N) {
            PySys_WriteStdout("%s", s.substr(i, N).c_str());
        }
        std::cout.flush();
    }

    virtual void writeendl() {
        PySys_WriteStdout("%s", "\n");
        std::cout.flush();
    }

    virtual void error(const std::string& msg) {
        std::string err = "raise Exception('''"+msg+"''')";
        PyRun_SimpleString(err.c_str());
    }
};

// Function for assigning elements of Array2D, since Cython has trouble
// with assigning to the reference returned by operator()
void CxxArray2D_set(Cantera::Array2D& array, size_t i, size_t j, double value)
{
    array(i,j) = value;
}

// Function which populates a 1D array
#define ARRAY_FUNC(PREFIX, CLASS_NAME, FUNC_NAME) \
    void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, double* data) \
    { object->FUNC_NAME(data); }

// function which takes a stride as the first argument and populates a 2D array
#define ARRAY_FUNC2(PREFIX, CLASS_NAME, FUNC_NAME) \
    void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, size_t dim, double* data) \
    { object->FUNC_NAME(dim, data); }


#define THERMO_1D(FUNC_NAME) ARRAY_FUNC(thermo, ThermoPhase, FUNC_NAME)
#define KIN_1D(FUNC_NAME) ARRAY_FUNC(kin, Kinetics, FUNC_NAME)
#define TRANSPORT_1D(FUNC_NAME) ARRAY_FUNC(tran, Transport, FUNC_NAME)
#define TRANSPORT_2D(FUNC_NAME) ARRAY_FUNC2(tran, Transport, FUNC_NAME)

THERMO_1D(getMassFractions)
THERMO_1D(setMassFractions)
THERMO_1D(getMoleFractions)
THERMO_1D(setMoleFractions)
THERMO_1D(getConcentrations)
THERMO_1D(setConcentrations)

THERMO_1D(getMolecularWeights)
THERMO_1D(getCharges)
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
THERMO_1D(getActivities)
THERMO_1D(getActivityCoefficients)

KIN_1D(getFwdRatesOfProgress)
KIN_1D(getRevRatesOfProgress)
KIN_1D(getNetRatesOfProgress)

KIN_1D(getEquilibriumConstants)
KIN_1D(getFwdRateConstants)
KIN_1D(getRevRateConstants)

KIN_1D(getDeltaEnthalpy)
KIN_1D(getDeltaGibbs)
KIN_1D(getDeltaEntropy)
KIN_1D(getDeltaSSEnthalpy)
KIN_1D(getDeltaSSGibbs)
KIN_1D(getDeltaSSEntropy)

KIN_1D(getCreationRates)
KIN_1D(getDestructionRates)
KIN_1D(getNetProductionRates)

TRANSPORT_1D(getMixDiffCoeffs)
TRANSPORT_1D(getMixDiffCoeffsMass)
TRANSPORT_1D(getMixDiffCoeffsMole)
TRANSPORT_1D(getThermalDiffCoeffs)
TRANSPORT_1D(getSpeciesViscosities)
TRANSPORT_1D(getMobilities)

TRANSPORT_2D(getMultiDiffCoeffs)
TRANSPORT_2D(getBinaryDiffCoeffs)
