// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/logger.h"
#include "cantera/numerics/eigen_sparse.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/kinetics/Kinetics.h"

#include "Python.h"

// Warning types supported by the Python C-API.
// See https://docs.python.org/3/c-api/exceptions.html#issuing-warnings
std::map<std::string, PyObject*> mapped_PyWarnings = {
    {"", PyExc_Warning},
    {"Bytes", PyExc_BytesWarning},
    {"Cantera", PyExc_UserWarning}, // pre-existing warning
    {"Deprecation", PyExc_DeprecationWarning},
    {"Future", PyExc_FutureWarning},
    {"Import", PyExc_ImportWarning},
    {"PendingDeprecation", PyExc_PendingDeprecationWarning},
    {"Resource", PyExc_ResourceWarning},
    {"Runtime", PyExc_RuntimeWarning},
    {"Syntax", PyExc_SyntaxWarning},
    {"Unicode", PyExc_UnicodeWarning},
    {"User", PyExc_UserWarning}
};

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

    virtual void warn(const std::string& warning, const std::string& msg) {
        if (mapped_PyWarnings.find(warning) != mapped_PyWarnings.end()) {
            PyErr_WarnEx(mapped_PyWarnings[warning], msg.c_str(), 1);
        } else {
            // issue generic warning
            PyErr_WarnEx(PyExc_Warning, msg.c_str(), 1);
        }
    }

    virtual void error(const std::string& msg) {
        PyErr_SetString(PyExc_RuntimeError, msg.c_str());
    }
};

// Function for assigning elements of Array2D, since Cython has trouble
// with assigning to the reference returned by operator()
void CxxArray2D_set(Cantera::Array2D& array, size_t i, size_t j, double value)
{
    array(i,j) = value;
}

// Service function to pass index/value triplets describing sparse matrix
size_t sparseTriplets(const Eigen::SparseMatrix<double>& mat,
    int* rows, int* cols, double* data, size_t length)
{
    size_t count = 0;
    for (int i = 0; i < mat.outerSize(); i++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
            if (count < length) {
                rows[count] = it.row();
                cols[count] = it.col();
                data[count] = it.value();
            }
            count++;
        }
    }
    if (count > length) {
        throw Cantera::CanteraError("sparseComponents",
            "Output arrays have insufficient length. Required size is {}, "
            "while provided length is {}.", count, length);
    }
    return count;
}

// Service function to pass CSC data describing sparse matrix
void sparseCscData(const Eigen::SparseMatrix<double>& mat,
    double* value, int* inner, int* outer)
{
    if (!mat.isCompressed()) {
        throw Cantera::CanteraError("sparseCscData",
            "Invalid input: Eigen matrix is not compressed.");
    }

    const double* valuePtr = mat.valuePtr();
    const int* innerPtr = mat.innerIndexPtr();
    for (size_t i = 0; i < mat.nonZeros(); ++i) {
        value[i] = valuePtr[i];
        inner[i] = innerPtr[i];
    }

    const int* outerPtr = mat.outerIndexPtr();
    for (size_t i = 0; i < mat.outerSize() + 1; ++i) {
        outer[i] = outerPtr[i];
    }
}

// Function which passes sparse matrix
#define SPARSE_MATRIX(PREFIX, CLASS_NAME, FUNC_NAME) \
    Eigen::SparseMatrix<double> PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object) \
    { return object->FUNC_NAME(); }

// Function which populates a 1D array
#define ARRAY_FUNC(PREFIX, CLASS_NAME, FUNC_NAME) \
    void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, double* data) \
    { object->FUNC_NAME(data); }

// function which takes a stride as the first argument and populates a 2D array
#define ARRAY_FUNC2(PREFIX, CLASS_NAME, FUNC_NAME) \
    void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, size_t dim, double* data) \
    { object->FUNC_NAME(dim, data); }

// Function which populates a 1D array, extra arguments
#define ARRAY_POLY(PREFIX, CLASS_NAME, FUNC_NAME) \
    void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, size_t i, double* data) \
    { object->FUNC_NAME(i, data); }

#define ARRAY_POLY_BINARY(PREFIX, CLASS_NAME, FUNC_NAME) \
    void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, size_t i, size_t j, double* data) \
    { object->FUNC_NAME(i, j, data); }

#define THERMO_1D(FUNC_NAME) ARRAY_FUNC(thermo, ThermoPhase, FUNC_NAME)
#define KIN_1D(FUNC_NAME) ARRAY_FUNC(kin, Kinetics, FUNC_NAME)
#define KIN_SPARSE_MATRIX(FUNC_NAME) SPARSE_MATRIX(kin, Kinetics, FUNC_NAME)
#define TRANSPORT_1D(FUNC_NAME) ARRAY_FUNC(tran, Transport, FUNC_NAME)
#define TRANSPORT_2D(FUNC_NAME) ARRAY_FUNC2(tran, Transport, FUNC_NAME)
#define TRANSPORT_POLY(FUNC_NAME) ARRAY_POLY(tran, Transport, FUNC_NAME)
#define TRANSPORT_POLY_BINARY(FUNC_NAME) ARRAY_POLY_BINARY(tran, Transport, FUNC_NAME)

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

KIN_SPARSE_MATRIX(reactantStoichCoeffs)
KIN_SPARSE_MATRIX(productStoichCoeffs)
KIN_SPARSE_MATRIX(revProductStoichCoeffs)

KIN_1D(getFwdRatesOfProgress)
KIN_1D(getRevRatesOfProgress)
KIN_1D(getNetRatesOfProgress)

KIN_SPARSE_MATRIX(fwdRatesOfProgress_ddX)
KIN_SPARSE_MATRIX(revRatesOfProgress_ddX)
KIN_SPARSE_MATRIX(netRatesOfProgress_ddX)

KIN_1D(getFwdRateConstants_ddT)
KIN_1D(getFwdRateConstants_ddP)
KIN_1D(getFwdRateConstants_ddC)

KIN_1D(getFwdRatesOfProgress_ddT)
KIN_1D(getRevRatesOfProgress_ddT)
KIN_1D(getNetRatesOfProgress_ddT)

KIN_1D(getFwdRatesOfProgress_ddP)
KIN_1D(getRevRatesOfProgress_ddP)
KIN_1D(getNetRatesOfProgress_ddP)

KIN_1D(getFwdRatesOfProgress_ddC)
KIN_1D(getRevRatesOfProgress_ddC)
KIN_1D(getNetRatesOfProgress_ddC)

KIN_1D(getEquilibriumConstants)
KIN_1D(getFwdRateConstants)
KIN_1D(getRevRateConstants)

KIN_1D(getDeltaEnthalpy)
KIN_1D(getDeltaGibbs)
KIN_1D(getDeltaEntropy)
KIN_1D(getDeltaSSEnthalpy)
KIN_1D(getDeltaSSGibbs)
KIN_1D(getDeltaSSEntropy)

KIN_1D(getThirdBodyConcentrations)

KIN_1D(getCreationRates)
KIN_1D(getDestructionRates)
KIN_1D(getNetProductionRates)

KIN_SPARSE_MATRIX(creationRates_ddX)
KIN_SPARSE_MATRIX(destructionRates_ddX)
KIN_SPARSE_MATRIX(netProductionRates_ddX)

KIN_1D(getCreationRates_ddT)
KIN_1D(getDestructionRates_ddT)
KIN_1D(getNetProductionRates_ddT)

KIN_1D(getCreationRates_ddP)
KIN_1D(getDestructionRates_ddP)
KIN_1D(getNetProductionRates_ddP)

KIN_1D(getCreationRates_ddC)
KIN_1D(getDestructionRates_ddC)
KIN_1D(getNetProductionRates_ddC)

TRANSPORT_1D(getMixDiffCoeffs)
TRANSPORT_1D(getMixDiffCoeffsMass)
TRANSPORT_1D(getMixDiffCoeffsMole)
TRANSPORT_1D(getThermalDiffCoeffs)
TRANSPORT_1D(getSpeciesViscosities)
TRANSPORT_1D(getMobilities)

TRANSPORT_2D(getMultiDiffCoeffs)
TRANSPORT_2D(getBinaryDiffCoeffs)

TRANSPORT_POLY(getViscosityPolynomial)
TRANSPORT_POLY(setViscosityPolynomial)
TRANSPORT_POLY(getConductivityPolynomial)
TRANSPORT_POLY(setConductivityPolynomial)
TRANSPORT_POLY_BINARY(getBinDiffusivityPolynomial)
TRANSPORT_POLY_BINARY(setBinDiffusivityPolynomial)
