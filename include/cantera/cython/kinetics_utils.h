// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef PY_KINETICS_UTILS_H
#define PY_KINETICS_UTILS_H

#include "cantera/kinetics/Kinetics.h"
#include "cantera/numerics/eigen_sparse.h"
#include "wrappers.h"

// Service function to pass index/value triplets describing sparse matrix
inline size_t sparseTriplets(const Eigen::SparseMatrix<double>& mat,
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
inline void sparseCscData(const Eigen::SparseMatrix<double>& mat,
    double* value, int* inner, int* outer)
{
    if (!mat.isCompressed()) {
        throw Cantera::CanteraError("sparseCscData",
            "Invalid input: Eigen matrix is not compressed.");
    }

    const double* valuePtr = mat.valuePtr();
    const int* innerPtr = mat.innerIndexPtr();
    for (int i = 0; i < mat.nonZeros(); ++i) {
        value[i] = valuePtr[i];
        inner[i] = innerPtr[i];
    }

    const int* outerPtr = mat.outerIndexPtr();
    for (int i = 0; i < mat.outerSize() + 1; ++i) {
        outer[i] = outerPtr[i];
    }
}

#define KIN_1D(FUNC_NAME) ARRAY_FUNC(kin, Kinetics, FUNC_NAME)

KIN_1D(getFwdRatesOfProgress)
KIN_1D(getRevRatesOfProgress)
KIN_1D(getNetRatesOfProgress)

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

KIN_1D(getCreationRates_ddT)
KIN_1D(getDestructionRates_ddT)
KIN_1D(getNetProductionRates_ddT)

KIN_1D(getCreationRates_ddP)
KIN_1D(getDestructionRates_ddP)
KIN_1D(getNetProductionRates_ddP)

KIN_1D(getCreationRates_ddC)
KIN_1D(getDestructionRates_ddC)
KIN_1D(getNetProductionRates_ddC)

#endif
