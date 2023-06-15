// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PY_TRANSPORT_UTILS_H
#define CT_PY_TRANSPORT_UTILS_H

#include "cantera/transport/Transport.h"
#include "wrappers.h"

#define TRANSPORT_1D(FUNC_NAME) ARRAY_FUNC(tran, Transport, FUNC_NAME)
#define TRANSPORT_2D(FUNC_NAME) ARRAY_FUNC2(tran, Transport, FUNC_NAME)
#define TRANSPORT_POLY(FUNC_NAME) ARRAY_POLY(tran, Transport, FUNC_NAME)
#define TRANSPORT_POLY_BINARY(FUNC_NAME) ARRAY_POLY_BINARY(tran, Transport, FUNC_NAME)

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

#endif
