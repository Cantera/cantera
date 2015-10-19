/**
 *   @file checkFinite.cpp Declarations for routines that check for the
 *       presence of NaNs in the code.
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "cantera/base/ct_defs.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera {

void checkFinite(const double tmp)
{
    if (!std::isfinite(tmp)) {
        if (std::isnan(tmp)) {
            throw CanteraError("checkFinite", "found NaN");
        } else if (tmp > 0) {
            throw CanteraError("checkFinite", "found +Inf");
        } else {
            throw CanteraError("checkFinite", "found -Inf");
        }
    }
}

void checkFinite(const std::string& name, double* values, size_t N)
{
    for (size_t i = 0; i < N; i++) {
        if (!std::isfinite(values[i])) {
            std::string message = name + " contains non-finite elements:\n\n";
            for (size_t j = 0; j < N; j++) {
                if (!std::isfinite(values[j])) {
                    message += name + "[" + int2str(j) + "] = " +
                               fp2str(values[j]) + "\n";
                }
            }
            throw CanteraError("checkFinite", message);
        }
    }
}

}
