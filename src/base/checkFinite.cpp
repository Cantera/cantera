/**
 *   @file checkFinite.cpp Declarations for routines that check for the
 *       presence of NaNs in the code.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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
                    message += fmt::format("{}[{}] = {}\n", name, j, values[j]);
                }
            }
            throw CanteraError("checkFinite", message);
        }
    }
}

}
