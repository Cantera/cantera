//! @file utils.h
#ifndef TPX_UTILS_H
#define TPX_UTILS_H

#include "cantera/tpx/Sub.h"

namespace tpx
{
//! @deprecated To be removed after %Cantera 3.0
Substance* GetSub(int isub);

//! Create a new Substance object corresponding to the specified name.
/**
 * Currently valid substances are:
 *
 * - water
 * - nitrogen
 * - methane
 * - hydrogen
 * - oxygen
 * - HFC-134a
 * - carbon-dioxide
 * - heptane
 */
Substance* newSubstance(const std::string& name);

}

#endif
