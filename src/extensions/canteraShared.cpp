//! @file canteraShared.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ct_defs.h"

namespace Cantera
{

bool usingSharedLibrary()
{
    // This implementation of usingSharedLibrary is compiled and embedded
    // only in the Cantera shared library
    return true;
}

}
