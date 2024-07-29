//! @file SundialsContext.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SUNDIALSCONTEXT_H
#define CT_SUNDIALSCONTEXT_H

#include "cantera/base/ct_defs.h"
#include "sundials/sundials_config.h"

#if SUNDIALS_VERSION_MAJOR >= 6
    #include "sundials/sundials_context.h"
#endif

namespace Cantera
{

//! A wrapper for managing a SUNContext object, need for Sundials >= 6.0
class SundialsContext
{
#if SUNDIALS_VERSION_MAJOR >= 6
public:
    SundialsContext() {
        #if SUNDIALS_VERSION_MAJOR >= 7
            SUNContext_Create(SUN_COMM_NULL, &m_context);

        #else
            SUNContext_Create(nullptr, &m_context);
        #endif
    }
    ~SundialsContext() {
        SUNContext_Free(&m_context);
    }
    SUNContext get() {
        return m_context;
    }

private:
    SUNContext m_context;
#endif
// For older Sundials versions, this is an empty class
};

}

#endif
