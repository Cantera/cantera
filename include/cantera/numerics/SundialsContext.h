//! @file SundialsContext.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SUNDIALSCONTEXT_H
#define CT_SUNDIALSCONTEXT_H

#include "cantera/base/ct_defs.h"
#include "sundials/sundials_config.h"
#include "sundials/sundials_context.h"

namespace Cantera
{

//! A wrapper for managing a SUNContext object
class SundialsContext
{
public:
    SundialsContext() {
        #if SUNDIALS_VERSION_MAJOR >= 7
            SUNContext_Create(SUN_COMM_NULL, &m_context);
            SUNContext_ClearErrHandlers(m_context);
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
};

}

#endif
