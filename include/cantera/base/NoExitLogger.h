//! @file NoExitLogger.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef NOEXITLOGGER_H
#define NOEXITLOGGER_H

#include <string>
#include "cantera/base/logger.h"

namespace Cantera {
/// Logger that doesn't exit when an error is thrown.
/// @ingroup textlogs
class NoExitLogger : public Logger {
public:
    NoExitLogger() {}
    virtual ~NoExitLogger() {}

    virtual void error(const std::string& msg)
    {
       std::cerr << msg << std::endl;
    }
};
}
#endif
