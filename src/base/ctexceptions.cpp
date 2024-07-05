//! @file ctexceptions.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctexceptions.h"
#include "application.h"
#include "cantera/base/global.h"

#define BOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED
#include <boost/stacktrace.hpp>

#include <sstream>

namespace Cantera
{

int CanteraError::traceDepth_ = 0;

CanteraError::CanteraError(const string& procedure) :
    procedure_(procedure)
{
    if (traceDepth_) {
        auto trace = boost::stacktrace::stacktrace(0, traceDepth_);
        traceback_ = boost::stacktrace::to_string(trace);
    }
}

const char* CanteraError::what() const throw()
{
    try {
        formattedMessage_ = "\n" + string(79, '*') + "\n";
        formattedMessage_ += getClass();
        if (procedure_.size()) {
            formattedMessage_ += " thrown by " + procedure_;
        }
        formattedMessage_ += ":\n" + getMessage();
        if (formattedMessage_.compare(formattedMessage_.size()-1, 1, "\n")) {
            formattedMessage_.append("\n");
        }
        if (traceDepth_) {
            formattedMessage_ += string(79, '-') + "\n" + traceback_;
        }
        formattedMessage_ += string(79, '*') + "\n";

    } catch (...) {
        // Something went terribly wrong and we couldn't even format the message.
    }
    return formattedMessage_.c_str();
}

string CanteraError::getMessage() const
{
    return msg_;
}

string CanteraError::getMethod() const
{
    return procedure_;
}

void CanteraError::setStackTraceDepth(int depth)
{
    traceDepth_ = depth;
}

string ArraySizeError::getMessage() const
{
    return fmt::format("Array size ({}) too small. Must be at least {}.",
                       sz_, reqd_);
}

string IndexError::getMessage() const
{
    return fmt::format("IndexError: {}[{}] outside valid range of 0 to {}.",
                       arrayName_, m_, mmax_);
}

} // namespace Cantera
