//! @file ctexceptions.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctexceptions.h"
#include "application.h"
#include "cantera/base/global.h"

#include <sstream>

namespace Cantera
{

// *** Exceptions ***

static const char* stars = "***********************************************************************\n";

CanteraError::CanteraError(const std::string& procedure) :
    procedure_(procedure)
{
}

const char* CanteraError::what() const throw()
{
    try {
        formattedMessage_ = "\n";
        formattedMessage_ += stars;
        formattedMessage_ += getClass();
        if (procedure_.size()) {
            formattedMessage_ += " thrown by " + procedure_;
        }
        formattedMessage_ += ":\n" + getMessage();
        if (formattedMessage_.compare(formattedMessage_.size()-1, 1, "\n")) {
            formattedMessage_.append("\n");
        }
        formattedMessage_ += stars;
    } catch (...) {
        // Something went terribly wrong and we couldn't even format the message.
    }
    return formattedMessage_.c_str();
}

std::string CanteraError::getMessage() const
{
    return msg_;
}

std::string ArraySizeError::getMessage() const
{
    return fmt::format("Array size ({}) too small. Must be at least {}.",
                       sz_, reqd_);
}

std::string IndexError::getMessage() const
{
    return fmt::format("IndexError: {}[{}] outside valid range of 0 to {}.",
                       arrayName_, m_, mmax_);
}

} // namespace Cantera
