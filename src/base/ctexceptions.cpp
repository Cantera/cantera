//! @file ctexceptions.cpp
#include "cantera/base/ctexceptions.h"
#include "application.h"

#ifdef HAVE_FENV_H
#include <fenv.h>
#endif
#include <sstream>

namespace Cantera
{

// *** Exceptions ***

static const char* stars = "***********************************************************************\n";

CanteraError::CanteraError(const std::string& procedure, const std::string& msg) :
    procedure_(procedure),
    msg_(msg),
    saved_(false)
{
    // Save the error in the global list of errors so that showError() can work
    save();
}

CanteraError::CanteraError(const std::string& procedure) :
    procedure_(procedure),
    saved_(false)
{
    // Save the error in the global list of errors so that showError() can work
    save();
}

void CanteraError::save()
{
    if (!saved_) {
        Application::Instance()->addError(procedure_, getMessage());
        saved_ = true;
    }
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
    std::stringstream ss;
    ss << "Array size (" << sz_ << ") too small. Must be at least " << reqd_ << ".";
    return ss.str();
}

std::string IndexError::getMessage() const
{
    std::stringstream ss;
    ss << "IndexError: " << arrayName_ << "[" << m_ << "]" <<
       " outside valid range of 0 to " << (mmax_) << ".";
    return ss.str();
}

bool check_FENV_OverUnder_Flow() {
#ifdef HAVE_FENV_H
     fexcept_t ff;
     fegetexceptflag(&ff, FE_OVERFLOW || FE_UNDERFLOW || FE_INVALID);
     if (ff) {
        return true;
     }
#endif
     return false;
};

void clear_FENV() {
#ifdef HAVE_FENV_H
     feclearexcept(FE_ALL_EXCEPT);
#endif
}

} // namespace Cantera
