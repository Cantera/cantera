/**
 * @file mlloger.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef MLLOGGER_H
#define MLLOGGER_H

#include "ctmatutils.h"
#include "cantera/base/logger.h"

namespace Cantera
{

class ML_Logger : public Logger
{
public:
    ML_Logger() {}
    virtual ~ML_Logger() {}

    virtual void write(const std::string& s) {
        mexPrintf("%s", s.c_str());
    }

    virtual void writeendl() {
        mexPrintf("\n");
    }

    virtual void error(const std::string& msg) {
        mexErrMsgTxt(msg.c_str());
    }
};

}

#endif
