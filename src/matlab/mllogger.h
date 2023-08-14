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
    ~ML_Logger() override {}

    void write(const string& s) override {
        mexPrintf("%s", s.c_str());
    }

    void writeendl() override {
        mexPrintf("\n");
    }

    void error(const string& msg) override {
        mexErrMsgTxt(msg.c_str());
    }
};

}

#endif
