/**
 * @file mlloger.h
 */
#ifndef MLLOGGER_H
#define MLLOGGER_H

#include "ctmatutils.h"
#include "cantera/base/logger.h"

#include <iostream>

static std::string ss = "disp(' ";

namespace Cantera
{

class ML_Logger : public Logger
{
public:
    ML_Logger() {}
    virtual ~ML_Logger() {}


    virtual void write(const std::string& s) {
        char ch = s[0];
        int n = 0;
        while (ch != '\0') {
            if (ch =='\n') {
                ss += "');";

                mexEvalString(ss.c_str());
                ss = "disp(' ";
            } else {
                ss += ch;
            }
            if (ch == '\'') {
                ss += ch;
            }
            n++;
            ch = s[n];
        }
    }

    virtual void writeendl(const std::string& msg) {
        mexPrintf("\n");
    }

    virtual void error(const std::string& msg) {
        std::string err = "error("+msg+");";
        //mexEvalString(err.c_str());
        mexErrMsgTxt(err.c_str());
    }

    virtual int env() {
        return 1;
    }
};

}

#endif
