#ifndef PYLOGGER_H
#define PYLOGGER_H

#include "Python.h"
#include <string>
#include "cantera/base/logger.h"

static std::string ss = "print \"\"\" ";

namespace Cantera
{

/// Logger for Python.
/// @ingroup textlogs
class Py_Logger : public Logger
{
public:
    Py_Logger() {}
    virtual ~Py_Logger() {}

    virtual void write(const std::string& s) {
        char ch = s[0];
        int n = 0;
        while (ch != '\0') {
            if (ch =='\n') {
                ss += "\"\"\"";
                PyRun_SimpleString((char*)ss.c_str());
                ss = "print \"\"\"";
            } else {
                ss += ch;
            }
            n++;
            ch = s[n];
        }
    }

    virtual void error(const std::string& msg) {
        std::string err = "raise \""+msg+"\"";
        PyRun_SimpleString((char*)err.c_str());
    }
};
}

#endif
