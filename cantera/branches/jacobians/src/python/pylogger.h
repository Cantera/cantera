#ifndef PYLOGGER_H
#define PYLOGGER_H

#include "Python.h"
#include <string>
#include "cantera/base/logger.h"

namespace Cantera
{

/// Logger for Python.
/// @ingroup textlogs
class Py_Logger : public Logger
{
public:
    Py_Logger() {
        PyRun_SimpleString("import sys");
    }
    virtual ~Py_Logger() {}

    virtual void write(const std::string& s) {
        std::string ss = "sys.stdout.write(\"\"\"";
        ss += s;
        ss += "\"\"\")";
        PyRun_SimpleString(ss.c_str());
        PyRun_SimpleString("sys.stdout.flush()");
    }

    virtual void error(const std::string& msg) {
        std::string err = "raise Exception(\"\"\""+msg+"\"\"\")";
        PyRun_SimpleString(err.c_str());
    }
};
}

#endif
