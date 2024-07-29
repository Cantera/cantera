// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PY_UTILS_UTILS_H
#define CT_PY_UTILS_UTILS_H

#include "cantera/base/logger.h"
#include "sundials/sundials_config.h"

#include "Python.h"

// Warning types supported by the Python C-API.
// See https://docs.python.org/3/c-api/exceptions.html#issuing-warnings
static std::map<std::string, PyObject*> mapped_PyWarnings = {
    {"", PyExc_Warning},
    {"Bytes", PyExc_BytesWarning},
    {"Cantera", PyExc_UserWarning}, // pre-existing warning
    {"Deprecation", PyExc_DeprecationWarning},
    {"Future", PyExc_FutureWarning},
    {"Import", PyExc_ImportWarning},
    {"PendingDeprecation", PyExc_PendingDeprecationWarning},
    {"Resource", PyExc_ResourceWarning},
    {"Runtime", PyExc_RuntimeWarning},
    {"Syntax", PyExc_SyntaxWarning},
    {"Unicode", PyExc_UnicodeWarning},
    {"User", PyExc_UserWarning}
};

// Cantera version for Python module compilation
inline std::string get_cantera_version_py()
{
    return CANTERA_VERSION;
}

// Git commit for Python module compilation
inline std::string get_cantera_git_commit_py()
{
#ifdef GIT_COMMIT
    return GIT_COMMIT;
#else
    return "unknown";
#endif
}

// Wrappers for preprocessor defines
inline std::string get_sundials_version()
{
    return SUNDIALS_VERSION;
}

class PythonLogger : public Cantera::Logger
{
public:
    void write(const std::string& s) override {
        // 1000 bytes is the maximum size permitted by PySys_WriteStdout
        static const size_t N = 999;
        for (size_t i = 0; i < s.size(); i+=N) {
            PySys_WriteStdout("%s", s.substr(i, N).c_str());
        }
        std::cout.flush();
    }

    void writeendl() override {
        PySys_WriteStdout("%s", "\n");
        std::cout.flush();
    }

    void warn(const std::string& warning, const std::string& msg) override {
        if (mapped_PyWarnings.find(warning) != mapped_PyWarnings.end()) {
            PyErr_WarnEx(mapped_PyWarnings[warning], msg.c_str(), 1);
        } else {
            // issue generic warning
            PyErr_WarnEx(PyExc_Warning, msg.c_str(), 1);
        }
    }

    void error(const std::string& msg) override {
        PyErr_SetString(PyExc_RuntimeError, msg.c_str());
    }
};

#endif
