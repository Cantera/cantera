//! @file pythonShim.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctexceptions.h"

#define BOOST_DLL_USE_STD_FS
#ifdef __MINGW32__
#define BOOST_DLL_FORCE_ALIAS_INSTANTIATION
#endif
#include <boost/dll/alias.hpp>

#include "Python.h"
#include <codecvt>

#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

namespace Cantera
{

namespace {

void checkPythonError(bool condition, const string& message) {
    if (condition) {
        if (PyErr_Occurred()) {
            PyErr_PrintEx(0);
        }
        throw CanteraError("loadCanteraPython", message);
    }
}

}

void loadCanteraPython()
{
    const char* venv_path = getenv("VIRTUAL_ENV");
    if (venv_path != nullptr) {
        PyConfig pyconf;
        PyConfig_InitPythonConfig(&pyconf);

        #ifdef _WIN32
            string suffix = "\\Scripts\\python.exe";
        #else
            string suffix = "/bin/python";
        #endif
        string path(venv_path);
        path += suffix;
        wstring wpath = wstring_convert<codecvt_utf8<wchar_t>>().from_bytes(path);
        PyStatus status = PyConfig_SetString(&pyconf, &pyconf.program_name,
                                                wpath.c_str());
        checkPythonError(PyStatus_Exception(status), "PyConfig_SetString failed");
        Py_InitializeFromConfig(&pyconf);
    } else {
        #if defined(CT_PYTHONHOME) && defined(_WIN32)
            const char* old_pythonhome = getenv("PYTHONHOME");
            if (old_pythonhome == nullptr || old_pythonhome[0] == '\0') {
                string pythonhome = "PYTHONHOME=";
                pythonhome += CT_PYTHONHOME;
                _putenv(pythonhome.c_str());
            }
        #endif
        Py_Initialize();
    }
    PyObject* pythonExt = PyImport_ImportModule("cantera");
    checkPythonError(pythonExt == nullptr, "cantera import failed");
    Py_DecRef(pythonExt);
}

// This creates and exports the name that is imported from Application::loadExtension
BOOST_DLL_ALIAS(loadCanteraPython,
                registerPythonExtensionManager);

}
