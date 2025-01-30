//! @file pythonShim.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"

#define BOOST_DLL_USE_STD_FS
#ifdef __MINGW32__
#define BOOST_DLL_FORCE_ALIAS_INSTANTIATION
#endif
#include <boost/dll/alias.hpp>

#include "Python.h"
#include <filesystem>

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

#ifdef _WIN32
void setPythonHome()
{
    string pythonhome = "PYTHONHOME=";
    const char* old_pythonhome = getenv("PYTHONHOME");
    if (old_pythonhome != nullptr && old_pythonhome[0] != '\0') {
        // Use existing setting of PYTHONHOME
        return;
    }

    const char* conda_prefix = getenv("CONDA_PREFIX");
    if (conda_prefix != nullptr && std::filesystem::exists(conda_prefix)) {
        pythonhome += conda_prefix;
        _putenv(pythonhome.c_str());
        return;
    }

    #ifdef CT_PYTHONHOME
        string build_pythonhome = stripnonprint(string(CT_PYTHONHOME));
        if (std::filesystem::exists(build_pythonhome)) {
            pythonhome += build_pythonhome;
            _putenv(pythonhome.c_str());
        }
    #endif
}
#endif

}

void loadCanteraPython()
{
    // Prevent output buffering managed by Python.
    // @todo: It may be better to avoid replacing the existing Logger instance
    //     with a PythonLogger in the case of an embedded Python interpreter.
    const char* venv_path = getenv("VIRTUAL_ENV");
    if (venv_path != nullptr) {
        PyConfig pyconf;
        // pyconf.buffered_stdio = 0;
        PyConfig_InitPythonConfig(&pyconf);

        #ifdef _WIN32
            string suffix = "\\Scripts\\python.exe";
        #else
            string suffix = "/bin/python";
        #endif
        string path(venv_path);
        path += suffix;
        wstring wpath = std::filesystem::path(path).wstring();
        PyStatus status = PyConfig_SetString(&pyconf, &pyconf.program_name,
                                                wpath.c_str());
        checkPythonError(PyStatus_Exception(status), "PyConfig_SetString failed");
        Py_InitializeFromConfig(&pyconf);
    } else {
        #ifdef _WIN32
            setPythonHome();
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
