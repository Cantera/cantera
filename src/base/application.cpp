//! @file application.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "application.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ExtensionManagerFactory.h"

#define BOOST_DLL_USE_STD_FS
#include <boost/dll/import.hpp>
#include <boost/algorithm/string.hpp>

#include <fstream>
#include <sstream>
#include <mutex>

namespace ba = boost::algorithm;

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef _MSC_VER
#pragma comment(lib, "advapi32")
#endif

namespace Cantera
{

//! Mutex for input directory access
static std::mutex dir_mutex;

//! Mutex for creating singletons within the application object
static std::mutex app_mutex;

void Application::Messages::addError(const string& r, const string& msg)
{
    if (msg.size() != 0) {
        errorMessage.push_back(
            "\n\n************************************************\n"
            "                Cantera Error!                  \n"
            "************************************************\n\n"
            "Procedure: " + r +
            "\nError:     " + msg + "\n");
    } else {
        errorMessage.push_back(r);
    }
}

int Application::Messages::getErrorCount()
{
    return static_cast<int>(errorMessage.size());
}

//! Mutex for access to string messages
static std::mutex msg_mutex;

Application::Messages* Application::ThreadMessages::operator ->()
{
    std::unique_lock<std::mutex> msgLock(msg_mutex);
    std::thread::id curId = std::this_thread::get_id();
    auto iter = m_threadMsgMap.find(curId);
    if (iter != m_threadMsgMap.end()) {
        return iter->second.get();
    }
    pMessages_t pMsgs(new Messages());
    m_threadMsgMap.insert({curId, pMsgs});
    return pMsgs.get();
}

void Application::ThreadMessages::removeThreadMessages()
{
    std::unique_lock<std::mutex> msgLock(msg_mutex);
    std::thread::id curId = std::this_thread::get_id();
    auto iter = m_threadMsgMap.find(curId);
    if (iter != m_threadMsgMap.end()) {
        m_threadMsgMap.erase(iter);
    }
}

Application::Application()
{
    // install a default logwriter that writes to standard output / standard error
    m_logwriter = make_unique<Logger>();
    setDefaultDirectories();
}

Application* Application::Instance()
{
    std::unique_lock<std::mutex> appLock(app_mutex);
    if (Application::s_app == 0) {
        Application::s_app = new Application();
    }
    return s_app;
}

void Application::ApplicationDestroy()
{
    std::unique_lock<std::mutex> appLock(app_mutex);
    if (Application::s_app != 0) {
        delete Application::s_app;
        Application::s_app = 0;
    }
}

void Application::warn_deprecated(const string& method, const string& extra)
{
    if (m_fatal_deprecation_warnings) {
        throw CanteraError(method, "Deprecated: " + extra);
    } else if (m_suppress_deprecation_warnings || warnings.count(method)) {
        return;
    }
    warnings.insert(method);
    warnlog("Deprecation", fmt::format("{}: {}", method, extra));
}

void Application::warn(const string& warning, const string& method, const string& extra)
{
    if (m_fatal_warnings) {
        throw CanteraError(method, extra);
    } else if (m_suppress_warnings) {
        return;
    }
    warnlog(warning, fmt::format("{}: {}", method, extra));
}

void Application::thread_complete()
{
    pMessenger.removeThreadMessages();
}

#ifdef _WIN32
long int Application::readStringRegistryKey(const string& keyName, const string& valueName,
        string& value, const string& defaultValue)
{
    HKEY key;
    long open_error = RegOpenKeyEx(HKEY_LOCAL_MACHINE, keyName.c_str(), 0, KEY_READ, &key);
    if (open_error != ERROR_SUCCESS) {
        return open_error;
    }
    value = defaultValue;
    CHAR buffer[1024];
    DWORD bufferSize = sizeof(buffer);
    ULONG error;
    error = RegQueryValueEx(key, valueName.c_str(), 0, NULL, (LPBYTE) buffer, &bufferSize);
    if (ERROR_SUCCESS == error) {
        value = buffer;
    }
    RegCloseKey(key);
    return error;
}
#endif

void Application::Messages::popError()
{
    if (!errorMessage.empty()) {
        errorMessage.pop_back();
    }
}

string Application::Messages::lastErrorMessage()
{
    if (!errorMessage.empty()) {
        return errorMessage.back();
    } else  {
        return "<no Cantera error>";
    }
}

void Application::Messages::getErrors(std::ostream& f)
{
    for (size_t j = 0; j < errorMessage.size(); j++) {
        f << errorMessage[j] << std::endl;
    }
    errorMessage.clear();
}

void Application::Messages::logErrors()
{
    for (size_t j = 0; j < errorMessage.size(); j++) {
        Cantera::writelog(errorMessage[j]);
        Cantera::writelogendl();
    }
    errorMessage.clear();
}

// Application methods

void Application::setLogger(unique_ptr<Logger> _logwriter)
{
    m_logwriter = std::move(_logwriter);
}

void Application::writelog(const string& msg)
{
    m_logwriter->write(msg);
}

void Application::writelogendl()
{
    m_logwriter->writeendl();
}

void Application::warnlog(const string& warning, const string& msg)
{
    m_logwriter->warn(warning, msg);
}

void Application::setDefaultDirectories()
{
    // always look in the local directory first
    inputDirs.push_back(".");

    // if environment variable CANTERA_DATA is defined, then add it to the
    // search path. CANTERA_DATA may include multiple directory, separated by
    // the OS-dependent path separator (in the same manner as the PATH
    // environment variable).
#ifdef _WIN32
    string pathsep = ";";
    string dirsep = "\\";
#else
    string pathsep = ":";
    string dirsep = "/";
#endif

    if (getenv("CANTERA_DATA") != nullptr) {
        string s(getenv("CANTERA_DATA"));
        size_t start = 0;
        size_t end = s.find(pathsep);
        while (end != npos) {
            inputDirs.push_back(s.substr(start, end-start));
            start = end + 1;
            end = s.find(pathsep, start);
        }
        inputDirs.push_back(s.substr(start,end));
    }

    // If a conda environment is active, add the location of the Cantera data directory
    // that may exist in that environment.
    if (getenv("CONDA_PREFIX") != nullptr) {
        string s(getenv("CONDA_PREFIX"));
        inputDirs.push_back(s + dirsep + "share" + dirsep + "cantera" + dirsep + "data");
    }

#ifdef _WIN32
    // Under Windows, the Cantera setup utility records the installation
    // directory in the registry. Data files are stored in the 'data'
    // subdirectory of the main installation directory.
    string installDir;
    readStringRegistryKey("SOFTWARE\\Cantera\\Cantera " CANTERA_SHORT_VERSION,
                          "InstallDir", installDir, "");
    if (installDir != "") {
        inputDirs.push_back(installDir + "data");

        // Scripts for converting mechanisms to YAML are installed in
        // the 'bin' subdirectory. Add that directory to the PYTHONPATH.
        const char* old_pythonpath = getenv("PYTHONPATH");
        string pythonpath = "PYTHONPATH=" + installDir + "\\bin";
        if (old_pythonpath) {
            pythonpath += ";";
            pythonpath.append(old_pythonpath);
        }
        _putenv(pythonpath.c_str());
    }

#endif

    // CANTERA_DATA is defined in file config.h. This file is written during the
    // build process (unix), and points to the directory specified by the
    // 'prefix' option to 'configure', or else to /usr/local/cantera.
#ifdef CANTERA_DATA
    string datadir = stripnonprint(string(CANTERA_DATA));
    inputDirs.push_back(datadir);
#endif
}

void Application::addDataDirectory(const string& dir)
{
    std::unique_lock<std::mutex> dirLock(dir_mutex);
    if (inputDirs.empty()) {
        setDefaultDirectories();
    }
    string d = stripnonprint(dir);

    // Expand "~/" to user's home directory, if possible
    if (d.find("~/") == 0 || d.find("~\\") == 0) {
        char* home = getenv("HOME"); // POSIX systems
        if (!home) {
            home = getenv("USERPROFILE"); // Windows systems
        }
        if (home) {
            d = home + d.substr(1, npos);
        }
    }

    // Remove any existing entry for this directory
    auto iter = std::find(inputDirs.begin(), inputDirs.end(), d);
    if (iter != inputDirs.end()) {
        inputDirs.erase(iter);
    }

    // Insert this directory at the beginning of the search path
    inputDirs.insert(inputDirs.begin(), d);
}

string Application::findInputFile(const string& name)
{
    std::unique_lock<std::mutex> dirLock(dir_mutex);
    string::size_type islash = name.find('/');
    string::size_type ibslash = name.find('\\');
    string::size_type icolon = name.find(':');
    vector<string>& dirs = inputDirs;

    // Expand "~/" to user's home directory, if possible
    if (name.find("~/") == 0 || name.find("~\\") == 0) {
        char* home = getenv("HOME"); // POSIX systems
        if (!home) {
            home = getenv("USERPROFILE"); // Windows systems
        }
        if (home) {
            string full_name = home + name.substr(1, npos);
            std::ifstream fin(full_name);
            if (fin) {
                return full_name;
            } else {
                throw CanteraError("Application::findInputFile",
                                   "Input file '{}' not found", name);
            }
        }
    }

    // If this is an absolute path, just look for the file there
    if (islash == 0 || ibslash == 0
        || (icolon == 1 && (ibslash == 2 || islash == 2)))
    {
        std::ifstream fin(name);
        if (fin) {
            return name;
        } else {
            throw CanteraError("Application::findInputFile",
                               "Input file '{}' not found", name);
        }
    }

    // Search the Cantera data directories for the input file, and return
    // the full path if a match is found
    size_t nd = dirs.size();
    for (size_t i = 0; i < nd; i++) {
        string full_name = dirs[i] + "/" + name;
        std::ifstream fin(full_name);
        if (fin) {
            return full_name;
        }
    }
    string msg = "\nInput file " + name + " not found in director";
    msg += (nd == 1 ? "y " : "ies ");
    for (size_t i = 0; i < nd; i++) {
        msg += "\n'" + dirs[i] + "'";
        if (i+1 < nd) {
            msg += ", ";
        }
    }
    msg += "\n\n";
    msg += "To fix this problem, either:\n";
    msg += "    a) move the missing files into the local directory;\n";
    msg += "    b) define environment variable CANTERA_DATA to\n";
    msg += "         point to the directory containing the file.";
    throw CanteraError("Application::findInputFile", msg);
}

void Application::loadExtension(const string& extType, const string& name)
{
    if (!usingSharedLibrary()) {
        throw CanteraError("Application::loadExtension",
            "Loading extensions requires linking to the Cantera shared library\n"
            "rather than the static library");
    }
    if (m_loaded_extensions.count({extType, name})) {
        return;
    }

    if (extType == "python" && !ExtensionManagerFactory::factory().exists("python")) {
        string errors;

        // type of imported symbol: void function with no arguments
        typedef void (loader_t)();

        // Only one Python module can be loaded at a time, and a handle needs to be held
        // to prevent it from being unloaded.
        static function<loader_t> loader;
        bool loaded = false;

        for (const auto& py_ver : m_pythonSearchVersions) {
            string py_ver_underscore = ba::replace_all_copy(py_ver, ".", "_");
            try {
                loader = boost::dll::import_alias<loader_t>(
                    "cantera_python" + py_ver_underscore, // library name
                    "registerPythonExtensionManager", // symbol to import
                    // append extensions and prefixes, search normal library path, and
                    // expose all loaded symbols (specifically, those from libpython)
                    boost::dll::load_mode::search_system_folders
                    | boost::dll::load_mode::append_decorations
                    | boost::dll::load_mode::rtld_global
                );
                loader();
                loaded = true;
                break;
            } catch (std::exception& err) {
                errors += fmt::format("\nPython {}: {}\n", py_ver, err.what());
            }
        }
        if (!loaded) {
            throw CanteraError("Application::loadExtension",
                "Error loading Python extension support. Tried the following:{}",
                errors);
        }
    }
    ExtensionManagerFactory::build(extType)->registerRateBuilders(name);
    m_loaded_extensions.insert({extType, name});
}

void Application::searchPythonVersions(const string& versions) {
    ba::split(m_pythonSearchVersions, versions, ba::is_any_of(","));
}

Application* Application::s_app = 0;

} // namespace Cantera
