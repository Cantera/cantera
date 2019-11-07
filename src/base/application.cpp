//! @file application.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "application.h"

#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"
#include "units.h"

#include <fstream>
#include <sstream>
#include <mutex>

using std::string;
using std::endl;

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
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

//! Mutex for controlling access to XML file storage
static std::mutex xml_mutex;

int get_modified_time(const std::string& path) {
#ifdef _WIN32
    HANDLE hFile = CreateFile(path.c_str(), NULL, NULL,
                              NULL, OPEN_EXISTING, 0, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        throw CanteraError("get_modified_time", "Couldn't open file:" + path);
    }
    FILETIME modified;
    GetFileTime(hFile, NULL, NULL, &modified);
    CloseHandle(hFile);
    return static_cast<int>(modified.dwLowDateTime);
#else
    struct stat attrib;
    stat(path.c_str(), &attrib);
    return static_cast<int>(attrib.st_mtime);
#endif
}

Application::Messages::Messages()
{
    // install a default logwriter that writes to standard
    // output / standard error
    logwriter.reset(new Logger());
}

void Application::Messages::addError(const std::string& r, const std::string& msg)
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

void Application::Messages::setLogger(Logger* _logwriter)
{
    logwriter.reset(_logwriter);
}

void Application::Messages::writelog(const std::string& msg)
{
    logwriter->write(msg);
}

void Application::Messages::writelogendl()
{
    logwriter->writeendl();
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

Application::Application() :
    m_suppress_deprecation_warnings(false),
    m_fatal_deprecation_warnings(false),
    m_suppress_thermo_warnings(false)
{
    // install a default logwriter that writes to standard
    // output / standard error
    setDefaultDirectories();
    Unit::units();
}

Application* Application::Instance()
{
    std::unique_lock<std::mutex> appLock(app_mutex);
    if (Application::s_app == 0) {
        Application::s_app = new Application();
    }
    return s_app;
}

Application::~Application()
{
    for (auto& f : xmlfiles) {
        f.second.first->unlock();
        delete f.second.first;
        f.second.first = 0;
    }
}

void Application::ApplicationDestroy()
{
    std::unique_lock<std::mutex> appLock(app_mutex);
    if (Application::s_app != 0) {
        delete Application::s_app;
        Application::s_app = 0;
    }
}

void Application::warn_deprecated(const std::string& method,
                                  const std::string& extra)
{
    if (m_fatal_deprecation_warnings) {
        throw CanteraError(method, "Deprecated: " + extra);
    } else if (m_suppress_deprecation_warnings || warnings.count(method)) {
        return;
    }
    warnings.insert(method);
    writelog(fmt::format("DeprecationWarning: {}: {}", method, extra));
    writelogendl();
}

void Application::warn_user(const std::string& method,
                            const std::string& extra)
{
    writelog(fmt::format("CanteraWarning: {}: {}", method, extra));
    writelogendl();
}

void Application::thread_complete()
{
    pMessenger.removeThreadMessages();
}

XML_Node* Application::get_XML_File(const std::string& file, int debug)
{
    std::unique_lock<std::mutex> xmlLock(xml_mutex);
    std::string path = findInputFile(file);
    int mtime = get_modified_time(path);

    if (xmlfiles.find(path) != xmlfiles.end()) {
        // Already have a parsed XML tree for this file cached. Check the
        // last-modified time.
        std::pair<XML_Node*, int> cache = xmlfiles[path];
        if (cache.second == mtime) {
            return cache.first;
        }
    }

    // Check whether or not the file is XML (based on the file extension). If
    // not, it will be first processed with the preprocessor.
    string::size_type idot = path.rfind('.');
    string ext;
    if (idot != string::npos) {
        ext = path.substr(idot, path.size());
    } else {
        ext = "";
    }
    XML_Node* x = new XML_Node("doc");
    if (ext != ".xml" && ext != ".ctml") {
        // Assume that we are trying to open a cti file. Do the conversion to XML.
        std::stringstream phase_xml(ct2ctml_string(path));
        x->build(phase_xml, path);
    } else {
        x->build(path);
    }
    x->lock();
    xmlfiles[path] = {x, mtime};
    return x;
}

XML_Node* Application::get_XML_from_string(const std::string& text)
{
    std::unique_lock<std::mutex> xmlLock(xml_mutex);
    std::pair<XML_Node*, int>& entry = xmlfiles[text];
    if (entry.first) {
        // Return existing cached XML tree
        return entry.first;
    }
    std::stringstream s;
    size_t start = text.find_first_not_of(" \t\r\n");
    if (text.substr(start,1) == "<") {
        s << text;
    } else {
        s << ct_string2ctml_string(text.substr(start));
    }
    entry.first = new XML_Node();
    entry.first->build(s, "[string]");
    return entry.first;
}

void Application::close_XML_File(const std::string& file)
{
    std::unique_lock<std::mutex> xmlLock(xml_mutex);
    if (file == "all") {
        for (const auto& f : xmlfiles) {
            f.second.first->unlock();
            delete f.second.first;
        }
        xmlfiles.clear();
    } else if (xmlfiles.find(file) != xmlfiles.end()) {
        xmlfiles[file].first->unlock();
        delete xmlfiles[file].first;
        xmlfiles.erase(file);
    }
}

#ifdef _WIN32
long int Application::readStringRegistryKey(const std::string& keyName, const std::string& valueName,
        std::string& value, const std::string& defaultValue)
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

std::string Application::Messages::lastErrorMessage()
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
        f << errorMessage[j] << endl;
    }
    errorMessage.clear();
}

void Application::Messages::logErrors()
{
    for (size_t j = 0; j < errorMessage.size(); j++) {
        writelog(errorMessage[j]);
        writelogendl();
    }
    errorMessage.clear();
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
    std::string pathsep = ";";
#else
    std::string pathsep = ":";
#endif

    if (getenv("CANTERA_DATA") != 0) {
        string s = string(getenv("CANTERA_DATA"));
        size_t start = 0;
        size_t end = s.find(pathsep);
        while(end != npos) {
            inputDirs.push_back(s.substr(start, end-start));
            start = end + 1;
            end = s.find(pathsep, start);
        }
        inputDirs.push_back(s.substr(start,end));
    }

#ifdef _WIN32
    // Under Windows, the Cantera setup utility records the installation
    // directory in the registry. Data files are stored in the 'data'
    // subdirectory of the main installation directory.
    std::string installDir;
    readStringRegistryKey("SOFTWARE\\Cantera\\Cantera " CANTERA_SHORT_VERSION,
                          "InstallDir", installDir, "");
    if (installDir != "") {
        inputDirs.push_back(installDir + "data");

        // Scripts for converting mechanisms to CTI and CMTL are installed in
        // the 'bin' subdirectory. Add that directory to the PYTHONPATH.
        const char* old_pythonpath = getenv("PYTHONPATH");
        std::string pythonpath = "PYTHONPATH=" + installDir + "\\bin";
        if (old_pythonpath) {
            pythonpath += ";";
            pythonpath.append(old_pythonpath);
        }
        _putenv(pythonpath.c_str());
    }

#endif

#ifdef DARWIN
    // add a default data location for Mac OS X
    inputDirs.push_back("/Applications/Cantera/data");
#endif

    // CANTERA_DATA is defined in file config.h. This file is written during the
    // build process (unix), and points to the directory specified by the
    // 'prefix' option to 'configure', or else to /usr/local/cantera.
#ifdef CANTERA_DATA
    string datadir = string(CANTERA_DATA);
    inputDirs.push_back(datadir);
#endif
}

void Application::addDataDirectory(const std::string& dir)
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

std::string Application::findInputFile(const std::string& name)
{
    std::unique_lock<std::mutex> dirLock(dir_mutex);
    string::size_type islash = name.find('/');
    string::size_type ibslash = name.find('\\');
    std::vector<string>& dirs = inputDirs;

    // Expand "~/" to user's home directory, if possible
    if (name.find("~/") == 0 || name.find("~\\") == 0) {
        char* home = getenv("HOME"); // POSIX systems
        if (!home) {
            home = getenv("USERPROFILE"); // Windows systems
        }
        if (home) {
            return home + name.substr(1, npos);
        }
    }

    if (islash == string::npos && ibslash == string::npos) {
        size_t nd = dirs.size();
        for (size_t i = 0; i < nd; i++) {
            string inname = dirs[i] + "/" + name;
            std::ifstream fin(inname);
            if (fin) {
                return inname;
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

    return name;
}

Application* Application::s_app = 0;

} // namespace Cantera
