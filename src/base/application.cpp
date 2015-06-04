//! @file application.cpp
#include "application.h"

#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"
#include "units.h"

#include <fstream>
#include <sstream>

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

// If running multiple threads in a cpp application, the Application class
// is the only internal object that is single instance with static data.

#ifdef THREAD_SAFE_CANTERA
cthreadId_t getThisThreadId()
{
#if defined(BOOST_HAS_WINTHREADS)
    return ::GetCurrentThreadId();
#elif defined(BOOST_HAS_PTHREADS)
    return pthread_self();
#endif
}

#endif

//! Mutex for input directory access
static mutex_t dir_mutex;

//! Mutex for creating singletons within the application object
static mutex_t app_mutex;

//! Mutex for controlling access to XML file storage
static mutex_t xml_mutex;

static int get_modified_time(const std::string& path) {
#ifdef _WIN32
    HANDLE hFile = CreateFile(path.c_str(), GENERIC_READ, FILE_SHARE_WRITE,
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

Application::Messages::Messages() :
    logwriter(0)
{
    // install a default logwriter that writes to standard
    // output / standard error
    logwriter = new Logger();
}

Application::Messages::Messages(const Messages& r) :
    errorMessage(r.errorMessage),
    errorRoutine(r.errorRoutine),
    logwriter(0)
{
    // install a default logwriter that writes to standard
    // output / standard error
    logwriter = new Logger(*(r.logwriter));
}

Application::Messages& Application::Messages::operator=(const Messages& r)
{
    if (this == &r) {
        return *this;
    }
    errorMessage = r.errorMessage;
    errorRoutine = r.errorRoutine;
    logwriter = new Logger(*(r.logwriter));
    return *this;
}

Application::Messages::~Messages()
{
    delete logwriter;
}

void Application::Messages::addError(const std::string& r, const std::string& msg)
{
    errorMessage.push_back(msg);
    errorRoutine.push_back(r);
}

int Application::Messages::getErrorCount()
{
    return static_cast<int>(errorMessage.size()) ;
}

void Application::Messages::setLogger(Logger* _logwriter)
{
    if (logwriter == _logwriter) {
        return ;
    }
    if (logwriter != 0) {
        delete logwriter;
        logwriter = 0 ;
    }
    logwriter = _logwriter;
}

void Application::Messages::logerror(const std::string& msg)
{
    Cantera::warn_deprecated("Application::Messages::logerror");
    logwriter->error(msg) ;
}

void Application::Messages::writelog(const std::string& msg)
{
    logwriter->write(msg);
}

void Application::Messages::writelogendl()
{
    logwriter->writeendl();
}

#ifdef THREAD_SAFE_CANTERA

//! Mutex for access to string messages
static mutex_t msg_mutex;

Application::Messages* Application::ThreadMessages::operator ->()
{
    ScopedLock msgLock(msg_mutex);
    cthreadId_t curId = getThisThreadId() ;
    threadMsgMap_t::iterator iter = m_threadMsgMap.find(curId) ;
    if (iter != m_threadMsgMap.end()) {
        return iter->second.get();
    }
    pMessages_t pMsgs(new Messages()) ;
    m_threadMsgMap.insert(std::pair< cthreadId_t, pMessages_t >(curId, pMsgs)) ;
    return pMsgs.get() ;
}

void Application::ThreadMessages::removeThreadMessages()
{
    ScopedLock msgLock(msg_mutex);
    cthreadId_t curId = getThisThreadId() ;
    threadMsgMap_t::iterator iter = m_threadMsgMap.find(curId) ;
    if (iter != m_threadMsgMap.end()) {
        m_threadMsgMap.erase(iter) ;
    }
}
#endif // THREAD_SAFE_CANTERA

Application::Application() :
    stop_on_error(false),
    m_suppress_deprecation_warnings(false)
{
#if !defined( THREAD_SAFE_CANTERA )
    pMessenger = std::auto_ptr<Messages>(new Messages());
#endif

    // install a default logwriter that writes to standard
    // output / standard error
    setDefaultDirectories();
#if defined(THREAD_SAFE_CANTERA)
    Unit::units() ;
#endif
}

Application* Application::Instance()
{
    ScopedLock appLock(app_mutex);
    if (Application::s_app == 0) {
        Application::s_app = new Application();
    }
    return s_app;
}

Application::~Application()
{
    std::map<std::string, std::pair<XML_Node*, int> >::iterator pos;
    for (pos = xmlfiles.begin(); pos != xmlfiles.end(); ++pos) {
        pos->second.first->unlock();
        delete pos->second.first;
        pos->second.first = 0;
    }
}

void Application::ApplicationDestroy()
{
    ScopedLock appLock(app_mutex);
    if (Application::s_app != 0) {
        delete Application::s_app;
        Application::s_app = 0;
    }
}

void Application::warn_deprecated(const std::string& method,
                                  const std::string& extra)
{
    if (m_suppress_deprecation_warnings || warnings.count(method)) {
        return;
    }
    warnings.insert(method);
    writelog("WARNING: '" + method + "' is deprecated. " + extra);
    writelogendl();
}

void Application::thread_complete()
{
#if defined(THREAD_SAFE_CANTERA)
    pMessenger.removeThreadMessages() ;
#endif
}

XML_Node* Application::get_XML_File(const std::string& file, int debug)
{
    ScopedLock xmlLock(xml_mutex);
    std::string path = "";
    path = findInputFile(file);
    int mtime = get_modified_time(path);

    if (xmlfiles.find(path) != xmlfiles.end()) {
        // Already have a parsed XML tree for this file cached. Check the
        // last-modified time.
        std::pair<XML_Node*, int> cache = xmlfiles[path];
        if (cache.second == mtime) {
            return cache.first;
        }
    }
    /*
     * Check whether or not the file is XML (based on the file extension). If
     * not, it will be first processed with the preprocessor.
     */
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
        x->build(phase_xml);
    } else {
        std::ifstream s(path.c_str());
        if (s) {
            x->build(s);
        } else {
            throw CanteraError("get_XML_File",
                "cannot open "+file+" for reading.\n"
                "Note, this error indicates a possible configuration problem.");
        }
    }
    x->lock();
    xmlfiles[path] = std::make_pair(x, mtime);
    return x;
}

XML_Node* Application::get_XML_from_string(const std::string& text)
{
    ScopedLock xmlLock(xml_mutex);
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
    entry.first->build(s);
    return entry.first;
}

void Application::close_XML_File(const std::string& file)
{
    ScopedLock xmlLock(xml_mutex);
    if (file == "all") {
        std::map<string, std::pair<XML_Node*, int> >::iterator
        b = xmlfiles.begin(),
        e = xmlfiles.end();
        for (; b != e; ++b) {
            b->second.first->unlock();
            delete b->second.first;
            xmlfiles.erase(b->first);
        }
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
        errorRoutine.pop_back() ;
        errorMessage.pop_back() ;
    }
}

std::string Application::Messages::lastErrorMessage()
{
    if (!errorMessage.empty()) {
        string head =
            "\n\n************************************************\n"
            "                Cantera Error!                  \n"
            "************************************************\n\n";
        return head+string("\nProcedure: ")+errorRoutine.back()
               +string("\nError:   ")+errorMessage.back();
    } else  {
        return "<no Cantera error>";
    }
}

void Application::Messages::getErrors(std::ostream& f)
{
    size_t i = errorMessage.size();
    if (i == 0) {
        return;
    }
    f << endl << endl;
    f << "************************************************" << endl;
    f << "                   Cantera Error!                  " << endl;
    f << "************************************************" << endl
      << endl;
    for (size_t j = 0; j < i; j++) {
        f << endl;
        f << "Procedure: " << errorRoutine[j] << endl;
        f << "Error:     " << errorMessage[j] << endl;
    }
    f << endl << endl;
    errorMessage.clear();
    errorRoutine.clear();
}

void Application::Messages::logErrors()
{
    size_t i = errorMessage.size();
    if (i == 0) {
        return;
    }
    writelog("\n\n");
    writelog("************************************************\n");
    writelog("                   Cantera Error!                  \n");
    writelog("************************************************\n\n");
    for (size_t j = 0; j < i; j++) {
        writelog("\n");
        writelog(string("Procedure: ")+ errorRoutine[j]+" \n");
        writelog(string("Error:     ")+ errorMessage[j]+" \n");
    }
    writelog("\n\n");
    errorMessage.clear();
    errorRoutine.clear();
}

void Application::setDefaultDirectories()
{
    std::vector<string>& dirs = inputDirs;

    // always look in the local directory first
    dirs.push_back(".");

#ifdef _WIN32
    // Under Windows, the Cantera setup utility records the installation
    // directory in the registry. Data files are stored in the 'data'
    // subdirectory of the main installation directory.

    std::string installDir;
    readStringRegistryKey("SOFTWARE\\Cantera\\Cantera 2.2",
                          "InstallDir", installDir, "");
    if (installDir != "") {
        dirs.push_back(installDir + "data");

        // Scripts for converting mechanisms to CTI and CMTL are installed in
        // the 'bin' subdirectory. Add that directory to the PYTHONPATH.
        const char* old_pythonpath = getenv("PYTHONPATH");
        std::string pythonpath = "PYTHONPATH=" + installDir + "\\bin";
        if (old_pythonpath) {
            pythonpath += ";";
            pythonpath.append(old_pythonpath);
        }
        putenv(pythonpath.c_str());
    }

#endif

#ifdef DARWIN
    //
    // add a default data location for Mac OS X
    //
    dirs.push_back("/Applications/Cantera/data");
#endif

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
            dirs.push_back(s.substr(start, end-start));
            start = end + 1;
            end = s.find(pathsep, start);
        }
        dirs.push_back(s.substr(start,end));
    }

    // CANTERA_DATA is defined in file config.h. This file is written
    // during the build process (unix), and points to the directory
    // specified by the 'prefix' option to 'configure', or else to
    // /usr/local/cantera.
#ifdef CANTERA_DATA
    string datadir = string(CANTERA_DATA);
    dirs.push_back(datadir);
#endif
}

void Application::addDataDirectory(const std::string& dir)
{
    ScopedLock dirLock(dir_mutex);
    if (inputDirs.empty()) {
        setDefaultDirectories();
    }
    string d = stripnonprint(dir);

    // Remove any existing entry for this directory
    std::vector<string>::iterator iter = std::find(inputDirs.begin(),
                                                   inputDirs.end(), d);
    if (iter != inputDirs.end()) {
        inputDirs.erase(iter);
    }

    // Insert this directory at the beginning of the search path
    inputDirs.insert(inputDirs.begin(), d);
}

std::string Application::findInputFile(const std::string& name)
{
    ScopedLock dirLock(dir_mutex);
    string::size_type islash = name.find('/');
    string::size_type ibslash = name.find('\\');
    string inname;
    std::vector<string>& dirs = inputDirs;

    // Expand "~/" to user's home directory, if possible
    if (name.find("~/") == 0) {
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
        inname = "";
        for (size_t i = 0; i < nd; i++) {
            inname = dirs[i] + "/" + name;
            std::ifstream fin(inname.c_str());
            if (fin) {
                fin.close();
                return inname;
            }
        }
        string msg;
        msg = "\nInput file " + name
              + " not found in director";
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
        throw CanteraError("findInputFile", msg);
    }

    return name;
}

Application* Application::s_app = 0;

} // namespace Cantera
