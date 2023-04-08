//! @file application.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BASE_APPLICATION_H
#define CT_BASE_APPLICATION_H

#include "cantera/base/config.h"
#include "cantera/base/logger.h"

#include <boost/algorithm/string/join.hpp>

#include <thread>

namespace Cantera
{

/*!
 * @defgroup globalData Global Data
 *
 * Global data are available anywhere. There are two kinds. Cantera has an
 * assortment of constant values for physical parameters. Also, Cantera
 * maintains a collection of global data which is specific to each process that
 * invokes Cantera functions. This process-specific data is stored in the class
 * Application.
 */


//!  Class to hold global data.
/*!
 * Class Application is the top-level class that stores data that should persist
 * for the duration of the process. The class should not be instantiated
 * directly; instead, it is instantiated as needed by the functions declared
 * here. At most one instance is created, and it is not destroyed until the
 * process terminates.
 *
 * @ingroup textlogs
 * @ingroup globalData
 */
class Application
{
protected:
    //! Class to carry out messages
    class Messages
    {
    public:
        Messages();

        Messages(const Messages& r) = delete;
        Messages& operator=(const Messages& r) = delete;

        //! Set an error condition in the application class without
        //! throwing an exception.
        /*!
         * This routine adds an error message to the end of the stack of errors
         * that Cantera accumulates in the Application class.
         * @param r    Procedure name which is generating the error condition
         * @param msg  Descriptive message of the error condition.
         *
         * If only one argument is specified, that string is used as the
         * entire message.
         * @ingroup errorhandling
         */
        void addError(const std::string& r, const std::string& msg="");

        //! Return the number of errors that have been encountered so far.
        /*!
         * @ingroup errorhandling
         */
        int getErrorCount();

        //! Discard the last error message
        /*!
         * %Cantera saves a stack of exceptions that it has caught in the
         * Application class. This routine eliminates the last exception to be
         * added to that stack.
         *
         * @ingroup errorhandling
         */
        void popError();

        //! Retrieve the last error message in a string
        /*!
         * This routine will retrieve the last error message and return it in
         * the return string.
         *
         * @ingroup errorhandling
         */
        std::string lastErrorMessage();

        //!  Prints all of the error messages to an ostream
        /*!
         * Write out all of the saved error messages to the ostream f using
         * the function Logger::writelog. Cantera saves a stack of exceptions
         * that it has caught in the Application class. This routine writes
         * out all of the error messages to the ostream and then clears them
         * from internal storage.
         *
         * @param f ostream which will receive the error messages
         *
         * @ingroup errorhandling
         */
        void getErrors(std::ostream& f);

        //!  Prints all of the error messages using writelog
        /*!
         * Print all of the error messages using function writelog. Cantera
         * saves a stack of exceptions that it has caught in the Application
         * class. This routine writes out all of the error messages and then
         * clears them from internal storage.
         *
         * @ingroup errorhandling
         */
        void logErrors();

        //!  Write a message to the screen.
        /*!
         * The string may be of any length, and may contain end-of-line
         * characters. This method is used throughout Cantera to write log
         * messages. It can also be called by user programs.  The advantage of
         * using writelog over writing directly to the standard output is that
         * messages written with writelog will display correctly even when
         * Cantera is used from MATLAB or other application that do not have a
         * standard output stream.
         *
         * @param msg  c++ string to be written to the screen
         * @ingroup textlogs
         */
        void writelog(const std::string& msg);

        //! Write an end of line character to the screen and flush output
        void writelogendl();

        //!  Write a warning message to the screen.
        /*!
         * @param warning  String specifying type of warning; see Logger::warn()
         * @param msg  String to be written to the screen
         * @ingroup textlogs
         */
        void warnlog(const std::string& warning, const std::string& msg);

        //! Install a logger.
        /*!
         * Called by the language interfaces to install an appropriate logger.
         * The logger is used for the writelog() function
         *
         * @param logwriter Pointer to a logger object
         * @see Logger.
         * @ingroup textlogs
         */
        void setLogger(Logger* logwriter);

    protected:
        //! Current list of error messages
        std::vector<std::string> errorMessage;

        //! Current pointer to the logwriter
        std::unique_ptr<Logger> logwriter;
    };

    //! Typedef for thread specific messages
    typedef shared_ptr<Messages> pMessages_t;

    //! Class that stores thread messages for each thread, and retrieves them
    //! based on the thread id.
    class ThreadMessages
    {
    public:
        //! Constructor
        ThreadMessages() {}

        //! Provide a pointer dereferencing overloaded operator
        /*!
         * @returns a pointer to Messages
         */
        Messages* operator->();

        //! Remove a local thread message
        void removeThreadMessages();

        //! Typedef for map between a thread and the message
        typedef std::map<std::thread::id, pMessages_t> threadMsgMap_t;

    private:
        //! Thread Msg Map
        threadMsgMap_t m_threadMsgMap;
    };

protected:
    //! Constructor for class sets up the initial conditions
    //! Protected ctor access thru static member function Instance
    Application();

public:
    //! Return a pointer to the one and only instance of class Application
    /*!
     * If the Application object has not yet been created, it is created
     */
    static Application* Instance();

    //! Destructor for class deletes global data
    virtual ~Application() {}

    //! Static function that destroys the application class's data
    static void ApplicationDestroy();

    //! @copydoc Messages::addError
    void addError(const std::string& r, const std::string& msg="") {
        pMessenger->addError(r, msg);
    }

    //! @copydoc Messages::getErrorCount
    int getErrorCount() {
        return pMessenger->getErrorCount();
    }

    //! @copydoc Messages::popError
    void popError() {
        pMessenger->popError();
    }

    //! @copydoc Messages::lastErrorMessage
    std::string lastErrorMessage() {
        return pMessenger->lastErrorMessage();
    }

    //! @copydoc Messages::getErrors
    void getErrors(std::ostream& f) {
        pMessenger->getErrors(f);
    }

    //! @copydoc Messages::logErrors
    void logErrors() {
        pMessenger->logErrors();
    }

    //!  Add a directory to the data file search path.
    /*!
     * @ingroup inputfiles
     *
     * @param dir  String name for the directory to be added to the search path
     */
    void addDataDirectory(const std::string& dir);

    //! Find an input file.
    /*!
     * This routine will search for a file in the default locations specified
     * for the application. See the routine setDefaultDirectories() listed
     * above. The first directory searched is usually the current working
     * directory.
     *
     * The default set of directories will not be searched if an absolute path
     * (for example, one starting with `/` or `C:\`) or a path relative to the
     * user's home directory (for example, starting with `~/`) is specified.
     *
     * The presence of the file is determined by whether the file can be
     * opened for reading by the current user.
     *
     * @param name Name of the input file to be searched for
     * @return  The absolute path name of the first matching file
     *
     * If the file is not found a CanteraError exception is thrown.
     *
     * @ingroup inputfiles
     */
    std::string findInputFile(const std::string& name);

    //! Get the Cantera data directories
    /*!
     * This routine returns a string including the names of all the
     * directories searched by Cantera for data files.
     *
     * @param sep Separator to use between directories in the string
     * @return A string of directories separated by the input sep
     *
     * @ingroup inputfiles
     */
    std::string getDataDirectories(const std::string& sep) {
        return boost::algorithm::join(inputDirs, sep);
    }

    //! Load an extension implementing user-defined models
    //! @param extType Specifies the interface / language of the extension, for example
    //!     "python"
    //! @param name Specifies the name of the extension. The meaning of this
    //!     parameter depends on the specific extension interface. For example, for
    //!     Python extensions, this is the name of the Python module containing the
    //!     models.
    //! @since New in Cantera 3.0
    void loadExtension(const std::string& extType, const std::string& name);

    //! Set the versions of Python to try when loading user-defined extensions,
    //! in order of preference. Separate multiple versions with commas, for example
    //! `"3.11,3.10"`.
    //! @since New in Cantera 3.0
    void searchPythonVersions(const string& versions);

#ifdef _WIN32
    long int readStringRegistryKey(const std::string& keyName, const std::string& valueName,
                                   std::string& value, const std::string& defaultValue);
#endif

    //! @copydoc Messages::writelog
    void writelog(const std::string& msg) {
        pMessenger->writelog(msg);
    }

    //! Write an endl to the screen and flush output
    void writelogendl() {
        pMessenger->writelogendl();
    }

    //! @copydoc Messages::warnlog
    void warnlog(const std::string& warning, const std::string& msg) {
        pMessenger->warnlog(warning, msg);
    }

    //! Print a warning indicating that *method* is deprecated. Additional
    //! information (removal version, alternatives) can be specified in
    //! *extra*. Deprecation warnings are printed once per method per
    //! invocation of the application.
    void warn_deprecated(const std::string& method, const std::string& extra="");

    //! Globally disable printing of deprecation warnings. Used primarily to
    //! prevent certain tests from failing.
    void suppress_deprecation_warnings() {
        m_suppress_deprecation_warnings = true;
        m_fatal_deprecation_warnings = false;
    }

    //! Turns deprecation warnings into exceptions. Activated within the test
    //! suite to make sure that no deprecated methods are being used.
    void make_deprecation_warnings_fatal() {
        m_fatal_deprecation_warnings = true;
    }

    //! Generate a general purpose warning; repeated warnings are not suppressed
    //! @param warning  Warning type; see Logger::warn()
    //! @param method  Name of method triggering the warning
    //! @param extra  Additional information printed for the warning
    void warn(const std::string& warning,
              const std::string& method, const std::string& extra="");

    //! Globally disable printing of (user) warnings. Used primarily to
    //! prevent certain tests from failing.
    void suppress_warnings() {
        m_suppress_warnings = true;
        m_fatal_warnings = false;
    }

    //! Returns `true` if warnings should be suppressed.
    bool warnings_suppressed() {
        return m_suppress_warnings;
    }

    //! Turns Cantera warnings into exceptions. Activated within the test
    //! suite to make sure that your warning message are being raised.
    void make_warnings_fatal() {
        m_fatal_warnings = true;
    }

    //! Globally disable printing of warnings about problematic thermo data,
    //! such as NASA polynomials with discontinuities at the midpoint temperature.
    void suppress_thermo_warnings(bool suppress=true) {
        m_suppress_thermo_warnings = suppress;
    }

    //! Returns `true` if thermo warnings should be suppressed.
    bool thermo_warnings_suppressed() {
        return m_suppress_thermo_warnings;
    }

    //! Set definition used for rate constant calculation.
    //! @see Kinetics::getFwdRateConstants()
    /*!
     * If set to 'false' (default value), rate constants of three-body reactions are
     * consistent with conventional definitions (for example Eq. 9.75 in Kee, Coltrin
     * and Glarborg, 'Chemically Reacting Flow', Wiley Interscience, 2003). If set to
     * 'true', output for rate constants of three-body reactions is multiplied by
     * third-body concentrations, consistent with Cantera's behavior prior to version 3.0.
     */
    void use_legacy_rate_constants(bool legacy=true) {
        m_use_legacy_rate_constants = legacy;
    }

    //! Returns `true` if legacy rate constant definition should be used
    bool legacy_rate_constants_used() {
        return m_use_legacy_rate_constants;
    }

    //! @copydoc Messages::setLogger
    void setLogger(Logger* logwriter) {
        pMessenger->setLogger(logwriter);
    }

    //! Delete and free memory allocated per thread in multithreaded applications
    /*!
     * Delete the memory allocated per thread by Cantera.  It should be called
     * from within the thread just before the thread terminates.  If your
     * version of Cantera has not been specifically compiled for thread safety
     * this function does nothing.
     */
    void thread_complete();

protected:
    //! Set the default directories for input files.
    /*!
     * %Cantera searches for input files along a path that includes platform-
     * specific default locations, and possibly user-specified locations.
     * This function installs the platform-specific directories on the search
     * path. It is invoked at startup by appinit(), and never should need to
     * be called by user programs.
     *
     * The current directory (".") is always searched first. Then, on Windows, the
     * registry is checked to find the Cantera installation directory, and the
     * 'data' subdirectory of the installation directory will be added to the search
     * path.
     *
     * On any platform, if environment variable CANTERA_DATA is set to a directory
     * name or a list of directory names separated with the OS-dependent path
     * separator (that is, ";" on Windows, ":" elsewhere), then these directories will
     * be added to the search path.
     *
     * Finally, the location where the data files were installed when
     * %Cantera was built is added to the search path.
     *
     * Additional directories may be added by calling function addDirectory.
     * @ingroup inputfiles
     */
    void setDefaultDirectories();

    //! Current vector of input directories to search for input files
    std::vector<std::string> inputDirs;

    //! Versions of Python to consider when attempting to load user extensions
    vector<string> m_pythonSearchVersions = {"3.11", "3.10", "3.9", "3.8"};

    //! Vector of deprecation warnings that have been emitted (to suppress
    //! duplicates)
    std::set<std::string> warnings;

    bool m_suppress_deprecation_warnings = false;
    bool m_fatal_deprecation_warnings = false;
    bool m_suppress_thermo_warnings = false;
    bool m_suppress_warnings = false;
    bool m_fatal_warnings = false;
    bool m_use_legacy_rate_constants = false;

    std::set<std::pair<std::string, std::string>> m_loaded_extensions;

    ThreadMessages pMessenger;

private:
    //! Pointer to the single Application instance
    static Application* s_app;
};

}

#endif
