//! @file application.h
#ifndef CT_BASE_APPLICATION_H
#define CT_BASE_APPLICATION_H

#include "cantera/base/config.h"
#include "cantera/base/ct_thread.h"
#include "cantera/base/logger.h"

#include <set>
#include <memory>

namespace Cantera
{

class XML_Node;

/*!
 * @defgroup globalData Global Data
 *
 * Global data are available anywhere. There are two kinds.
 * Cantera has an assortment of constant values for physical parameters.
 * Also, Cantera maintains a collection of global data which is specific
 * to each process that invokes Cantera functions. This process-specific
 * data is stored in the class Application.
 */


//!  Class to hold global data.
/*!
 * Class Application is the top-level
 * class that stores data that should persist for the duration of
 * the process. The class should not be instantiated directly;
 * instead, it is instantiated as needed by the functions declared
 * here. At most one instance is created, and it is not destroyed
 * until the process terminates.
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
        //! Constructor for the Messages class
        /*! Constructor for the Messages class which is a subclass
         *  of the Application class.
         */
        Messages();

        Messages(const Messages& r);
        Messages& operator=(const Messages& r);
        ~Messages();

        //! Set an error condition in the application class without
        //! throwing an exception.
        /*!
         * This routine adds an error message to the end of the stack of errors
         * that Cantera accumulates in the Application class.
         * @param r    Procedure name which is generating the error condition
         * @param msg  Descriptive message of the error condition.
         * @ingroup errorhandling
         */
        void addError(const std::string& r, const std::string& msg);

        //! Return the number of errors that have been encountered so far.
        /*!
         * @ingroup errorhandling
         */
        int getErrorCount() ;

        //! Discard the last error message
        /*!
         * %Cantera saves a stack of exceptions that it has caught in the
         * Application class. This routine eliminates the last exception to be
         * added to that stack.
         *
         * @ingroup errorhandling
         */
        void popError() ;

        //! Retrieve the last error message in a string
        /*!
         * This routine will retrieve the last error message and return it in
         * the return string.
         *
         * @ingroup errorhandling
         */
        std::string lastErrorMessage() ;

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
        void getErrors(std::ostream& f) ;

        //!  Prints all of the error messages using writelog
        /*!
         * Print all of the error messages using function writelog.
         * Cantera saves a stack of exceptions that it
         * has caught in the Application class. This routine writes
         * out all of the error messages
         * and then clears them from internal storage.
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

        //! Write an error message and quit.
        /*!
         *  The default behavior is to write to the standard error stream, and
         *  then call exit(). Note that no end-of-line character is appended
         *  to the message, and so if one is desired it must be included in
         *  the string. Note that this default behavior will terminate the
         *  application Cantera is invoked from (MATLAB, Excel, etc.) If this
         *  is not desired, then derive a class and reimplement this method.
         *
         * @param msg    Error message to be written to cerr.
         * @deprecated To be removed after Cantera 2.2
         */
        void logerror(const std::string& msg) ;

        //! Install a logger.
        /*!
         *  Called by the language interfaces to install an appropriate logger.
         *  The logger is used for the writelog() function
         *
         * @param logwriter Pointer to a logger object
         * @see Logger.
         * @ingroup textlogs
         */
        void setLogger(Logger* logwriter) ;

    protected:
        //! Current list of error messages
        std::vector<std::string> errorMessage;

        //! Current error Routine
        std::vector<std::string> errorRoutine;

        //! Current pointer to the logwriter
        Logger* logwriter;
    } ;

#ifdef THREAD_SAFE_CANTERA
    //! Typedef for thread specific messages
    typedef boost::shared_ptr< Messages >   pMessages_t ;

    //! Class that stores thread messages for each thread, and retrieves them
    //! based on the thread id.
    class ThreadMessages
    {
    public:
        //! Constructor
        ThreadMessages() {}

        //! Provide a pointer dereferencing overloaded operator
        /*!
         * @return  returns a pointer to Messages
         */
        Messages* operator->();

        //! Remove a local thread message
        void removeThreadMessages();

        //! Typedef for map between a thread and the message
        typedef std::map< cthreadId_t, pMessages_t > threadMsgMap_t ;

    private:
        //! Thread Msg Map
        threadMsgMap_t   m_threadMsgMap ;
    } ;
#endif


protected:
    //! Constructor for class sets up the initial conditions
    //! Protected ctor access thru static member function Instance
    Application();

public:
    //! Return a pointer to the one and only instance of class Application
    /*!
     * If the an Application object has not yet been created it is created
     */
    static Application* Instance();

    //! Destructor for class deletes global data
    /*!
     * Delete any open XML trees, the logwriter, and
     * the XML log, if any.
     */
    virtual ~Application();

    //! Static function that destroys the application class's data
    static void ApplicationDestroy();

    //! @copydoc Messages::addError
    void addError(const std::string& r, const std::string& msg) {
        pMessenger->addError(r, msg) ;
    }

    //! @copydoc Messages::getErrorCount
    int getErrorCount() {
        return pMessenger->getErrorCount() ;
    }

    //! @copydoc Messages::popError
    void popError() {
        pMessenger->popError() ;
    }

    //! @copydoc Messages::lastErrorMessage
    std::string lastErrorMessage() {
        return pMessenger->lastErrorMessage() ;
    }

    //! @copydoc Messages::getErrors
    void getErrors(std::ostream& f) {
        pMessenger->getErrors(f) ;
    }

    //! @copydoc Messages::logErrors
    void logErrors() {
        pMessenger->logErrors() ;
    }

    //!  Add a directory to the data file search path.
    /*!
     * @ingroup inputfiles
     *
     * @param dir  String name for the directory to be added to the search path
     */
    void addDataDirectory(const std::string& dir) ;

    //! Find an input file.
    /*!
     *  This routine will search for a file in the default locations specified
     *  for the application. See the routine setDefaultDirectories() listed
     *  above.
     *
     *  The default set of directories specified for the application will be
     *  searched if a '/' or an '\\' is found in the name. If either is found
     *  then a relative path name is presumed, and the default directories are
     *  not searched.
     *
     *  The presence of the file is determined by whether the file can be
     *  opened for reading by the current user.
     *
     *  @param name Name of the input file to be searched for
     *
     *  @return  The absolute path name of the first matching file is
     *      returned. If a relative path name is indicated, the relative path
     *      name is returned.
     *
     *  If the file is not found, a message is written to stdout and a
     *  CanteraError exception is thrown.
     *
     * @ingroup inputfiles
     */
    std::string findInputFile(const std::string& name) ;

    //! Return a pointer to the XML tree for a Cantera input file.
    /*!
     *  This routine will find the file and read the XML file into an
     *  XML tree structure. Then, a pointer will be returned. If the
     *  file has already been processed, then just the pointer will
     *  be returned.
     *
     * @param file String containing the relative or absolute file name
     * @param debug Debug flag
     */
    XML_Node* get_XML_File(const std::string& file, int debug=0) ;

    //! Read a CTI or CTML string and fill up an XML tree.
    /*!
     *  Return a pointer to the XML tree corresponding to the specified
     *  CTI or XML string. If the given string has been processed before,
     *  the cached XML tree will be returned. Otherwise, the XML tree
     *  will be generated and stored in the cache.
     *  @param text    CTI or CTML string
     *  @return        Root of the corresponding XML tree
     */
    XML_Node* get_XML_from_string(const std::string& text);

    //! Close an XML File
    /*!
     * Close a file that is opened by this application object
     *
     * @param file String containing the relative or absolute file name
     */
    void close_XML_File(const std::string& file) ;

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

     //! @copydoc Messages::logerror
    void logerror(const std::string& msg) {
        pMessenger->logerror(msg);
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
    void thread_complete() ;

protected:
    //! Set the default directories for input files.
    /*!
     * %Cantera searches for input files along a path that includes platform-
     * specific default locations, and possibly user-specified locations.
     * This function installs the platform-specific directories on the search
     * path. It is invoked at startup by appinit(), and never should need to
     * be called by user programs.
     *
     * The current directory (".") is always searched first. Then, on Windows
     * platforms, if environment variable COMMONPROGRAMFILES is set (which it
     * should be on Win XP or Win 2000), then directories under this one will
     * be added to the search path. The %Cantera Windows installer installs
     * data files to this location.
     *
     * On the Mac, directory '/Applications/Cantera/data' is added to the
     * search path.
     *
     * On any platform, if environment variable CANTERA_DATA is set to a
     * directory name, then this directory is added to the search path.
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
    //! Current list of error messages
    //vector<string> errorMessage;
    //! Current list of warning messages
    //vector<string> warning;
    //! Current error Routine
    //vector<string> errorRoutine;
    //! Last error message
    //string msglog;
    //! Current line length
    // size_t linelen;
    //! Current value of stop_on_error
    bool stop_on_error;
    //! Current map of options
    std::map<std::string, std::string>     options;
    //! Current value of tmp_dir
    std::string tmp_dir;
    //! Current vector of XML file trees that have been previously parsed
    //! The second element of the value is used to store the last-modified time
    //! for the file, to enable change detection.
    std::map<std::string, std::pair<XML_Node*, int> > xmlfiles;
    //! Vector of deprecation warnings that have been emitted (to suppress duplicates)
    std::set<std::string> warnings;

    bool m_suppress_deprecation_warnings;

    //! Current pointer to the logwriter
    //Logger* logwriter;

#if defined(THREAD_SAFE_CANTERA)
    ThreadMessages   pMessenger ;
#else
    std::auto_ptr<Messages> pMessenger ;
#endif

private:
    //! Pointer to the single Application instance
    static Application* s_app ;
};

}

#endif
