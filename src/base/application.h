#ifndef CT_BASE_APPLICATION_H
#define CT_BASE_APPLICATION_H

#include "cantera/base/config.h"
#include "cantera/base/ct_thread.h"
#include "cantera/base/logger.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

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
 * @ingroup HTML_logs
 * @ingroup textlogs
 * @ingroup globalData
 */
class Application
{
protected:
    //! Class to carry out messages
    /*!
     * @ingroup HTML_logs
     */
    class Messages
    {
    public:
        //! Constructor for the Messages class
        /*! Constructor for the Messages class which is a subclass
         *  of the Application class.
         */
        Messages();

        //! Copy Constructor for the Messages class
        /*! Constructor for the Messages class which is a subclass
         * of the Application class.
         *  @param r   Message to be copied
         */
        Messages(const Messages& r);

        //! Assignment operator
        /*!
         *  @param r   Message to be copied
         */
        Messages& operator=(const Messages& r);

        //! Destructor for the Messages class
        ~Messages();

        //! Set an error condition in the application class without
        //! throwing an exception.
        /*!
         * This routine adds an error message to the end of the stack
         * of errors that Cantera accumulates in the Application
         * class.
         * @param r   location
         * @param msg  Description of the error
         * @ingroup errorhandling
         */
        void addError(std::string r, std::string msg);

        //! Return the number of errors that have been encountered so far.
        /*!
         * @ingroup errorhandling
         */
        int getErrorCount() ;

        //! Discard the last error message
        /*!
         * %Cantera saves a stack of exceptions that it
         * has caught in the Application class. This routine eliminates
         * the last exception to be added to that stack.
         *
         * @ingroup errorhandling
         */
        void popError() ;

        //! Retrieve the last error message in a string
        /*!
         * This routine will retrieve the last error message and return
         * it in the return string.
         *
         * @ingroup errorhandling
         */
        std::string lastErrorMessage() ;

        //!  Prints all of the error messages to an ostream
        /*!
         * Print all of the error messages using function writelog.
         * Write out all of the saved error messages to the ostream f
         * Cantera saves a stack of exceptions that it
         * has caught in the Application class. This routine writes
         * out all of the error messages to the ostream
         * and then clears them from internal storage.
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
         * The string may be of any
         * length, and may contain end-of-line characters. This method is
         * used throughout %Cantera to write log messages.
         *
         * @param msg  c++ string to be written to the screen
         * @ingroup textlogs
         */
        void writelog(const std::string& msg);

        //! Write an end of line and flush output
        void writelogendl();

        //!  Write a message to the screen.
        /*!
         * The string may be of any
         * length, and may contain end-of-line characters. This method is
         * used throughout %Cantera to write log messages.
         *
         * @param pszmsg  c character string to be written to the screen
         * @ingroup textlogs
         */
        void writelog(const char* pszmsg) ;

        //! Write an error message and quit.
        /*!
         *  The default behavior is
         *  to write to the standard error stream, and then call
         *  exit(). Note that no end-of-line character is appended to
         *  the message, and so if one is desired it must be included
         *  in the string. Note that this default behavior will
         *  terminate the application Cantera is invoked from (MATLAB,
         *  Excel, etc.) If this is not desired, then derive a class
         *  and reimplement this method.
         *
         * @param msg    Error message to be written to cerr.
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

#ifdef WITH_HTML_LOGS

        //!Create a new group for log messages.
        /*!
         *  Usually this is called
         *  upon entering the function, with the title parameter equal to
         *  the name of the function or method. Subsequent messages
         *  written with addLogEntry will appear grouped under this
         *  heading, until endLogGroup() is called.
         *
         *  @param title String name of the LogGroup
         *  @param loglevel loglevel of the group.
         *  @ingroup HTML_logs
         */
        void beginLogGroup(std::string title, int loglevel) ;

        //! Add an entry to an HTML log file.
        /*!
         *  Entries appear in the form "tag:value".
         *
         * @param tag      tag
         * @param value    double value
         *
         * @ingroup HTML_logs
         */
        void addLogEntry(std::string tag, std::string value) ;

        //! Add an entry to an HTML log file.
        /*!
         *  Entries appear in the form "tag:value".
         *
         * @param tag      tag
         * @param value    double value
         *
         * @ingroup HTML_logs
         */
        void addLogEntry(std::string tag, doublereal value) ;

        //! Add an entry to an HTML log file.
        /*!
         *  Entries appear in the form "tag:value".
         *
         * @param tag      tag
         * @param value    double value
         *
         * @ingroup HTML_logs
         */
        void addLogEntry(std::string tag, int value) ;

        //! Add an entry to an HTML log file.
        /*!
         *  Entries appear in the form "msg".
         *
         * @param msg  Message to be added
         *
         * @ingroup HTML_logs
         */
        void addLogEntry(std::string msg) ;

        //! Close the current group of log messages.
        /*!
         *  This is typically
         *  called just before leaving a function or method, to close the
         *  group of messages that were output from this
         *  function. Subsequent messages written with addLogEntry() will
         *  appear at the next-higher level in the outline, unless
         *  beginLogGroup() is called first to create a new group.
         *
         * @param title Name of the log group. It defaults to the most recent
         *              log group created.
         */
        void endLogGroup(std::string title) ;

        //! Write the HTML log file.
        /*!
         *  Log entries are stored in memory in
         *  an XML tree until this function is called, which writes the
         *  tree to a file and clears the entries stored in memory.  The
         *  output file will have the name specified in the 'file'
         *  argument.  If this argument has no extension, the extension
         *  '.html' will be appended. Also, if the file already exists, an
         *  integer will be appended to the name so that no existing log
         *  file will be overwritten.
         *  WITH_HTML_LOGS must be defined.
         *
         *  @param  file Name of the file to be written
         */
        void write_logfile(std::string file);
#endif

    protected:
        //! Current list of error messages
        std::vector<std::string> errorMessage;

        //! Current error Routine
        std::vector<std::string> errorRoutine;

        //! Current pointer to the logwriter
        Logger* logwriter;
#ifdef WITH_HTML_LOGS
        //! Current pointer to the top of the XML_Node tree for the current HTML log
        XML_Node* xmllog;

        //! Pointer to the last current position in the XML_Node tree for the current HTML log
        XML_Node* current;

        //! Current value of the loglevel
        int loglevel;

        //! Vector of loglevels for loggroups that are open
        std::vector<int> loglevels;

        //! Current vector of loggroups that are open
        std::vector<std::string> loggroups;
#endif
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

    //! Set an error condition in the application class without
    //! throwing an exception.
    /*!
     * This routine adds an error message to the end of the stack
     * of errors that Cantera accumulates in the Application
     * class.
     * @param r   location
     * @param msg  Description of the error
     * @ingroup errorhandling
     */
    void addError(std::string r, std::string msg) {
        pMessenger->addError(r, msg) ;
    }

    //! Return the number of errors that have been encountered so far.
    /*!
     * @ingroup error handling
     */
    int getErrorCount() {
        return pMessenger->getErrorCount() ;
    }

    //! Discard the last error message
    /*!
     * %Cantera saves a stack of exceptions that it
     * has caught in the Application class. This routine eliminates
     * the last exception to be added to that stack.
     *
     * @ingroup errorhandling
     */
    void popError() {
        pMessenger->popError() ;
    }

    //! Retrieve the last error message in a string
    /*!
     * This routine will retrieve the last error message and return
     * it in the return string.
     *
     * @ingroup errorhandling
     */
    std::string lastErrorMessage() {
        return pMessenger->lastErrorMessage() ;
    }

    //!  Prints all of the error messages to an ostream
    /*!
     * Print all of the error messages using function writelog.
     * Write out all of the saved error messages to the ostream f
     * Cantera saves a stack of exceptions that it
     * has caught in the Application class. This routine writes
     * out all of the error messages to the ostream
     * and then clears them from internal storage.
     *
     * @param f ostream which will receive the error messages
     *
     * @ingroup errorhandling
     */
    void getErrors(std::ostream& f) {
        pMessenger->getErrors(f) ;
    }

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
    void logErrors() {
        pMessenger->logErrors() ;
    }

    //!  Add a directory to the data file search path.
    /*!
     * @ingroup inputfiles
     *
     * @param dir  String name for the directory to be added to the search path
     */
    void addDataDirectory(std::string dir) ;

    //! Find an input file.
    /*!
     *    This routine will search for a file in the default
     *    locations specified for the application.
     *    See the routine setDefaultDirectories() listed above.
     *
     *    The default set of directories specified for the application
     *    will be searched if a '/' or an '\\' is found in the
     *    name. If either is found then a relative path name is
     *    presumed, and the default directories are not searched.
     *
     *    The presence of the file is determined by whether the file
     *    can be opened for reading by the current user.
     *
     *  @param name Name of the input file to be searched for
     *
     *    @return
     *
     *      The absolute path name of the first matching
     *      file is returned. If a relative path name
     *      is indicated, the relative path name is returned.
     *
     *      If the file is not found, a message is written to
     *      stdout and  a CanteraError exception is thrown.
     *
     * @ingroup inputfiles
     */
    std::string findInputFile(std::string name) ;

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
    XML_Node* get_XML_File(std::string file, int debug=0) ;

    //! Close an XML File
    /*!
     * Close a file that is opened by this application object
     *
     * @param file String containing the relative or absolute file name
     */
    void close_XML_File(std::string file) ;

    //!  Write a message to the screen.
    /*!
     * The string may be of any
     * length, and may contain end-of-line characters. This method is
     * used throughout %Cantera to write log messages.
     *
     * @param msg  c++ string to be written to the screen
     * @ingroup textlogs
     */
#ifdef _WIN32
    long int readStringRegistryKey(const std::string& keyName, const std::string& valueName,
                                   std::string& value, const std::string& defaultValue);
#endif

    void writelog(const std::string& msg) {
        pMessenger->writelog(msg);
    }


    //! Write an endl to the screen and flush output
    /*!
     * @ingroup textlogs
     */
    void writelogendl() {
        pMessenger->writelogendl();
    }

    //!  Write a message to the screen.
    /*!
     * The string may be of any
     * length, and may contain end-of-line characters. This method is
     * used throughout %Cantera to write log messages.
     *
     * @param pszmsg  c null terminated string to be written to the screen
     * @ingroup textlogs
     */
    void writelog(const char* pszmsg) {
        pMessenger->writelog(pszmsg);
    }

    //! Write an error message and quit.
    /*!
     *  The default behavior is
     *  to write to the standard error stream, and then call
     *  exit(). Note that no end-of-line character is appended to
     *  the message, and so if one is desired it must be included
     *  in the string. Note that this default behavior will
     *  terminate the application Cantera is invoked from (MATLAB,
     *  Excel, etc.) If this is not desired, then derive a class
     *  and reimplement this method.
     *
     * @param msg    Error message to be written to cerr.
     */
    void logerror(const std::string& msg) {
        pMessenger->logerror(msg);
    }

    //!  Install a logger -  Called by the language interfaces to install an
    //!  appropriate logger.
    /*!
     * @param logwriter Pointer to a logger object
     *  @see Logger.
     *  @ingroup textlogs
     */
    void setLogger(Logger* logwriter) {
        pMessenger->setLogger(logwriter);
    }

    //! Delete Messenger object allocated per thread.
    void thread_complete() ;

#ifdef WITH_HTML_LOGS
    //!Create a new group for log messages.
    /*!
     *  Usually this is called
     *  upon entering the function, with the title parameter equal to
     *  the name of the function or method. Subsequent messages
     *  written with addLogEntry will appear grouped under this
     *  heading, until endLogGroup() is called.
     *
     *  @param title String name of the LogGroup
     *  @param loglevel loglevel of the group.
     *  @ingroup HTML_logs
     */
    void beginLogGroup(std::string title, int loglevel) {
        pMessenger->beginLogGroup(title,loglevel);
    }

    //! Add an entry to an HTML log file.
    /*!
     *  Entries appear in the form "tag:value".
     *
     * @param tag      tag
     * @param value    double value
     *
     * @ingroup HTML_logs
     */
    void addLogEntry(std::string tag, std::string value) {
        pMessenger->addLogEntry(tag, value);
    }

    //! Add an entry to an HTML log file.
    /*!
     *  Entries appear in the form "tag:value".
     *
     * @param tag      tag
     * @param value    double value
     *
     * @ingroup HTML_logs
     */
    void addLogEntry(std::string tag, doublereal value) {
        pMessenger->addLogEntry(tag, value);
    }

    //! Add an entry to an HTML log file.
    /*!
     *  Entries appear in the form "tag:value".
     *
     * @param tag      tag
     * @param value    double value
     *
     * @ingroup HTML_logs
     */
    void addLogEntry(std::string tag, int value) {
        pMessenger->addLogEntry(tag, value);
    }

    //! Add an entry to an HTML log file.
    /*!
     *  Entries appear in the form "msg".
     *
     * @param msg      Message to be added to file
     *
     * @ingroup HTML_logs
     */
    void addLogEntry(std::string msg) {
        pMessenger->addLogEntry(msg);
    }

    //! Close the current group of log messages.
    /*!
     *  This is typically
     *  called just before leaving a function or method, to close the
     *  group of messages that were output from this
     *  function. Subsequent messages written with addLogEntry() will
     *  appear at the next-higher level in the outline, unless
     *  beginLogGroup() is called first to create a new group.
     *
     * @param title Name of the log group. It defaults to the most recent
     *              log group created.
     * @ingroup HTML_logs
     */
    void endLogGroup(std::string title) {
        pMessenger->endLogGroup(title) ;
    }

    //! Write the HTML log file.
    /*!
     *  Log entries are stored in memory in
     *  an XML tree until this function is called, which writes the
     *  tree to a file and clears the entries stored in memory.  The
     *  output file will have the name specified in the 'file'
     *  argument.  If this argument has no extension, the extension
     *  '.html' will be appended. Also, if the file already exists, an
     *  integer will be appended to the name so that no existing log
     *  file will be overwritten.
     *  WITH_HTML_LOGS must be defined.
     *
     *  @param  file Name of the file to be written
     */
    void write_logfile(std::string file) {
        pMessenger->write_logfile(file) ;
    }
#endif

protected:
    //! Set the default directories for input files.
    /*!
     * %Cantera searches
     * for input files along a path that includes platform-specific
     * default locations, and possibly user-specified locations.  This
     * function installs the platform-specific directories on the
     * search path. It is invoked at startup by appinit(), and never
     * should need to be called by user programs.
     *
     * The current directory (".") is always searched first. Then, on
     * Windows platforms, if environment variable COMMONPROGRAMFILES
     * is set (which it should be on Win XP or Win 2000), then
     * directories under this one will be added to the search
     * path. The %Cantera Windows installer installs data files to this
     * location.
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
    //! Current vector of xml file trees that have been previously parsed
    std::map<std::string, XML_Node*> xmlfiles;
    //! Current pointer to the logwriter
    //Logger* logwriter;
#ifdef WITH_HTML_LOGS
    //! Current pointer to the top of the XML_Node tree for the current HTML log
    //XML_Node *xmllog;
    //! Pointer to the last current position in the XML_Node tree for the current HTML log
    //XML_Node *current;
    //! Current value of loglevel
    //int loglevel;
    //! Vector of loglevels for loggroups that are open
    //vector<int> loglevels;
    //! Current vector of loggroups that are open
    //vector<string> loggroups;
#endif

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
