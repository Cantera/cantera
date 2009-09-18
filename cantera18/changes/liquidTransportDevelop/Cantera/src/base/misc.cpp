/**
 *  @file misc.cpp
 * This file contains a miscellaneous collection of global data
 * functions and text for the globalData module (see \ref globalData), the
 * logs module (see \ref logs), and the HTML_logs module (see \ref HTML_logs).
 *
 * These modules are defined here:
 *      globalData 
 *      logs
 *      HTML_logs
 */
/*
 *  $Id: misc.cpp,v 1.16 2009/07/13 16:47:40 hkmoffa Exp $
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#pragma warning(disable:4996)
#endif

// RFB:TODO May need OS specifc include to use varargs.h instead
#include <stdarg.h>

#include "global.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "units.h"
#include "xml.h"
#include "ctml.h"
#include "FactoryBase.h"
#include "logger.h"


#include <fstream>
#include <memory>

#ifdef WIN32
#include <algorithm>
#include <functional>
#include <new>
#endif

using namespace std;

// RFB :  If running multiple threads in a cpp application,
// the Application class is the only internal object
// that is single instance with static data.  Synchronize access to those data
// structures.
// Using macros to avoid polluting code with alot of ifdef's

#ifdef THREAD_SAFE_CANTERA

#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
static boost::mutex  dir_mutex;  // For input directory access
static boost::mutex  msg_mutex;  // For access to string messages
static boost::mutex  app_mutex;  // Application state including creating singleton
//static boost::mutex  log_mutex;  // Logger pointer
static boost::mutex  xml_mutex;  // XML file storage

//! Macro for locking input directory access
#define DIR_LOCK() boost::mutex::scoped_lock   d_lock(dir_mutex)

//! Macro for locking access to string messages
#define MSG_LOCK() boost::mutex::scoped_lock   m_lock(msg_mutex)

//! Macro for locking creating singletons in the application state
#define APP_LOCK() boost::mutex::scoped_lock   a_lock(app_mutex)

//! Macro for locking XML file writing
#define XML_LOCK() boost::mutex::scoped_lock   x_lock(xml_mutex)

   #ifdef WITH_HTML_LOGS
      //static boost::mutex  html_mutex; // html logs
      //#define HTML_LOCK() boost::mutex::scoped_lock   h_lock(html_mutex)
   #endif

   #if defined(BOOST_HAS_WINTHREADS)
      #include <windows.h>
      typedef unsigned int cthreadId_t ;
      class thread_equal 
      {
      public:
         bool operator()(cthreadId_t L, cthreadId_t R)
         {
            return L == R ;
         }
      } ;
      cthreadId_t getThisThreadId()
      {
         return ::GetCurrentThreadId() ;
      }
   #elif defined(BOOST_HAS_PTHREADS)

      typedef pthread_t cthreadId_t ;
      class thread_equal 
      {
      public:
         bool operator()(cthreadId_t L, cthreadId_t R)
         {
            return pthread_equal(L, R) ;
         }
      } ;
      cthreadId_t getThisThreadId()
      {
         return pthread_self() ;
      }
   #elif defined(BOOST_HAS_MPTASKS)
      typedef MPTaskID  cthreadId_t ;
      class thread_equal 
      {
      public:
         bool operator()(cthreadId_t L, cthreadId_t R)
         {
            return L == R ;
         }
      } ;
      cthreadId_t getThisThreadId()
      {
         return MPCurrentTaskID() ;
      }
   #endif

#else
   #define DIR_LOCK() 
   #define MSG_LOCK() 
   #define APP_LOCK() 
   //#define LOG_LOCK() 
   #define XML_LOCK() 

   #ifdef WITH_HTML_LOGS
      //#define HTML_LOCK() 
   #endif

#endif

namespace Cantera {

  /*!
   * @defgroup globalData Global Data
   *
   * Global data are available anywhere. There are two kinds.
   * Cantera has an assortment of constant values for physical parameters.
   * Also, Cantera maintains a collection of global data which is specific
   * to each process that invokes Cantera functions. This process-specific
   * data is storred in the class Application.
   */
  //@{
  //@}
    
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
  class Application {

  protected:
    //! Class to carry out messages
    /*!
     * @ingroup HTML_logs
     */
    class Messages {
    public:
      //! Constructor for the Messages class
      /*!
       * constructor for the Messages class which is a subclass
       * of the Application class.
       */
      Messages() {
	// install a default logwriter that writes to standard
	// output / standard error
	logwriter = new Logger();
#ifdef WITH_HTML_LOGS
	xmllog = 0; 
	current = 0;
	loglevel = 0;
#endif
      }

      //! Destructor for the Messages class
      ~Messages() {
	delete logwriter;
#ifdef WITH_HTML_LOGS
	if (xmllog) {
	  write_logfile("orphan");
	}
#endif
      }
      
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
      void addError(std::string r, std::string msg) ;
    
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
       *  to write to the standard eror stream, and then call
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

      //! Returns an integer specifying the application environment
      int getUserEnv() ;

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
      vector<string> errorMessage;
     
      //! Current error Routine
      vector<string> errorRoutine;
   
      //! Current pointer to the logwriter
      Logger* logwriter;
#ifdef WITH_HTML_LOGS
      //! Current pointer to the top of the XML_Node tree for the current HTML log
      XML_Node *xmllog;

      //! Pointer to the last current position in the XML_Node tree for the current HTML log
      XML_Node *current;

      //! Current value of the loglevel
      int loglevel;

      //! Vector of loglevels for loggroups that are open
      vector<int> loglevels;

      //! Current vector of loggroups that are open
      vector<string> loggroups;
#endif
    } ;
    
#ifdef THREAD_SAFE_CANTERA
      //! Typedef for thread specific messages
      typedef boost::shared_ptr< Messages >   pMessages_t ;

      //! Typedef for map between a thread and the message
      typedef std::map< cthreadId_t, pMessages_t > threadMsgMap_t ;

      class ThreadMessages
      {
      public:
	//! Constructor
	ThreadMessages()
	{
	}
	
	//! Provide a pointer deferencing overloaded operator
	/*!
	 * @return  returns a pointer to Messages
	 */
	Messages* operator->()
	{
	  MSG_LOCK() ;
	  cthreadId_t curId = getThisThreadId() ;
	  threadMsgMap_t::iterator iter = m_threadMsgMap.find( curId ) ;
	  if ( iter != m_threadMsgMap.end() )
            {
	      return (iter->second.get()) ;
            }
	  pMessages_t pMsgs( new Messages() ) ;
	  m_threadMsgMap.insert( std::pair< cthreadId_t, pMessages_t >( curId, pMsgs ) ) ;
	  return pMsgs.get() ;
	}

	//! Remove a local thread message
	void removeThreadMessages() {
	  MSG_LOCK() ;
	  cthreadId_t curId = getThisThreadId() ;
	  threadMsgMap_t::iterator iter = m_threadMsgMap.find( curId ) ;
	  if ( iter != m_threadMsgMap.end() ) {
	    m_threadMsgMap.erase( iter ) ;
	  }
	}

      private:
	//! Thread Msg Map
         threadMsgMap_t   m_threadMsgMap ;
      } ;
    #endif


  protected:   //RFB Protected ctor access thru static member function Instance
    //! Constructor for class sets up the initial conditions
    Application() : /*linelen(0),*/ stop_on_error(false),
                                    tmp_dir("."), m_sleep("1")              
    {
#if !defined( THREAD_SAFE_CANTERA )
     pMessenger = std::auto_ptr<Messages>(new Messages());
#endif
      // if TMP or TEMP is set, use it for the temporary
      // directory
      char* ctmpdir = getenv("CANTERA_TMPDIR");
      if (ctmpdir != 0) {
        tmp_dir = string(ctmpdir);
      } else {
        char* tmpdir = getenv("TMP");
        if (tmpdir == 0) 
	  tmpdir = getenv("TEMP");
        if (tmpdir != 0)
 	  tmp_dir = string(tmpdir);
      }

      // if SLEEP is set, use it as the sleep time
      char* sleepstr = getenv("SLEEP");
      if (sleepstr != 0) {
	m_sleep = string(sleepstr);
      }

      // install a default logwriter that writes to standard
      // output / standard error
      //      logwriter = new Logger();
      //#ifdef WITH_HTML_LOGS
      //      // HTML log files
      //      xmllog = 0; 
      //      current = 0;
      //      loglevel = 0;
      //#endif
      setDefaultDirectories();
#if defined(THREAD_SAFE_CANTERA)
      Unit::units() ;
#endif
    }

  public:
    //! Return a pointer to the one and only instance of class Application
    /*
     * If the an Application object has not yet been created it is created
     */
    static Application* Instance() {
      APP_LOCK();
      if (Application::s_app == 0) {
	Application::s_app = new Application();
      }
      return s_app;
    }
    
    //! Destructor for class deletes global data
    /*!
     * Delete any open XML trees, the logwriter, and
     * the XML log, if any.
     */
    virtual ~Application() {
      std::map<std::string, XML_Node*>::iterator pos;
      for (pos = xmlfiles.begin(); pos != xmlfiles.end(); ++pos) {
	pos->second->unlock();
	delete pos->second;
	pos->second = 0;
      }
    }

    //! Static function that destroys the application class's data
    static void ApplicationDestroy() {
       APP_LOCK() ;
       if (Application::s_app != 0) {
	 delete Application::s_app;
	 Application::s_app = 0;
       }
    }

    //! Set the directory where temporary files may be written.
    /*!
     * @param tmp Name of the directory
     */
    void setTmpDir(std::string tmp) ;

    //! get the name of a directory where temporary files may be written.
    std::string getTmpDir() ;

    //! Pause the application for a certain amount
    /*!
     * This has been needed in the past for some network file systems.
     */
    std::string sleep() ;

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
    void addError(std::string r, std::string msg) { pMessenger->addError(r, msg) ; }

    //! Return the number of errors that have been encountered so far.
    /*!
     * @ingroup error handling
     */
    int getErrorCount() { return pMessenger->getErrorCount() ; }

    //! Discard the last error message
    /*!
     * %Cantera saves a stack of exceptions that it
     * has caught in the Application class. This routine eliminates
     * the last exception to be added to that stack.
     *
     * @ingroup errorhandling
     */
    void popError() { pMessenger->popError() ; }

    //! Retrieve the last error message in a string
    /*!
     * This routine will retrieve the last error message and return
     * it in the return string.
     *
     * @ingroup errorhandling
     */
    std::string lastErrorMessage() { return pMessenger->lastErrorMessage() ; }

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
    void getErrors( std::ostream& f) { pMessenger->getErrors(f) ; }

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
    void logErrors() { pMessenger->logErrors() ; }

    //!  Add a directory to the data file search path.
    /*!
     * @ingroup inputfiles
     *
     * @param dir  String name for the directory to be added to the search path
     */
    void addDataDirectory( std::string dir ) ;
    
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
    void writelog(const std::string& msg) { pMessenger->writelog(msg); }


    //! Write an endl to the screen and flush output
    /*!
     * @ingroup textlogs
     */
    void writelogendl() { pMessenger->writelogendl(); }

    //!  Write a message to the screen.
    /*!
     * The string may be of any
     * length, and may contain end-of-line characters. This method is
     * used throughout %Cantera to write log messages.
     *
     * @param pszmsg  c null terminated string to be written to the screen
     * @ingroup textlogs
     */
    void writelog(const char* pszmsg) { pMessenger->writelog(pszmsg); }

    //! Write an error message and quit.
    /*!
     *  The default behavior is
     *  to write to the standard eror stream, and then call
     *  exit(). Note that no end-of-line character is appended to
     *  the message, and so if one is desired it must be included
     *  in the string. Note that this default behavior will
     *  terminate the application Cantera is invoked from (MATLAB,
     *  Excel, etc.) If this is not desired, then derive a class
     *  and reimplement this method.
     *
     * @param msg    Error message to be written to cerr.
     */
    void logerror(const std::string& msg) {pMessenger->logerror(msg); }

    //! Returns an integer specifying the application environment.
    int getUserEnv() { return pMessenger->getUserEnv() ; }

    //!  Install a logger -  Called by the language interfaces to install an
    //!  appropriate logger. 
    /*!
     * @param logwriter Pointer to a logger object
     *  @see Logger.
     *  @ingroup textlogs
     */
    void setLogger(Logger* logwriter) {pMessenger->setLogger(logwriter);}
   
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
    void beginLogGroup(std::string title, int loglevel) { pMessenger->beginLogGroup(title,loglevel);} 

    //! Add an entry to an HTML log file.
    /*!
     *  Entries appear in the form "tag:value".
     *
     * @param tag      tag
     * @param value    double value
     *
     * @ingroup HTML_logs
     */
    void addLogEntry(std::string tag, std::string value) { pMessenger->addLogEntry(tag, value);}

    //! Add an entry to an HTML log file.
    /*!
     *  Entries appear in the form "tag:value".
     *
     * @param tag      tag
     * @param value    double value
     *
     * @ingroup HTML_logs
     */
    void addLogEntry(std::string tag, doublereal value) { pMessenger->addLogEntry(tag, value);}

    //! Add an entry to an HTML log file.
    /*!
     *  Entries appear in the form "tag:value".
     *
     * @param tag      tag
     * @param value    double value
     *
     * @ingroup HTML_logs
     */
    void addLogEntry(std::string tag, int value) { pMessenger->addLogEntry(tag, value);}

    //! Add an entry to an HTML log file.
    /*!
     *  Entries appear in the form "msg".
     *
     * @param msg      Message to be added to file
     *
     * @ingroup HTML_logs
     */
    void addLogEntry(std::string msg) { pMessenger->addLogEntry(msg); }

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
    void endLogGroup(std::string title) { pMessenger->endLogGroup(title) ;}

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
    void write_logfile(std::string file) { pMessenger->write_logfile(file) ; }
#endif

protected:
    void setDefaultDirectories();

    //! Current vector of input directories to search for input files
    vector<string> inputDirs;
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
    map<string, string>     options;
    //! Current value of tmp_dir
    string tmp_dir;
    //! Current vector of xml file trees that have been previously parsed
    map<string, XML_Node*> xmlfiles;
    //! Current sleep command.
    string m_sleep;
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
      auto_ptr< Messages >   pMessenger ;
#endif
    
  private:
     static Application* s_app ;
  };
            
  //! Return a pointer to the application object
  inline Application* app() {
    return Application::Instance() ;
  }


  /// Pointer to the single Application instance
  Application* Application::s_app = 0;

    /**
     * Definition of the static member of the Unit class.
     */
    Unit* Unit::s_u = 0;
      #if defined(THREAD_SAFE_CANTERA)
        boost::mutex Unit::units_mutex ;
      #endif


    /*
     * Delete all global data.  It should be called at the end of the
     * application if leak checking is to be done.
     */
    void appdelete() {
        Application::ApplicationDestroy() ;
        FactoryBase::deleteFactories();
        Unit::deleteUnit();
    }

    /*
     * Delete Messenger object allocated per thread.
     */
    void thread_complete() {
       app()->thread_complete() ;
    }
  
    void Application::thread_complete() {
#if defined(THREAD_SAFE_CANTERA)
         pMessenger.removeThreadMessages() ;
#endif
    }

    XML_Node* get_XML_File(std::string file, int debug) {
       XML_Node* xtmp = app()->get_XML_File(file, debug) ;
       //writelog("get_XML_File: returned from app:get_XML_FILE " + int2str(int(xtmp)) + "\n");
       return xtmp;
    }

    XML_Node* Application::get_XML_File(std::string file, int debug) {
        XML_LOCK();
		std::string path = "";
        /*
        try {
            path = findInputFile(file);
        }
        catch (CanteraError) {
            string::size_type idot = file.rfind('.');
            string ext = "";
            if (idot != string::npos) {
                ext = file.substr(idot, file.size());
                string ctifile = file.substr(0,idot)+".cti";
                try {
                    path = findInputFile(ctifile);
                }
                catch (CanteraError) {
                    path = findInputFile(file);
                }
            }
            else
                path = findInputFile(file);
        }
        */
        // The code above will try to process a cti file if the xml
        // file is not found. But I (dgg) don't think it makes much sense,
        // so it is replaced by:
        path = findInputFile(file);
        // 
#ifdef WIN32
        // RFB: For Windows make the path POSIX compliant so code looking for directory 
        // separators is simpler.  Just look for '/' not both '/' and '\\'
        replace_if( path.begin(), path.end(), bind2nd( equal_to<char>(), '\\'), '/' ) ; 
#endif

        string ff = path;
        if (xmlfiles.find(path) 
            == xmlfiles.end()) {
            /*
             * Check whether or not the file is XML. If not, it will
             * be first processed with the preprocessor. We determine
             * whether it is an XML file by looking at the file extension.
             */
            string::size_type idot = path.rfind('.');
            string ext;
            if (idot != string::npos) {
                ext = path.substr(idot, path.size());
            } else {
                ext = "";
                idot = path.size();
            }
            if (ext != ".xml" && ext != ".ctml") {
                /*
                 * We will assume that we are trying to open a cti file.
                 * First, determine the name of the xml file, ff, derived from
                 * the cti file.
                 * In all cases, we will write the xml file to the current
                 * directory.
                 */
                string::size_type islash = path.rfind('/');
                if (islash != string::npos) 
                    ff = string("./")+path.substr(islash+1,idot-islash - 1) + ".xml";
                else {
                    ff = string("./")+path.substr(0,idot) + ".xml";
                }
#ifdef DEBUG_PATHS
                cout << "get_XML_File(): Expected location of xml file = "
                     << ff << endl;
#endif
                /*
                 * Do a search of the existing XML trees to determine if we have
                 * already processed this file. If we have, return a pointer to
                 * the processed xml tree.
                 */
                if (xmlfiles.find(ff) != xmlfiles.end()) {
#ifdef DEBUG_PATHS
                    cout << "get_XML_File(): File, " << ff << ", was previously read."
                         << " Retrieving the storred xml tree." << endl;
#endif
                    return xmlfiles[ff];      
                }
                /*
                 * Ok, we didn't find the processed XML tree. Do the conversion
                 * to xml, possibly overwriting the file, ff, in the process.
                 */
                ctml::ct2ctml(path.c_str(),debug);
            }
            else {
                ff = path;
            }
            /*
             * Take the XML file ff, open it, and process it, creating an
             * XML tree, and then adding an entry in the map. We will store
             * the absolute pathname as the key for this map.
             */
            ifstream s(ff.c_str());

            XML_Node* x = new XML_Node("doc");
            if (s) {
                x->build(s);
                x->lock();
                xmlfiles[ff] = x;
            }
            else {
                string estring = "cannot open "+ff+" for reading.";
                estring += "Note, this error indicates a possible configuration problem."; 
                throw CanteraError("get_XML_File", estring);
            }
        }

        /*
         * Return the XML node pointer. At this point, we are sure that the
         * lookup operation in the return statement will return a valid
         * pointer. 
         */
        return xmlfiles[ff];
    }

  // Close an XML File 
  /*
   * Close a file that is openned by the application object
   *
   * @param file string name of the file
   */
  void close_XML_File(std::string file) {
    app()->close_XML_File(file) ;
  }

  void Application::close_XML_File(std::string file) {
    XML_LOCK();
    if (file == "all") {
      map<string, XML_Node*>::iterator 
	b = xmlfiles.begin(), 
	e = xmlfiles.end();
      for(; b != e; ++b) {
	b->second->unlock();
	delete b->second;
	xmlfiles.erase(b->first);
      }
    }
    else if (xmlfiles.find(file) != xmlfiles.end()) {
      xmlfiles[file]->unlock();
      delete xmlfiles[file];
      xmlfiles.erase(file);
    }
  }

    void setTmpDir(std::string tmp) { app()->setTmpDir(tmp); }
    void Application::setTmpDir(std::string tmp) {
       APP_LOCK();
       tmp_dir = tmp ;
    }
    
    string tmpDir() { return app()->getTmpDir(); }
    std::string Application::getTmpDir() {
       APP_LOCK();
       return tmp_dir ;
    }
    
    string sleep() { return app()->sleep(); }
    std::string Application::sleep() {
       APP_LOCK();
       return m_sleep ;
    }


  /*
   * Return the number of errors that have been encountered so far.
   * \ingroup errorhandling
   */
  int nErrors() {
    return app()->getErrorCount();
  }
  
  // Return the number of errors encountered so far
  int Application::Messages::getErrorCount() {
    return static_cast<int>(errorMessage.size()) ;
  }

  /*
   * popError eliminates the last error message that Cantera
   * has saved. Cantera saves a stack of exceptions that it
   * has caught in the Application class. This routine eliminates
   * the last exception to be added to that stack.
   * \ingroup errorhandling
   */
  void popError() {
    app()->popError() ;
  }

  // Eliminate the last error message
  void Application::Messages::popError() {
    if ( static_cast<int>(errorMessage.size()) > 0) {
      errorRoutine.pop_back() ;
      errorMessage.pop_back() ;
    }
  }

  /*
   * Retrieve the last error message.
   * This routine will retrieve the last error message and return
   * it in the return string.
   * \ingroup errorhandling
   */
  string lastErrorMessage() {
    return app()->lastErrorMessage() ;
  }

  // Retrieve the last error message
  std::string Application::Messages::lastErrorMessage() {
    if ( static_cast<int>(errorMessage.size()) > 0) {
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

    /*
     * Prints all of the error messages to stream f.
     * Write out to ostream, f, all of the saved error messages.
     * Cantera saves a stack of exceptions that it
     * has caught in the Application class. This routine writes
     * out all of the error messages to ostream f, and then
     * clears them from internal storage.
     * \ingroup errorhandling
     */
    void showErrors(std::ostream& f) {
        app()->getErrors(f) ;
    }

  // Prints all error messages to stream f
  void Application::Messages::getErrors( std::ostream& f) {
    int i = static_cast<int>(errorMessage.size());
    if (i == 0) return;
    f << endl << endl;
    f << "************************************************" << endl;
    f << "                   Cantera Error!                  " << endl;
    f << "************************************************" << endl
      << endl;
    int j;
    for (j = 0; j < i; j++) {
      f << endl;
      f << "Procedure: " << errorRoutine[j] << endl;
      f << "Error:     " << errorMessage[j] << endl;
    } 
    f << endl << endl;
    errorMessage.clear();
    errorRoutine.clear();
  }
  
    /*
     * Print all of the error messages using function writelog.
     * Write out all of the saved error messages to the log device.
     * Cantera saves a stack of exceptions that it
     * has caught in the Application class. This routine writes
     * out all of the error messages to the log, usually stdout,
     *  and then clears them from internal storage.
     * \ingroup errorhandling
     */
    void showErrors() {
        app()->logErrors() ;
    }

  // Print out all of the error messages
  void Application::Messages::logErrors() {
    int i = static_cast<int>(errorMessage.size());
    if (i == 0) return;
    writelog("\n\n");
    writelog("************************************************\n");
    writelog("                   Cantera Error!                  \n");
    writelog("************************************************\n\n");
    int j;
    for (j = 0; j < i; j++) {
      writelog("\n");
            writelog(string("Procedure: ")+ errorRoutine[j]+" \n");
            writelog(string("Error:     ")+ errorMessage[j]+" \n");
    } 
    writelog("\n\n");
    errorMessage.clear();
    errorRoutine.clear();
  }

  /*
   * Set an error condition in the application class without 
   * throwing an exception.
   * This routine adds an error message to the end of the stack
   * of errors that Cantera accumulates in the Application
   * class.
   * \ingroup errorhandling
   */
  void setError(std::string r, std::string msg) {
    app()->addError(r, msg) ;
  }

  // Set an error condition in the application class without throwing an exception
  void Application::Messages::addError(std::string r, std::string msg) {
    errorMessage.push_back(msg);
    errorRoutine.push_back(r);
  }
  
   
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
    void Application::setDefaultDirectories() {
        vector<string>& dirs = inputDirs;

        // always look in the local directory first
        dirs.push_back(".");


#ifdef WIN32
        //
        // Under Windows, the Cantera setup utility puts data files in
        // a directory 'Cantera\data' below the one the environment
        // variable COMMONPROGRAMFILES points to. (This is usually
        // C:\Program Files\Common Files.) If this environment
        // variable is defined, then this directory is assumed to
        // exist and is added to the search path.
        //
        const char* comfiles = getenv("COMMONPROGRAMFILES");
        if (comfiles != 0) {
            string cfiles = string(comfiles);

            // remove quotes if necessary
            if (cfiles[0] == '\'') 
                cfiles = cfiles.substr(1,1000);
            if (cfiles[cfiles.size()-1] == '\'') cfiles[cfiles.size()-1] = '\n';

            string datadir = string(comfiles) + "/Cantera/data";
            string tmpldir = string(comfiles) + "/Cantera/templates";
            dirs.push_back(datadir);
            dirs.push_back(tmpldir);
        }
#endif

#ifdef DARWIN
        //
        // add a default data location for Mac OS X
        //
        dirs.push_back("/Applications/Cantera/data");
#endif

        //
        // if environment variable CANTERA_DATA is defined, then add
        // it to the search path
        //
        if (getenv("CANTERA_DATA") != 0) {
            string datadir = string(getenv("CANTERA_DATA"));
            dirs.push_back(datadir);
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


  // Add a directory to the input file search path.
  // @ingroup inputfiles
  void addDirectory(std::string dir) {
    app()->addDataDirectory( dir ) ;
  }
  void Application::addDataDirectory( std::string dir ) {
    DIR_LOCK() ;
    if (inputDirs.size() == 0) setDefaultDirectories();
    string d = stripnonprint(dir);
    size_t m, n = inputDirs.size();
    
    // don't add if already present
    for (m = 0; m < n; m++) {
      if (d == inputDirs[m]) return;
    }
    
    inputDirs.push_back(d);
  }
  
  /*    
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
   *    \return
   *    
   *      The absolute path name of the first matching
   *      file is returned. If a relative path name
   *      is indicated, the relative path name is returned.
   *  
   *      If the file is not found, a message is written to 
   *      stdout and  a CanteraError exception is thrown.
   */
  std::string findInputFile(std::string name) {
    return app()->findInputFile(name) ;
  }
  
  std::string Application::findInputFile(std::string name) {
    DIR_LOCK() ;
    string::size_type islash = name.find('/');
    string::size_type ibslash = name.find('\\');
    string inname;
    vector<string>& dirs = inputDirs;
 
    int nd;
    if (islash == string::npos && ibslash == string::npos) {
      nd = static_cast<int>(dirs.size());
      int i;
      inname = "";
      for (i = 0; i < nd; i++) {
	inname = dirs[i] + "/" + name;
	ifstream fin(inname.c_str());
	if (fin) {
	  fin.close();
	  return inname;
	}
      }
      string msg;
      msg = "\nInput file " + name 
	+ " not found in director";
      msg += (nd == 1 ? "y " : "ies ");
      for (i = 0; i < nd; i++) {
	msg += "\n'" + dirs[i] + "'";
	if (i < nd-1) msg += ", ";
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

    doublereal toSI(std::string unit) {
        doublereal f = Unit::units()->toSI(unit);
        if (f) return f;
        else throw CanteraError("toSI","unknown unit string: "+unit);
        //return 1.0;
    }

    doublereal actEnergyToSI(std::string unit) {
        doublereal f = Unit::units()->actEnergyToSI(unit);
        if (f) return f;
        else return 1.0;
    }

    string canteraRoot() {
        char* ctroot = 0;
        ctroot = getenv("CANTERA_ROOT");
        if (ctroot != 0) { return string(ctroot); }
        else {
#ifdef CANTERA_ROOT
            return string(CANTERA_ROOT);
#else
            return "";
#endif
        }
    }

    // exceptions

    CanteraError::CanteraError(std::string proc, std::string msg) {
        app()->addError(proc, msg);
    }
    
    ArraySizeError::ArraySizeError(std::string proc, int sz, int reqd) :
        CanteraError(proc, "Array size ("+int2str(sz)+
            ") too small. Must be at least "+int2str(reqd)) {}

    ElementRangeError::ElementRangeError(std::string func, int m, int mmax) :
        CanteraError(func, "Element index " + int2str(m) + 
            " outside valid range of 0 to " + int2str(mmax-1)) {}



    ///////////////////////////////////////////////////////////
    //
    //  Warnings 
    //
    //////////////////////////////////////////////////////////

    // Print a warning when a deprecated method is called.
    // @param classnm Class the method belongs to
    // @param oldnm Name of the deprecated method
    // @param newnm Name of the method users should use instead
    void deprecatedMethod(std::string classnm, std::string oldnm, std::string newnm) {
        writelog(">>>> WARNING: method "+oldnm+" of class "+classnm
            +" is deprecated.\n");
        writelog("         Use method "+newnm+" instead.\n");
        writelog("         (If you want to rescue this method from deprecated\n");
        writelog("         status, see http://www.cantera.org/deprecated.html)");
    }


    void removeAtVersion(std::string func, std::string version) {
      //if (version >= "CANTERA_VERSION") {
      writelog("Removed procedure: "+func+"\n");
      writelog("Removed in version: "+version+"\n");
      throw CanteraError("removeAtVersion: "+ func,"procedure has been removed.");
      // }
    }


  // defgroup logs Diagnostic Output (global.h)
  //
  // Writing diagnostic information to the screen or to a file.
  // It is often useful to be able to write diagnostic messages to
  // the screen or to a file. Cantera provides two sets of
  // procedures for this purpose. The first set is designed to
  // write text messages to the screen to document the progress of
  // a complex calculation, such as a flame simulation.The second
  // set writes nested lists in HTML format. This is useful to
  // print debugging output for a complex calculation that calls
  // many different procedures.
  
  
  // defgroup textlogs Writing messages to the screen (global.h)
  // ingroup logs
  //
  // Write a message to the screen. The string may be of any
  // length, and may contain end-of-line characters. This method is
  // used throughout Cantera to write log messages. It can also be
  // called by user programs.  The advantage of using writelog over
  // writing directly to the standard output is that messages
  // written with writelog will display correctly even when Cantera
  // is used from MATLAB or other application that do not have a
  // standard output stream. @ingroup textlogs

  // Write a message to the screen
   void writelog(const std::string& msg) {
     app()->writelog(msg);
   }

  void writelogAM(const std::string& msg) {
     app()->writelog(msg);
   }

  // Write a message to the screen
  void Application::Messages::writelog(const std::string& msg) {
    logwriter->write(msg);
  }

  // Write a message to the screen.
  void writelog(const char* msg) {
    app()->writelog(msg);
  }
  void writelogAM(const char* msg) {
    app()->writelog(msg);
  }

  // Write a message to the screen
  void Application::Messages::writelog(const char* pszmsg) {
    logwriter->write( pszmsg ) ;
  }

  // Write an endl to the screen and flush output
  void Application::Messages::writelogendl() {
    logwriter->writeendl();
  }
 
  // Write a message to the screen using printf format
  void writelogf(const char* fmt,...) {
    enum { BUFSIZE = 2048 } ;
    char sbuf[BUFSIZE] ;
       
    va_list args ;
       
    va_start( args, fmt ) ;
       
#if defined(WIN32) && defined(MSC_VER)
    _vsnprintf( sbuf, BUFSIZE, fmt, args ) ; 
#else
    vsprintf( sbuf, fmt, args ) ;
#endif
       
    writelog( sbuf ) ;
       
    va_end(args) ;
  }

  void writelogendl() {
    app()->writelogendl();
  }

  // Write an error message and terminate execution. test.
  // @ingroup textlogs
  void error(const std::string& msg) {
    app()->logerror(msg);
  }

  // Write an error message and terminate execution
  void Application::Messages::logerror(const std::string& msg) {
    logwriter->error(msg) ;
  }

    // @ingroup textlogs
   int userInterface() {
       return app()->getUserEnv();
   }

  // Return an integer specifying the application environment.
  int Application::Messages::getUserEnv() {
    return logwriter->env() ;
  }

  // Install a logger. Called by the language interfaces to install an
  // appropriate logger. 
  // @see Logger.
  // @ingroup textlogs
  void setLogger(Logger* logwriter) {
		try {
          app()->setLogger(logwriter) ;
		}
		catch (std::bad_alloc) {
          logwriter->error("bad alloc thrown by app()");
		}
  }

    void Application::Messages::setLogger(Logger* _logwriter) {
       if ( logwriter == _logwriter ) return ;
       if (logwriter != 0 ) {
         delete logwriter;
         logwriter = 0 ;
       }
       logwriter = _logwriter;
    }


#ifdef WITH_HTML_LOGS

  ////////////////////////////////////////////////////////////////
  // 
  // defgroup HTML_logs Writing HTML Logfiles
  // ingroup logs
  // 
  //  These functions are designed to allow writing HTML diagnostic
  //  messages in a manner that allows users to control how much
  //  diagnostic output to print. It works like this: Suppose you
  //  have function A that invokes function B that invokes function
  //  C. You want to be able to print diagnostic messages just from
  //  function A, or from A and B, or from A, B, and C, or to turn
  //  off printing diagnostic messages altogether. All you need to
  //  do is call 'beginLogGroup' within function A, and specify a
  //  loglevel value. Then in B, call beginLogGroup again, but
  //  without an explicit value for loglevel. By default, the
  //  current level is decremented by one in beginLogGroup. If it
  //  is <= 0, no log messages are written. Thus, if each function
  //  begins with beginLogGroup and calls endLogGroup before
  //  returning, then setting loglevel = 3 will cause messages from
  //  A, B, and C to be written (in nested HTML lists), loglevel =
  //  2 results in messages only being written from A and B, etc.
  //
  /////////////////////////////////////////////////////////////////


  // Create a new group for log messages.  Usually this is called
  // upon entering the function, with the title parameter equal to
  // the name of the function or method. Subsequent messages
  // written with addLogEntry will appear grouped under this
  // heading, until endLogGroup() is called.
  // @ingroup HTML_logs
  void beginLogGroup(std::string title, int loglevel) {
    app()->beginLogGroup(title, loglevel) ;
  }

  void Application::Messages::beginLogGroup(std::string title, int _loglevel /*=-99*/) {
    // loglevel is a member of the Messages class.
    if (_loglevel != -99) {
      loglevel = _loglevel;
    } else {
      loglevel--;
    }
    if (loglevel <= 0) return;
    // Add the current loglevel to the vector of loglevels
    loglevels.push_back(loglevel);
    // Add the title of the current logLevel to the vector of titles
    loggroups.push_back(title);
    // If we haven't started an XML tree for the log file, do so here
    if (xmllog == 0) {
      // The top of this tree will have a zero pointer.
      xmllog = new XML_Node("html");
      current = &xmllog->addChild("ul");
    }
    // Add two children to the XML tree.
    current = &current->addChild("li","<b>"+title+"</b>");
    current = &current->addChild("ul");
  }

  // Add an entry to the log file. Entries appear in the form "tag:
  // value".
  // @ingroup HTML_logs
  void addLogEntry(std::string tag, std::string value) {
    app()->addLogEntry(tag, value) ;
  }

  void Application::Messages::addLogEntry(std::string tag, std::string value) {
    if (loglevel > 0 && current) 
      current->addChild("li",tag+": "+value);
  }

  // Add an entry to the log file. Entries appear in the form "tag:
  // value".
  // @ingroup HTML_logs
  void addLogEntry(std::string tag, doublereal value) {
    app()->addLogEntry( tag, value ) ;
  }

  void Application::Messages::addLogEntry(std::string tag, doublereal value) {
    if (loglevel > 0 && current) 
      current->addChild("li",tag+": "+fp2str(value));
  }

  // Add an entry to the log file. Entries appear in the form "tag:
  // value".
  // @ingroup HTML_logs
  void addLogEntry(std::string tag, int value) {
    app()->addLogEntry( tag, value ) ;
  }

  void Application::Messages::addLogEntry(std::string tag, int value) {
    if (loglevel > 0 && current) 
      current->addChild("li",tag+": "+int2str(value));
  }

  // Add an entry to the log file.
  // @ingroup HTML_logs
  void addLogEntry(std::string msg) {
    app()->addLogEntry(msg) ;
  }

  void Application::Messages::addLogEntry(std::string msg) {
    if (loglevel > 0 && current)
      current->addChild("li",msg);
  }

  // Close the current group of log messages. This is typically
  // called just before leaving a function or method, to close the
  // group of messages that were output from this
  // function. Subsequent messages written with addLogEntry will
  // appear at the next-higher level in the outline, unless
  // beginLogGroup is called first to create a new group.  
  // @ingroup HTML_logs
  void endLogGroup(std::string title) {
    app()->endLogGroup(title) ;
  }

  void Application::Messages::endLogGroup(std::string title) {
    if (loglevel <= 0) return;
    AssertThrowMsg(current, "Application::Messages::endLogGroup",
		   "Error while ending a LogGroup. This is probably due to an unmatched"
		   " beginnning and ending group");
    current = current->parent();
    AssertThrowMsg(current, "Application::Messages::endLogGroup",
		   "Error while ending a LogGroup. This is probably due to an unmatched"
		   " beginnning and ending group");
    current = current->parent();
    // Get the loglevel of the previous level and get rid of 
    // vector entry in loglevels.
    loglevel = loglevels.back();
    loglevels.pop_back();
    if (title != "" && title != loggroups.back()) {
      writelog("Logfile error."
	       "\n   beginLogGroup: "+ loggroups.back()+
	       "\n   endLogGroup:   "+title+"\n");
      write_logfile("logerror"); 
    } else if (loggroups.size() == 1) {
      write_logfile(loggroups.back()+"_log"); 
      loggroups.clear();
      loglevels.clear();
    } else {
      loggroups.pop_back();
    }
  }


  // Write the HTML log file.
  /*
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
   *  @ingroup HTML_logs
   */
  void write_logfile(std::string file) {
    app()->write_logfile(file) ;
  }

  // Write the HTML log file.
  void Application::Messages::write_logfile(std::string file) {
    if (!xmllog) {
      return;
    }
    string::size_type idot = file.rfind('.');
    string ext = "";
    string nm = file;
    if (idot != string::npos) {
      ext = file.substr(idot, file.size());
      nm = file.substr(0,idot);
    }
    else {
      ext = ".html";
      nm = file;
    }

    // see if file exists. If it does, find an integer that
    // can be appended to the name to create the name of a file
    // that does not exist.
    string fname = nm + ext;
    ifstream f(fname.c_str());
    if (f) {
      int n = 0;
      while (1 > 0) {
	n++;
	fname = nm + int2str(n) + ext;
	ifstream f(fname.c_str());
	if (!f) break;
      }
    }

    // Now we have a file name that does not correspond to any 
    // existing file. Open it as an output stream, and dump the 
    // XML (HTML) tree to it.

    if (xmllog) {
      ofstream f(fname.c_str());
      // go to the top of the tree, and write it all.
      xmllog->root().write(f);
      f.close();
      writelog("Log file " + fname + " written.\n");
      delete xmllog;
      xmllog = 0;
      current = 0;
    }
  }

#endif // WITH_HTML_LOGS


  /// split a string at a '#' sign. Used to separate a file name
  /// from an id string.
  static void split_at_pound(const std::string& src, std::string& file, std::string& id) { 
    string::size_type ipound = src.find('#');
    if (ipound != string::npos) {
      id = src.substr(ipound+1,src.size());
      file = src.substr(0,ipound);
    }
    else {
      id = "";
      file = src;
    }
  }
  /*
   * This routine will locate an XML node in either the input
   * XML tree or in another input file specified by the file
   * part of the file_ID string. Searches are based on the
   * ID attribute of the XML element only.
   *
   * param file_ID This is a concatenation of two strings seperated
   *                by the "#" character. The string before the
   *                pound character is the file name of an xml
   *                file to carry out the search. The string after
   *                the # character is the ID attribute 
   *                of the xml element to search for. 
   *                The string is interpreted as a file string if
   *                no # character is in the string.
   *
   * param root    If the file string is empty, searches for the
   *                xml element with matching ID attribute are
   *                carried out from this XML node.
   */
  XML_Node* get_XML_Node(const std::string& file_ID, XML_Node* root) {
    std::string fname, idstr;
    XML_Node *db, *doc;
    split_at_pound(file_ID, fname, idstr);
    if (fname == "") {
      if (!root) throw CanteraError("get_XML_Node",
				    "no file name given. file_ID = "+file_ID);
      db = root->findID(idstr, 3);
    } 
    else {
      doc = get_XML_File(fname);
      if (!doc) throw CanteraError("get_XML_Node", 
				   "get_XML_File failed trying to open "+fname);
      db = doc->findID(idstr, 3);
    }
    if (!db) {
      throw CanteraError("get_XML_Node", 
			 "id tag '"+idstr+"' not found.");
    }
    return db;
  }


  /*
   * This routine will locate an XML node in either the input
   * XML tree or in another input file specified by the file
   * part of the file_ID string. Searches are based on the
   * XML element name and the ID attribute of the XML element.
   * An exact match of both is usually required. However, the
   * ID attribute may be set to "", in which case the first
   * xml element with the correct element name will be returned.
   *
   * @param nameTarget This is the XML element name to look for.
   *                   
   * @param file_ID This is a concatenation of two strings seperated
   *                by the "#" character. The string before the
   *                pound character is the file name of an xml
   *                file to carry out the search. The string after
   *                the # character is the ID attribute 
   *                of the xml element to search for. 
   *                The string is interpreted as a file string if
   *                no # character is in the string.
   *
   * @param root    If the file string is empty, searches for the
   *                xml element with matching ID attribute are
   *                carried out from this XML node.
   */
  XML_Node* get_XML_NameID(const std::string& nameTarget,
			   const std::string& file_ID, 
			   XML_Node* root) {
    string fname, idTarget;
    XML_Node *db, *doc;
    split_at_pound(file_ID, fname, idTarget);
    if (fname == "") {
      if (!root) return 0;
      db = root->findNameID(nameTarget, idTarget);
    } else {
      doc = get_XML_File(fname);
      if (!doc) return 0;
      db = doc->findNameID(nameTarget, idTarget);
    }
    return db;
  }

  //! Assigned storage for the static member of the FactoryBase class
  std::vector<FactoryBase*> FactoryBase::s_vFactoryRegistry;
}


