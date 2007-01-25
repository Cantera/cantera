/**
 * @file global.h
 *
 * These functions handle various utility functions, and store 
 * some parameters in 
 * global storage that are accessible at all times. 
 *  
 */

/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_GLOBAL_H
#define CT_GLOBAL_H

#include <string>
#include "ct_defs.h"

namespace Cantera {

    class XML_Node;
    class Logger;

  //! Return the number of errors that have been encountered so far
  /*!
   * @ingroup errorhandling
   */
  int nErrors();

  //! Returns the last error message
  /*!
   * @return String containing the description of the last error 
   *         message.
   *
   * @ingroup errorhandling
   */
  std::string lastErrorMessage();

  //! Set an error condition in the application class without throwing an exception.
  /*!
   * This routine adds an error message to the end of the stack
   * of errors that Cantera accumulates in the Application
   * class.
   * @param r Procedure name which is generating the error condition
   * @param msg Descriptive message of the error condition.
   *
   * \ingroup errorhandling
   */
  void setError(std::string r, std::string msg);

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
   * \ingroup errorhandling
   */
   void showErrors(std::ostream& f);

  //! Print all of the error messages using function writelog.
  /*!
   * Print all of the error messages using function writelog.
   * Write out all of the saved error messages to the log device.
   * Cantera saves a stack of exceptions that it
   * has caught in the Application class. This routine writes
   * out all of the error messages to the log, usually stdout,
   *  and then clears them from internal storage.
   * \ingroup errorhandling
   */
  void showErrors();

  //! Discard the last error message
  /*!
   * Cantera saves a stack of exceptions that it
   * has caught in the Application class. This routine eliminates
   * the last exception to be added to that stack.
   *
   * \ingroup errorhandling
   */
  void popError();

    /// Find an input file.
    std::string findInputFile(std::string name);

    void addDirectory(std::string dir);

    void appdelete();

    /// The root directory where Cantera is installed
    std::string canteraRoot();

    /// Set the temporary file directory. The default is to use the 
    /// directory specified by enviroment variable TMP or TEMP. If neither 
    /// of these are defined, then the current working directory will be 
    /// used for temporary files. Call this function to specify some other
    /// place to put temporary files.
    void setTmpDir(std::string tmp);

    /// The directory where temporary files may be created
    std::string tmpDir();

    /// Delay time in seconds.
    std::string sleep();

    void writelog(const std::string& msg);

    void writelog(const char* msg);

    void error(const std::string& msg);

    // returns 1 for MATLAB, 2 for Python, and 0 for C++ or Fortran.
    int userInterface();

    void setLogger(Logger* logwriter);

    /// Return the conversion factor to convert unit std::string 'unit' to
    /// SI units.
    doublereal toSI(std::string unit);

    /// Return the conversion factor to convert activation energy unit
    /// std::string 'unit' to Kelvin.
    doublereal actEnergyToSI(std::string unit);

    /// Return a pointer to the XML tree for a Cantera input file. 
    XML_Node* get_XML_File(std::string file);

    /// Close a Cantera input file.
    void close_XML_File(std::string file);

#ifdef WITH_HTML_LOGS
    void beginLogGroup(std::string title, int loglevel=-99);
    void addLogEntry(std::string tag, std::string value);
    void addLogEntry(std::string tag, doublereal value);
    void addLogEntry(std::string tag, int value);
    void addLogEntry(std::string msg);
    void endLogGroup(std::string title="");
    void write_logfile(std::string file = "log.html");
#else
    inline void beginLogGroup(std::string title, int loglevel=-99) {}
    inline void addLogEntry(std::string tag, std::string value) {}
    inline void addLogEntry(std::string tag, doublereal value) {}
    inline void addLogEntry(std::string tag, int value) {}
    inline void addLogEntry(std::string msg) {}
    inline void endLogGroup(std::string title="") {}
    inline void write_logfile(std::string file = "log.html") {}
#endif

}

#endif
