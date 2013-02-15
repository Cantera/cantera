/**
 * @file global.h
 * This file contains definitions for utility functions and text for modules,
 * inputfiles, logs, textlogs, HTML_logs (see \ref inputfiles, \ref logs,
 * \ref textlogs and \ref HTML_logs).
 *
 * @ingroup utils
 *
 * These functions store
 * some parameters in  global storage that are accessible at all times
 * from the calling application.
 * Contains module definitions for
 *     -  inputfiles  (see \ref inputfiles)
 *     -  logs        (see \ref logs)
 *     -  textlogs    (see \ref textlogs)
 *     -  HTML_logs   (see \ref HTML_logs)
 */
// Copyright 2001  California Institute of Technology
#ifndef CT_GLOBAL_H
#define CT_GLOBAL_H

#include "ct_defs.h"
#include "IndependentVars.h"

#include <vector>

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
 * @ingroup errorhandling
 */
void setError(const std::string& r, const std::string& msg);

//!  Prints all of the error messages to an ostream
/*!
 * Write out all of the saved error messages to the ostream f
 * using the member function writelog of class logger.
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

//! Print all of the error messages using function writelog of class logger.
/*!
 * Print all of the error messages
 * using the member function writelog of class logger.
 * Write out all of the saved error messages to the log device.
 * Cantera saves a stack of exceptions that it
 * has caught in the Application class. This routine writes
 * out all of the error messages to the log, usually stdout,
 * and then clears them from internal storage.
 *
 * \ingroup errorhandling
 */
void showErrors();

//! Discard the last error message
/*!
 * %Cantera saves a stack of exceptions that it
 * has caught in the Application class. This routine eliminates
 * the last exception to be added to that stack.
 *
 * \ingroup errorhandling
 */
void popError();

/*!
 * @defgroup inputfiles Input File Handling
 *
 * The properties of phases and interfaces are specified in
 * text files. These procedures handle various aspects of reading
 * these files.
 *
 * For input files not specified by an absolute pathname,
 * %Cantera searches
 * for input files along a path that includes platform-specific
 * default locations, and possibly user-specified locations.
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
 *
 * There are two different types of input files within %Cantera:
 *    ctml: This is an xml file laid out in such a way that %Cantera can
 *          interpret the contents.
 *    cti:  A human-readable ascii format for information that %Cantera
 *          will read.
 *
 *   %Cantera can take its input from both types of files. However, given
 *   a file in cti format, the initial operation that %Cantera will perform
 *   is to translate the cti file into a ctml file.
 *   The translation is carried out via a system call to a python interpreter
 *   program that actually carries out the translation. In general, a new
 *   ctml file is created by the translation that is written to the current
 *   local directory.
 *   The ctml file is then read back into %Cantera as the input.
 *
 * Other input routines in other modules:
 *   @see importKinetics()
 *
 * @{
 */

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
std::string findInputFile(const std::string& name);

//!  Add a directory to the input file search path.
/*!
 * @ingroup inputfiles
 *
 * @param dir  String name for the directory to be added to the search path
 */
void addDirectory(const std::string& dir);

//@}

//! Delete and free all memory associated with the application
/*!
 * Delete all global data.  It should be called at the end of the
 * application if leak checking is to be done.
 */
void appdelete();

//! Delete and free memory allocated per thread in multithreaded applications
/*!
 * Delete the memory allocated per thread by Cantera.  It should be called from
 * within the thread just before the thread terminates.  If your version of Cantera has not
 * been specifically compiled for thread safety this function does nothing.
 */
void thread_complete();

//! Returns root directory where %Cantera where installed
/*!
 * @return
 *  Returns a string containing the name of the base directory where %Cantera is installed.
 *  If the environmental variable CANTERA_ROOT is defined, this function will
 *  return its value, preferentially.
 *
 * @ingroup inputfiles
 */
std::string canteraRoot();

/*!
 * @defgroup logs Diagnostic Output
 *
 * Writing diagnostic information to the screen or to a file.
 * It is often useful to be able to write diagnostic messages to
 * the screen or to a file. Cantera provides two sets of
 * procedures for this purpose. The first set is designed to
 * write text messages to the screen to document the progress of
 * a complex calculation, such as a flame simulation.The second
 * set writes nested lists in HTML format. This is useful to
 * print debugging output for a complex calculation that calls
 * many different procedures.
 */

/*!
 * @defgroup textlogs Writing messages to the screen
 * @ingroup logs
 */

//!  Write a message to the screen.
/*!
 * The string may be of any
 * length, and may contain end-of-line characters. This method is
 * used throughout Cantera to write log messages. It can also be
 * called by user programs.  The advantage of using writelog over
 * writing directly to the standard output is that messages
 * written with writelog will display correctly even when Cantera
 * is used from MATLAB or other application that do not have a
 * standard output stream.
 *
 *  This routine is part of the interface suite whose behavior changes
 *  with the interface. The interface suite has been moved to the
 *  class logger and inherited classes of logger.
 *
 * @param msg  String message to be written to the screen
 * @ingroup textlogs
 */
void writelog(const std::string& msg);

inline void writelog(const std::string& msg, int loglevel)
{
    if (loglevel > 0) {
        writelog(msg);
    }
}

inline void writelog1(const std::string& msg, int loglevel)
{
    if (loglevel > 0) {
        writelog(msg);
    }
}

//!  Write a message to the screen.
/*!
 * The string may be of any
 * length, and may contain end-of-line characters. This method is
 * used throughout %Cantera to write log messages.
 *
 *  This routine is part of the interface suite whose behavior changes
 *  with the interface. The interface suite has been moved to the
 *  class logger and inherited classes of logger.
 *
 * @param msg  c character string to be written to the screen
 * @ingroup textlogs
 */
void writelog(const char* msg);

//! Write a formatted message to the screen
/*!
 * Using the printf formatting of C, write a message to the screen
 * with variable values.
 *
 * Here, we format an internal string with the correct values
 * and then feed it into writelog().
 *
 * @param fmt  c format string for the following arguments
 * @ingroup textlogs
 */
void writelogf(const char* fmt, ...);

//! Write an end of line character to the screen and flush output
/*!
 * Some implementations differentiate between \n and endl in
 * terms of when the output is flushed.
 */
void writelogendl();

//! Write an error message and terminate execution.
/*!
 *  This routine is part of the interface suite whose behavior changes
 *  with the interface. The interface suite has been moved to the
 *  class logger and inherited classes of logger.
 *
 *  @param msg Error message to be written to the screen.
 *  @ingroup textlogs
 */
void error(const std::string& msg);

//! Install a logger.
/*!
 *  Called by the language interfaces to install an appropriate logger.
 *  The logger is used for the writelog() function
 *
 * @param logwriter Pointer to a logger object
 * @see Logger.
 * @ingroup textlogs
 */
void setLogger(Logger* logwriter);

//! Set the default form of the independent variables
/*!
 *   Set the default form of the independent variables for the application. This is used to
 *   construct a default form of the independent variables for ThermoPhase, Kinetics, and transport
 *   objects.
 *
 *   The jacobian wrt the Independent variables is calculated for calls to property evaluators.
 *
 *   @param indvars                INDVAR_FORM enum.   Set the default form of the independent variables to
 *                                 an enum given by this value. The enum is defined in the file
 *                                 IndependentVars.h
 *   @param hasVoltage             Boolean indicating we should add a voltage variable.
 *   @param hasSurfaceTension      Boolean indicating we should add a surface tension variable.
 */
void setDefaultIndependentVars(INDVAR_FORM indvars, bool hasVoltage = false, bool hasSurfaceTension = false);

//! Return the conversion factor to convert unit std::string 'unit'
//! to SI units.
/*!
 * @param unit String containing the units
 */
doublereal toSI(const std::string& unit);

/// Return the conversion factor to convert activation energy unit
/// std::string 'unit' to Kelvin.
/*!
 * @param unit  String containing the activation energy units
 */
doublereal actEnergyToSI(const std::string& unit);

/// Return a pointer to the XML tree for a Cantera input file.
/*!
 *  This routine will find the file and read the XML file into an
 *  XML tree structure. Then, a pointer will be returned. If the
 *  file has already been processed, then just the pointer will
 *  be returned.
 *
 * @param file String containing the relative or absolute file name
 * @param debug Debug flag
 */
XML_Node* get_XML_File(const std::string& file, int debug = 0);

/// Close a Cantera input file.
/*!
 * @param file String containing the relative or absolute file name
 */
void close_XML_File(const std::string& file);

#ifdef WITH_HTML_LOGS

/*!
 * @defgroup HTML_logs Writing HTML Logfiles
 * @ingroup logs
 *
 *  These functions are designed to allow writing HTML diagnostic
 *  messages in a manner that allows users to control how much
 *  diagnostic output to print. It works like this: Suppose you
 *  have function A that invokes function B that invokes function
 *  C. You want to be able to print diagnostic messages just from
 *  function A, or from A and B, or from A, B, and C, or to turn
 *  off printing diagnostic messages altogether. All you need to
 *  do is call 'beginLogGroup' within function A, and specify a
 *  loglevel value. Then in B, call beginLogGroup again, but
 *  without an explicit value for loglevel. By default, the
 *  current level is decremented by one in beginLogGroup. If it
 *  is <= 0, no log messages are written. Thus, if each function
 *  begins with beginLogGroup and calls endLogGroup before
 *  returning, then setting loglevel = 3 will cause messages from
 *  A, B, and C to be written (in nested HTML lists), loglevel =
 *  2 results in messages only being written from A and B, etc.
 */

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
void beginLogGroup(const std::string& title, int loglevel = -99);

//! Add an entry to an HTML log file.
/*!
 *  Entries appear in the form "tag:value".
 *
 * @param tag      tag
 * @param value    string value
 *
 * @ingroup HTML_logs
 */
void addLogEntry(const std::string& tag, const std::string& value);

//! Add an entry to an HTML log file.
/*!
 *  Entries appear in the form "tag:value".
 *
 * @param tag      tag
 * @param value    double value
 *
 * @ingroup HTML_logs
 */
void addLogEntry(const std::string& tag, doublereal value);

//! Add an entry to an HTML log file.
/*!
 *  Entries appear in the form "tag:value".
 *
 * @param tag      tag
 * @param value    int value
 *
 * @ingroup HTML_logs
 */
void addLogEntry(const std::string& tag, int value);

//! Add an entry msg string to an HTML log file.
/*!
 * Add a message string to the HTML log file
 *
 * @param msg      string mesg
 *
 * @ingroup HTML_logs
 */
void addLogEntry(const std::string& msg);

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
void endLogGroup(const std::string& title = "");

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
 *  @ingroup HTML_logs
 */
void write_logfile(const std::string& file = "log.html");

#else
inline void beginLogGroup(const std::string& title, int loglevel=-99) {}
inline void addLogEntry(const std::string& tag, const std::string& value) {}
inline void addLogEntry(const std::string& tag, doublereal value) {}
inline void addLogEntry(const std::string& tag, int value) {}
inline void addLogEntry(const std::string& msg) {}
inline void endLogGroup(const std::string& title="") {}
inline void write_logfile(const std::string& file = "log.html") {}
#endif

//! This routine will locate an XML node in either the input
//! XML tree or in another input file specified by the file
//! part of the file_ID string.
/*!
 *    Searches are based on the
 *    ID attribute of the XML element only.
 *
 * @param file_ID This is a concatenation of two strings separated
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
 *
 * @return
 *    Returns the XML_Node, if found. Returns null if not found.
 */
XML_Node* get_XML_Node(const std::string& file_ID, XML_Node* root);

//! This routine will locate an XML node in either the input
//! XML tree or in another input file specified by the file
//! part of the file_ID string.
/*!
 * Searches are based on the
 * XML element name and the ID attribute of the XML element.
 * An exact match of both is usually required. However, the
 * ID attribute may be set to "", in which case the first
 * xml element with the correct element name will be returned.
 *
 * @param nameTarget This is the XML element name to look for.
 *
 * @param file_ID This is a concatenation of two strings separated
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
 *
 * @return
 *    Returns the XML_Node, if found. Returns null if not found.
 */
XML_Node* get_XML_NameID(const std::string& nameTarget, const std::string& file_ID, XML_Node* root);

//! Clip *value* such that lower <= value <= upper
template<class T>
inline T clip(const T& value, const T& lower, const T& upper)
{
    return std::max(lower, std::min(upper, value));
}

void assignVectorFadToDouble(std::vector<doublereal>& left, const std::vector<Cantera::doubleFAD>& right);


template<typename T1, typename T2>
void assignVectorVectorClass(std::vector<std::vector<T1> >& left, const std::vector<std::vector<T2> >& right)
{
    size_t nn = (right).size();
    left.resize(nn);
    for (size_t k = 0; k < nn; k++) {
        const std::vector<T2> & pfrom = (right)[k];
        std::vector<T1> & pto = left[k];
        size_t mm = pfrom.size();
        pto.resize(mm);
        for (size_t l = 0; l < mm; l++) {
            pto[l] = pfrom[l];
        }
    }
};
}

#endif

