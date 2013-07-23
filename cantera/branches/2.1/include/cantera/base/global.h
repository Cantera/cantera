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

namespace Cantera
{

class XML_Node;
class Logger;

//! Declaration for whether the Debug mode is turned on within Cantera
/*!
 *  Turn on the mode by using the following compile time syntax
 *
 *       scons debug_verbose=True build
 */
extern const int g_DEBUG_MODE;

//! Return the number of errors that have been encountered so far
/*!
 * @ingroup errorhandling
 */
int nErrors();

//! @copydoc Application::Messages::lastErrorMessage
std::string lastErrorMessage();

//! @copydoc Application::Messages::addError
void setError(const std::string& r, const std::string& msg);

//! @copydoc Application::Messages::getErrors
void showErrors(std::ostream& f);

//! @copydoc Application::Messages::logErrors
void showErrors();

//! @copydoc Application::Messages::popError
void popError();

/*!
 * @defgroup inputfiles Input File Handling
 *
 * The properties of phases and interfaces are specified in text files. These
 * procedures handle various aspects of reading these files.
 *
 * For input files not specified by an absolute pathname, %Cantera searches
 * for input files along a path that includes platform-specific default
 * locations, and possibly user-specified locations.
 *
 * The current directory (".") is always searched first. Then, on Windows
 * platforms, if environment variable COMMONPROGRAMFILES is set (which it
 * should be on Win XP or Win 2000), then directories under this one will be
 * added to the search path. The %Cantera Windows installer installs data
 * files to this location.
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
 *  - ctml: This is an xml file laid out in such a way that %Cantera can
 *          interpret the contents.
 *  - cti:  A human-readable ascii format for information that %Cantera
 *          will read.
 *
 * %Cantera can take its input from both types of files. However, given a file
 * in cti format, the initial operation that %Cantera will perform is to
 * translate the cti file into a ctml file. The translation is carried out via
 * a system call to a python interpreter program that actually carries out the
 * translation. In general, a new ctml file is created by the translation that
 * is written to the current local directory. The ctml file is then read back
 * into %Cantera as the input.
 *
 * Other input routines in other modules:
 *   @see importKinetics()
 * @{
 */

//! @copydoc Application::findInputFile
std::string findInputFile(const std::string& name);

//! @copydoc Application::addDataDirectory
void addDirectory(const std::string& dir);
//@}

//! Delete and free all memory associated with the application
/*!
 * Delete all global data.  It should be called at the end of the
 * application if leak checking is to be done.
 */
void appdelete();

//! @copydoc Application::thread_complete
void thread_complete() ;

//! Returns root directory where %Cantera is installed
/*!
 * @return Returns a string containing the name of the base directory where
 *     %Cantera is installed. If the environmental variable CANTERA_ROOT is
 *     defined, this function will return its value, preferentially.
 *
 * @ingroup inputfiles
 */
std::string canteraRoot();

/*!
 * @defgroup logs Diagnostic Output
 *
 * Writing diagnostic information to the screen or to a file. It is often
 * useful to be able to write diagnostic messages to the screen or to a file.
 * Cantera provides two sets of procedures for this purpose. The first set is
 * designed to write text messages to the screen to document the progress of a
 * complex calculation, such as a flame simulation.The second set writes
 * nested lists in HTML format. This is useful to print debugging output for a
 * complex calculation that calls many different procedures.
 */

/*!
 * @defgroup textlogs Writing messages to the screen
 * @ingroup logs
 */

//! @copydoc Application::Messages::writelog(const std::string&)
void writelog(const std::string& msg);

//! Write a message to the log only if loglevel > 0
inline void writelog(const std::string& msg, int loglevel)
{
    if (loglevel > 0) {
        writelog(msg);
    }
}

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
void writelogf(const char* fmt,...);

//! Write an end of line character to the screen and flush output
void writelogendl();

//! @copydoc Application::Messages::logerror
void error(const std::string& msg);

//! @copydoc Application::warn_deprecated
void warn_deprecated(const std::string& method, const std::string& extra="");

//! @copydoc Application::suppress_deprecation_warnings
void suppress_deprecation_warnings();

//! @copydoc Application::Messages::setLogger
void setLogger(Logger* logwriter);

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

//! @copydoc Application::get_XML_File
XML_Node* get_XML_File(const std::string& file, int debug = 0);

//! @copydoc Application::close_XML_File
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

//! @copydoc Application::Messages::beginLogGroup
void beginLogGroup(const std::string& title, int loglevel=-99);

//! @copydoc Application::Messages::addLogEntry(const std::string&, const std::string&)
void addLogEntry(const std::string& tag, const std::string& value);

//! @copydoc Application::Messages::addLogEntry(const std::string&, doublereal)
void addLogEntry(const std::string& tag, doublereal value);

//! @copydoc Application::Messages::addLogEntry(const std::string&, int)
void addLogEntry(const std::string& tag, int value);

//! @copydoc Application::Messages::addLogEntry(const std::string&)
void addLogEntry(const std::string& msg);

//! @copydoc Application::Messages::endLogGroup
void endLogGroup(const std::string& title="");

//! @copydoc Application::Messages::write_logfile
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
XML_Node* get_XML_NameID(const std::string& nameTarget,
                         const std::string& file_ID,
                         XML_Node* root);

//! Clip *value* such that lower <= value <= upper
template <class T>
inline T clip(const T& value, const T& lower, const T& upper)
{
    return std::max(lower, std::min(upper, value));
}

}

#endif
