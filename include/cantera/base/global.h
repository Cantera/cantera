/**
 * @file global.h
 * This file contains definitions for utility functions and text for modules,
 * inputfiles, logs, textlogs, (see \ref inputfiles, \ref logs, and
 * \ref textlogs).
 *
 * @ingroup utils
 *
 * These functions store some parameters in global storage that are accessible
 * at all times from the calling application. Contains module definitions for
 *     -  inputfiles  (see \ref inputfiles)
 *     -  logs        (see \ref logs)
 *     -  textlogs    (see \ref textlogs)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_GLOBAL_H
#define CT_GLOBAL_H

#include "ct_defs.h"
#include "cantera/base/fmt.h"

namespace Cantera
{

class XML_Node;
class Logger;

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
 * The current directory (".") is always searched first. Then, on Windows, the
 * registry is checked to find the Cantera installation directory, and the
 * 'data' subdirectory of the installation directory will be added to the search
 * path.
 *
 * On the Mac, directory '/Applications/Cantera/data' is added to the
 * search path.
 *
 * On any platform, if environment variable CANTERA_DATA is set to a directory
 * name or a list of directory names separated with the OS-dependent path
 * separator (i.e. ";" on Windows, ":" elsewhere), then these directories will
 * be added to the search path.
 *
 * Finally, the location where the data files were installed when
 * %Cantera was built is added to the search path.
 *
 * Additional directories may be added by calling function addDirectory.
 *
 * There are currently three different types of input files within %Cantera. The
 * YAML format is new in Cantera 2.5, and replaces both the CTI and CTML (XML)
 * formats, which are deprecated and will be removed in Cantera 3.0. The scripts
 * `cti2yaml.py` and `ctml2yaml.py` can be used to convert legacy input files to
 * the YAML format.
 *
 *  - YAML: A human-readable input file written using YAML syntax which
 *    defines species, phases, and reactions, and contains thermodynamic,
 *    chemical kinetic, and transport data needed by %Cantera.
 *
 *  - CTI: A human-readable input file written using Python syntax. Some options
 *    for non-ideal equations of state available in the CTML format have not
 *    been implemented for the CTI format.
 *
 *  - CTML: This is an XML file laid out in such a way that %Cantera can
 *    interpret the contents directly. Given a file in CTI format, %Cantera will
 *    convert the CTI file into the CTML format on-the-fly using a Python script
 *    (ctml_writer). This process is done in-memory without writing any new
 *    files to disk. Explicit use of the CTML format is not recommended.
 *
 * %Cantera provides converters (`ck2yaml` and `ckcti`) for converting
 * Chemkin-format mechanisms to the YAML and CTI formats, respectively.
 *
 * Other input routines in other modules:
 *   @see importKinetics()
 *
 * @{
 */

//! @copydoc Application::findInputFile
std::string findInputFile(const std::string& name);

//! @copydoc Application::addDataDirectory
void addDirectory(const std::string& dir);

//! @copydoc Application::getDataDirectories
std::string getDataDirectories(const std::string& sep);
//@}

//! Delete and free all memory associated with the application
/*!
 * Delete all global data.  It should be called at the end of the
 * application if leak checking is to be done.
 */
void appdelete();

//! @copydoc Application::thread_complete
void thread_complete();

//! Returns the hash of the git commit from which Cantera was compiled, if known
std::string gitCommit();

//! Returns root directory where %Cantera is installed
/*!
 * @returns a string containing the name of the base directory where %Cantera is
 *     installed. If the environmental variable CANTERA_ROOT is defined, this
 *     function will return its value, preferentially.
 *
 * @ingroup inputfiles
 */
std::string canteraRoot();

/*!
 * @defgroup logs Diagnostic Output
 *
 * Writing diagnostic information to the screen or to a file. It is often
 * useful to be able to write diagnostic messages to the screen or to a file.
 * Cantera a set of procedures for this purpose designed to write text messages
 * to the screen to document the progress of a complex calculation, such as a
 * flame simulation.
 */

/*!
 * @defgroup textlogs Writing messages to the screen
 * @ingroup logs
 */

//! @copydoc Application::Messages::writelog(const std::string&)
void writelog_direct(const std::string& msg);

//! Write a message to the log only if loglevel > 0
inline void debuglog(const std::string& msg, int loglevel)
{
    if (loglevel > 0) {
        writelog_direct(msg);
    }
}

//! Write a formatted message to the screen.
//!
//! This function passes its arguments to the fmt library 'format' function to
//! generate a formatted string from a Python-style (curly braces) format
//! string. This method is used throughout Cantera to write log messages. It can
//! also be called by user programs. The advantage of using writelog over
//! writing directly to the standard output is that messages written with
//! writelog will display correctly even when Cantera is used from MATLAB or
//! other application that do not have a standard output stream.
//! @ingroup textlogs
template <typename... Args>
void writelog(const std::string& fmt, const Args&... args) {
    if (sizeof...(args) == 0) {
        writelog_direct(fmt);
    } else {
        writelog_direct(fmt::format(fmt, args...));
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
 * @param args arguments used to interpolate the format string
 * @ingroup textlogs
 */
template <typename... Args>
void writelogf(const char* fmt, const Args& ... args) {
    writelog_direct(fmt::sprintf(fmt, args...));
}

//! Write an end of line character to the screen and flush output
void writelogendl();

void writeline(char repeat, size_t count,
               bool endl_after=true, bool endl_before=false);

//! @copydoc Application::warn_deprecated
void warn_deprecated(const std::string& method, const std::string& extra="");

//! @copydoc Application::suppress_deprecation_warnings
void suppress_deprecation_warnings();

//! helper function passing user warning to global handler
void _warn_user(const std::string& method, const std::string& extra);

/*!
 * Print a user warning raised from *method*.
 *
 * @param method  method name
 * @param msg  Python-style format string with the following arguments
 * @param args  arguments for the format string
 */
template <typename... Args>
void warn_user(const std::string& method, const std::string& msg,
               const Args&... args) {
    if (sizeof...(args) == 0) {
        _warn_user(method, msg);
    } else {
        _warn_user(method, fmt::format(msg, args...));
    }
}

//! @copydoc Application::make_deprecation_warnings_fatal
void make_deprecation_warnings_fatal();

//! @copydoc Application::suppress_thermo_warnings
void suppress_thermo_warnings(bool suppress=true);

//! @copydoc Application::thermo_warnings_suppressed
bool thermo_warnings_suppressed();

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

//! @copydoc Application::get_XML_from_string
XML_Node* get_XML_from_string(const std::string& text);

//! @copydoc Application::close_XML_File
void close_XML_File(const std::string& file);

//! This routine will locate an XML node in either the input
//! XML tree or in another input file specified by the file
//! part of the file_ID string.
/*!
 * Searches are based on the ID attribute of the XML element only.
 *
 * @param file_ID This is a concatenation of two strings separated by the "#"
 *                character. The string before the pound character is the file
 *                name of an XML file to carry out the search. The string after
 *                the # character is the ID attribute of the XML element to
 *                search for. The string is interpreted as a file string if no #
 *                character is in the string.
 * @param root    If the file string is empty, searches for the XML element with
 *                matching ID attribute are carried out from this XML node.
 * @returns the XML_Node, if found. Returns null if not found.
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
 */
XML_Node* get_XML_Node(const std::string& file_ID, XML_Node* root);


//! This routine will locate an XML node in either the input XML tree or in
//! another input file specified by the file part of the file_ID string.
/*!
 * Searches are based on the XML element name and the ID attribute of the XML
 * element. An exact match of both is usually required. However, the ID
 * attribute may be set to "", in which case the first XML element with the
 * correct element name will be returned.
 *
 * @param nameTarget This is the XML element name to look for.
 * @param file_ID This is a concatenation of two strings separated by the "#"
 *                character. The string before the pound character is the file
 *                name of an XML file to carry out the search. The string after
 *                the # character is the ID attribute of the XML element to
 *                search for. The string is interpreted as a file string if no #
 *                character is in the string.
 * @param root    If the file string is empty, searches for the XML element with
 *                matching ID attribute are carried out from this XML node.
 * @returns the XML_Node, if found. Returns null if not found.
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
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

//! Sign of a number. Returns -1 if x < 0, 1 if x > 0 and 0 if x == 0.
template <typename T> int sign(T x) {
    return (T(0) < x) - (x < T(0));
}

}

#endif
