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

class Logger;
class AnyMap;

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
 * separator (that is, ";" on Windows, ":" elsewhere), then these directories will
 * be added to the search path.
 *
 * Finally, the location where the data files were installed when
 * %Cantera was built is added to the search path.
 *
 * Additional directories may be added by calling function addDirectory.
 *
 * %Cantera input files are written using YAML syntax. For more information on using
 * YAML files in Cantera, see the
 * <a href="https://cantera.org/tutorials/yaml/defining-phases.html">YAML Users' Guide</a>
 * or the <a href="../../../../sphinx/html/yaml/index.html">YAML Input File API Reference</a>.
 *
 * %Cantera provides the `ck2yaml` tool for converting Chemkin-format mechanisms to the
 * YAML format. The scripts `cti2yaml.py` and `ctml2yaml.py` can be used to convert
 * legacy CTI and XML input files (from Cantera 2.6 and earlier) to the YAML format.
 *
 * @{
 */

//! @copydoc Application::findInputFile
std::string findInputFile(const std::string& name);

//! @copydoc Application::addDataDirectory
void addDirectory(const std::string& dir);

//! @copydoc Application::getDataDirectories
std::string getDataDirectories(const std::string& sep);
//! @}

//! @copydoc Application::loadExtension
void loadExtension(const std::string& extType, const std::string& name);

//! Load extensions providing user-defined models from the `extensions` section of the
//! given node. @see Application::loadExtension
//!
//! @since New in Cantera 3.0
void loadExtensions(const AnyMap& node);

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

//! Returns true if Cantera was compiled in debug mode. Used for handling some cases
//! where behavior tested in the test suite changes depending on whether the `NDEBUG`
//! preprocessor macro is defined.
bool debugModeEnabled();

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

//! helper function passing deprecation warning to global handler
void _warn_deprecated(const std::string& method, const std::string& extra="");

//! Print a deprecation warning raised from *method*.
/*!
*  @see Application::warn_deprecated
 * @param method  method name
 * @param msg  Python-style format string with the following arguments
 * @param args  arguments for the format string
 */
template <typename... Args>
void warn_deprecated(const std::string& method, const std::string& msg,
                     const Args&... args) {
    if (sizeof...(args) == 0) {
        _warn_deprecated(method, msg);
    } else {
        _warn_deprecated(method, fmt::format(msg, args...));
    }
}

//! @copydoc Application::suppress_deprecation_warnings
void suppress_deprecation_warnings();

//! helper function passing generic warning to global handler
void _warn(const std::string& warning,
           const std::string& method, const std::string& extra);

//! Print a generic warning raised from *method*.
/*!
 * @see Application::warn
 * @param warning  type of warning; See Logger::warn
 * @param method  method name
 * @param msg  Python-style format string with the following arguments
 * @param args  arguments for the format string
 */
template <typename... Args>
void warn(const std::string& warning, const std::string& method,
          const std::string& msg, const Args&... args) {
    if (sizeof...(args) == 0) {
        _warn(warning, method, msg);
    } else {
        _warn(warning, method, fmt::format(msg, args...));
    }
}

//! Print a user warning raised from *method* as `CanteraWarning`.
/*!
 * @param method  method name
 * @param msg  Python-style format string with the following arguments
 * @param args  arguments for the format string
 */
template <typename... Args>
void warn_user(const std::string& method, const std::string& msg,
               const Args&... args) {
    if (sizeof...(args) == 0) {
        _warn("Cantera", method, msg);
    } else {
        _warn("Cantera", method, fmt::format(msg, args...));
    }
}

//! @copydoc Application::make_deprecation_warnings_fatal
void make_deprecation_warnings_fatal();

//! @copydoc Application::make_warnings_fatal
void make_warnings_fatal();

//! @copydoc Application::suppress_thermo_warnings
void suppress_thermo_warnings(bool suppress=true);

//! @copydoc Application::thermo_warnings_suppressed
bool thermo_warnings_suppressed();

//! @copydoc Application::suppress_warnings
void suppress_warnings();

//! @copydoc Application::warnings_suppressed
bool warnings_suppressed();

//! @copydoc Application::use_legacy_rate_constants
void use_legacy_rate_constants(bool legacy=true);

//! @copydoc Application::legacy_rate_constants_used
bool legacy_rate_constants_used();

//! @copydoc Application::Messages::setLogger
void setLogger(Logger* logwriter);

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

//! Convert a type name to a human readable string, using `boost::core::demangle` if
//! available. Also has a set of short names for some common types.
//! @note Mainly for internal use by AnyMap and Delegator
std::string demangle(const std::type_info& type);

}

#endif
