/**
 * @file global.h
 * This file contains definitions for utility functions and text for modules,
 * inputfiles and logging, (see @ref inputGroup, and @ref logGroup).
 *
 * These functions store some parameters in global storage that are accessible
 * at all times from the calling application. Contains module definitions for
 *     -  inputfiles  (see @ref inputGroup)
 *     -  text logs  (see @ref logGroup)
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

//! @defgroup ioGroup File Input/Output
//! @details Classes and functions used for reading and writing of %Cantera input files.

/**
 * @defgroup inputGroup Input File Handling
 * @brief Handling of %Cantera input files.
 *
 * The properties of phases and interfaces are specified in text files. These
 * procedures handle various aspects of reading these files. Input files should be read
 * using the newSolution() or newInterface() service functions (see @ref solnGroup).
 *
 * For input files not specified by an absolute pathname, %Cantera searches
 * for input files along a path that includes platform-specific default
 * locations, and possibly user-specified locations:
 *
 *    - The current directory @c "." is always searched first. Then, on Windows, the
 *      registry is checked to find the %Cantera installation directory, and the
 *      @c data subdirectory of the installation directory will be added to the search
 *      path.
 *    - On any platform, if environment variable @c CANTERA_DATA is set to a directory
 *      name or a list of directory names separated with the OS-dependent path
 *      separator (that is, @c ";" on Windows, @c ":" elsewhere), then these directories
 *      will be added to the search path.
 *    - Finally, the location where the data files were installed when %Cantera was
 *      built is added to the search path.
 *    - Additional directories may be added dynamically by calling the addDirectory()
 *      function.
 *
 * %Cantera input files are written using YAML syntax. For more information on using
 * YAML files in %Cantera, see the [YAML Users' Guide](../userguide/input-tutorial.html)
 * or the [YAML Input File API Reference](../yaml/index.html). %Cantera provides the
 * [`ck2yaml`](../userguide/ck2yaml-tutorial.html) tool for converting Chemkin-format
 * mechanisms to the YAML format. The utilities [`cti2yaml`](../yaml/cti2yaml.html) and
 * [`ctml2yaml`](../yaml/ctml2yaml.html) should be used to convert legacy CTI and XML
 * input files (from %Cantera 2.6 and earlier) to the YAML format.
 *
 * @ingroup ioGroup
 * @{
 */

//! @copydoc Application::findInputFile
string findInputFile(const string& name);

//! @copydoc Application::addDataDirectory
void addDirectory(const string& dir);

//! @copydoc Application::getDataDirectories
string getDataDirectories(const string& sep);
//! @}

//! @copydoc Application::loadExtension
void loadExtension(const string& extType, const string& name);

//! Load extensions providing user-defined models from the `extensions` section of the
//! given node. @see Application::loadExtension
//!
//! @since New in %Cantera 3.0
void loadExtensions(const AnyMap& node);

//! @copydoc Application::searchPythonVersions
void searchPythonVersions(const string& versions);

//! Delete and free all memory associated with the application
/*!
 * Delete all global data.  It should be called at the end of the
 * application if leak checking is to be done.
 */
void appdelete();

//! @copydoc Application::thread_complete
void thread_complete();

//! @defgroup globalSettings  Global Cantera Settings
//! @brief Functions for accessing global %Cantera settings.
//! @ingroup globalData

//! @addtogroup globalSettings
//! @{

//! Returns `true` if %Cantera was loaded as a shared library in the current
//! application. Returns `false` if it was statically linked.
//! @since New in %Cantera 3.0
bool usingSharedLibrary();

//! @name %Cantera Version Information
//! @{

//! Returns the %Cantera version. This function is used to access the version from a
//! library, rather than the @c CANTERA_VERSION macro that is available at compile time.
//! @since New in %Cantera 3.0
string version();

//! Returns the hash of the git commit from which %Cantera was compiled, if known
string gitCommit();

//! @}

//! Returns true if %Cantera was compiled in debug mode. Used for handling some cases
//! where behavior tested in the test suite changes depending on whether the `NDEBUG`
//! preprocessor macro is defined.
bool debugModeEnabled();

//! Returns true if %Cantera was compiled with C++ @c HDF5 support.
//! @since New in %Cantera 3.0.
bool usesHDF5();

//! @}

/**
 * @defgroup logGroup Logging
 * @brief Logging and generation of diagnostic output.
 *
 * Writing diagnostic information to the screen or to a file. It is often
 * useful to be able to write diagnostic messages to the screen or to a file.
 * %Cantera defines a set of procedures for this purpose designed to write text messages
 * to the screen to document the progress of a complex calculation, such as a
 * flame simulation.
 * @ingroup debugGroup
 */

//! @addtogroup logGroup
//! @{

//! @copydoc Application::Messages::writelog(const string&)
void writelog_direct(const string& msg);

//! Write a message to the log only if loglevel > 0
inline void debuglog(const string& msg, int loglevel)
{
    if (loglevel > 0) {
        writelog_direct(msg);
    }
}

//! Write a formatted message to the screen.
//!
//! This function passes its arguments to the fmt library 'format' function to
//! generate a formatted string from a Python-style (curly braces) format
//! string. This method is used throughout %Cantera to write log messages. It can
//! also be called by user programs. The advantage of using writelog over
//! writing directly to the standard output is that messages written with
//! writelog will display correctly even when %Cantera is used from MATLAB or
//! other application that do not have a standard output stream.
template <typename... Args>
void writelog(const string& fmt, const Args&... args) {
    if (sizeof...(args) == 0) {
        writelog_direct(fmt);
    } else {
        writelog_direct(fmt::format(fmt::runtime(fmt), args...));
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
 */
template <typename... Args>
void writelogf(const char* fmt, const Args& ... args) {
    writelog_direct(fmt::sprintf(fmt, args...));
}

//! Write an end of line character to the screen and flush output
void writelogendl();

void writeline(char repeat, size_t count,
               bool endl_after=true, bool endl_before=false);

//! @} End of logging group

//! @defgroup warnGroup  Warnings
//! @ingroup debugGroup
//! @brief Warnings raised by %Cantera
//! @{

//! @cond

//! helper function passing deprecation warning to global handler
void _warn_deprecated(const string& method, const string& extra="");

//! @endcond

//! Print a deprecation warning raised from *method*.
/*!
*  @see Application::warn_deprecated
 * @param method  method name
 * @param msg  Python-style format string with the following arguments
 * @param args  arguments for the format string
 */
template <typename... Args>
void warn_deprecated(const string& method, const string& msg, const Args&... args) {
    if (sizeof...(args) == 0) {
        _warn_deprecated(method, msg);
    } else {
        _warn_deprecated(method, fmt::format(fmt::runtime(msg), args...));
    }
}

//! @cond

//! helper function passing generic warning to global handler
void _warn(const string& warning, const string& method, const string& extra);

//! @endcond

//! Print a generic warning raised from *method*.
/*!
 * @see Application::warn
 * @param warning  type of warning; See Logger::warn
 * @param method  method name
 * @param msg  Python-style format string with the following arguments
 * @param args  arguments for the format string
 */
template <typename... Args>
void warn(const string& warning, const string& method,
          const string& msg, const Args&... args) {
    if (sizeof...(args) == 0) {
        _warn(warning, method, msg);
    } else {
        _warn(warning, method, fmt::format(fmt::runtime(msg), args...));
    }
}

//! Print a user warning raised from *method* as `CanteraWarning`.
/*!
 * @param method  method name
 * @param msg  Python-style format string with the following arguments
 * @param args  arguments for the format string
 */
template <typename... Args>
void warn_user(const string& method, const string& msg, const Args&... args) {
    if (sizeof...(args) == 0) {
        _warn("Cantera", method, msg);
    } else {
        _warn("Cantera", method, fmt::format(fmt::runtime(msg), args...));
    }
}

//! @} End of warnings group

//! @addtogroup globalSettings
//! @{

//! @name Global Warning Settings
//! @{

//! @copydoc Application::suppress_deprecation_warnings
void suppress_deprecation_warnings();

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

//! @} End of warning settings

//! @copydoc Application::use_legacy_rate_constants
void use_legacy_rate_constants(bool legacy=true);

//! @copydoc Application::legacy_rate_constants_used
bool legacy_rate_constants_used();

// @} End of globalSettings group

//! @copydoc Application::Messages::setLogger
//! @ingroup logGroup
void setLogger(Logger* logwriter);

//! Enables printing a stacktrace to `std::err` if a segfault occurs. The Boost
//! documentation says doing this from an error handler is not safe on all platforms
//! and risks deadlocking. However, it can be useful for debugging and is therefore
//! enabled when running tests.
//! @since New in %Cantera 3.0
//! @ingroup globalSettings
void printStackTraceOnSegfault();

//! Clip *value* such that lower <= value <= upper
//! @ingroup mathTemplates
template <class T>
inline T clip(const T& value, const T& lower, const T& upper)
{
    return std::max(lower, std::min(upper, value));
}

//! Sign of a number. Returns -1 if x < 0, 1 if x > 0 and 0 if x == 0.
//! @ingroup mathTemplates
template <typename T> int sign(T x) {
    return (T(0) < x) - (x < T(0));
}

//! Convert a type name to a human readable string, using `boost::core::demangle` if
//! available. Also has a set of short names for some common types.
//! @note Mainly for internal use by AnyMap and Delegator
string demangle(const std::type_info& type);

}

#endif
