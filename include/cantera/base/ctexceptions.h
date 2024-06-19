/**
 * @file ctexceptions.h
 *   Definitions for the classes that are
 *   thrown when %Cantera experiences an error condition
 *   (also contains errorhandling module text - see @ref errGroup).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CTEXCEPTIONS_H
#define CT_CTEXCEPTIONS_H

#include "ct_defs.h"
#include "cantera/base/fmt.h"
#include <exception>

namespace Cantera
{

/**
 * @defgroup debugGroup Errors and Diagnostics
*/

/**
 * @defgroup errGroup Errors
 *
 * @brief Handling of errors and unknown events within %Cantera.
 *
 * The general idea is that exceptions are thrown using the common base class
 * called CanteraError. Derived types of CanteraError characterize what type of
 * error is thrown. A list of all of the thrown errors is kept in the
 * Application class.
 *
 * Any exceptions which are not caught cause a fatal error exit from the
 * program.
 *
 * A group of defines may be used during debugging to assert conditions which
 * should be true. These are named AssertTrace(), AssertThrow(), and
 * AssertThrowMsg(). Examples of their usage is given below.
 *
 * @code
 * AssertTrace(p == OneAtm);
 * AssertThrow(p == OneAtm, "Kinetics::update");
 * AssertThrowMsg(p == OneAtm, "Kinetics::update",
 *              "Algorithm limited to atmospheric pressure");
 * @endcode
 *
 * Their first argument is a boolean. If the boolean is not true, a CanteraError
 * is thrown, with descriptive information indicating where the error occurred.
 * The Assert* checks are skipped if the NDEBUG preprocessor symbol is defined,
 * for example with the compiler option -DNDEBUG.
 * @ingroup debugGroup
 */


//! Base class for exceptions thrown by %Cantera classes.
/*!
 * This class is the base class for exceptions thrown by Cantera. It inherits
 * from std::exception so that normal error handling operations from
 * applications may automatically handle the errors in their own way.
 *
 * @ingroup errGroup
 */
class CanteraError : public std::exception
{
public:
    //! Normal Constructor for the CanteraError base class
    /*!
     * @param procedure Name of the function within which the error was
     *             generated. For member functions, this should be written as
     *             `ClassName::functionName`. For constructors, this should be
     *             `ClassName::ClassName`. Arguments can be specified to
     *             disambiguate overloaded functions, such as
     *             `ClassName::functionName(int, int)`.
     * @param msg  Descriptive string describing the type of error message. This
     *     can be a fmt-style format string (that is, using curly braces to indicate
     *     fields), which will be used with additional arguments to generate a
     *     formatted error message
     * @param args Arguments which will be used to interpolate the format string
     */
    template <typename... Args>
    CanteraError(const string& procedure, const string& msg, const Args&... args)
        : CanteraError(procedure)
    {
        if (sizeof...(args) == 0) {
            msg_ = msg;
        } else {
            msg_ = fmt::format(fmt::runtime(msg), args...);
        }
    }

    //! Destructor for base class does nothing
    virtual ~CanteraError() throw() {};

    //! Get a description of the error
    const char* what() const throw() override;

    //! Method overridden by derived classes to format the error message
    virtual string getMessage() const;

    //! Get the name of the method that threw the exception
    virtual string getMethod() const;

    //! Method overridden by derived classes to indicate their type
    virtual string getClass() const {
        return "CanteraError";
    }

    //! Set the number of stack frames to include when a CanteraError is displayed. By
    //! default, or if the depth is set to 0, no stack information will be shown.
    static void setStackTraceDepth(int depth);

protected:
    //! Protected default constructor discourages throwing errors containing no
    //! information.
    CanteraError() {};

    //! Constructor used by derived classes that override getMessage()
    explicit CanteraError(const string& procedure);

    //! The name of the procedure where the exception occurred
    string procedure_;
    mutable string formattedMessage_; //!< Formatted message returned by what()

private:
    string msg_; //!< Message associated with the exception

    string traceback_; //!< Stack trace to location where exception was thrown
    static int traceDepth_; //!< Number of stack frames to show. 0 to disable.
};


//! Array size error.
/*!
 * This error is thrown if a supplied length to a vector supplied to %Cantera is
 * too small.
 *
 * @ingroup errGroup
 */
class ArraySizeError : public CanteraError
{
public:
    //! Constructor
    /*!
     * The length needed is supplied by the argument, reqd, and the
     * length supplied is given by the argument sz.
     *
     * @param procedure String name for the function within which the error was
     *             generated.
     * @param sz   This is the length supplied to Cantera.
     * @param reqd This is the required length needed by Cantera
     */
    ArraySizeError(const string& procedure, size_t sz, size_t reqd) :
        CanteraError(procedure), sz_(sz), reqd_(reqd) {}

    string getMessage() const override;
    string getClass() const override {
        return "ArraySizeError";
    }

private:
    size_t sz_, reqd_;
};


//! An array index is out of range.
/*!
 *  @ingroup errGroup
 */
class IndexError : public CanteraError
{
public:
    //! Constructor
    /*!
     * This class indicates an out-of-bounds array index.
     *
     * @param func String name for the function within which the error was
     *             generated.
     * @param arrayName name of the corresponding array
     * @param m   This is the value of the out-of-bounds index.
     * @param mmax This is the maximum allowed value of the index. The
     *             minimum allowed value is assumed to be 0.
     */
    IndexError(const string& func, const string& arrayName, size_t m, size_t mmax) :
        CanteraError(func), arrayName_(arrayName), m_(m), mmax_(mmax) {}

    ~IndexError() throw() override {};
    string getMessage() const override;
    string getClass() const override {
        return "IndexError";
    }

private:
    string arrayName_;
    size_t m_, mmax_;
};

//! An error indicating that an unimplemented function has been called
//! @ingroup errGroup
class NotImplementedError : public CanteraError
{
public:
    //! @param func Name of the unimplemented function, such as
    //!     `ClassName::functionName`
    NotImplementedError(const string& func) :
        CanteraError(func, "Not implemented.") {}

    //! Alternative constructor taking same arguments as CanteraError
    template <typename... Args>
    NotImplementedError(const string& func, const string& msg, const Args&... args) :
        CanteraError(func, msg, args...) {}

    string getClass() const override {
        return "NotImplementedError";
    }
};

//! Provides a line number
#define XSTR_TRACE_LINE(s) STR_TRACE_LINE(s)

//! Provides a line number
#define STR_TRACE_LINE(s) #s

//! Provides a string variable containing the file and line number
/*!
 * This is a std:string containing the file name and the line number
 */
#define STR_TRACE (string(__FILE__) + ":" + XSTR_TRACE_LINE(__LINE__))

#ifdef NDEBUG
#ifndef AssertTrace
#  define AssertTrace(expr)                        ((void) (0))
#endif
#ifndef AssertThrow
#  define AssertThrow(expr, procedure)             ((void) (0))
#endif
#ifndef AssertThrowMsg
#  define AssertThrowMsg(expr,procedure, ...)  ((void) (0))
#endif
#else

//! Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a CanteraError is thrown. A diagnostic string
 * containing the file and line number, indicating where the error occurred is
 * added to the thrown object.
 *
 * @param expr  Boolean expression that must be true
 *
 * @ingroup errGroup
 */
#ifndef AssertTrace
#  define AssertTrace(expr)  ((expr) ? (void) 0 : throw CanteraError(STR_TRACE, string("failed assert: ") + #expr))
#endif

//! Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a CanteraError is thrown. A diagnostic string
 * indicating where the error occurred is added to the thrown object.
 *
 * @param expr  Boolean expression that must be true
 * @param procedure  Character string or std:string expression indicating the
 *     procedure where the assertion failed
 * @ingroup errGroup
 */
#ifndef AssertThrow
#  define AssertThrow(expr, procedure)   ((expr) ? (void) 0 : throw CanteraError(procedure, string("failed assert: ") + #expr))
#endif

//! Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a CanteraError is thrown. A diagnostic string
 * indicating where the error occurred is added to the thrown object.
 *
 * @param expr  Boolean expression that must be true
 * @param procedure  Character string or std:string expression indicating
 *                   the procedure where the assertion failed
 * Additional arguments are passed on to the constructor for CanteraError to
 * generate a formatted error message.
 *
 * @ingroup errGroup
 */
#ifndef AssertThrowMsg
#  define AssertThrowMsg(expr, procedure, ...)  ((expr) ? (void) 0 : throw CanteraError(procedure + string(":\nfailed assert: \"") + string(#expr) + string("\""), __VA_ARGS__))
#endif

#endif

//! Throw an exception if the specified exception is not a finite number.
//! @ingroup errGroup
#ifndef AssertFinite
#  define AssertFinite(expr, procedure, ...) AssertThrowMsg(expr < BigNumber && expr > -BigNumber, procedure, __VA_ARGS__)
#endif

}

#endif
