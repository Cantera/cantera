/**
 * @file ctexceptions.h
 *   Definitions for the classes that are
 *   thrown when %Cantera experiences an error condition
 *   (also contains errorhandling module text - see \ref errorhandling).
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_CTEXCEPTIONS_H
#define CT_CTEXCEPTIONS_H

#include <exception>
#include <string>

namespace Cantera
{

/*!
 * @defgroup errorhandling Error Handling
 *
 * \brief These classes and related functions are used to handle errors and
 *        unknown events within Cantera.
 *
 *  The general idea is that exceptions are thrown using the common
 *  base class called CanteraError. Derived types of CanteraError
 *  characterize what type of error is thrown. A list of all
 *  of the thrown errors is kept in the Application class.
 *
 *  Any exceptions which are not caught cause a fatal error exit
 *  from the program.
 *
 *  Below is an example of how to catch errors that throw the CanteraError class.
 *  In general, all Cantera C++ programs will have this basic structure.
 *
 *  \include demo1a.cpp
 *
 *  The function showErrors() will print out the fatal error
 *  condition to standard output.
 *
 *  A group of defines may be used during debugging to assert
 *  conditions which should be true. These are named AssertTrace(),
 *  AssertThrow(), and AssertThrowMsg().  Examples of their usage is
 *  given below.
 *
 * @code
 *       AssertTrace(p == OneAtm);
 *       AssertThrow(p == OneAtm, "Kinetics::update");
 *       AssertThrowMsg(p == OneAtm, "Kinetics::update",
 *                    "Algorithm limited to atmospheric pressure");
 * @endcode
 *
 *  Their first argument is a boolean. If the boolean is not true, a
 *  CanteraError is thrown, with descriptive information indicating
 *  where the error occurred. The Assert* checks are skipped if the NDEBUG
 *  preprocessor symbol is defined, e.g. with the compiler option -DNDEBUG.
 */

//! Enum containing Cantera's behavior for situations where overflow or underflow of real variables
//! may occur.
/*!
 * Note this frequently occurs when taking exponentials of delta Gibbs energies of reactions
 * or when taking the exponentials of logs of activity coefficients.
 */
enum CT_RealNumber_Range_Behavior {

   //! For this specification of range behavior, nothing is done. This is the fastest
   //! behavior when all calculations are believed to be ranged well. For situations
   //! where there are range errors, NaN's or INF's will be introduced.
   DONOTHING_CTRB = -1,

   //! For this specification of range behavior, the overflow or underflow calculation is changed.
   //! Cantera will proceed by bounding the real number to maintain its viability, silently
   //! changing the actual answer.
   CHANGE_OVERFLOW_CTRB,

   //! When an overflow or underflow occurs, Cantera will throw an error
   THROWON_OVERFLOW_CTRB,

   //! Cantera will use the fenv check capability introduced in C99 to check for
   //! overflow and underflow conditions at crucial points.
   //! It will throw an error if these conditions occur.
   FENV_CHECK_CTRB,

   //!  Cantera will throw an error in debug mode but will not in production mode.
   //! (default)
   THROWON_OVERFLOW_DEBUGMODEONLY_CTRB
};


//! Base class for exceptions thrown by Cantera classes.
/*!
 * This class is the base class for exceptions thrown by Cantera.
 * It inherits from std::exception so that normal error handling
 * operations from applications may automatically handle the
 * errors in their own way.
 *
 * @ingroup errorhandling
 */
class CanteraError : public std::exception
{
public:
    //! Normal Constructor for the CanteraError base class
    /*!
     * In the constructor, a call to the Application class is made to store
     * the strings associated with the generated error condition.
     *
     * @param procedure String name for the function within which the error was
     *             generated.
     * @param msg  Descriptive string describing the type of error message.
     */
    CanteraError(const std::string& procedure, const std::string& msg);

    //! Destructor for base class does nothing
    virtual ~CanteraError() throw() {};

    //! Get a description of the error
    const char* what() const throw();

    //! Function to put this error onto Cantera's error stack
    void save();

    //! Method overridden by derived classes to format the error message
    virtual std::string getMessage() const;

    //! Method overridden by derived classes to indicate their type
    virtual std::string getClass() const {
        return "CanteraError";
    }

protected:
    //! Protected default constructor discourages throwing errors containing no information.
    CanteraError() : saved_(false) {};

    //! Constructor used by derived classes that override getMessage()
    explicit CanteraError(const std::string& procedure);

    //! The name of the procedure where the exception occurred
    std::string procedure_;
    mutable std::string formattedMessage_; //!< Formatted message returned by what()

private:
    std::string msg_; //!< Message associated with the exception
    bool saved_; //!< Exception has already been saved to Cantera's error stack
};


//! Array size error.
/*!
 *  This error is thrown if a supplied length to a vector supplied
 *  to Cantera is too small.
 *
 *  @ingroup errorhandling
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
    ArraySizeError(const std::string& procedure, size_t sz, size_t reqd) :
        CanteraError(procedure), sz_(sz), reqd_(reqd) {}

    virtual std::string getMessage() const;
    virtual std::string getClass() const {
        return "ArraySizeError";
    }

private:
    size_t sz_, reqd_;
};


//! An array index is out of range.
/*!
 *  @ingroup errorhandling
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
    IndexError(const std::string& func, const std::string& arrayName, size_t m, size_t mmax) :
        CanteraError(func), arrayName_(arrayName), m_(m), mmax_(mmax) {}

    virtual ~IndexError() throw() {};
    virtual std::string getMessage() const;
    virtual std::string getClass() const {
        return "IndexError";
    }

private:
    std::string arrayName_;
    size_t m_, mmax_;
};

//! An error indicating that an unimplemented function has been called
class NotImplementedError : public CanteraError
{
public:
    NotImplementedError(const std::string& func) :
        CanteraError(func, "Not implemented.") {}

    virtual std::string getClass() const {
        return "NotImplementedError";
    }
};

//! Quick check on whether there has been an underflow or overflow condition in the floating point unit
/*!
 *    @return Returns true if there has been such a condition and it has not been cleared. returns false
 *            if there hasn't been an overflow, underflow or invalid condition.
 */
bool check_FENV_OverUnder_Flow();

//! Clear all the flags for floating-point exceptions
void clear_FENV();


//! Provides a line number
#define XSTR_TRACE_LINE(s) STR_TRACE_LINE(s)

//! Provides a line number
#define STR_TRACE_LINE(s) #s

//! Provides a std::string variable containing the file and line number
/*!
 *   This is a std:string containing the file name and the line number
 */
#define STR_TRACE   (std::string(__FILE__) +  ":" + XSTR_TRACE_LINE(__LINE__))

#ifdef NDEBUG
#ifndef AssertTrace
#  define AssertTrace(expr)                        ((void) (0))
#endif
#ifndef AssertThrow
#  define AssertThrow(expr, procedure)             ((void) (0))
#endif
#ifndef AssertThrowMsg
#  define AssertThrowMsg(expr,procedure, message)  ((void) (0))
#endif
#else

//! Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a CanteraError is thrown. A diagnostic string containing the
 * file and line number,  indicating where the error
 * occurred is added to the thrown object.
 *
 * @param expr  Boolean expression that must be true
 *
 * @ingroup errorhandling
 */
#ifndef AssertTrace
#  define AssertTrace(expr)  ((expr) ? (void) 0 : throw Cantera::CanteraError(STR_TRACE, std::string("failed assert: ") + #expr))
#endif

//!  Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a CanteraError is thrown. A diagnostic string indicating where the error
 * occurred is added to the thrown object.
 *
 * @param expr  Boolean expression that must be true
 * @param procedure  Character string or std:string expression indicating the procedure where the assertion failed
 * @ingroup errorhandling
 */
#ifndef AssertThrow
#  define AssertThrow(expr, procedure)   ((expr) ? (void) 0 : throw Cantera::CanteraError(procedure, std::string("failed assert: ") + #expr))
#endif

//!  Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a CanteraError is thrown. A
 * diagnostic string indicating where the error occurred is added
 * to the thrown object.
 *
 * @param expr  Boolean expression that must be true
 * @param procedure  Character string or std:string expression indicating
 *                   the procedure where the assertion failed
 * @param message  Character string or std:string expression containing
 *    a descriptive message is added to the thrown error condition.
 *
 * @ingroup errorhandling
 */
#ifndef AssertThrowMsg
#  define AssertThrowMsg(expr, procedure, message)  ((expr) ? (void) 0 : throw Cantera::CanteraError(procedure + std::string(": at failed assert: \"") + std::string(#expr) + std::string("\""), message))
#endif

#endif

//! Throw an exception if the specified exception is not a finite number.
#ifndef AssertFinite
#  define AssertFinite(expr, procedure, message) AssertThrowMsg(expr < BigNumber && expr > -BigNumber, procedure, message)
#endif

}

#endif
