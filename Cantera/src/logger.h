/// @file logger.h
/// Class Logger
///

#ifndef CT_LOGGER_H
#define CT_LOGGER_H

#include <iostream>
using namespace std;

namespace Cantera {

    ///
    /// Base class for 'loggers' that write text messages to log files.
    ///
    /// This class is used to direct log messages to application- or
    /// environment-specific output. The default is to simply print
    /// the messages to the standard output stream or standard error
    /// stream, but classes may be derived from Logger that implement
    /// other output options. This is important when Cantera is used
    /// in applications that do not display the standard output, such
    /// as MATLAB.  The Cantera MATLAB interface derives a class from
    /// Logger that implements these methods with MATLAB-specific
    /// procedures, insuring that the messages will be passed through
    /// to the user. It would also be possible to derive a class that
    /// displayed the messages in a pop-up window, or redirected them
    /// to a file, etc.
    ///  
    /// To install a logger, call function setLogger (global.h / misc.cpp).
    ///
    /// See the files Cantera/python/src/pylogger.h and
    /// Cantera/matlab/cantera/private/mllogger.h for examples of
    /// deriving logger classes.
    /// @ingroup textlogs
    ///
    class Logger {
    public:

      Logger() {}
      virtual ~Logger() {}

      /// Write a log message. The default behavior is to write to
      /// the standard output. Note that no end-of-line character is
      /// appended to the message, and so if one is desired it must
      /// be included in the string.
	  
	virtual void write(const string& msg) {
	  cout << msg;
        }

        /// Write an error message and quit. The default behavior is
        /// to write to the standard eror stream, and then call
        /// exit(). Note that no end-of-line character is appended to
        /// the message, and so if one is desired it must be included
        /// in the string. Note that this default behavior will
        /// terminate the application Cantera is invoked from (MATLAB,
        /// Excel, etc.) If this is not desired, then derive a class
        /// and reimplement this method.
	    virtual void error(const string& msg) {
            cerr << msg << endl;
            exit(-1);
        }

        /// Return an integer specifying the application environment.
        virtual int env() { return 0; }
    };

}
#endif
