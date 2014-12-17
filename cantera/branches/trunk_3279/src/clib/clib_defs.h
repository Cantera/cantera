/**
 * @file clib_defs.h
 */
#ifndef CTC_DEFS_H
#define CTC_DEFS_H

#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"
#include <iostream>

#ifdef _WIN32
// Windows (MSVC or MinGW)
# ifdef CANTERA_USE_INTERNAL
#  define CANTERA_CAPI extern __declspec(dllexport)
# else
#  define CANTERA_CAPI extern __declspec(dllimport)
# endif
#else
// Non-Windows platform
# ifdef CANTERA_USE_INTERNAL
#  define CANTERA_CAPI extern
# else
#  define CANTERA_CAPI
# endif
#endif

// Values returned for error conditions
#ifndef ERR
# define ERR -999
#endif

#ifndef DERR
# define DERR -999.999
#endif

namespace Cantera
{

//! Exception handler used at language interface boundaries.
/*!
 * When called from a "catch (...)" block, this function will attempt to save
 * an error message in global error stack and return a value indicating the
 * type of exception caught.
 *
 * @param ctErrorCode Value to return if a CanteraError is caught
 * @param otherErrorCode Value to return if a different exception is caught
 */
template <typename T>
T handleAllExceptions(T ctErrorCode, T otherErrorCode)
{
    // Rethrow the previous exception, then catch a more
    // specific exception type if possible.
    try {
        throw;
    } catch (CanteraError& cterr) {
        cterr.save();
        return ctErrorCode;
    } catch (std::exception& err) {
        std::cerr << "Cantera: caught an instance of "
                  << err.what() << std::endl;
        setError("handleAllExceptions", err.what());
        return otherErrorCode;
    } catch (...) {
        std::cerr << "Cantera: caught an instance of "
                  "an unknown exception type" << std::endl;
        setError("handleAllExceptions", "unknown exception");
        return otherErrorCode;
    }
}

}

#endif
