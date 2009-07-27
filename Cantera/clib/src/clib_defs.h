/**
 * @file clib_defs.h
 */
/*
 *      $Id: clib_defs.h,v 1.5 2009/07/11 17:16:09 hkmoffa Exp $
 */


#ifndef CTC_DEFS_H
#define CTC_DEFS_H


#ifdef WIN32
//   Either build as a DLL under Windows or  not. 
//   the decision relies upon whether the NO_DLL_BUILD define is 
//   set or not.
#ifdef NO_DLL_BUILD
#define DLL_EXPORT
#define DLL_IMPORT
#else
#define DLL_IMPORT __declspec(dllimport)
#define DLL_EXPORT __declspec(dllexport)
#endif

#pragma warning(disable:4786)
#pragma warning(disable:4267)
#pragma warning(disable:4503)
#else
//   On other platforms, we turn off the DLL macros.
#define DLL_EXPORT
#define DLL_IMPORT
#endif

#ifdef CANTERA_USE_INTERNAL
#define DLL_CPREFIX DLL_EXPORT
#define EEXXTT extern
#else
#define DLL_CPREFIX DLL_IMPORT
#define EEXXTT  
#endif

// Values returned for error conditions
#ifndef ERR
#define ERR -999
#endif
#ifndef DERR
#define DERR -999.999
#endif

namespace Cantera {}
//using namespace Cantera;
namespace std {}
//using namespace std;

#endif
