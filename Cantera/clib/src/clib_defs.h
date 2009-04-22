#ifndef CTC_DEFS_H
#define CTC_DEFS_H

// Build as a DLL under Windows
#ifdef WIN32
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
#define DLL_EXPORT
#define DLL_IMPORT
#endif

// Values returned for error conditions
#define ERR -999
#define DERR -999.999

namespace Cantera {}
using namespace Cantera;
namespace std {}
using namespace std;

#endif
