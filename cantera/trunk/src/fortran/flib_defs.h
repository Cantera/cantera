#ifndef FCTC_DEFS_H
#define FCTC_DEFS_H

// Build as a DLL under Windows
#ifdef _WIN32
#define DLL_IMPORT __declspec(dllimport)
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#define DLL_IMPORT
#endif

// Values returned for error conditions
#define ERR -999
#define DERR -999.999

#include "cantera/base/config.h"

typedef integer status_t;

namespace Cantera {}
namespace std {}

#endif
