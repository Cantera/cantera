//! @file fmt.h Wrapper for either system-installed or local headers for fmt
#include "ct_defs.h"

//! Do not use the fmt macro from fmtlib because it shadows a function of
//! the same name in kinetics/Group.h
#define FMT_NO_FMT_STRING_ALIAS

#if CT_USE_SYSTEM_FMT
  #include "fmt/format.h"
  #if defined(FMT_VERSION) && FMT_VERSION >= 40000
    #include "fmt/printf.h"
  #endif
  #include "fmt/ostream.h"
#else
  #include "cantera/ext/fmt/format.h"
  #if defined(FMT_VERSION) && FMT_VERSION >= 40000
    #include "cantera/ext/fmt/printf.h"
  #endif
  #include "cantera/ext/fmt/ostream.h"
#endif
