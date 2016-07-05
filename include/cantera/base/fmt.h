//! @file fmt.h Wrapper for either system-installed or local headers for fmt
#include "ct_defs.h"

#if CT_USE_SYSTEM_FMT
  #include "fmt/format.h"
  #include "fmt/ostream.h"
#else
  #include "cantera/ext/fmt/format.h"
  #include "cantera/ext/fmt/ostream.h"
#endif
