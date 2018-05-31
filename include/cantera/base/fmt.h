//! @file fmt.h Wrapper for either system-installed or local headers for fmt
#ifndef CT_FMT_H
#define CT_FMT_H
#include "ct_defs.h"

//! Do not use the fmt macro from fmtlib because it shadows a function of
//! the same name in kinetics/Group.h
#define FMT_NO_FMT_STRING_ALIAS

//! Use header-only library to avoid relocation issues with linking to the
//! static libfmt.a
#define FMT_HEADER_ONLY

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

#if !defined(FMT_VERSION) || FMT_VERSION < 50000
namespace fmt {
using memory_buffer = MemoryWriter;
}
template <typename... Args>
void format_to(fmt::memory_buffer& b, Args... args) {
    b.write(args...);
}
inline std::string to_string(fmt::memory_buffer& b) {
    return b.str();
}

#endif

#endif
