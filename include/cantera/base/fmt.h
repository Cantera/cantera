//! @file fmt.h Wrapper for either system-installed or local headers for fmt
#ifndef CT_FMT_H
#define CT_FMT_H
#include "ct_defs.h"

//! Do not use the fmt macro from fmtlib because it shadows a function of
//! the same name in kinetics/Group.h
#define FMT_NO_FMT_STRING_ALIAS

//! Versions 6.2.0 and 6.2.1 of fmtlib do not include this define before they
//! include windows.h, breaking builds on Windows. Fixed in fmtlib 7.0.0 and
//! newer. https://github.com/fmtlib/fmt/pull/1616
#if defined(_WIN32) && !defined(NOMINMAX)
#define NOMINMAX
#endif
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

#if !defined(FMT_VERSION) ||  FMT_VERSION < 80000
template <typename... Args>
void fmt_append(fmt::memory_buffer& b, Args... args) {
    format_to(b, args...);
}
#else
template <typename... Args>
void fmt_append(fmt::memory_buffer& b, Args... args) {
    format_to(fmt::appender(b), args...);
}
#endif

#endif
