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
  #include <fmt/format.h>
  #include <fmt/printf.h>
  #include <fmt/ostream.h>
#else
  #include "cantera/ext/fmt/format.h"
  #include "cantera/ext/fmt/printf.h"
  #include "cantera/ext/fmt/ostream.h"
#endif

#if FMT_VERSION < 80000
template <typename... Args>
void fmt_append(fmt::memory_buffer& b, const std::string& tmpl, Args... args) {
    format_to(b, tmpl, args...);
}
namespace fmt {
template <typename T>
T runtime(T arg) {
    return arg;
}
}
#else
template <typename... Args>
void fmt_append(fmt::memory_buffer& b, const std::string& tmpl, Args... args) {
    format_to(fmt::appender(b), fmt::runtime(tmpl), args...);
}
#endif

#if FMT_VERSION > 100000
  #if CT_USE_SYSTEM_FMT
    #include <fmt/ranges.h>
  #else
    #include <fmt/join.h>
  #endif
#endif

#endif
