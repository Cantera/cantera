//! @file fmt.h Wrapper for either system-installed or local headers for fmt
#ifndef CT_FMT_H
#define CT_FMT_H
#include "ct_defs.h"

//! Do not use the fmt macro from fmtlib because it shadows a function of
//! the same name in kinetics/Group.h
#define FMT_NO_FMT_STRING_ALIAS

#if CT_USE_SYSTEM_FMT
  #include <fmt/format.h>
  #include <fmt/printf.h>
  #include <fmt/ostream.h>
#else
  #include "cantera/ext/fmt/format.h"
  #include "cantera/ext/fmt/printf.h"
  #include "cantera/ext/fmt/ostream.h"
#endif

template <typename... Args>
void fmt_append(fmt::memory_buffer& b, const std::string& tmpl, Args... args) {
    format_to(fmt::appender(b), fmt::runtime(tmpl), args...);
}

#if FMT_VERSION > 100000
  #if CT_USE_SYSTEM_FMT
    #include <fmt/ranges.h>
  #else
    #include <fmt/join.h>
  #endif
#endif

#endif
