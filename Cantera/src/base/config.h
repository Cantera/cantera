#ifndef SRC_CONFIG_H
#define SRC_CONFIG_H
#ifdef _WIN32
#undef WIN32
#define WIN32
#endif
#ifdef WIN32
#ifdef CANTERA_APP
#include "../winconfig.h"
#else
#include "../../../winconfig.h"
#endif
#else
#ifdef CANTERA_APP
#include "../config.h"
#else
#include "../../../config.h"
#endif
#endif
#endif
