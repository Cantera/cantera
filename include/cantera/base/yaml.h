#include "cantera/base/ct_defs.h"
#if CT_USE_SYSTEM_YAMLCPP
#include "yaml-cpp/yaml.h"
#else
#ifndef YAML_CPP_STATIC_DEFINE
#define YAML_CPP_STATIC_DEFINE
#endif
#include "cantera/ext/yaml-cpp/yaml.h"
#endif
