/**
 * @file ctxml.h
 */
#ifndef CTC_XML_H
#define CTC_XML_H

#include "clib_defs.h"

#ifdef CANTERA_USE_INTERNAL
#include "cantera/base/config.h"
#else
#include "cantera/base/config.h"
#endif

extern "C" {
    CANTERA_CAPI int xml_new(const char* name);
    CANTERA_CAPI int xml_get_XML_File(const char* file, int debug = 0);
    CANTERA_CAPI int xml_del(int i);
    CANTERA_CAPI int xml_clear();
    CANTERA_CAPI int xml_copy(int i);
    CANTERA_CAPI int xml_assign(int i, int j);
    CANTERA_CAPI int xml_build(int i, const char* file);
    CANTERA_CAPI int xml_preprocess_and_build(int i, const char* file, int debug);
    CANTERA_CAPI int xml_attrib(int i, const char* key, char* value);
    CANTERA_CAPI int xml_addAttrib(int i, const char* key, const char* value);
    CANTERA_CAPI int xml_addComment(int i, const char* comment);
    CANTERA_CAPI int xml_value(int i, char* value);
    CANTERA_CAPI int xml_tag(int i, char* tag);
    CANTERA_CAPI int xml_child(int i, const char* loc);
    CANTERA_CAPI int xml_child_bynumber(int i, int m);
    CANTERA_CAPI int xml_findID(int i, const char* id);
    CANTERA_CAPI int xml_findByName(int i, const char* nm);
    CANTERA_CAPI int xml_nChildren(int i);
    CANTERA_CAPI int xml_addChild(int i, const char* name, const char* value);
    CANTERA_CAPI int xml_addChildNode(int i, int j);
    CANTERA_CAPI int xml_write(int i, const char* file);
    CANTERA_CAPI int xml_removeChild(int i, int j);
    CANTERA_CAPI int ctml_getFloatArray(int i, size_t n, double* data, int iconvert=0);
}

#endif

