/**
 * @file ctxml.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_XML_H
#define CTC_XML_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int xml_new(const char* name);
    CANTERA_CAPI int xml_get_XML_File(const char* file, int debug);
    CANTERA_CAPI int xml_del(int i);
    CANTERA_CAPI int xml_copy(int i);
    CANTERA_CAPI int xml_build(int i, const char* file);
    CANTERA_CAPI int xml_attrib(int i, const char* key, size_t lenval, char* value);
    CANTERA_CAPI int xml_addAttrib(int i, const char* key, const char* value);
    CANTERA_CAPI int xml_addComment(int i, const char* comment);
    CANTERA_CAPI int xml_value(int i, size_t lenval, char* value);
    CANTERA_CAPI int xml_tag(int i, size_t lentag, char* tag);
    CANTERA_CAPI int xml_child(int i, const char* loc);
    CANTERA_CAPI int xml_child_bynumber(int i, int m);
    CANTERA_CAPI int xml_findID(int i, const char* id);
    CANTERA_CAPI int xml_findByName(int i, const char* nm);
    CANTERA_CAPI int xml_nChildren(int i);
    CANTERA_CAPI int xml_addChild(int i, const char* name, const char* value);
    CANTERA_CAPI int xml_addChildNode(int i, int j);
    CANTERA_CAPI int xml_write(int i, const char* file);
    CANTERA_CAPI int xml_removeChild(int i, int j);
    CANTERA_CAPI int ct_clearXML();

#ifdef __cplusplus
}
#endif

#endif
