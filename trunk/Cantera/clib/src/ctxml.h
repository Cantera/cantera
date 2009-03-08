#ifndef CTC_XML_H
#define CTC_XML_H

#include "clib_defs.h"

extern "C" {  

    int DLL_IMPORT xml_new(const char* name);
    int DLL_IMPORT xml_get_XML_File(const char* file);
    int DLL_IMPORT xml_del(int i);
    int DLL_IMPORT xml_clear();
    int DLL_IMPORT xml_copy(int i);
    int DLL_IMPORT xml_assign(int i, int j);
    int DLL_IMPORT xml_build(int i, const char* file);
    int DLL_IMPORT xml_preprocess_and_build(int i, const char* file);
    int DLL_IMPORT xml_attrib(int i, const char* key, char* value);
    int DLL_IMPORT xml_addAttrib(int i, const char* key, const char* value);
    int DLL_IMPORT xml_addComment(int i, const char* comment);
    int DLL_IMPORT xml_value(int i, char* value);
    int DLL_IMPORT xml_tag(int i, char* tag);
    int DLL_IMPORT xml_child(int i, const char* loc);
    int DLL_IMPORT xml_child_bynumber(int i, int m);
    int DLL_IMPORT xml_findID(int i, const char* id);
    int DLL_IMPORT xml_findByName(int i, const char* nm);
    int DLL_IMPORT xml_nChildren(int i);
    int DLL_IMPORT xml_addChild(int i, const char* name, const char* value);
    int DLL_IMPORT xml_addChildNode(int i, int j);
    int DLL_IMPORT xml_write(int i, const char* file);
    int DLL_IMPORT xml_removeChild(int i, int j);
    int DLL_IMPORT ctml_getFloatArray(int i, int n, double* data, int iconvert=0);
}

#endif

