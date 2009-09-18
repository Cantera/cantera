/**
 * @file ctxml.h
 */
/*
 *      $Id: ctxml.h,v 1.6 2009/07/11 17:16:09 hkmoffa Exp $
 */

#ifndef CTC_XML_H
#define CTC_XML_H

#include "clib_defs.h"

#ifdef CANTERA_USE_INTERNAL
#include "config.h"
#else
#include "cantera/config.h"
#endif

extern "C" {  

    EEXXTT int DLL_CPREFIX xml_new(const char* name);
    EEXXTT int DLL_CPREFIX xml_get_XML_File(const char* file, int debug = 0);
    EEXXTT int DLL_CPREFIX xml_del(int i);
    EEXXTT int DLL_CPREFIX xml_clear();
    EEXXTT int DLL_CPREFIX xml_copy(int i);
    EEXXTT int DLL_CPREFIX xml_assign(int i, int j);
    EEXXTT int DLL_CPREFIX xml_build(int i, const char* file);
    EEXXTT int DLL_CPREFIX xml_preprocess_and_build(int i, const char* file, int debug);
    EEXXTT int DLL_CPREFIX xml_attrib(int i, const char* key, char* value);
    EEXXTT int DLL_CPREFIX xml_addAttrib(int i, const char* key, const char* value);
    EEXXTT int DLL_CPREFIX xml_addComment(int i, const char* comment);
    EEXXTT int DLL_CPREFIX xml_value(int i, char* value);
    EEXXTT int DLL_CPREFIX xml_tag(int i, char* tag);
    EEXXTT int DLL_CPREFIX xml_child(int i, const char* loc);
    EEXXTT int DLL_CPREFIX xml_child_bynumber(int i, int m);
    EEXXTT int DLL_CPREFIX xml_findID(int i, const char* id);
    EEXXTT int DLL_CPREFIX xml_findByName(int i, const char* nm);
    EEXXTT int DLL_CPREFIX xml_nChildren(int i);
    EEXXTT int DLL_CPREFIX xml_addChild(int i, const char* name, const char* value);
    EEXXTT int DLL_CPREFIX xml_addChildNode(int i, int j);
    EEXXTT int DLL_CPREFIX xml_write(int i, const char* file);
    EEXXTT int DLL_CPREFIX xml_removeChild(int i, int j);
    EEXXTT int DLL_CPREFIX ctml_getFloatArray(int i, int n, double* data, int iconvert=0);
}

#endif

