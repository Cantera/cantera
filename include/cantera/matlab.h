/**
 * @file ctxml.h
 */
#ifndef CTC_XML_H
#define CTC_XML_H

#include <stddef.h>

int xml_new(const char* name);
int xml_get_XML_File(const char* file, int debug);
int xml_del(int i);
int xml_clear();
int xml_copy(int i);
int xml_build(int i, const char* file);
int xml_preprocess_and_build(int i, const char* file, int debug);
int xml_attrib(int i, const char* key, char* value);
int xml_addAttrib(int i, const char* key, const char* value);
int xml_addComment(int i, const char* comment);
int xml_value(int i, char* value);
int xml_tag(int i, char* tag);
int xml_child(int i, const char* loc);
int xml_child_bynumber(int i, int m);
int xml_findID(int i, const char* id);
int xml_findByName(int i, const char* nm);
int xml_nChildren(int i);
int xml_addChild(int i, const char* name, const char* value);
int xml_addChildNode(int i, int j);
int xml_write(int i, const char* file);
int xml_removeChild(int i, int j);
int ctml_getFloatArray(int i, size_t n, double* data, int iconvert);

#endif

