/**
 * @file ctxml.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#define CANTERA_USE_INTERNAL
#include "cantera/clib/ctxml.h"

// Cantera includes
#include "cantera/base/ctml.h"
#include "Cabinet.h"

#include <string.h>
#include <fstream>

using namespace std;
using namespace Cantera;

typedef Cabinet<XML_Node, false> XmlCabinet;
template<> XmlCabinet* XmlCabinet::s_storage = 0;

extern "C" {

    int xml_new(const char* name = 0)
    {
        try {
            XML_Node* x;
            if (!name) {
                x = new XML_Node;
            } else {
                x = new XML_Node(name);
            }
            return XmlCabinet::add(x);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_get_XML_File(const char* file, int debug)
    {
        try {
            XML_Node* x = get_XML_File(file, debug);
            return XmlCabinet::add(x);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int ct_clearXML()
    {
        try {
            XmlCabinet::clear();
            close_XML_File("all");
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_del(int i)
    {
        try {
            XmlCabinet::del(i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_removeChild(int i, int j)
    {
        try {
            XmlCabinet::item(i).removeChild(&XmlCabinet::item(j));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_copy(int i)
    {
        try {
            return XmlCabinet::newCopy(i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_build(int i, const char* file)
    {
        try {
            writelog("WARNING: xml_build called. Use get_XML_File instead.");
            string path = findInputFile(file);
            XmlCabinet::item(i).build(path);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_attrib(int i, const char* key, size_t lenval, char* value)
    {
        try {
            XML_Node& node = XmlCabinet::item(i);
            if (node.hasAttrib(key)) {
                return static_cast<int>(copyString(node[key], value, lenval));
            } else {
                throw CanteraError("xml_attrib","node "
                                   " has no attribute '"+string(key)+"'");
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_addAttrib(int i, const char* key, const char* value)
    {
        try {
            XmlCabinet::item(i).addAttribute(key, value);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int xml_addComment(int i, const char* comment)
    {
        try {
            XmlCabinet::item(i).addComment(comment);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int xml_tag(int i, size_t lentag, char* tag)
    {
        try {
            return static_cast<int>(copyString(XmlCabinet::item(i).name(), tag, lentag));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_value(int i, size_t lenval, char* value)
    {
        try {
            return static_cast<int>(copyString(XmlCabinet::item(i).value(), value, lenval));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_child(int i, const char* loc)
    {
        try {
            XML_Node& c = XmlCabinet::item(i).child(string(loc));
            return XmlCabinet::add(&c);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int xml_child_bynumber(int i, int m)
    {
        try {
            XML_Node& c = XmlCabinet::item(i).child(m);
            return XmlCabinet::add(&c);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int xml_findID(int i, const char* id)
    {
        try {
            XML_Node* c = XmlCabinet::item(i).findID(id);
            if (c) {
                return XmlCabinet::add(c);
            } else {
                throw CanteraError("xml_findID", "id not found: '{}'", id);
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int xml_findByName(int i, const char* nm)
    {
        try {
            XML_Node* c = XmlCabinet::item(i).findByName(nm);
            if (c) {
                return XmlCabinet::add(c);
            } else {
                throw CanteraError("xml_findByName", "name '{}' not found", nm);
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int xml_nChildren(int i)
    {
        try {
            return (int) XmlCabinet::item(i).nChildren();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_addChild(int i, const char* name, const char* value)
    {
        try {
            XML_Node& c = XmlCabinet::item(i).addChild(name, value);
            return XmlCabinet::add(&c);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int xml_addChildNode(int i, int j)
    {
        try {
            XML_Node& c = XmlCabinet::item(i).addChild(XmlCabinet::item(j));
            return XmlCabinet::add(&c);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int xml_write(int i, const char* file)
    {
        try {
            ofstream f(file);
            if (f) {
                XmlCabinet::item(i).write(f);
            } else {
                throw CanteraError("xml_write",
                                   "file "+string(file)+" not found.");
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

}
