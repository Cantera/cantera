/**
 * @file ctxml.cpp
 */
#define CANTERA_USE_INTERNAL
#include "ctxml.h"

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

    int xml_clear()
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

    int xml_assign(int i, int j)
    {
        try {
            return XmlCabinet::assign(i,j);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_build(int i, const char* file)
    {
        try {
            writelog("WARNING: xml_build called. Use get_XML_File instead.");
            string path = findInputFile(file);
            ifstream f(path.c_str());
            if (!f) {
                throw CanteraError("xml_build",
                                   "file "+string(file)+" not found.");
            }
            XmlCabinet::item(i).build(f);
            f.close();
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_preprocess_and_build(int i, const char* file, int debug)
    {
        try {
            XmlCabinet::item(i) = *get_XML_File(file);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    int xml_attrib(int i, const char* key, char* value)
    {
        try {
            XML_Node& node = XmlCabinet::item(i);
            if (node.hasAttrib(key)) {
                string v = node[key];
                strncpy(value, v.c_str(), 80);
            } else
                throw CanteraError("xml_attrib","node "
                                   " has no attribute '"+string(key)+"'");
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
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

    int xml_tag(int i, char* tag)
    {
        try {
            string v = XmlCabinet::item(i).name();
            strncpy(tag, v.c_str(), 80);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    int xml_value(int i, char* value)
    {
        try {
            string v = XmlCabinet::item(i).value();
            strncpy(value, v.c_str(), 80);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
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
                throw CanteraError("xml_find_id","id not found: "+string(id));
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
            } else
                throw CanteraError("xml_findByName","name "+string(nm)
                                   +" not found");
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

    int ctml_getFloatArray(int i, size_t n, doublereal* data, int iconvert)
    {
        try {
            XML_Node& node = XmlCabinet::item(i);
            vector_fp v;
            bool conv = false;
            if (iconvert > 0) {
                conv = true;
            }
            getFloatArray(node, v, conv);
            size_t nv = v.size();

            // array not big enough
            if (n < nv) {
                throw CanteraError("ctml_getFloatArray",
                                   "array must be dimensioned at least "+int2str(nv));
            }

            for (size_t i = 0; i < nv; i++) {
                data[i] = v[i];
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

}
