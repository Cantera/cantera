/**
 * @file ctxml.cpp
 */
/*
 *      $Id: ctxml.cpp,v 1.13 2009/07/22 01:21:35 hkmoffa Exp $
 */


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#define CANTERA_USE_INTERNAL
#include "ctxml.h"


// Cantera includes
#include "ctml.h"
#include "Cabinet.h"
#include "Storage.h"

#include <string.h>

using namespace std;
using namespace Cantera;


// Assign storage for the static member of the Templated Cabinet class
// class Cabinet<XML_Node>;
template<> Cabinet<XML_Node>*   Cabinet<XML_Node>::__storage = 0;

inline XML_Node* _xml(int i) {
    return Cabinet<XML_Node>::cabinet(false)->item(i);
}

extern "C" {  

    int DLL_EXPORT xml_new(const char* name = 0) {
        XML_Node* x;
        if (!name) 
            x = new XML_Node;
        else 
            x = new XML_Node(name);
        return Cabinet<XML_Node>::cabinet(true)->add(x);
    }

    int DLL_EXPORT xml_get_XML_File(const char* file, int debug) {
        try {
            XML_Node* x = get_XML_File(std::string(file), debug);
            return Cabinet<XML_Node>::cabinet(false)->add(x);
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT xml_clear() {
        try {
            Cabinet<XML_Node>::cabinet(false)->clear();
            close_XML_File("all");
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT xml_del(int i) {
        Cabinet<XML_Node>::cabinet(false)->del(i);
        return 0;
    }

    int DLL_EXPORT xml_removeChild(int i, int j) {
        _xml(i)->removeChild(_xml(j));
        return 0;
    }

    int DLL_EXPORT xml_copy(int i) {
        return Cabinet<XML_Node>::cabinet(false)->newCopy(i);
    }

    int DLL_EXPORT xml_assign(int i, int j) {
        return Cabinet<XML_Node>::cabinet(false)->assign(i,j);
    }

    int DLL_EXPORT xml_build(int i, const char* file) {
        try {
            writelog("WARNING: xml_build called. Use get_XML_File instead.");
            string path = findInputFile(string(file));
            ifstream f(path.c_str());
            if (!f) {
                throw CanteraError("xml_build",
                    "file "+string(file)+" not found.");
            }
            _xml(i)->build(f);
            f.close();
            return 0;
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT xml_preprocess_and_build(int i, const char* file, int debug) {
        try {
            get_CTML_Tree(_xml(i), string(file), debug);
            return 0;
        }
        catch (CanteraError) { return -1; }
    }



    int DLL_EXPORT xml_attrib(int i, const char* key, char* value) {
        try {
            string ky = string(key);
            XML_Node& node = *_xml(i);
            if (node.hasAttrib(ky)) {
                string v = node[ky];
                strncpy(value, v.c_str(), 80);
            }
            else 
                throw CanteraError("xml_attrib","node "
                    " has no attribute '"+ky+"'");
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT xml_addAttrib(int i, const char* key, const char* value) {
        try {
            string ky = string(key);
            string val = string(value);
            XML_Node& node = *_xml(i);
            node.addAttribute(ky, val);
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT xml_addComment(int i, const char* comment) {
        try {
            string c = string(comment);
            XML_Node& node = *_xml(i);
            node.addComment(c);
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT xml_tag(int i, char* tag) {
        try {
            XML_Node& node = *_xml(i);
            const string v = node.name();
            strncpy(tag, v.c_str(), 80);
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT xml_value(int i, char* value) {
        try {
            XML_Node& node = *_xml(i);
            const string v = node.value();
            strncpy(value, v.c_str(), 80);
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT xml_child(int i, const char* loc) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& c = node.child(string(loc));
            return Cabinet<XML_Node>::cabinet()->add(&c);
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT xml_child_bynumber(int i, int m) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& c = node.child(m);
            return Cabinet<XML_Node>::cabinet()->add(&c);
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT xml_findID(int i, const char* id) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node* c = node.findID(string(id));
            if (c) {
                return Cabinet<XML_Node>::cabinet()->add(c);
            }
            else 
                throw CanteraError("xml_find_id","id not found: "+string(id));
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT xml_findByName(int i, const char* nm) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node* c = node.findByName(string(nm));
            if (c) {
                return Cabinet<XML_Node>::cabinet()->add(c);
            }
            else 
                throw CanteraError("xml_findByName","name "+string(nm)
                    +" not found");
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT xml_nChildren(int i) {
        try {
            XML_Node& node = *_xml(i);
            return node.nChildren();
        }
        catch (CanteraError) { return -1; }
    }

    int DLL_EXPORT xml_addChild(int i, const char* name, const char* value) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& c = node.addChild(string(name),string(value));
            return Cabinet<XML_Node>::cabinet()->add(&c);
        }
        catch (CanteraError) { showErrors(cout); return -1; }
        return 0;
    }

    int DLL_EXPORT xml_addChildNode(int i, int j) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& chld = *_xml(j);
            XML_Node& c = node.addChild(chld);
            return Cabinet<XML_Node>::cabinet()->add(&c);
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT xml_write(int i, const char* file) {
        try {
            ofstream f(file);
            if (f) {
                XML_Node& node = *_xml(i);
                node.write(f);
            }
            else {
                throw CanteraError("xml_write",
                    "file "+string(file)+" not found.");
            }
            return 0;
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

    int DLL_EXPORT ctml_getFloatArray(int i, int n, doublereal* data, int iconvert) {
        try {
            XML_Node& node = *_xml(i);
            vector_fp v;
            bool conv = false;
            if (iconvert > 0) conv = true;
            getFloatArray(node, v, conv);
            int nv = v.size();

            // array not big enough
            if (n < nv) {
                throw CanteraError("ctml_getFloatArray",
                    "array must be dimensioned at least "+int2str(nv));
            }
            
            for (int i = 0; i < nv; i++) {
                data[i] = v[i];
            }
            n = nv;
        }
        catch (CanteraError) { return -1; }
        return 0;
    }

}
