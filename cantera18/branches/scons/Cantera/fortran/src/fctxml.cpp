/**
 * @file fctxml.cpp
 *
 */

/*
 * $Revision: 1.9 $
 * $Date: 2009/07/23 17:03:04 $
 */

// Copyright 2001  California Institute of Technology


#include "flib_defs.h"

#include "ctml.h"

#include <cstring>

using namespace ctml;
using namespace std;

#include "../../clib/src/Cabinet.h"

// Assign storage for the templated classes static member
template<> Cabinet<XML_Node> * Cabinet<XML_Node>::__storage = 0;

inline XML_Node* _xml(const integer* i) {
    return Cabinet<XML_Node>::cabinet(false)->item(*i);
}

static void handleError() {
    error(lastErrorMessage());
}

std::string f2string(const char* s, ftnlen n);

extern "C" {  

    integer DLL_EXPORT fxml_new_(const char* name, ftnlen namelen) {
        XML_Node* x;
        if (!name) 
            x = new XML_Node;
        else 
            x = new XML_Node(f2string(name, namelen), 0);
        return Cabinet<XML_Node>::cabinet(true)->add(x);
    }

    status_t DLL_EXPORT fxml_get_xml_file_(const char* file, ftnlen filelen) {
        try {
            XML_Node* x = get_XML_File(f2string(file, filelen));
            int ix = Cabinet<XML_Node>::cabinet(false)->add(x);
            return ix;
        }
        catch (CanteraError) {
            handleError();
            return -1; 
        }
    }

    status_t DLL_EXPORT fxml_clear_() {
        try {
            Cabinet<XML_Node>::cabinet(false)->clear();
            close_XML_File("all");
            return 0;
        }
        catch (CanteraError) { handleError(); return -1;}
    }

    status_t DLL_EXPORT fxml_del_(const integer* i) {
        Cabinet<XML_Node>::cabinet(false)->del(*i);
        return 0;
    }

    status_t DLL_EXPORT fxml_removechild_(const integer* i, const integer* j) {
        _xml(i)->removeChild(_xml(j));
        return 0;
    }

    status_t DLL_EXPORT fxml_copy_(const integer* i) {
        return Cabinet<XML_Node>::cabinet(false)->newCopy(*i);
    }

    status_t DLL_EXPORT fxml_assign_(const integer* i, const integer* j) {
        return Cabinet<XML_Node>::cabinet(false)->assign(*i,*j);
    }

    status_t DLL_EXPORT fxml_attrib_(const integer* i, const char* key, 
        char* value, ftnlen keylen, ftnlen valuelen) {
        try {
            std::string ky = f2string(key, keylen);
            XML_Node& node = *_xml(i);
            if (node.hasAttrib(ky)) {
                std::string v = node[ky];
                strncpy(value, v.c_str(), valuelen);
            }
            else 
                throw CanteraError("fxml_attrib","node "
                    " has no attribute '"+ky+"'");
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_addattrib_(const integer* i, 
        const char* key, const char* value, ftnlen keylen, ftnlen valuelen) {
        try {
            std::string ky = f2string(key, keylen);
            std::string val = f2string(value, valuelen);
            XML_Node& node = *_xml(i);
            node.addAttribute(ky, val);
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_addcomment_(const integer* i, const char* comment, 
        ftnlen commentlen) {
        try {
            std::string c = f2string(comment, commentlen);
            XML_Node& node = *_xml(i);
            node.addComment(c);
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_tag_(const integer* i, char* tag, ftnlen taglen) {
        try {
            XML_Node& node = *_xml(i);
            const std::string v = node.name();
            strncpy(tag, v.c_str(), taglen);
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_value_(const integer* i, char* value, ftnlen valuelen) {
        try {
            XML_Node& node = *_xml(i);
            const std::string v = node.value();
            strncpy(value, v.c_str(), valuelen);
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_child_(const integer* i, const char* loc, ftnlen loclen) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& c = node.child(f2string(loc, loclen));
            return Cabinet<XML_Node>::cabinet()->add(&c);
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_child_bynumber_(const integer* i, const integer* m) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& c = node.child(*m);
            return Cabinet<XML_Node>::cabinet()->add(&c);
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_findid_(const integer* i, const char* id, ftnlen idlen) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node* c = node.findID(f2string(id, idlen));
            if (c) {
                return Cabinet<XML_Node>::cabinet()->add(c);
            }
            else 
                throw CanteraError("fxml_find_id","id not found: "+f2string(id, idlen));
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_findbyname_(const integer* i, const char* nm, ftnlen nmlen) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node* c = node.findByName(f2string(nm, nmlen));
            if (c) {
                return Cabinet<XML_Node>::cabinet()->add(c);
            }
            else 
                throw CanteraError("fxml_findByName","name "+f2string(nm, nmlen)
                    +" not found");
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    integer DLL_EXPORT fxml_nchildren_(const integer* i) {
        try {
            XML_Node& node = *_xml(i);
            return node.nChildren();
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_addchild_(const integer* i, const char* name, 
        const char* value, ftnlen namelen, ftnlen valuelen) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& c = node.addChild(f2string(name, namelen),
                f2string(value,valuelen));
            return Cabinet<XML_Node>::cabinet()->add(&c);
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_addchildnode_(const integer* i, const integer* j) {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& chld = *_xml(j);
            XML_Node& c = node.addChild(chld);
            return Cabinet<XML_Node>::cabinet()->add(&c);
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT fxml_write_(const integer* i, const char* file, ftnlen filelen) {
        try {
            std::string ff(file, filelen);
            ofstream f(ff.c_str());
            if (f) {
                XML_Node& node = *_xml(i);
                node.write(f);
            }
            else {
                throw CanteraError("fxml_write",
                    "file "+f2string(file, filelen)+" not found.");
            }
            return 0;
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

    status_t DLL_EXPORT ctml_getfloatarray_(const integer* i, const integer* n, 
        doublereal* data, const integer* iconvert) {
        try {
            XML_Node& node = *_xml(i);
            vector_fp v;
            bool conv = false;
            if (*iconvert > 0) conv = true;
            getFloatArray(node, v, conv);
            int nv = v.size();

            // array not big enough
            if (*n < nv) {
                throw CanteraError("ctml_getfloatarray",
                    "array must be dimensioned at least "+int2str(nv));
            }
            
            for (int i = 0; i < nv; i++) {
                data[i] = v[i];
            }
            //n = nv;
        }
        catch (CanteraError) { handleError(); }
        return 0;
    }

}
