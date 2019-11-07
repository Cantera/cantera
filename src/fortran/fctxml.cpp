/**
 * @file fctxml.cpp
 *
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "clib/clib_utils.h"
#include "cantera/base/ctml.h"

#include <cstring>
#include <fstream>

using namespace std;
using Cantera::XML_Node;
using Cantera::CanteraError;
using Cantera::handleAllExceptions;

#include "clib/Cabinet.h"

typedef Cabinet<XML_Node, false> XmlCabinet;
template<> XmlCabinet* XmlCabinet::s_storage = 0;

typedef integer status_t;

namespace {

XML_Node* _xml(const integer* i)
{
    return &XmlCabinet::item(*i);
}

} // unnamed namespace

std::string f2string(const char* s, ftnlen n);

extern "C" {

    integer fxml_new_(const char* name, ftnlen namelen)
    {
        try {
            XML_Node* x;
            if (!name) {
                x = new XML_Node;
            } else {
                x = new XML_Node(f2string(name, namelen), 0);
            }
            return XmlCabinet::add(x);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t fxml_get_xml_file_(const char* file, ftnlen filelen)
    {
        try {
            XML_Node* x = Cantera::get_XML_File(f2string(file, filelen));
            return XmlCabinet::add(x);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t fxml_clear_()
    {
        try {
            XmlCabinet::clear();
            Cantera::close_XML_File("all");
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t fxml_del_(const integer* i)
    {
        try {
            XmlCabinet::del(*i);
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t fxml_removechild_(const integer* i, const integer* j)
    {
        try {
            _xml(i)->removeChild(_xml(j));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t fxml_copy_(const integer* i)
    {
        try {
            return XmlCabinet::newCopy(*i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t fxml_attrib_(const integer* i, const char* key, char* value,
                          ftnlen keylen, ftnlen valuelen)
    {
        try {
            std::string ky = f2string(key, keylen);
            XML_Node& node = *_xml(i);
            if (node.hasAttrib(ky)) {
                std::string v = node[ky];
                strncpy(value, v.c_str(), valuelen);
            } else {
                throw CanteraError("fxml_attrib","node "
                                   " has no attribute '"+ky+"'");
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t fxml_addattrib_(const integer* i, const char* key,
                             const char* value, ftnlen keylen, ftnlen valuelen)
    {
        try {
            std::string ky = f2string(key, keylen);
            std::string val = f2string(value, valuelen);
            XML_Node& node = *_xml(i);
            node.addAttribute(ky, val);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t fxml_addcomment_(const integer* i, const char* comment,
                              ftnlen commentlen)
    {
        try {
            std::string c = f2string(comment, commentlen);
            XML_Node& node = *_xml(i);
            node.addComment(c);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t fxml_tag_(const integer* i, char* tag, ftnlen taglen)
    {
        try {
            XML_Node& node = *_xml(i);
            const std::string v = node.name();
            strncpy(tag, v.c_str(), taglen);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t fxml_value_(const integer* i, char* value, ftnlen valuelen)
    {
        try {
            XML_Node& node = *_xml(i);
            const std::string v = node.value();
            strncpy(value, v.c_str(), valuelen);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t fxml_child_(const integer* i, const char* loc, ftnlen loclen)
    {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& c = node.child(f2string(loc, loclen));
            return XmlCabinet::add(&c);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t fxml_child_bynumber_(const integer* i, const integer* m)
    {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& c = node.child(*m);
            return XmlCabinet::add(&c);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t fxml_findid_(const integer* i, const char* id, ftnlen idlen)
    {
        try {
            XML_Node& node = *_xml(i);
            XML_Node* c = node.findID(f2string(id, idlen));
            if (c) {
                return XmlCabinet::add(c);
            } else {
                throw CanteraError("fxml_findid", "id not found: '{}'",
                                   f2string(id, idlen));
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t fxml_findbyname_(const integer* i, const char* nm, ftnlen nmlen)
    {
        try {
            XML_Node& node = *_xml(i);
            XML_Node* c = node.findByName(f2string(nm, nmlen));
            if (c) {
                return XmlCabinet::add(c);
            } else {
                throw CanteraError("fxml_findByName", "name '{}' not found",
                                   f2string(nm, nmlen));
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    integer fxml_nchildren_(const integer* i)
    {
        try {
            XML_Node& node = *_xml(i);
            return node.nChildren();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t fxml_addchild_(const integer* i, const char* name,
                            const char* value, ftnlen namelen, ftnlen valuelen)
    {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& c = node.addChild(f2string(name, namelen),
                                        f2string(value,valuelen));
            return XmlCabinet::add(&c);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t fxml_addchildnode_(const integer* i, const integer* j)
    {
        try {
            XML_Node& node = *_xml(i);
            XML_Node& chld = *_xml(j);
            XML_Node& c = node.addChild(chld);
            return XmlCabinet::add(&c);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t fxml_write_(const integer* i, const char* file, ftnlen filelen)
    {
        try {
            std::string ff(file, filelen);
            ofstream f(ff.c_str());
            if (f) {
                XML_Node& node = *_xml(i);
                node.write(f);
            } else {
                throw CanteraError("fxml_write",
                                   "file "+f2string(file, filelen)+" not found.");
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t ctml_getfloatarray_(const integer* i, const integer* n,
                                 doublereal* data, const integer* iconvert)
    {
        try {
            XML_Node& node = *_xml(i);
            Cantera::vector_fp v;
            bool conv = false;
            if (*iconvert > 0) {
                conv = true;
            }
            getFloatArray(node, v, conv);
            int nv = v.size();

            // array not big enough
            if (*n < nv) {
                throw CanteraError("ctml_getfloatarray",
                                   "array must be dimensioned at least {}", nv);
            }

            for (int i = 0; i < nv; i++) {
                data[i] = v[i];
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

}
