/**
 * @file ctml.cpp
 *
 * Functions to read and write CTML.
 *
 */

/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2002  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ctml.h"

#define CTML_VERSION_1_4_1

namespace ctml {

    static doublereal fpValue(string val) {
        return atof(stripws(val).c_str());
    }

    void addBool(XML_Node& node, string title, bool val) {
        string v = (val ? "true" : "false");
        XML_Node& f = node.addChild("bool",v);
        f.addAttribute("title",title);
    }

    void addInteger(XML_Node& node, string title, int val, 
        string units, string type) {
        XML_Node& f = node.addChild("integer",val);
        f.addAttribute("title",title);
        if (type != "") f.addAttribute("type",type);
        if (units != "") f.addAttribute("units",units);
    } 

    void addIntegerArray(XML_Node& node, string title, int n, 
        const int* vals, string units, string type,
        doublereal minval, doublereal maxval) {
        string fmt = "%8d";
        int i;
        string v = "";
        for (i = 0; i < n; i++) {
            v += int2str(vals[i],fmt);
            if (i == n-1) v += "\n"; 
            else if (i > 0 && (i+1) % 3 == 0) v += ",\n";
            else v += ", ";
        }
        XML_Node& f = node.addChild("intArray",v);
        f.addAttribute("title",title);
        if (type != "") f.addAttribute("type",type);        
        f.addAttribute("size",n);
        if (units != "") f.addAttribute("units",units);
        if (minval != Undef) f.addAttribute("min",minval);
        if (maxval != Undef) f.addAttribute("max",maxval);
    }

    void addFloat(XML_Node& node, 
        string title, 
        doublereal val, 
        string units, 
        string type, 
        doublereal minval,
        doublereal maxval) {
        string fmt = "%17.9E";
#ifdef CTML_VERSION_1_4
        XML_Node& f = node.addChild("float",val,fmt);
        f.addAttribute("title",title);
#endif
#ifdef CTML_VERSION_1_4_1
        XML_Node& f = node.addChild(title,val,fmt); 
#endif
        if (type != "") f.addAttribute("type",type);
        if (units != "") f.addAttribute("units",units);
        if (minval != Undef) f.addAttribute("min",minval);
        if (maxval != Undef) f.addAttribute("max",maxval);
    }


    void addFloatArray(XML_Node& node, string title, int n, 
        const double* vals, string units, string type,
        doublereal minval, doublereal maxval) {
        string fmt = "%17.9E";
        int i;
        string v = "";
        for (i = 0; i < n; i++) {
            v += fp2str(vals[i],fmt);
            if (i == n-1) v += "\n"; 
            else if (i > 0 && (i+1) % 3 == 0) v += ",\n";
            else v += ", ";
        }
        XML_Node& f = node.addChild("floatArray",v);
        f.addAttribute("title",title);
        if (type != "") f.addAttribute("type",type);        
        f.addAttribute("size",n);
        if (units != "") f.addAttribute("units",units);
        if (minval != Undef) f.addAttribute("min",minval);
        if (maxval != Undef) f.addAttribute("max",maxval);
    }

    void addString(XML_Node& node, string title, string val, 
        string type) {
        XML_Node& f = node.addChild("string",val);
        f.addAttribute("title",title);
        if (type != "") f.addAttribute("type",type);        
    }

    XML_Node* getByTitle(XML_Node& node, string title) {
        XML_Node* s = node.findByAttr("title",title);
        if (s->parent() == &node) return s;
        else return 0;
    }
        
    string getString(XML_Node& parent, string name) {
        if (!parent.hasChild(name)) return "";
        return parent(name);
    }

    void getString(XML_Node& node, string title, string& val, 
        string& type) {
        val = "";
        type = "";
        XML_Node* s = getByTitle(node, title);
        if (s) 
            if (s->name() == "string") {
                val = (*s).value();
                type = (*s)["type"];
                return;
            }
    }

    void getIntegers(XML_Node& node, map<string,int>& v) {
        vector<XML_Node*> f;
        node.getChildren("integer",f);
        int n = f.size();
        integer x, x0, x1;
        string typ, title, vmin, vmax;
        for (int i = 0; i < n; i++) {
            XML_Node& fi = *(f[i]);
            x = atoi(fi().c_str());
            title = fi["title"];
            vmin = fi["min"];
            vmax = fi["max"];
            if (vmin != "")
                x0 = atoi(vmin.c_str());
            if (fi["max"] != "")
                x1 = atoi(vmax.c_str());
            v[title] = x;
        }
    }


    void getStrings(XML_Node& node, map<string,string>& v) {
        vector<XML_Node*> f;
        node.getChildren("string",f);
        int n = f.size();
        string typ, title;
        for (int i = 0; i < n; i++) {
            XML_Node& fi = *(f[i]);
            title = fi["title"];
            v[title] = fi();
        }
    }


    void getFloats(XML_Node& node, map<string,double>& v, bool convert) {
        vector<XML_Node*> f;
        node.getChildren("float",f);
        int n = f.size();
        doublereal x, x0, x1, fctr;
        string typ, title, units, vmin, vmax;
        for (int i = 0; i < n; i++) {
            XML_Node& fi = *(f[i]);
            x = atof(fi().c_str());
            x0 = Undef;
            x1 = Undef;
            typ = fi["type"];
            title = fi["title"];
            units = fi["units"];
            vmin = fi["min"];
            vmax = fi["max"];
            if (vmin != "") {
                x0 = atof(vmin.c_str());
                if (x < x0 - Tiny) {
                write("\nWarning: value "+fi()+" is below lower limit of "
                    +vmin+".\n");
                }
            }
            if (fi["max"] != "") {
                x1 = atof(vmax.c_str());
                if (x > x1 + Tiny) {
                write("\nWarning: value "+fi()+" is above upper limit of "
                    +vmax+".\n");
                }
            }
            fctr = (convert ? toSI(units) : 1.0); // toSI(typ,units);
            v[title] = fctr*x;
        }
    }


    /**
     * Get a floating-point value from a child element.  Returns a
     * double value for the child named 'name' of element 'parent'. If
     * 'type' is supplied and matches a known unit type, unit
     * conversion to SI will be done if the child element has an attribute
     * 'units'.
     */
    doublereal getFloat(XML_Node& parent, string name, string type) {
        if (!parent.hasChild(name)) 
            throw CanteraError("getFloat (called from XML Node \"" +
			       parent.name() + "\"): ",
			       "no child XML element named " + name);
        XML_Node& node = parent.child(name);
        doublereal x, x0, x1, fctr = 1.0;
        string units, vmin, vmax;
        x = atof(node().c_str());
        x0 = Undef;
        x1 = Undef;
        units = node["units"];
        vmin = node["min"];
        vmax = node["max"];
        if (vmin != "") {
            x0 = atof(vmin.c_str());
            if (x < x0 - Tiny) {
                write("\nWarning: value "+node()+" is below lower limit of "
                    +vmin+".\n");
                }
            }
        if (node["max"] != "") {
            x1 = atof(vmax.c_str());
            if (x > x1 + Tiny) {
                write("\nWarning: value "+node()+" is above upper limit of "
                    +vmax+".\n");
            }
        }
        if (type == "actEnergy" && units != "") 
            fctr = actEnergyToSI(units);
        else if (type != "" && units != "")
            fctr = toSI(units);
        return fctr*x;
    }


    void getFloatArray(XML_Node& node, vector_fp& v, bool convert) {
        int icom;
        string numstr;
        if (node.name() != "floatArray") 
            throw CanteraError("getFloatArray","wrong element type: "
                +node.name());

        v.clear();
        doublereal vmin = Undef, vmax = Undef;

        doublereal funit = 1.0;
        if (node["units"] != "" && convert) {
            funit = toSI(node["units"]);
        }

        if (node["min"] != "") 
            vmin = atof(node["min"].c_str());
        if (node["max"] != "") 
            vmax = atof(node["max"].c_str());

        doublereal vv;
        string val = node.value();
        while (1 > 0) {
            icom = val.find(',');
            if (icom >= 0) {
                numstr = val.substr(0,icom);
                val = val.substr(icom+1,val.size());
                v.push_back(atof(numstr.c_str()));
            }
            else {
                v.push_back(atof(val.c_str()));
                break;
            }
            vv = v.back();
            if (vmin != Undef && vv < vmin - Tiny) {
                write("\nWarning: value "+fp2str(vv)+
                    " is below lower limit of " +fp2str(vmin)+".\n");
            }
            if (vmax != Undef && vv > vmax + Tiny) {
                write("\nWarning: value "+fp2str(vv)+
                    " is above upper limit of " +fp2str(vmin)+".\n");
            }
        }
        int nv = v.size();
        for (int n = 0; n < nv; n++) {
            v[n] *= funit;
        }
    }

    void getMap(XML_Node& node, map<string, string>& m) {
        vector<string> v;
        getStringArray(node, v);
        string key, val;
        int n = v.size();
        int icolon;
        for (int i = 0; i < n; i++) {
            icolon = v[i].find(":");
            if (icolon < 0) {
                throw CanteraError("getMap","missing colon in map entry ("
                    +v[i]+")");
            }
            key = v[i].substr(0,icolon);
            val = v[i].substr(icolon+1, v.size());
            m[key] = val;
        }
    }

    void getPairs(XML_Node& node, vector<string>& key, vector<string>& val) {
        vector<string> v;
        getStringArray(node, v);
        int n = v.size();
        int icolon;
        for (int i = 0; i < n; i++) {
            icolon = v[i].find(":");
            if (icolon < 0) {
                throw CanteraError("getMap","missing colon in map entry ("
                    +v[i]+")");
            }
            key.push_back(v[i].substr(0,icolon));
            val.push_back(v[i].substr(icolon+1, v.size()));
        }
    }

    void getStringArray(XML_Node& node, vector<string>& v) {
        int ibegin, iend;

        v.clear();

        string val = node.value();
        while (1 > 0) {
            ibegin = val.find_first_not_of(" \n\t");
            if (ibegin >= 0) {
                val = val.substr(ibegin,val.size());
                iend = val.find_first_of(" \n\t");
                if (iend > 0) {
                    v.push_back(val.substr(0,iend));
                    val = val.substr(iend+1,val.size());
                }
                else {
                    v.push_back(val.substr(0,iend));
                    break;
                }
            }
            else
                break;
        }
    }

    void getFunction(XML_Node& node, string& type, doublereal& xmin,
        doublereal& xmax, vector_fp& coeffs) {
        XML_Node& c = node.child("floatArray");
        coeffs.clear();
        getFloatArray(c,coeffs);
        xmin = Undef;
        if (node["min"] != "") xmin = fpValue(node["min"]);
        if (node["max"] != "") xmax = fpValue(node["max"]);
        type = node["type"];
    }
}
