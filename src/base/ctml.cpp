/**
 * @file ctml.cpp
 * Definitions for functions to read and write CTML.
 */
// Copyright 2002  California Institute of Technology

#include "cantera/base/ctml.h"

//@{
#define CTML_VERSION_1_4_1
//@}

#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

using namespace std;
using namespace Cantera;

namespace ctml
{
std::string FP_Format = "%23.15E";
std::string INT_Format = "%8d";

void addInteger(Cantera::XML_Node& node, const std::string& title, const int val,
                const std::string& units, const std::string& type)
{
#ifdef CTML_VERSION_1_4
    XML_Node& f = node.addChild("integer", val);
    f.addAttribute("title", title);
#else
    XML_Node& f = node.addChild(title, val);
#endif
    f.addAttribute("vtype", "integer");
    if (type != "") {
        f.addAttribute("type",type);
    }
    if (units != "") {
        f.addAttribute("units",units);
    }
}

void addFloat(Cantera::XML_Node& node, const std::string& title,
              const doublereal val, const std::string& units,
              const std::string& type, const doublereal minval,
              const doublereal maxval)
{
#ifdef CTML_VERSION_1_4
    XML_Node& f = node.addChild("float", val, ctml::FP_Format);
    f.addAttribute("title", title);
#else
    XML_Node& f = node.addChild(title, val, ctml::FP_Format);
#endif
    if (type != "") {
        f.addAttribute("type",type);
    }
    if (units != "") {
        f.addAttribute("units",units);
    }
    f.addAttribute("vtype", "float");
    if (minval != Undef) {
        f.addAttribute("min",minval);
    }
    if (maxval != Undef) {
        f.addAttribute("max",maxval);
    }
}

void addFloatArray(Cantera::XML_Node& node, const std::string& title, const size_t n,
                   const doublereal* const vals, const std::string& units,
                   const std::string& type,
                   const doublereal minval, const doublereal maxval)
{
    size_t i;
    std::string v = "";
    for (i = 0; i < n; i++) {
        v += fp2str(vals[i],FP_Format);
        if (i == n-1) {
            v += "\n";
        } else if (i > 0 && (i+1) % 3 == 0) {
            v += ",\n";
        } else {
            v += ", ";
        }
    }
    XML_Node& f = node.addChild("floatArray",v);
    f.addAttribute("title",title);
    if (type != "") {
        f.addAttribute("type",type);
    }
    f.addAttribute("size", double(n));
    if (units != "") {
        f.addAttribute("units",units);
    }
    if (minval != Undef) {
        f.addAttribute("min",minval);
    }
    if (maxval != Undef) {
        f.addAttribute("max",maxval);
    }
}

void addNamedFloatArray(Cantera::XML_Node& node, const std::string& name, const int n,
                        const doublereal* const vals, const std::string units,
                        const std::string type, const doublereal minval,
                        const doublereal maxval)
{
    int i;
    std::string v = "";
    for (i = 0; i < n; i++) {
        v += fp2str(vals[i],FP_Format);
        if (i == n-1) {
            v += "\n";
        } else if (i > 0 && (i+1) % 3 == 0) {
            v += ",\n";
        } else {
            v += ", ";
        }
    }
    XML_Node& f = node.addChild(name, v);
    if (type != "") {
        f.addAttribute("type",type);
    }
    /*
     *  Add vtype, which indicates the type of the value. Here we specify it as a list of floats separated
     *  by commas, with a length given by size attribute.
     */
    f.addAttribute("vtype", "floatArray");

    f.addAttribute("size", n);
    if (units != "") {
        f.addAttribute("units", units);
    }
    if (minval != Undef) {
        f.addAttribute("min", minval);
    }
    if (maxval != Undef) {
        f.addAttribute("max", maxval);
    }
}

void addString(Cantera::XML_Node& node, const std::string& titleString,
               const std::string& valueString,
               const std::string& typeString)
{
    XML_Node& f = node.addChild("string", valueString);
    f.addAttribute("title", titleString);
    if (typeString != "") {
        f.addAttribute("type", typeString);
    }
}

XML_Node* getByTitle(const Cantera::XML_Node& node, const std::string& title)
{
    XML_Node* s = node.findByAttr("title", title);
    if (!s) {
        return 0;
    }
    if (s->parent() == &node) {
        return s;
    }
    return 0;
}

std::string getChildValue(const Cantera::XML_Node& parent, const std::string& nameString)
{
    if (!parent.hasChild(nameString)) {
        return "";
    }
    return parent(nameString);
}

void getString(const Cantera::XML_Node& node, const std::string& titleString, std::string& valueString,
               std::string& typeString)
{
    valueString = "";
    typeString = "";
    XML_Node* s = getByTitle(node, titleString);
    if (s)
        if (s->name() == "string") {
            valueString = (*s).value();
            typeString = (*s)["type"];
            return;
        }
}

void getNamedStringValue(const Cantera::XML_Node& node, const std::string& nameString, std::string& valueString,
                         std::string& typeString)
{
    valueString = "";
    typeString = "";
    if (node.hasChild(nameString)) {
        XML_Node& xc = node.child(nameString);
        valueString = xc.value();
        typeString = xc["type"];
    } else {
        XML_Node* s = getByTitle(node, nameString);
        if (s) {
            if (s->name() == "string") {
                valueString = (*s).value();
                typeString = (*s)["type"];
                return;
            }
        }
    }
}

void getIntegers(const Cantera::XML_Node& node,
                 std::map<std::string, int>& v)
{
    std::vector<XML_Node*> f;
    node.getChildren("integer",f);
    int n = static_cast<int>(f.size());
    integer x;
    std::string typ, title, vmin, vmax;
    for (int i = 0; i < n; i++) {
        const XML_Node& fi = *(f[i]);
        x = atoi(fi().c_str());
        title = fi["title"];
        vmin = fi["min"];
        vmax = fi["max"];
        if (vmin != "")
            if (fi["max"] != "") {
                v[title] = x;
            }
    }
}

void getFloats(const Cantera::XML_Node& node, std::map<std::string, double>& v,
               const bool convert)
{
    warn_deprecated("ctml::getFloats",
                    "To be removed in Cantera 2.2.");
    std::vector<XML_Node*> f;
    node.getChildren("float",f);
    int n = static_cast<int>(f.size());
    doublereal x, x0, x1, fctr;
    std::string typ, title, units, vmin, vmax;
    for (int i = 0; i < n; i++) {
        const XML_Node& fi = *(f[i]);
        x = fpValue(fi());
        x0 = Undef;
        x1 = Undef;
        typ = fi["type"];
        title = fi["title"];
        units = fi["units"];
        vmin = fi["min"];
        vmax = fi["max"];
        if (vmin != "") {
            x0 = fpValue(vmin);
            if (x < x0 - Tiny) {
                writelog("\nWarning: value "+fi()+" is below lower limit of "
                         +vmin+".\n");
            }
        }
        if (fi["max"] != "") {
            x1 = fpValue(vmax);
            if (x > x1 + Tiny) {
                writelog("\nWarning: value "+fi()+" is above upper limit of "
                         +vmax+".\n");
            }
        }
        fctr = (convert ? toSI(units) : 1.0); // toSI(typ,units);
        v[title] = fctr*x;
    }
}

doublereal getFloat(const Cantera::XML_Node& parent,
                    const std::string& name,
                    const std::string& type)
{
    if (!parent.hasChild(name))
        throw CanteraError("getFloat (called from XML Node \"" +
                           parent.name() + "\"): ",
                           "no child XML element named \"" + name + "\" exists");
    const XML_Node& node = parent.child(name);
    return getFloatCurrent(node, type);
}

doublereal getFloatCurrent(const Cantera::XML_Node& node,
                           const std::string& type)
{
    doublereal x, x0, x1, fctr = 1.0;
    string units, vmin, vmax;
    x = fpValue(node());
    x0 = Undef;
    x1 = Undef;
    units = node["units"];
    vmin = node["min"];
    vmax = node["max"];
    if (vmin != "") {
        x0 = fpValue(vmin);
        if (x < x0 - Tiny) {
            writelog("\nWarning: value "+node()+" is below lower limit of "
                     +vmin+".\n");
        }
    }
    if (node["max"] != "") {
        x1 = fpValue(vmax);
        if (x > x1 + Tiny) {
            writelog("\nWarning: value "+node()+" is above upper limit of "
                     +vmax+".\n");
        }
    }
    // Note, most types of converters default to toSI() type atm.
    // This may change and become more specific in the future.
    if (type == "actEnergy" && units != "") {
        fctr = actEnergyToSI(units);
    } else if (type == "toSI" && units != "") {
        fctr = toSI(units);
    } else if (type == "temperature" && units != "") {
        fctr = toSI(units);
    } else if (type == "density" && units != "") {
        fctr = toSI(units);
    } else if (type == "pressure" && units != "") {
        fctr = toSI(units);
    } else if (type != "" && units != "") {
        fctr = toSI(units);
#ifdef DEBUG_MODE
        writelog("\nWarning: conversion toSI() was done on node value "  + node.name() +
                 "but wasn't explicitly requested. Type was \"" + type + "\"\n");
#endif
    }
    // Note, below currently produces a lot of output due to transport blocks.
    // This needs to be addressed.
#ifdef DEBUG_MODE_MORE
    else if (type == "" && units != "") {
        writelog("\nWarning: XML node "  + node.name() +
                 "has a units attribute, \"" + units + "\","
                 "but no conversion was done because the getFloat() command didn't have a type\n");
    }
#endif
    return fctr*x;
}

bool getOptionalFloat(const Cantera::XML_Node& parent,
                      const std::string& name,
                      doublereal& fltRtn,
                      const std::string& type)
{
    if (parent.hasChild(name)) {
        fltRtn= getFloat(parent, name, type);
        return true;
    }
    return false;
}

doublereal getFloatDefaultUnits(const Cantera::XML_Node& parent,
                                const std::string& name,
                                const std::string& defaultUnits,
                                const std::string& type)
{
    doublereal fctr = 1.0;
    if (defaultUnits == "") {
        throw CanteraError("getFloatDefaultUnits",
                           "need to supply an actual value of defaultUnits");
    }
    if (type == "actEnergy") {
        fctr = actEnergyToSI(defaultUnits);
    } else if (type == "toSI") {
        fctr = toSI(defaultUnits);
    } else if (defaultUnits == "temperature") {
        fctr = toSI(defaultUnits);
    } else if (type == "density") {
        fctr = toSI(defaultUnits);
    } else if (type == "pressure") {
        fctr = toSI(defaultUnits);
    } else {
        throw CanteraError("getFloatDefaultUnits",
                           "type of units must be supplied and understood");
    }
    doublereal val = getFloat(parent, name, type);
    val /= fctr;
    return val;
}

bool getOptionalModel(const Cantera::XML_Node& parent, const std::string& nodeName,
                      std::string& modelName)
{
    if (parent.hasChild(nodeName)) {
        const XML_Node& node = parent.child(nodeName);
        modelName = node["model"];
        return true;
    }
    return false;
}

int getInteger(const Cantera::XML_Node& parent, const std::string& name)
{
    if (!parent.hasChild(name)) {
        throw CanteraError("getInteger (called from XML Node \"" +
                           parent.name() + "\"): ",
                           "no child XML element named " + name);
    }
    const XML_Node& node = parent.child(name);
    int x, x0, x1;
    string units, vmin, vmax;
    x = atoi(node().c_str());
    x0 = -9999999;
    x1 =  9999999;
    vmin = node["min"];
    vmax = node["max"];
    if (vmin != "") {
        x0 = atoi(vmin.c_str());
        if (x < x0) {
            writelog("\nWarning: value "+node()+" is below lower limit of "
                     +vmin+".\n");
        }
    }
    if (node["max"] != "") {
        x1 = atoi(vmax.c_str());
        if (x > x1) {
            writelog("\nWarning: value "+node()+" is above upper limit of "
                     +vmax+".\n");
        }
    }
    return x;
}

size_t getFloatArray(const Cantera::XML_Node& node, std::vector<doublereal> & v,
                     const bool convert, const std::string& unitsString,
                     const std::string& nodeName)
{
    std::string::size_type icom;
    string numstr;
    doublereal dtmp;
    string nn = node.name();
    const Cantera::XML_Node* readNode = &node;
    if (nn != nodeName) {
        vector<Cantera::XML_Node*> ll;
        node.getChildren(nodeName, ll);
        if (ll.size() == 0) {
            throw CanteraError("getFloatArray",
                               "wrong xml element type/name: was expecting "
                               + nodeName + "but accessed " + node.name());
        } else {
            readNode = ll[0];
            ll.clear();
            readNode->getChildren("floatArray", ll);
            if (ll.size() > 0) {
                readNode = ll[0];
            }
        }
    }

    v.clear();
    doublereal vmin = Undef, vmax = Undef;

    doublereal funit = 1.0;
    /*
     * Get the attributes field, units, from the XML node
     */
    std::string units = (*readNode)["units"];
    if (units != "" && convert) {
        if (unitsString == "actEnergy" && units != "") {
            funit = actEnergyToSI(units);
        } else if (unitsString != "" && units != "") {
            funit = toSI(units);
        }
    }

    if ((*readNode)["min"] != "") {
        vmin = fpValueCheck((*readNode)["min"]);
    }
    if ((*readNode)["max"] != "") {
        vmax = fpValueCheck((*readNode)["max"]);
    }

    doublereal vv;
    std::string val = readNode->value();
    while (1 > 0) {
        icom = val.find(',');
        if (icom != string::npos) {
            numstr = val.substr(0,icom);
            val = val.substr(icom+1,val.size());
            dtmp = fpValueCheck(numstr);
            v.push_back(dtmp);
        } else {
            /*
             * This little bit of code is to allow for the
             * possibility of a comma being the last
             * item in the value text. This was allowed in
             * previous versions of Cantera, even though it
             * would appear to be odd. So, we keep the
             * possibility in for backwards compatibility.
             */
            if (!val.empty()) {
                dtmp = fpValueCheck(val);
                v.push_back(dtmp);
            }
            break;
        }
        vv = v.back();
        if (vmin != Undef && vv < vmin - Tiny) {
            writelog("\nWarning: value "+fp2str(vv)+
                     " is below lower limit of " +fp2str(vmin)+".\n");
        }
        if (vmax != Undef && vv > vmax + Tiny) {
            writelog("\nWarning: value "+fp2str(vv)+
                     " is above upper limit of " +fp2str(vmin)+".\n");
        }
    }
    for (size_t n = 0; n < v.size(); n++) {
        v[n] *= funit;
    }
    return v.size();
}

void getMap(const Cantera::XML_Node& node, std::map<std::string, std::string>& m)
{
    std::vector<std::string> v;
    getStringArray(node, v);
    std::string key, val;
    int n = static_cast<int>(v.size());
    string::size_type icolon;
    for (int i = 0; i < n; i++) {
        icolon = v[i].find(":");
        if (icolon == string::npos) {
            throw CanteraError("getMap","missing colon in map entry ("
                               +v[i]+")");
        }
        key = v[i].substr(0,icolon);
        val = v[i].substr(icolon+1, v[i].size());
        m[key] = val;
    }
}

int getPairs(const Cantera::XML_Node& node, std::vector<std::string>& key,
             std::vector<std::string>& val)
{
    vector<string> v;
    getStringArray(node, v);
    int n = static_cast<int>(v.size());
    string::size_type icolon;
    for (int i = 0; i < n; i++) {
        icolon = v[i].find(":");
        if (icolon == string::npos) {
            throw CanteraError("getPairs","Missing a colon in the Pair entry ("
                               +v[i]+")");
        }
        key.push_back(v[i].substr(0,icolon));
        val.push_back(v[i].substr(icolon+1, v[i].size()));
        //cout << "getPairs: " << key.back() << " " << val.back() << endl;
    }
    return n;
}

void getMatrixValues(const Cantera::XML_Node& node,
                     const std::vector<std::string>& keyStringRow,
                     const std::vector<std::string>& keyStringCol,
                     Cantera::Array2D& retnValues, const bool convert,
                     const bool matrixSymmetric)
{
    size_t szKey1 = keyStringRow.size();
    size_t szKey2 = keyStringCol.size();
    size_t nrow   = retnValues.nRows();
    size_t ncol   = retnValues.nColumns();
    if (szKey1 > nrow) {
        throw CanteraError("getMatrixValues",
                           "size of key1 greater than numrows");
    }
    if (szKey2 > ncol) {
        throw CanteraError("getMatrixValues",
                           "size of key2 greater than num cols");
    }
    if (matrixSymmetric) {
        if (nrow != ncol) {
            throw CanteraError("getMatrixValues",
                               "nrow != ncol for a symmetric matrix");
        }
    }

    /*
     * Get the attributes field, units, from the XML node
     * and determine the conversion factor, funit.
     */
    doublereal funit = 1.0;
    string units = node["units"];
    if (units != "" && convert) {
        funit = toSI(units);
    }

    string key1;
    string key2;
    string rmm;
    string val;
    vector<string> v;
    getStringArray(node, v);
    string::size_type icolon;
    for (size_t i = 0; i < v.size(); i++) {
        icolon = v[i].find(":");
        if (icolon == string::npos) {
            throw CanteraError("getMatrixValues","Missing two colons ("
                               +v[i]+")");
        }
        key1 = v[i].substr(0,icolon);
        rmm = v[i].substr(icolon+1, v[i].size());
        icolon = rmm.find(":");
        if (icolon == string::npos) {
            throw CanteraError("getMatrixValues","Missing one colon ("
                               +v[i]+")");
        }
        key2 = rmm.substr(0,icolon);
        val = rmm.substr(icolon+1, rmm.size());

        size_t icol = npos;
        size_t irow = npos;
        for (size_t j = 0; j < szKey1; j++) {
            if (key1 == keyStringRow[j]) {
                irow = j;
                break;
            }
        }
        if (irow == npos) {
            throw CanteraError("getMatrixValues","Row not matched by string: "
                               + key1);
        }
        for (size_t j = 0; j < szKey2; j++) {
            if (key2 == keyStringCol[j]) {
                icol = j;
                break;
            }
        }
        if (icol == npos) {
            throw CanteraError("getMatrixValues","Col not matched by string: "
                               + key2);
        }
        double dval = fpValueCheck(val);
        dval *= funit;
        /*
         * Finally, insert the value;
         */
        retnValues(irow, icol) = dval;
        if (matrixSymmetric) {
            retnValues(icol, irow) = dval;
        }
    }
}

void getStringArray(const Cantera::XML_Node& node, std::vector<std::string>& v)
{
    std::string val = node.value();
    tokenizeString(val, v);
}

}
