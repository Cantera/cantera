/**
 * @file ctml.cpp
 * Definitions for functions to read and write CTML.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/Array.h"

using namespace std;

namespace Cantera
{
std::string FP_Format = "%23.15E";

void addFloat(XML_Node& node, const std::string& title,
              const doublereal val, const std::string& units,
              const std::string& type, const doublereal minval,
              const doublereal maxval)
{
    XML_Node& f = node.addChild(title, val, FP_Format);
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

void addFloatArray(XML_Node& node, const std::string& title, const size_t n,
                   const doublereal* const vals, const std::string& units,
                   const std::string& type,
                   const doublereal minval, const doublereal maxval)
{
    std::string v = "";
    for (size_t i = 0; i < n; i++) {
        v += fmt::sprintf(FP_Format, vals[i]);
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

void addNamedFloatArray(XML_Node& node, const std::string& name, const size_t n,
                        const doublereal* const vals, const std::string units,
                        const std::string type, const doublereal minval,
                        const doublereal maxval)
{
    std::string v = "";
    for (size_t i = 0; i < n; i++) {
        v += fmt::sprintf(FP_Format, vals[i]);
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

    // Add vtype, which indicates the type of the value. Here we specify it as a
    // list of floats separated by commas, with a length given by size
    // attribute.
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

void addString(XML_Node& node, const std::string& titleString,
               const std::string& valueString,
               const std::string& typeString)
{
    XML_Node& f = node.addChild("string", valueString);
    f.addAttribute("title", titleString);
    if (typeString != "") {
        f.addAttribute("type", typeString);
    }
}

XML_Node* getByTitle(const XML_Node& node, const std::string& title)
{
    XML_Node* s = node.findByAttr("title", title);
    if (s && s->parent() == &node) {
        return s;
    }
    return 0;
}

std::string getChildValue(const XML_Node& parent, const std::string& nameString)
{
    if (!parent.hasChild(nameString)) {
        return "";
    }
    return parent(nameString);
}

void getString(const XML_Node& node, const std::string& titleString, std::string& valueString,
               std::string& typeString)
{
    XML_Node* s = getByTitle(node, titleString);
    if (s && s->name() == "string") {
        valueString = s->value();
        typeString = s->attrib("type");
    } else {
        valueString = "";
        typeString = "";
    }
}

void getIntegers(const XML_Node& node,
                 std::map<std::string, int>& v)
{
    std::vector<XML_Node*> f = node.getChildren("integer");
    for (size_t i = 0; i < f.size(); i++) {
        const XML_Node& fi = *f[i];
        if (fi["min"] != "" && fi["max"] != "") {
            v[fi["title"]] = fi.int_value();
        }
    }
}

doublereal getFloat(const XML_Node& parent,
                    const std::string& name,
                    const std::string& type)
{
    if (!parent.hasChild(name)) {
        throw CanteraError("getFloat (called from XML Node \"" +
                           parent.name() + "\"): ",
                           "no child XML element named \"" + name + "\" exists");
    }
    const XML_Node& node = parent.child(name);
    return getFloatCurrent(node, type);
}

doublereal getFloatCurrent(const XML_Node& node, const std::string& type)
{
    doublereal fctr = 1.0;
    doublereal x = node.fp_value();
    const string& units = node["units"];
    const string& vmin = node["min"];
    const string& vmax = node["max"];
    if (vmin != "" && x < fpValue(vmin) - Tiny) {
        writelog("\nWarning: value "+node.value()+" is below lower limit of "
                 +vmin+".\n");
    }
    if (node["max"] != "" && x > fpValue(vmax) + Tiny) {
        writelog("\nWarning: value "+node.value()+" is above upper limit of "
                 +vmax+".\n");
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
        writelog("\nWarning: conversion toSI() was done on node value " + node.name() +
                 "but wasn't explicitly requested. Type was \"" + type + "\"\n");
    }
    return fctr*x;
}

bool getOptionalFloat(const XML_Node& parent,
                      const std::string& name,
                      doublereal& fltRtn,
                      const std::string& type)
{
    if (parent.hasChild(name)) {
        fltRtn = getFloat(parent, name, type);
        return true;
    }
    return false;
}

bool getOptionalModel(const XML_Node& parent, const std::string& nodeName,
                      std::string& modelName)
{
    if (parent.hasChild(nodeName)) {
        modelName = parent.child(nodeName)["model"];
        return true;
    }
    return false;
}

int getInteger(const XML_Node& parent, const std::string& name)
{
    if (!parent.hasChild(name)) {
        throw CanteraError("getInteger (called from XML Node \"" +
                           parent.name() + "\"): ",
                           "no child XML element named " + name);
    }
    const XML_Node& node = parent.child(name);
    int x = node.int_value();
    const string& vmin = node["min"];
    const string& vmax = node["max"];
    if (vmin != "" && x < intValue(vmin)) {
        writelog("\nWarning: value "+node.value()+" is below lower limit of "
                 +vmin+".\n");
    }
    if (node["max"] != "" && x > intValue(vmax)) {
        writelog("\nWarning: value "+node.value()+" is above upper limit of "
                 +vmax+".\n");
    }
    return x;
}

size_t getFloatArray(const XML_Node& node, vector_fp & v,
                     const bool convert, const std::string& unitsString,
                     const std::string& nodeName)
{
    const XML_Node* readNode = &node;
    if (node.name() != nodeName) {
        vector<XML_Node*> ll = node.getChildren(nodeName);
        if (ll.size() == 0) {
            throw CanteraError("getFloatArray",
                               "wrong XML element type/name: was expecting "
                               + nodeName + "but accessed " + node.name());
        } else {
            readNode = ll[0];
            ll = readNode->getChildren("floatArray");
            if (ll.size() > 0) {
                readNode = ll[0];
            }
        }
    }

    v.clear();
    doublereal vmin = Undef, vmax = Undef;
    doublereal funit = 1.0;

    // Get the attributes field, units, from the XML node
    std::string units = readNode->attrib("units");
    if (units != "" && convert) {
        if (unitsString == "actEnergy" && units != "") {
            funit = actEnergyToSI(units);
        } else if (unitsString != "" && units != "") {
            funit = toSI(units);
        }
    }

    if (readNode->attrib("min") != "") {
        vmin = fpValueCheck(readNode->attrib("min"));
    }
    if (readNode->attrib("max") != "") {
        vmax = fpValueCheck(readNode->attrib("max"));
    }

    std::string val = readNode->value();
    while (true) {
        size_t icom = val.find(',');
        if (icom != string::npos) {
            string numstr = val.substr(0,icom);
            val = val.substr(icom+1,val.size());
            v.push_back(fpValueCheck(numstr));
        } else {
            // This little bit of code is to allow for the possibility of a
            // comma being the last item in the value text. This was allowed in
            // previous versions of Cantera, even though it would appear to be
            // odd. So, we keep the possibility in for backwards compatibility.
            if (!val.empty()) {
                v.push_back(fpValueCheck(val));
            }
            break;
        }
        doublereal vv = v.back();
        if (vmin != Undef && vv < vmin - Tiny) {
            writelog("\nWarning: value {} is below lower limit of {}.\n",
                     vv, vmin);
        }
        if (vmax != Undef && vv > vmax + Tiny) {
            writelog("\nWarning: value {} is above upper limit of {}.\n",
                     vv, vmax);
        }
    }
    for (size_t n = 0; n < v.size(); n++) {
        v[n] *= funit;
    }
    return v.size();
}

void getMap(const XML_Node& node, std::map<std::string, std::string>& m)
{
    std::vector<std::string> v;
    getStringArray(node, v);
    for (size_t i = 0; i < v.size(); i++) {
        size_t icolon = v[i].find(":");
        if (icolon == string::npos) {
            throw CanteraError("getMap","missing colon in map entry ("
                               +v[i]+")");
        }
        m[v[i].substr(0,icolon)] = v[i].substr(icolon+1, v[i].size());
    }
}

int getPairs(const XML_Node& node, std::vector<std::string>& key,
             std::vector<std::string>& val)
{
    vector<string> v;
    getStringArray(node, v);
    int n = static_cast<int>(v.size());
    for (int i = 0; i < n; i++) {
        size_t icolon = v[i].find(":");
        if (icolon == string::npos) {
            throw CanteraError("getPairs","Missing a colon in the Pair entry ("
                               +v[i]+")");
        }
        key.push_back(v[i].substr(0,icolon));
        val.push_back(v[i].substr(icolon+1, v[i].size()));
    }
    return n;
}

void getMatrixValues(const XML_Node& node,
                     const std::vector<std::string>& keyStringRow,
                     const std::vector<std::string>& keyStringCol,
                     Array2D& retnValues, const bool convert,
                     const bool matrixSymmetric)
{
    if (keyStringRow.size() > retnValues.nRows()) {
        throw CanteraError("getMatrixValues",
                           "size of key1 greater than numrows");
    } else if (keyStringCol.size() > retnValues.nColumns()) {
        throw CanteraError("getMatrixValues",
                           "size of key2 greater than num cols");
    } else if (matrixSymmetric && retnValues.nRows() != retnValues.nColumns()) {
        throw CanteraError("getMatrixValues",
                           "nrow != ncol for a symmetric matrix");
    }

    // Get the attributes field, units, from the XML node and determine the
    // conversion factor, funit.
    doublereal funit = 1.0;
    if (convert && node["units"] != "") {
        funit = toSI(node["units"]);
    }

    vector<string> v;
    getStringArray(node, v);
    for (size_t i = 0; i < v.size(); i++) {
        size_t icolon = v[i].find(":");
        if (icolon == string::npos) {
            throw CanteraError("getMatrixValues","Missing two colons ("
                               +v[i]+")");
        }
        string key1 = v[i].substr(0,icolon);
        string rmm = v[i].substr(icolon+1, v[i].size());

        icolon = rmm.find(":");
        if (icolon == string::npos) {
            throw CanteraError("getMatrixValues","Missing one colon ("
                               +v[i]+")");
        }

        size_t irow = find(keyStringRow.begin(), keyStringRow.end(), key1)
                      - keyStringRow.begin();
        if (irow == keyStringRow.size()) {
            throw CanteraError("getMatrixValues","Row not matched by string: "
                               + key1);
        }

        string key2 = rmm.substr(0,icolon);
        size_t icol = find(keyStringCol.begin(), keyStringCol.end(), key2)
                      - keyStringCol.begin();
        if (icol == keyStringCol.size()) {
            throw CanteraError("getMatrixValues","Col not matched by string: "
                               + key2);
        }
        double dval = fpValueCheck(rmm.substr(icolon+1, rmm.size())) * funit;

        // Finally, insert the value;
        retnValues(irow, icol) = dval;
        if (matrixSymmetric) {
            retnValues(icol, irow) = dval;
        }
    }
}

void getStringArray(const XML_Node& node, std::vector<std::string>& v)
{
    tokenizeString(node.value(), v);
}

}
