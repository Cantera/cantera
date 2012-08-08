/**
 * @file ctml.cpp
 * Definitions for functions to read and write CTML.
 *
 */
// Copyright 2002  California Institute of Technology

#include "cantera/base/ctml.h"

//@{
#define CTML_VERSION_1_4_1
//@}

#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

#include <cctype>
#include <cstring>
#include <cstdlib>


using namespace std;
using namespace Cantera;

namespace ctml
{

std::string FP_Format = "%23.15E";
std::string INT_Format = "%8d";

//====================================================================================================================
//! Convert a floating point value from a string to a double
/*!
 *   @param val String value input
 *
 *   @return Returns a double
 */
static doublereal fpValue(std::string val)
{
    return atof(stripws(val).c_str());
}

//====================================================================================================================
//  This function adds a child node with the name, "integer", with a value
//  consisting of a single integer
/*
 *   This function will add a child node to the current XML node, with the
 *   name "integer". It will have a title attribute, and the body
 *   of the XML node will be filled out with a single integer
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
 const XML_Node &node;
 std::string titleString = "maxIterations";
 int  value = 1000;
 std::string typeString = "optional";
 std::string units = "";
 addInteger(node, titleString, value, typeString, units);
 @endverbatim
 *
 *  Creates the following the snippet in the XML file:
 *  @verbatim
 <parentNode>
 <integer title="maxIterations" type="optional">
 100
 <\integer>
 <\parentNode>
 @endverbatim
 *
 *   @param node          reference to the XML_Node object of the parent XML element
 *   @param titleString   String name of the title attribute
 *   @param value         Value - single integer
 *   @param unitsString   String name of the Units attribute. The default is to
 *                        have an empty string.
 *   @param typeString    String type. This is an optional parameter. The default
 *                        is to have an empty string.
 */
void addInteger(Cantera::XML_Node& node, const std::string& title, const int val,
                const std::string units, const std::string type)
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

//====================================================================================================================
//  This function adds a child node with the name, "float", with a value
//  consisting of a single floating point number
/*
 *   This function will add a child node to the current XML node, with the
 *   name "float". It will have a title attribute, and the body
 *   of the XML node will be filled out with a single float
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
 const XML_Node &node;
 std::string titleString = "activationEnergy";
 double  value = 50.3;
 doublereal maxval = 1.0E3;
 doublereal minval = 0.0;
 std::string typeString = "optional";
 std::string unitsString = "kcal/gmol";
 addFloat(node, titleString, value, unitsString, typeString, minval, maxval);
 @endverbatim
 *
 *  Creates the following the snippet in the XML file:
 *  @verbatim
 <parentNode>
 <float title="activationEnergy" type="optional" units="kcal/gmol" min="0.0" max="1.0E3">
 50.3
 <\float>
 <\parentNode>
 @endverbatim
 *
 *   @param node          reference to the XML_Node object of the parent XML element
 *   @param titleString   String name of the title attribute
 *   @param value         Value - single integer
 *   @param unitsString   String name of the Units attribute. The default is to
 *                        have an empty string.
 *   @param typeString    String type. This is an optional parameter. The default
 *                        is to have an empty string.
 *
 * @todo minval and maxval should be codified. typeString should be codified
 *       as to its usage.
 */
void addFloat(Cantera::XML_Node& node, const std::string& title,
              const doublereal val, const std::string units,
              const std::string type, const doublereal minval,
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
//====================================================================================================================
//  This function adds a child node with the name, "floatArray", with a value
//  consisting of a comma separated list of floats
/*
 *   This function will add a child node to the current XML node, with the
 *   name "floatArray". It will have a title attribute, and the body
 *   of the XML node will be filled out with a comma separated list of
 *   integers
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
   const XML_Node &node;
   std::string titleString = "additionalTemperatures";
   int  n = 3;
   int Tcases[3] = [273.15, 298.15, 373.15];
   std::string typeString = "optional";
   std::string units = "Kelvin";
   addFloatArray(node, titleString, n, &cases[0], typeString, units);
   @endverbatim
 *
 *  Creates the following the snippet in the XML file:
 *  @verbatim
   <parentNode>
     <floatArray title="additionalTemperatures" type="optional" units="Kelvin">
        273.15, 298.15, 373.15
     <\floatArray>
   <\parentNode>
 @endverbatim
 *
 *   @param node          reference to the XML_Node object of the parent XML element
 *   @param titleString   String name of the title attribute
 *   @param n             Length of the doubles vector.
 *   @param values        Pointer to a vector of doubles
 *   @param unitsString   String name of the Units attribute. This is an optional
 *                        parameter. The default is to
 *                        have an empty string.
 *   @param type          String type. This is an optional parameter. The default
 *                        is to have an empty string.
 *   @param minval        Minimum allowed value of the int. This is an optional
 *                        parameter. The default is the
 *                        special double, Cantera::Undef, which means to ignore the
 *                        entry.
 *   @param maxval        Maximum allowed value of the int. This is an optional
 *                        parameter. The default is the
 *                        special double, Cantera::Undef, which means to ignore the
 *                        entry.
 *
 * @todo typeString should be codified as to its usage.
 */
void addFloatArray(Cantera::XML_Node& node, const std::string& title, const size_t n,
                   const doublereal* const vals, const std::string units,
                   const std::string type,
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

//====================================================================================================================
//  This function adds a child node with the name string with a string value
//  to the current node
/*
 *   This function will add a child node to the current XML node, with the
 *   name "string". It will have a title attribute, and the body
 *   of the XML node will be filled out with the valueString argument verbatim.
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
 const XML_Node &node;
 addString(XML_Node& node, std::string titleString, std::string valueString,
 std::string typeString);
 @endverbatim
 *
 *  Creates the following the snippet in the XML file:
 *  @verbatim
 <string title="titleString" type="typeString">
 valueString
 <\string>
 @endverbatim
 *
 *   @param node          reference to the XML_Node object of the parent XML element
 *   @param valueString   Value string to be used in the new XML node.
 *   @param titleString   String name of the title attribute
 *   @param typeString    String type. This is an optional parameter.
 */
void addString(Cantera::XML_Node& node, const std::string& titleString,
               const std::string& valueString,
               const std::string typeString)
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
//====================================================================================================================
//  This function reads a child node with the name string and returns
//  its xml value as the return string
/*
 *   If the child XML_node named "name" doesn't exist, the empty string is returned.
 *
 * Code snippet:
 *       @verbatim
 const XML_Node &parent;
 string name = "vacency_species";
 string valueString = getChildValue(parent, name
 std::string typeString);
 @endverbatim
 *
 *  returns valueString = "O(V)"
 *
 *  from the following the snippet in the XML file:
 *
 *  @verbatim
 <vacencySpecies>
 O(V)
 <\vancencySpecies>
 @endverbatim
 *
 *   @param parent   parent reference to the XML_Node object of the parent XML element
 *   @param name     Name of the child XML_Node to read the value from.
 *
 *   @return         String value of the child XML_Node
 */
std::string getChildValue(const Cantera::XML_Node& parent, const std::string& nameString)
{
    if (!parent.hasChild(nameString)) {
        return "";
    }
    return parent(nameString);
}

//====================================================================================================================
//  Get a vector of integer values from a child element.
/*
 *  Returns a std::map containing a keyed values for child XML_Nodes
 *  of the current node with the name, "integer".
 *  In the keyed mapping there will be a list of titles vs. values
 *  for all of the XML nodes.
 *  The integer XML_nodes are expected to be in a particular form created
 *  by the function addInteger(). One value per XML_node is expected.
 *
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
   const XML_Node &State_XMLNode;
   std::map<std::string, integer> v;
   getinteger(State_XMLNode, v);
 @endverbatim
 *
 *  reads the corresponding XML file:
 *
 *  @verbatim
 <state>
   <integer title="i1">   1  <\integer>
   <integer title="i2">   2  <\integer>
   <integer title="i3">   3  <\integer>
 <\state>
 @endverbatim
 *
 *  Will produce the mapping:
 *
 *         v["i1"] = 1
 *         v["i2"] = 2
 *         v["i3"] = 3
 *
 *
 *   @param node     Current XML node to get the values from
 *   @param v        Output map of the results.
 */
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

//====================================================================================================================
//  Get a vector of floating-point values from a child element.
/*
 *  Returns a std::map containing a keyed values for child XML_Nodes
 *  of the current node with the name, "float".
 *  In the keyed mapping there will be a list of titles vs. values
 *  for all of the XML nodes.
 *  The float XML_nodes are expected to be in a particular form created
 *  by the function addFloat(). One value per XML_node is expected.
 *
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
   const XML_Node &State_XMLNode;
   std::map<std::string,double> v;
   bool convert = true;
   getFloats(State_XMLNode, v, convert);
 @endverbatim
 *
 *  reads the corresponding XML file:
 *
 *  @verbatim
 <state>
   <float title="a1" units="m3">   32.4 <\float>
   <float title="a2" units="cm3">   1.  <\float>
   <float title="a3">             100.  <\float>
 <\state>
 @endverbatim
 *
 *  Will produce the mapping:
 *
 *         v["a1"] = 32.4
 *         v["a2"] = 1.0E-6
 *         v["a3"] = 100.
 *
 *
 *   @param node     Current XML node to get the values from
 *   @param v        Output map of the results.
 *   @param convert  Turn on conversion to SI units
 */
void getFloats(const Cantera::XML_Node& node, std::map<std::string, double>& v,
               const bool convert)
{
    std::vector<XML_Node*> f;
    node.getChildren("float",f);
    int n = static_cast<int>(f.size());
    doublereal x, x0, x1, fctr;
    std::string typ, title, units, vmin, vmax;
    for (int i = 0; i < n; i++) {
        const XML_Node& fi = *(f[i]);
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
                writelog("\nWarning: value "+fi()+" is below lower limit of "
                         +vmin+".\n");
            }
        }
        if (fi["max"] != "") {
            x1 = atof(vmax.c_str());
            if (x > x1 + Tiny) {
                writelog("\nWarning: value "+fi()+" is above upper limit of "
                         +vmax+".\n");
            }
        }
        fctr = (convert ? toSI(units) : 1.0); // toSI(typ,units);
        v[title] = fctr*x;
    }
}

//====================================================================================================================
//  Get a floating-point value from a child element.
/*
 *  Returns a double value for the child named 'name' of element 'parent'. If
 *  'type' is supplied and matches a known unit type, unit
 *  conversion to SI will be done if the child element has an attribute
 *  'units'.
 *
 *  Note, it's an error for the child element not to exist.
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
 const XML_Node &State_XMLNode;
 doublereal pres = OneAtm;
 if (state_XMLNode.hasChild("pressure")) {
 pres = getFloat(State_XMLNode, "pressure", "toSI");
 }
 @endverbatim
 *
 *  reads the corresponding XML file:
 *  @verbatim
 <state>
 <pressure units="Pa"> 101325.0 </pressure>
 <\state>
 @endverbatim
 *
 *   @param parent reference to the XML_Node object of the parent XML element
 *   @param name   Name of the XML child element
 *   @param type   String type. Currently known types are "toSI" and "actEnergy",
 *                 and "" , for no conversion. The default value is ""
 *                 which implies that no conversion is allowed.
 */
doublereal getFloat(const Cantera::XML_Node& parent,
                    const std::string& name,
                    const std::string type)
{
    if (!parent.hasChild(name))
        throw CanteraError("getFloat (called from XML Node \"" +
                           parent.name() + "\"): ",
                           "no child XML element named \"" + name + "\" exists");
    const XML_Node& node = parent.child(name);
    return getFloatCurrent(node, type);
}

//====================================================================================================================
//  Get a floating-point value from the current XML element
/*
 *  Returns a doublereal value from the current element. If
 *  'type' is supplied and matches a known unit type, unit
 *  conversion to SI will be done if the child element has an attribute  'units'.
 *
 *  Note, it's an error for the child element not to exist.
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
 const XML_Node &State_XMLNode;
 doublereal pres = OneAtm;
 if (state_XMLNode.hasChild("pressure")) {
   XML_Node *pres_XMLNode = State_XMLNode.getChild("pressure");
   pres = getFloatCurrent(pres_XMLNode, "toSI");
 }
 @endverbatim
 *
 *  Rreads the corresponding XML file:
 *  @verbatim
 <state>
   <pressure units="Pa"> 101325.0 </pressure>
 <\state>
 @endverbatim
 *
 *   @param currXML reference to the current XML_Node object
 *   @param type   String type. Currently known types are "toSI" and "actEnergy",
 *                 and "" , for no conversion. The default value is "",
 *                 which implies that no conversion is allowed.
 */
doublereal getFloatCurrent(const Cantera::XML_Node& node,
                           const std::string type)
{
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
            writelog("\nWarning: value "+node()+" is below lower limit of "
                     +vmin+".\n");
        }
    }
    if (node["max"] != "") {
        x1 = atof(vmax.c_str());
        if (x > x1 + Tiny) {
            writelog("\nWarning: value "+node()+" is above upper limit of "
                     +vmax+".\n");
        }
    }
    // Note, most type's of converters default to toSI() type atm.
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
                 "but wasn't explicity requested. Type was \"" + type + "\"\n");
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

//====================================================================================================================
bool getOptionalFloat(const Cantera::XML_Node& parent,
                      const std::string& name,
                      doublereal& fltRtn,
                      const std::string type)
{
    if (parent.hasChild(name)) {
        fltRtn= getFloat(parent, name, type);
        return true;
    }
    return false;
}
//====================================================================================================================
//  Get an optional floating-point value from a child element.
/*
 *  Returns a doublereal value for the child named 'name' of element 'parent'. If
 *  'type' is supplied and matches a known unit type, unit
 *  conversion to SI will be done if the child element has an attribute
 *  'units'.
 *
 *
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
 const XML_Node &State_XMLNode;
 doublereal pres = OneAtm;
 bool exists = getOptionalFloat(State_XMLNode, "pressure", pres, "toSI");
 @endverbatim
 *
 *  reads the corresponding XML file:
 *  @verbatim
 <state>
   <pressure units="Pa"> 101325.0 </pressure>
 <\state>
 @endverbatim
 *
 *   @param parent reference to the XML_Node object of the parent XML element
 *   @param name   Name of the XML child element
 *   @param fltRtn Float Return. It will be overridden if the XML
 *                 element exists.
 *   @param type   String type. Currently known types are "toSI"
 *                 and "actEnergy",
 *                 and "" , for no conversion. The default value is "",
 *                 which implies that no conversion is allowed.
 *
 * @return returns true if the child element named "name" exists
 */
doublereal getFloatDefaultUnits(const Cantera::XML_Node& parent, std::string name,
                                std::string defaultUnits, std::string type)
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
//====================================================================================================================
//  Get an optional model name from a named child node.
/*
 *  Returns the model name attribute for the child named 'nodeName' of element 'parent'.
 *  Note, it's optional for the child node to exist
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
 std::string modelName = "";
 bool exists = getOptionalModel(transportNode, "compositionDependence",
                  modelName);
 @endverbatim
 *
 *  Reads the corresponding XML file:
 *
 *  @verbatim
  <transport>
    <compositionDependence model="Solvent_Only"/>
  </transport>
 @endverbatim
 *
 *   On return modelName is set to "Solvent_Only".
 *
 *   @param parent     Reference to the XML_Node object of the parent XML element
 *   @param nodeName   Name of the XML child element
 *   @param modelName  On return this contains the contents of the model attribute
 *
 *   @return True if the nodeName XML node exists. False otherwise.
 */
bool getOptionalModel(const Cantera::XML_Node& parent, const std::string nodeName,
                      std::string& modelName)
{
    if (parent.hasChild(nodeName)) {
        const XML_Node& node = parent.child(nodeName);
        modelName = node["model"];
        return true;
    }
    return false;
}
//====================================================================================================================
//  Get an integer value from a child element.
/*
 *  Returns an integer value for the child named 'name' of element 'parent'.
 *
 *  Note, it's an error for the child element not to exist.
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
 const XML_Node &State_XMLNode;
 int number = 1;
 if (state_XMLNode.hasChild("NumProcs")) {
 number = getInteger(State_XMLNode, "numProcs");
 }
 @endverbatim
 *
 *  reads the corresponding XML file:
 *  @verbatim
 <state>
 <numProcs> 10 <\numProcs>
 <\state>
 @endverbatim
 *
 *   @param parent reference to the XML_Node object of the parent XML element
 *   @param name   Name of the XML child element
 */
int getInteger(const Cantera::XML_Node& parent, std::string name)
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
//====================================================================================================================
//  This function reads the current node or a  child node of the current node
//  with the default name, "floatArray", with a value field
//  consisting of a comma separated list of floats
/*
 *   This function will read either the current XML node or a  child node
 *   to the current XML node, with the
 *   name "floatArray". It will have a title attribute, and the body
 *   of the XML node will be filled out with a comma separated list of
 *   doublereals.
 *     Get an array of floats from the XML Node. The argument field
 *   is assumed to consist of an arbitrary number of comma
 *   separated floats, with an arbitrary amount of white space
 *   separating each field.
 *      If the node array has an units attribute field, then
 *   the units are used to convert the floats, iff convert is true.
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
   const XML_Node &State_XMLNode;
   vector_fp v;
   bool convert = true;
   unitsString = "";
   nodeName="floatArray";
   getFloatArray(State_XMLNode, v, convert, unitsString, nodeName);
 @endverbatim
 *
 *  reads the corresponding XML file:
 *
 *  @verbatim
 <state>
   <floatArray  units="m3">   32.4, 1, 100. <\floatArray>
 <\state>
 @endverbatim
 *
 *  Will produce the vector
 *
 *         v[0] = 32.4
 *         v[1] = 1.0
 *         v[2] = 100.
 *
 *
 *   @param  node         XML parent node of the floatArray
 *   @param  v            Output vector of floats containing the floatArray information.
 *   @param  convert      Conversion to SI is carried out if this boolean is
 *                        True. The default is true.
 *   @param  unitsString  String name of the type attribute. This is an optional
 *                        parameter. The default is to have an empty string.
 *                        The only string that is recognized is actEnergy.
 *                        Anything else has no effect. This affects what
 *                        units converter is used.
 *   @param  nodeName     XML Name of the XML node to read.
 *                        The default value for the node name is floatArray
 *
 *   @return              Returns the number of floats read
 *
 *  @note change the v to a std::vector to eliminate a doxygen error. No idea why doxygen needs this.
 */
size_t getFloatArray(const Cantera::XML_Node& node, std::vector<doublereal> & v,
                     const bool convert, const std::string unitsString,
                     const std::string nodeName)
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
        vmin = atofCheck((*readNode)["min"].c_str());
    }
    if ((*readNode)["max"] != "") {
        vmax = atofCheck((*readNode)["max"].c_str());
    }

    doublereal vv;
    std::string val = readNode->value();
    while (1 > 0) {
        icom = val.find(',');
        if (icom != string::npos) {
            numstr = val.substr(0,icom);
            val = val.substr(icom+1,val.size());
            dtmp = atofCheck(numstr.c_str());
            v.push_back(dtmp);
        } else {
            /*
             * This little bit of code is to allow for the
             * possibility of a comma being the last
             * item in the value text. This was allowed in
             * previous versions of Cantera, even though it
             * would appear to be odd. So, we keep the
             * possibilty in for backwards compatibility.
             */
            if (!val.empty()) {
                dtmp = atofCheck(val.c_str());
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

//====================================================================================================================
// This routine is used to interpret the value portions of XML
// elements that contain colon separated pairs.
/*
 *  These are used, for example, in describing the element
 *  composition of species.
 *
 *       <atomArray> H:4 C:1 <atomArray\>
 *
 *  The string is first separated into a string vector according
 * to the location of white space. Then each string is again
 * separated into two parts according to the location of a
 * colon in the string. The first part of the string is
 * used as the key, while the second part of the string is
 * used as the value, in the return map.
 * It is an error to not find a colon in each string pair.
 *
 *  @param node Current node
 *  @param m    Output Map containing the pairs of values found
 *              in the XML Node
 */
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
//====================================================================================================================
// This function interprets the value portion of an XML element
// as a series of "Pairs" separated by white space.
/*
 * Each pair consists of nonwhite-space characters.
 * The first ":" found in the pair string is used to separate
 * the string into two parts. The first part is called the "key"
 * The second part is called the "val".
 * String vectors of key[i] and val[i] are returned in the
 * argument list.
 * Warning: No spaces are allowed in each pair. Quotes get included as part
 *          of the string.
 *   Example: @verbatim
     <xmlNode>
        red:112    blue:34
        green:banana
     </xmlNode>
            @endverbatim
 *
 * Returns:
 *          key       val
 *     0:   "red"     "112"
 *     1:   "blue"    "34"
 *     2:   "green"   "banana"
 *
 *
 *  @param node             XML Node
 *  @param key              Vector of keys for each entry
 *  @param val              Vector of values for each entry
 *
 *  @return Returns the number of pairs found
 */
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
//====================================================================================================================
// This function interprets the value portion of an XML element
// as a series of "Matrix ids and entries" separated by white space.
/*
 * Each pair consists of nonwhite-space characters.
 * The first two ":" found in the pair string is used to separate
 * the string into three parts. The first part is called the first
 * key. The second part is the second key. Both parts must match
 * an entry in the keyString1 and keyString2, respectively,
 * in order to provide a location to
 * place the object in the matrix.
 * The third part is called the value. It is expected to be
 * a double. It is translated into a double and placed into the
 * correct location in the matrix.
 *
 * Warning: No spaces are allowed in each triplet. Quotes are part
 *          of the string.
 *   Example
 *         keyString = red, blue, black, green
 *    <xmlNode>
 *        red:green:112
 *        blue:black:3.3E-23
 *
 *    </xmlNode>
 *
 * Returns:
 *     retnValues(0, 3) = 112
 *     retnValues(1, 2) = 3.3E-23
 *
 *
 *  @param node          XML Node containing the information for the matrix
 *  @param keyStringRow  Key string for the row
 *  @param keyStringCol  Key string for the column entries
 *  @param returnValues  Return Matrix.
 *  @param convert       If this is true, and if the node has a units
 *                       attribute, then conversion to si units is carried
 *                       out. Default is true.
 *  @param matrixSymmetric  If true entries are made so that the matrix
 *                       is always symmetric. Default is false.
 */
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
        double dval = atofCheck(val.c_str());
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
//====================================================================================================================
// This function interprets the value portion of an XML element
// as a string. It then separates the string up into tokens
// according to the location of white space.
/*
 * The separate tokens are returned in the string vector
 *
 * @param node   Node to get the value from
 * @param v      Output vector containing the string tokens
 */
void getStringArray(const Cantera::XML_Node& node, std::vector<std::string>& v)
{
    std::string val = node.value();
    tokenizeString(val, v);
}

}
