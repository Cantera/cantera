/**
 * @file ctml.h
 * CTML ("Cantera Markup Language") is the variant of XML that Cantera uses
 * to store data. These functions read and write it.
 * (see \ref inputfiles and importCTML, ck2ctml)
 */
// Copyright 2002  California Institute of Technology

#ifndef CT_CTML_H
#define CT_CTML_H

#include "ct_defs.h"
#include "xml.h"
#include "Array.h"

//! The ctml namespace adds functionality to the XML object, by providing
//! standard functions that read, write, and interpret XML files and
//! object trees.
/*!
 *  Standardization of reads and write from Cantera files occur here.
 */
namespace ctml
{

//! const Specifying the CTML version number
/*!
 * @todo Codify what the CTML_Version number means.
 */
const std::string CTML_Version = "1.4.1";

extern std::string FP_Format;

extern std::string INT_Format;

//!  This function adds a child node with the name, "integer", with a value
//!  consisting of a single integer
/*!
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
 *
 * @todo I don't think this is used. Figure out what is used for writing integers,
 *       and codify that. unitsString shouldn't be here, since it's an int.
 *       typeString should be codified as to its usage.
 */
void addInteger(Cantera::XML_Node& node, const std::string& titleString,
                const int value, const std::string unitsString="",
                const std::string typeString="");

//!  This function adds a child node with the name, "float", with a value
//!  consisting of a single floating point number
/*!
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
   doublereal  value = 50.3;
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
 *   @param minval        Minimum allowed value of the float. The default is the
 *                        special double, Cantera::Undef, which means to ignore the
 *                        entry.
 *   @param maxval        Maximum allowed value of the float. The default is the
 *                        special double, Cantera::Undef, which means to ignore the
 *                        entry.
 *
 * @todo I don't think this is used. Figure out what is used for writing floats,
 *       and codify that. minval and maxval should be codified.
 *       typeString should be codified as to its usage.
 */
void addFloat(Cantera::XML_Node& node, const std::string& titleString,
              const doublereal value, const std::string unitsString="",
              const std::string typeString="", const doublereal minval = Cantera::Undef,
              const doublereal maxval = Cantera::Undef);

//!  This function adds a child node with the name, "floatArray", with a value
//!  consisting of a comma separated list of floats
/*!
 *   This function will add a child node to the current XML node, with the
 *   name "floatArray". It will have a title attribute, and the body
 *   of the XML node will be filled out with a comma separated list of
 *   doublereals.
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
 *   @param typeString    String type. This is an optional parameter. The default
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
 * @todo I don't think this is used. Figure out what is used for writing integers,
 *       and codify that. unitsString shouldn't be here, since it's an int.
 *       typeString should be codified as to its usage.
 */
void addFloatArray(Cantera::XML_Node& node,  const std::string& titleString,
                   const size_t n,  const doublereal* const values,
                   const std::string unitsString="", const std::string typeString="",
                   const doublereal minval = Cantera::Undef,
                   const doublereal maxval = Cantera::Undef);

//!  This function adds a child node with the name string with a string value
//!  to the current node
/*!
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
void addString(Cantera::XML_Node& node,  const std::string& titleString,
               const std::string& valueString, const std::string typeString="");


//!  This function reads the current node or a  child node of the current node
//!  with the default name, "floatArray", with a value field
//!  consisting of a comma separated list of floats
/*!
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
 *   @return              Returns the number of floats read into v.
 */
size_t getFloatArray(const Cantera::XML_Node& node, std::vector<doublereal> & v,
                     const bool convert=true, const std::string unitsString="",
                     const std::string nodeName = "floatArray");

//! This function interprets the value portion of an XML element
//! as a string. It then separates the string up into tokens
//! according to the location of white space.
/*!
 * The separate tokens are returned in the string vector
 *
 * @param node   Node to get the value from
 * @param v      Output vector containing the string tokens
 */
void getStringArray(const Cantera::XML_Node& node, std::vector<std::string>& v);


//! This routine is used to interpret the value portions of XML
//! elements that contain colon separated pairs.
/*!
 *  These are used, for example, in describing the element
 *  composition of species.
 * @verbatim
        <atomArray> H:4 C:1 <atomArray\>
   @endverbatim
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
void getMap(const Cantera::XML_Node& node, std::map<std::string, std::string>& m);

//! This function interprets the value portion of an XML element
//! as a series of "Pairs" separated by white space.
/*!
 * Each pair consists of non-whitespace characters.
 * The first ":" found in the pair string is used to separate
 * the string into two parts. The first part is called the "key"
 * The second part is called the "val".
 * String vectors of key[i] and val[i] are returned in the
 * argument list.
 * Warning: No spaces are allowed in each pair. Quotes get included as part
 *          of the string.
 *   Example:
 *       @verbatim
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
             std::vector<std::string>& val);

//! This function interprets the value portion of an XML element
//! as a series of "Matrix ids and entries" separated by white space.
/*!
 * Each pair consists of non-whitespace characters.
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
 *
 * @verbatim
     <xmlNode>
         red:green:112
         blue:black:3.3E-23

     </xmlNode>
 @endverbatim
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
 *                       attribute, then conversion to SI units is carried
 *                       out. Default is true.
 *  @param matrixSymmetric  If true entries are made so that the matrix
 *                       is always symmetric. Default is false.
 */
void getMatrixValues(const Cantera::XML_Node& node,
                     const std::vector<std::string>& keyStringRow,
                     const std::vector<std::string>& keyStringCol,
                     Cantera::Array2D& returnValues, const bool convert = true,
                     const bool matrixSymmetric = false);


//!  Get a vector of integer values from a child element.
/*!
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
void getIntegers(const Cantera::XML_Node& node, std::map<std::string,int>& v);


//!  Get a floating-point value from a child element.
/*!
 *  Returns a doublereal value for the child named 'name' of element 'parent'. If
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
 *                 and "" , for no conversion. The default value is "",
 *                 which implies that no conversion is allowed.
 */
doublereal getFloat(const Cantera::XML_Node& parent, const std::string& name,
                    const std::string type="");

//!  Get a floating-point value from the current XML element
/*!
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
 *  Reads the corresponding XML file:
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
doublereal getFloatCurrent(const Cantera::XML_Node& currXML, const std::string type="");

//!  Get an optional floating-point value from a child element.
/*!
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
bool getOptionalFloat(const Cantera::XML_Node& parent, const std::string& name,
                      doublereal& fltRtn, const std::string type="");


//!  Get a vector of floating-point values from a child element.
/*!
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
               const bool convert=true);

//!  Get an integer value from a child element.
/*!
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
   <numProcs> 10 <numProcs/>
 <\state>
 @endverbatim
 *
 *   @param parent reference to the XML_Node object of the parent XML element
 *   @param name   Name of the XML child element
 */
int getInteger(const Cantera::XML_Node& parent, const std::string& name);

//!  Get a floating-point value from a child element with a defined units field
/*!
 *  Returns a doublereal value for the child named 'name' of element 'parent'.
 *  'type' must be supplied and match a known unit type.
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
 pres = getFloatDefaultUnits(State_XMLNode, "pressure", "Pa", "toSI");
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
 *   @param defaultUnits Default units string to be found in the units attribute.
 *                 If the units string in the XML field is equal to defaultUnits,
 *                 no units conversion will be carried out.
 *   @param type   String type. Currently known types are "toSI" and "actEnergy",
 *                 and "" , for no conversion. The default value is "",
 *                 which implies that no conversion is allowed.
 */
doublereal getFloatDefaultUnits(const Cantera::XML_Node& parent,
                                const std::string& name,
                                const std::string& defaultUnits,
                                const std::string& type="toSI");

//!  Get an optional model name from a named child node.
/*!
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
 *  reads the corresponding XML file:
 *  @verbatim
  <transport model="Simple">
    <compositionDependence model="Solvent_Only"/>
  </transport>
 @endverbatim
 *
 *   On return modelName is set to "Solvent_Only".
 *
 *   @param parent reference to the XML_Node object of the parent XML element
 *   @param nodeName   Name of the XML child element
 *   @param modelName  On return this contains the contents of the model attribute
 *
 *   @return True if the nodeName XML node exists. False otherwise.
 */
bool getOptionalModel(const Cantera::XML_Node& parent, const std::string nodeName,
                      std::string& modelName);

//! Search the child nodes of the current node for an XML Node with a Title
//! attribute of a given name.
/*!
 *   @param node   Current node from which to conduct the search
 *   @param title  Name of the title attribute
 *
 *   @return  Returns a pointer to the matched child node. Returns 0 if no node is
 *            found.
 */
Cantera::XML_Node* getByTitle(const Cantera::XML_Node& node, const std::string& title);

//!  This function attempts to read a named child node and returns with the contents in the value string.
//!  title attribute named "titleString"
/*!
 *   This function will read a child node to the current XML node, with the
 *   name "string". It must have a title attribute, named titleString, and the body
 *   of the XML node will be read into the valueString output argument.
 *
 *   If the child node is not found then the empty string is returned.
 *
 *  Example:
 *
 * Code snippet:
 *       @verbatim
   const XML_Node &node;
   std::string valueString;
   std::string typeString;
   std::string nameString = "timeIncrement";
   getString(XML_Node& node, nameString, valueString, valueString, typeString);
 @endverbatim
 *
 *  Reads the following the snippet in the XML file:
 *
 *  *  @verbatim
 <nameString type="typeString">
   valueString
 <\nameString>
 @endverbatim
 *
 *  or alternatively as a retrofit and special case, it also reads the following case
 *
 *  @verbatim
 <string title="nameString" type="typeString">
   valueString
 <\string>
 @endverbatim
 *
 *   @param node          Reference to the XML_Node object of the parent XML element
 *   @param nameString    Name of the XML Node                               input  variable
 *   @param valueString   Value string that is found in the child node.      output variable
 *   @param typeString    String type. This is an optional output variable. It is filled
 *                        with the attribute "type" of the XML entry.         output variable
 */
void getNamedStringValue(const Cantera::XML_Node& node, const std::string& nameString, std::string& valueString,
                         std::string& typeString);


//!  This function reads a child node with the name, nameString, and returns
//!  its xml value as the return string
/*!
 *   If the child XML_node named "name" doesn't exist, the empty string is returned.
 *
 * Code snippet:
 *       @verbatim
 const XML_Node &parent;
 string nameString = "vacancy_species";
 string valueString = getChildValue(parent, nameString
 std::string typeString);
 @endverbatim
 *
 *  returns valueString = "O(V)"
 *
 *  from the following the snippet in the XML file:
 *
 *  @verbatim
 <vacancySpecies>
   O(V)
 <\vacancySpecies>
 @endverbatim
 *
 *   @param parent     parent reference to the XML_Node object of the parent XML element
 *   @param nameString Name of the child XML_Node to read the value from.
 *
 *   @return           String value of the child XML_Node
 */
std::string getChildValue(const Cantera::XML_Node& parent,
                          const std::string& nameString);

//! Read an ctml file from a file and fill up an XML tree
/*!
 *  This is the main routine that reads a ctml file and puts it into
 *  an XML_Node tree
 *
 *  @param node    Root of the tree
 *  @param file    Name of the file
 *  @param debug   Turn on debugging printing
 */
void get_CTML_Tree(Cantera::XML_Node* node, const std::string file,
                   const int debug = 0);

//! Read an ctml file from a file and fill up an XML tree.
//!   @param file    Name of the file
//!   @return        Root of the tree
Cantera::XML_Node getCtmlTree(const std::string file);

//! Convert a cti file into a ctml file
/*!
 *
 *  @param   file    Pointer to the file
 *  @param   debug   Turn on debug printing
 *
 *  @ingroup inputfiles
 */
void ct2ctml(const char* file, const int debug = 0);

//! Convert a Chemkin-format mechanism into a CTI file.
/*!
 * @param in_file         input file containing species and reactions
 * @param thermo_file     optional input file containing thermo data
 * @param transport_file  optional input file containing transport parameters
 * @param id_tag          id of the phase
 */
void ck2cti(const std::string& in_file, const std::string& thermo_file="",
            const std::string& transport_file="",
            const std::string& id_tag="gas");
}

#endif
