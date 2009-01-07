/**
 * @file ctml.cpp
 * Definitions for functions to read and write CTML.
 *
 */

/*
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

#include "global.h"
#include <cctype>
#include <cstring>


using namespace std;
using namespace Cantera;

namespace ctml {

    static doublereal fpValue(std::string val) {
        return atof(stripws(val).c_str());
    }

    void addBool(XML_Node& node, string title, bool val) {
        string v = (val ? "true" : "false");
        XML_Node& f = node.addChild("bool",v);
        f.addAttribute("title",title);
    }

  //  This function adds a child node with the name, "integer", with a value
  //  consisting of a single integer
  /*
   *   This function will add a child node to the current XML node, with the
   *   name "integer". It will have a title attribute, and the body
   *   of the XML node will be filled out with a single integer
   *
   *  Example:
   *
   * Code snipet:
   *       @verbatum
         const XML_Node &node;
	 std::string titleString = "maxIterations";
	 int  value = 1000;
	 std::string typeString = "optional";
	 std::string units = "";
         addInteger(node, titleString, value, typeString, units);
   @endverbatum
   *
   *  Creates the following the snippet in the XML file:
   *  @verbatum
     <parentNode>
       <integer title="maxIterations" type="optional">
          100
       <\integer>
     <\parentNode>
   @endverbatum
   *
   *   @param node          reference to the XML_Node object of the parent XML element
   *   @param titleString   String name of the title attribute
   *   @param value         Value - single integer
   *   @param unitsString   String name of the Units attribute. The default is to
   *                        have an empty string.
   *   @param typeString    String type. This is an optional parameter. The default
   *                        is to have an empty string.
   */
    void addInteger(XML_Node& node, const std::string &title, const int val, 
                    const std::string units, const std::string type) {
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
#else
    XML_Node& f = node.addChild(title,val,fmt); 
#endif
    if (type != "") f.addAttribute("type",type);
    if (units != "") f.addAttribute("units",units);
    if (minval != Undef) f.addAttribute("min",minval);
    if (maxval != Undef) f.addAttribute("max",maxval);
  }

    /**
     *  Add a floatArray XML type to the xml file.
     *  This is a generic XML entry containing a vector
     *  of doubles as its values and containing a set
     *  of attributes that describes the length of the
     *  vector, and optionally the units of the vector.
     *
     *  Note, a comma is not put after the last double
     *  entry anymore.
     */
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

  //  This function adds a child node with the name string with a string value
  //  to the current node
  /* 
   *   This function will add a child node to the current XML node, with the
   *   name "string". It will have a title attribute, and the body
   *   of the XML node will be filled out with the valueString argument verbatim.
   *
   *  Example:  
   *
   * Code snipet:
   *       @verbatum
         const XML_Node &node;
	 addString(XML_Node& node, std::string titleString, std::string valueString, 
	           std::string typeString);
   @endverbatum
   *
   *  Creates the following the snippet in the XML file:
   *  @verbatum
     <string title="titleString" type="typeString">
        valueString
     <\string>
   @endverbatum
   *
   *   @param node          reference to the XML_Node object of the parent XML element
   *   @param valueString   Value string to be used in the new XML node.
   *   @param titleString   String name of the title attribute
   *   @param typeString    String type. This is an optional parameter.
   */
  void addString(XML_Node& node, std::string titleString, std::string valueString, 
		 std::string typeString) {
    XML_Node& f = node.addChild("string", valueString);
    f.addAttribute("title", titleString);
    if (typeString != "") f.addAttribute("type", typeString);        
  }

    XML_Node* getByTitle(XML_Node& node, string title) {
        XML_Node* s = node.findByAttr("title",title);
        if (!s) return 0;
        if (s->parent() == &node) return s;
        else return 0;
    }

  //  This function reads a child node with the name string and returns
  //  its xml value as the return string
  /*
   *   If the child XML_node named "name" doesn't exist, the empty string is returned.
   *  
   * Code snipet:
   *       @verbatum
         const XML_Node &parent;
	 string name = "vacency_species";
	 string valueString = getChildValue(parent, name
	           std::string typeString);
   @endverbatum
   *
   *  returns valueString = "O(V)"
   * 
   *  from the following the snippet in the XML file:
   *
   *  @verbatum
     <vacencySpecies>
        O(V)
     <\vancencySpecies>
   @endverbatum
   *
   *   @param parent   parent reference to the XML_Node object of the parent XML element
   *   @param name     Name of the childe XML_Node to read the value from.
   *
   *   @return         String value of the child XML_Node
   */
  std::string getChildValue(const XML_Node& parent, const std::string &nameString) {
    if (!parent.hasChild(nameString)) return "";
    return parent(nameString);
  }

  //  This function reads a child node with the name, "string", with a specific
  //  title attribute named "titleString"
  /* 
   *   This function will read a child node to the current XML node, with the
   *   name "string". It must have a title attribute, named titleString, and the body
   *   of the XML node will be read into the valueString output argument.
   *
   *  Example:  
   *
   * Code snipet:
   *       @verbatum
         const XML_Node &node;
	 getString(XML_Node& node, std::string titleString, std::string valueString, 
	           std::string typeString);
   @endverbatum
   *
   *  Reads the following the snippet in the XML file:
   *  @verbatum
     <string title="titleString" type="typeString">
        valueString
     <\string>
   @endverbatum
   *
   *   @param node          reference to the XML_Node object of the parent XML element
   *   @param titleString   String name of the title attribute of the child node
   *   @param valueString   Value string that is found in the child node. output variable
   *   @param typeString    String type. This is an optional output variable
   */
  void getString(XML_Node& node, const std::string &titleString, std::string& valueString, 
		 std::string& typeString) {
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
  
  void getIntegers(const XML_Node& node, map<string,int>& v) {
        vector<XML_Node*> f;
        node.getChildren("integer",f);
        int n = static_cast<int>(f.size());
        integer x, x0, x1;
        string typ, title, vmin, vmax;
        for (int i = 0; i < n; i++) {
            const XML_Node& fi = *(f[i]);
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


  

    void getFloats(const XML_Node& node, map<string,double>& v, bool convert) {
        vector<XML_Node*> f;
        node.getChildren("float",f);
        int n = static_cast<int>(f.size());
        doublereal x, x0, x1, fctr;
        string typ, title, units, vmin, vmax;
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
   * Code snipet:
   *       @verbatum
      const XML_Node &State_XMLNode;
      doublereal pres = OneAtm;
      if (state_XMLNode.hasChild("pressure")) {
        pres = getFloat(State_XMLNode, "pressure", "toSI");
      }
          @endverbatum
   *
   *  reads the corresponding XML file:
   *  @verbatum
      <state>
        <pressure units="Pa"> 101325.0 </pressure>
      <\state>
      @endverbatum
   *
   *   @param parent reference to the XML_Node object of the parent XML element
   *   @param name   Name of the XML child element
   *   @param type   String type. Currently known types are "toSI" and "actEnergy",
   *                 and "" , for no conversion. The default value is ""
   *                 which implies that no conversion is allowed.
   */
  doublereal getFloat(const XML_Node& parent, std::string name, std::string type) {
    if (!parent.hasChild(name)) 
      throw CanteraError("getFloat (called from XML Node \"" +
			 parent.name() + "\"): ",
			 "no child XML element named \"" + name + "\" exists");
    const XML_Node& node = parent.child(name);
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

  
  doublereal getFloatDefaultUnits(const Cantera::XML_Node& parent, std::string name,
				  std::string defaultUnits, std::string type) {

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


  //  Get an integer value from a child element. 
  /* 
   *  Returns an integer value for the child named 'name' of element 'parent'.
   *
   *  Note, it's an error for the child element not to exist.
   *
   *  Example:  
   *
   * Code snipet:
   *       @verbatum
         const XML_Node &State_XMLNode;
         int number = 1;
         if (state_XMLNode.hasChild("NumProcs")) {
           number = getInteger(State_XMLNode, "numProcs");
         }
   @endverbatum
   *
   *  reads the corresponding XML file:
   *  @verbatum
     <state>
        <numProcs> 10 <numProcs/>
     <\state>
   @endverbatum
   *
   *   @param parent reference to the XML_Node object of the parent XML element
   *   @param name   Name of the XML child element
   */
  int getInteger(const XML_Node& parent, string name) {
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

    /*
     * getFloatArray():
     *
     * Get an array of floats from the XML Node. The argument field 
     * is assumed to consist of an arbitrary number of comma 
     * separated floats, with an arbitrary amount of white space
     * separating each field.
     *    If the node array has an units attribute field, then
     * the units are used to convert the floats, iff convert is true.
     *
     *  nodeName = The default value for the node name is floatArray
     */
    void getFloatArray(const XML_Node& node, vector_fp& v, bool convert,
                       string type, string nodeName) {
	string::size_type icom;
        string numstr;
	doublereal dtmp;
        string nn = node.name();
        if (nn != nodeName) 
            throw CanteraError("getFloatArray",
                "wrong xml element type/name: was expecting "
                 + nodeName + "but accessed " + node.name());

        v.clear();
        doublereal vmin = Undef, vmax = Undef;

        doublereal funit = 1.0;
        /*
         * Get the attributes field, units, from the XML node
         */
        string units = node["units"];
        if (units != "" && convert) {
          if (type == "actEnergy" && units != "") {
            funit = actEnergyToSI(units);
          } else if (type != "" && units != "") {
            funit = toSI(units);
	  }
	}

	if (node["min"] != "") 
	    vmin = atofCheck(node["min"].c_str());
	if (node["max"] != "") 
	    vmax = atofCheck(node["max"].c_str());

	doublereal vv;
        string val = node.value();
        while (1 > 0) {
            icom = val.find(',');
            if (icom != string::npos) {
                numstr = val.substr(0,icom);
                val = val.substr(icom+1,val.size());
		dtmp = atofCheck(numstr.c_str());
                v.push_back(dtmp);
            }
            else {
	      /*
	       * This little bit of code is to allow for the
	       * possibility of a comma being the last 
	       * item in the value text. This was allowed in
	       * previous versions of Cantera, even though it
	       * would appear to be odd. So, we keep the
	       * possibilty in for backwards compatibility.
	       */
	      int nlen = strlen(val.c_str());
	      if (nlen > 0) {
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
        int nv = v.size();
        for (int n = 0; n < nv; n++) {
            v[n] *= funit;
        }
    }

    /**
     * This routine is used to interpret the value portions of XML
     * elements that contain colon separated pairs. These are used,
     * for example, in describing the element composition of species.
     *       <atomArray> H:4 C:1 <atomArray\>
     * The string is first separated into a string vector according
     * to the location of white space. Then each string is again
     * separated into two parts according to the location of a
     * colon in the string. The first part of the string is 
     * used as the key, while the second part of the string is
     * used as the value, in the return map.
     * It is an error to not find a colon in each string pair.
     */
    void getMap(const XML_Node& node, map<string, string>& m) {
        vector<string> v;
        getStringArray(node, v);
        string key, val;
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

  
  //! This function interprets the value portion of an XML element
  //! as a series of "Pairs" separated by white space.
  /*!
   * Each pair consists of nonwhite-space characters.
   * The first ":" found in the pair string is used to separate
   * the string into two parts. The first part is called the "key"
   * The second part is called the "val".
   * String vectors of key[i] and val[i] are returned in the
   * argument list.
   * Warning: No spaces are allowed in each pair. Quotes are part
   *          of the string.
   *   Example: @verbatum
   *    <xmlNode> 
           red:112    blue:34
           green:banana
       </xmlNode>      @endverbatum
   * 
   * Returns:
   *          key       val
   *     0:   "red"     "112"
   *     1:   "blue"    "34"
   *     2:   "green"   "banana"
   */
  void getPairs(const XML_Node& node, vector<string>& key, 
		vector<string>& val) {
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
  }

    /**
     * This function interprets the value portion of an XML element
     * as a series of "Matrix ids and entries" separated by white space.
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
     */
    void getMatrixValues(const XML_Node& node,
			 const vector<string>& keyString1,
			 const vector<string>& keyString2,
			 Array2D &retnValues, bool convert,
                         bool matrixSymmetric) {
	int szKey1 = keyString1.size();
	int szKey2 = keyString2.size();
	int nrow   = retnValues.nRows();
	int ncol   = retnValues.nColumns();
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
	int icol, irow;
        int n = static_cast<int>(v.size());
	string::size_type icolon;
        for (int i = 0; i < n; i++) {
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
	    icol = -1;
	    irow = -1;
	    for (int j = 0; j < szKey1; j++) {
	      if (key1 == keyString1[j]) {
		irow = j;
		break;
	      }
	    }
	    if (irow == -1) {
	      throw CanteraError("getMatrixValues","Row not matched by string: "
				 + key1);		
	    }
	    for (int j = 0; j < szKey2; j++) {
	      if (key2 == keyString2[j]) {
		icol = j;
		break;
	      }
	    }
	    if (icol == -1) {
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

  static string::size_type findFirstWS(const string& val) {
    string::size_type ibegin = string::npos;
    int j = 0;
    std::string::const_iterator i = val.begin();
    for ( ; i != val.end(); i++) {
      char ch = *i;
      int ll = (int) ch;
      if (isspace(ll)) {
	ibegin = (string::size_type) j;
	break;
      }
      j++;
    }
    return ibegin;  
  }

  static string::size_type findFirstNotOfWS(const string& val) {
    string::size_type ibegin = string::npos;
    int j = 0;
    std::string::const_iterator i = val.begin();
    for ( ; i != val.end(); i++) {
      char ch = *i;
      int ll = (int) ch;
      if (!isspace(ll)) {
	ibegin = (string::size_type) j;
	break;
      }
      j++;
    }
    return ibegin;  
  }

    /**
     * This function interprets the value portion of an XML element
     * as a string. It then separates the string up into tokens
     * according to the location of white space.
     * The separate tokens are returned in the string vector,
     * v.
     */
  void getStringArray(const XML_Node& node, vector<string>& v) {
    string val = node.value();
    getStringArray(val, v);
  }

    /**
     * This function interprets the value portion of an XML element
     * as a string. It then separates the string up into tokens
     * according to the location of white space.
     * The separate tokens are returned in the string vector,
     * v.
     */
  void getStringArray(const std::string& oval, vector<string>& v) {
    std::string val(oval);
    string::size_type ibegin, iend;
    v.clear();
    while (1 > 0) {
      ibegin = findFirstNotOfWS(val);
      //val.find_first_not_of(" \n\t");
      if (ibegin != string::npos) {
	val = val.substr(ibegin,val.size());
	//iend = val.find_first_of(" \n\t");
	iend = findFirstWS(val);
	if (iend == string::npos) {
	  v.push_back(val);
	  break;
	} else {
	  v.push_back(val.substr(0,iend));
	  val = val.substr(iend+1,val.size());
	}
      }
      else {
	break;
      }
    }
  }

  void getFunction(const XML_Node& node, string& type, doublereal& xmin,
		   doublereal& xmax, vector_fp& coeffs) {
    const XML_Node& c = node.child("floatArray");
    coeffs.clear();
    getFloatArray(c,coeffs);
    xmin = Undef;
    if (node["min"] != "") xmin = fpValue(node["min"]);
    if (node["max"] != "") xmax = fpValue(node["max"]);
    type = node["type"];
  }
}
