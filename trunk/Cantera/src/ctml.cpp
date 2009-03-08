/**
 * @file ctml.cpp
 *
 * Functions to read and write CTML.
 *
 */

/* $Author: hkmoffa $
 * $Revision: 1.19 $
 * $Date: 2006/06/13 17:04:22 $
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

    void addString(XML_Node& node, string title, string val, 
        string type) {
        XML_Node& f = node.addChild("string",val);
        f.addAttribute("title",title);
        if (type != "") f.addAttribute("type",type);        
    }

    XML_Node* getByTitle(XML_Node& node, string title) {
        XML_Node* s = node.findByAttr("title",title);
        if (!s) return 0;
        if (s->parent() == &node) return s;
        else return 0;
    }
        
    string getString(const XML_Node& parent, string name) {
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


    void getStrings(const XML_Node& node, map<string,string>& v) {
        vector<XML_Node*> f;
        node.getChildren("string",f);
        int n = static_cast<int>(f.size());
        string typ, title;
        for (int i = 0; i < n; i++) {
            const XML_Node& fi = *(f[i]);
            title = fi["title"];
            v[title] = fi();
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


    /**
     * Get a floating-point value from a child element.  Returns a
     * double value for the child named 'name' of element 'parent'. If
     * 'type' is supplied and matches a known unit type, unit
     * conversion to SI will be done if the child element has an attribute
     * 'units'.
     */
    doublereal getFloat(const XML_Node& parent, string name, string type) {
        if (!parent.hasChild(name)) 
            throw CanteraError("getFloat (called from XML Node \"" +
			       parent.name() + "\"): ",
			       "no child XML element named " + name);
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
        if (type == "actEnergy" && units != "") 
            fctr = actEnergyToSI(units);
        else if (type != "" && units != "")
            fctr = toSI(units);
        return fctr*x;
    }


    /**
     * Get an integer value from a child element.  Returns an
     * integer value for the child named 'name' of element 'parent'.
     */
    int getInteger(const XML_Node& parent, string name) {
        if (!parent.hasChild(name)) 
            throw CanteraError("getInteger (called from XML Node \"" +
			       parent.name() + "\"): ",
			       "no child XML element named " + name);
        const XML_Node& node = parent.child(name);
        int x, x0, x1;
        string units, vmin, vmax;
        x = atoi(node().c_str());
        x0 = Undefined;
        x1 = Undefined;
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

   /**
     * This function interprets the value portion of an XML element
     * as a series of "Pairs" separated by white space.
     * Each pair consists of nonwhite-space characters.
     * The first ":" found in the pair string is used to separate
     * the string into two parts. The first part is called the "key"
     * The second part is called the "val".
     * String vectors of key[i] and val[i] are returned in the
     * argument list.
     * Warning: No spaces are allowed in each pair. Quotes are part
     *          of the string.
     *   Example
     *    <xmlNode> 
     *        red:112    blue:34
     *        green:banana
     *    </xmlNode>
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


    /**
     * This function interprets the value portion of an XML element
     * as a string. It then separates the string up into tokens
     * according to the location of white space.
     * The separate tokens are returned in the string vector,
     * v.
     */
    void getStringArray(const XML_Node& node, vector<string>& v) {
	string::size_type ibegin, iend;
        v.clear();
        string val = node.value();
        while (1 > 0) {
            ibegin = val.find_first_not_of(" \n\t");
            if (ibegin != string::npos) {
                val = val.substr(ibegin,val.size());
                iend = val.find_first_of(" \n\t");
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
