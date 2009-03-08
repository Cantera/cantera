/**
 * @file ctml.h
 *
 * CTML ("Cantera Markup Language") is the variant of XML that Cantera uses
 * to store data. These functions read and write it. 
 * 
 * see also: importCTML, ck2ctml.
 */


/* $Author: hkmoffa $
 * $Revision: 1.9 $
 * $Date: 2006/06/13 17:04:21 $
 */

// Copyright 2002  California Institute of Technology



#ifndef CT_CTML_H
#define CT_CTML_H

#include "ct_defs.h"
#include "xml.h"
#include "Array.h"
using namespace Cantera;

namespace ctml {

    const string CTML_Version = "1.4.1";

    bool isBuiltin(string nm);

    void addBool(XML_Node& node, 
        string title, 
        bool val);

    void addInteger(XML_Node& node, 
        string title, 
        int val, 
        string units="", 
        string type="");

    void addFloat(XML_Node& node, 
        string title, 
        doublereal val, 
        string units="", 
        string type="", 
        doublereal minval = Undef,
        doublereal maxval = Undef);

    void addIntegerArray(XML_Node& node, 
        string title, 
        int n, 
        const int* vals, 
        string units="", 
        string type="",
        doublereal minval=Undef, 
        doublereal maxval=Undef);

    void addFloatArray(XML_Node& node, 
        string title, 
        int n, 
        const double* vals, 
        string units="", 
        string type="",
        doublereal minval = Undef,
        doublereal maxval = Undef);

    void addString(XML_Node& node, 
        string title, 
        string val, 
        string type="");

    void getFloatArray(const XML_Node& node, vector_fp& v, 
                       bool convert=true, string type="",
		       string nodeName = "floatArray");

    void getStringArray(const XML_Node& node, vector<string>& v);
    void getMap(const XML_Node& node, map<string, string>& m);
    void getPairs(const XML_Node& node, vector<string>& key, 
		  vector<string>& val);
    void getMatrixValues(const XML_Node& node, 
			 const vector<string>& keyString1,
                         const vector<string>& keyString2,
                         Array2D &returnValues, bool convert = true,
                         bool matrixSymmetric = false);

    void getIntegers(const XML_Node& node, map<string,int>& v);
    void getFloats(const XML_Node& node, map<string,double>& v,
		   bool convert=true);
    doublereal getFloat(const XML_Node& parent, string name, string type="");
    int getInteger(const XML_Node& parent, string name);
    
    void getStrings(const XML_Node& node, map<string,string>& v);
    void getFunction(const XML_Node& node, string& type, doublereal& xmin,
        doublereal& xmax, vector_fp& coeffs);
    XML_Node* getByTitle(XML_Node& node, string title);
    void getString(XML_Node& node, string title, string& val, 
        string& type);

    string getString(const XML_Node& parent, string name);

    // these are defined in ct2ctml.cpp
    void get_CTML_Tree(XML_Node* node, string file);
    void ct2ctml(const char* file);
}

#endif
