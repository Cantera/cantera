/**
 * @file ctml.h
 *
 * CTML ("Cantera Markup Language") is the variant of XML that Cantera uses
 * to store data. These functions read and write it. 
 * 
 * see also: importCTML, ck2ctml.
 */


/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2002  California Institute of Technology



#ifndef CT_CTML_H
#define CT_CTML_H

#include "ct_defs.h"
#include "xml.h"
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

    void getFloatArray(XML_Node& node, 
        vector_fp& v, bool convert=true);

    void getStringArray(XML_Node& node, vector<string>& v);
    void getMap(XML_Node& node, map<string, string>& m);
    void getPairs(XML_Node& node, vector<string>& key, vector<string>& val);

    void getIntegers(XML_Node& node, map<string,int>& v);
    void getFloats(XML_Node& node, map<string,double>& v, bool convert=true);
    doublereal getFloat(XML_Node& parent, string name, string type="");
    void getStrings(XML_Node& node, map<string,string>& v);
    void getFunction(XML_Node& node, string& type, doublereal& xmin,
        doublereal& xmax, vector_fp& coeffs);
    XML_Node* getByTitle(XML_Node& node, string title);
    void getString(XML_Node& node, string title, string& val, 
        string& type);

    string getString(XML_Node& parent, string name);

    void get_CTML_Tree(XML_Node*, string file);
    void ct2ctml(const char* file);
}

#endif
