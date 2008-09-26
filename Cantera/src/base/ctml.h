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
#include "Array.h"

namespace ctml {

    const std::string CTML_Version = "1.4.1";

    bool isBuiltin(std::string nm);

    void addBool(Cantera::XML_Node& node, 
        std::string title, 
        bool val);

    void addInteger(Cantera::XML_Node& node, 
        std::string title, 
        int val, 
        std::string units="", 
        std::string type="");

    void addFloat(Cantera::XML_Node& node, 
        std::string title, 
        doublereal val, 
        std::string units="", 
        std::string type="", 
        doublereal minval = Cantera::Undef,
        doublereal maxval = Cantera::Undef);

    void addIntegerArray(Cantera::XML_Node& node, 
        std::string title, 
        int n, 
        const int* vals, 
        std::string units="", 
        std::string type="",
        doublereal minval=Cantera::Undef, 
        doublereal maxval=Cantera::Undef);

    void addFloatArray(Cantera::XML_Node& node, 
        std::string title, 
        int n, 
        const double* vals, 
        std::string units="", 
        std::string type="",
        doublereal minval = Cantera::Undef,
        doublereal maxval = Cantera::Undef);

    void addString(Cantera::XML_Node& node, 
        std::string title, 
        std::string val, 
        std::string type="");

    void getFloatArray(const Cantera::XML_Node& node, Cantera::vector_fp& v, 
        bool convert=true, std::string type="",
        std::string nodeName = "floatArray");

    void getStringArray(const Cantera::XML_Node& node, std::vector<std::string>& v);
    void getStringArray(const std::string& val, std::vector<std::string>& v);
    void getMap(const Cantera::XML_Node& node, std::map<std::string, std::string>& m);
    void getPairs(const Cantera::XML_Node& node, std::vector<std::string>& key, 
        std::vector<std::string>& val);
    void getMatrixValues(const Cantera::XML_Node& node, 
        const std::vector<std::string>& keyString1,
        const std::vector<std::string>& keyString2,
        Cantera::Array2D &returnValues, bool convert = true,
        bool matrixSymmetric = false);

    void getIntegers(const Cantera::XML_Node& node, std::map<std::string,int>& v);
    void getFloats(const Cantera::XML_Node& node, std::map<std::string,double>& v,
		   bool convert=true);
    doublereal getFloat(const Cantera::XML_Node& parent, std::string name,
        std::string type="");
    int getInteger(const Cantera::XML_Node& parent, std::string name);
    
    void getStrings(const Cantera::XML_Node& node, std::map<std::string,
        std::string>& v);
    void getFunction(const Cantera::XML_Node& node, std::string& type, 
        doublereal& xmin, doublereal& xmax, Cantera::vector_fp& coeffs);
    Cantera::XML_Node* getByTitle(Cantera::XML_Node& node, std::string title);
    void getString(Cantera::XML_Node& node, std::string title, 
        std::string& val, std::string& type);

    std::string getString(const Cantera::XML_Node& parent, std::string name);

    // these are defined in ct2ctml.cpp
    void get_CTML_Tree(Cantera::XML_Node* node, std::string file, int debug = 0);

    //! Convert a cti file into a ctml file
    /*!
     *  
     *  @ingroup inputfiles
     */
    void ct2ctml(const char* file, int debug = 0);
}

#endif
