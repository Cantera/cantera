#ifndef CT_CK2CTML_H
#define CT_CK2CTML_H

#include <iostream>
#include <string>

#include "ctml.h"

namespace ckr{
  class CKReader;
}

namespace ctml {


//     void addNASA(XML_Node& node,  
//         const vector_fp& low, const vector_fp& high,  
//         doublereal minx=-999.0, doublereal midx=-999.0, 
//         doublereal maxx=-999.0);

//     void addShomate(XML_Node& node, 
//         const vector_fp& low, const vector_fp& high,  
//         doublereal minx=-999.0, doublereal midx=-999.0, 
//         doublereal maxx=-999.0);

//     void addArrhenius(XML_Node& node,  
//         doublereal A, doublereal b, doublereal E, int order, 
//         string unitsys, string E_units);

//     void addRateCoeff(XML_Node& node, string title, string type,  
//         string direction, int order, doublereal A, doublereal b, 
//         doublereal E, string unitsys, string E_units);

//     void addTroeFalloff(XML_Node& node, const vector_fp& params);

//     void addSRIFalloff(XML_Node& node, const vector_fp& params);

//     void addElement(XML_Node& node, string idtag, const ckr::Element& el);
    
//     void addSpecies(XML_Node& node, string idtag, const ckr::Species& sp);

//     void addReaction(XML_Node& node, string idtag, int i, 
//         const ckr::Reaction& rxn, const ckr::ReactionUnits& runits,
//         doublereal version);

    //    void addTransport(istream& s, XML_Node& node);
 
    void ck2ctml(string idtag, ckr::CKReader& r,
                 Cantera::XML_Node& root);

    int convert_ck(const char * const in_file, const char * const db_file,
        const char * const tr_file, const char * const out_file, const char * const id_tag);

}

#endif

