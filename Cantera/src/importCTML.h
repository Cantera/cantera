/**
 *  @file importCTML.h
 *
 * $Author$
 * $Revision$
 * $Date$
 *
 */

// Copyright 2002  California Institute of Technology


#ifndef CT_IMPORTCTML_H
#define CT_IMPORTCTML_H

#include <string>
using namespace std;

#include "ThermoPhase.h"
//#include "Kinetics.h"
//#include "transport/TransportBase.h"

namespace Cantera {

    class Kinetics;
    //class ThermoPhase;
    class XML_Node;

    bool isCTMLFile(string infile);
    XML_Node* get_XML_Node(const string& src, XML_Node* root);
    bool importPhase(XML_Node& phase, ThermoPhase* th);
    bool importKinetics(const XML_Node& phase, vector<ThermoPhase*> th, 
        Kinetics* kin);
    bool installReactionArrays(const XML_Node& parent, Kinetics& kin, 
        string default_phase, bool check_for_duplicates = false);
    ThermoPhase* newPhase(XML_Node& phase);
    bool buildSolutionFromXML(XML_Node& root, string id, string nm, 
        ThermoPhase* th, Kinetics* k);

}
 
#endif

