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

#include "Kinetics.h"
#include "transport/TransportBase.h"

namespace Cantera {

    bool isCTMLFile(string infile);
    bool importPhase(XML_Node& phase, thermophase_t* th);
    bool importKinetics(XML_Node& phase, vector<thermophase_t*> th, 
        Kinetics* kin);
    bool installReactionArrays(XML_Node& parent, Kinetics& kin, 
        string default_phase, bool check_for_duplicates = false);
    ThermoPhase* newPhase(XML_Node& phase);
    bool buildSolutionFromXML(XML_Node& root, string id, string nm, 
        ThermoPhase* th, Kinetics* k);

}
 
#endif

