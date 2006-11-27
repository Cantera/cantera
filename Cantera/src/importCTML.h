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
//using namespace std;

#include "ThermoPhase.h"
#include "Kinetics.h"

namespace Cantera {

    class Kinetics;
    class SpeciesThermoFactory;
    //class ThermoPhase;
    class XML_Node;

    bool isCTMLFile(std::string infile);

    XML_Node* get_XML_Node(const std::string& src, XML_Node* root);

    XML_Node* get_XML_NameID(const std::string& nameTarget,
			     const std::string& file_ID, XML_Node* root);

    bool installSpecies(int k, const XML_Node& s, thermo_t& p, 
			SpeciesThermo& spthermo, int rule,
			SpeciesThermoFactory* spfactory);

    bool importPhase(XML_Node& phase, ThermoPhase* th, 
        SpeciesThermoFactory* spfactory = 0);

    doublereal isDuplicateReaction(std::map<int, doublereal>& r1,
                                   std::map<int, doublereal>& r2);

    bool importKinetics(const XML_Node& phase, std::vector<ThermoPhase*> th, 
        Kinetics* kin);
    bool installReactionArrays(const XML_Node& parent, Kinetics& kin, 
        std::string default_phase, bool check_for_duplicates = false);

    /***
     * This function will check a specific reaction to see if it the
     * elements balance.
     */
    void checkRxnElementBalance(Kinetics& kin, 
        const ReactionData &rdata, doublereal errorTolerance = 1.0e-3);

    void getRateCoefficient(const XML_Node& kf, kinetics_t& kin, 
			    ReactionData& rdata, int negA);

    ThermoPhase* newPhase(XML_Node& phase);
    ThermoPhase* newPhase(std::string file, std::string id);
    bool buildSolutionFromXML(XML_Node& root, std::string id, std::string nm, 
        ThermoPhase* th, Kinetics* k);

    const XML_Node *speciesXML_Node(std::string kname,
				    const XML_Node *phaseSpecies);

}
 
#endif

