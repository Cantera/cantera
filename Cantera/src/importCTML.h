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

namespace Cantera {

    class Kinetics;
    //class ThermoPhase;
    class XML_Node;

    bool isCTMLFile(string infile);

    /**
     * This routine will locate an XML node in either the input
     * XML tree or in another input file specified by the file
     * part of the file_ID string. Searches are based on the
     * ID attribute of the XML element only.
     *
     * @param file_ID This is a concatenation of two strings seperated
     *                by the "#" character. The string before the
     *                pound character is the file name of an xml
     *                file to carry out the search. The string after
     *                the # character is the ID attribute 
     *                of the xml element to search for. 
     *                The string is interpreted as a file string if
     *                no # character is in the string.
     *
     * @param root    If the file string is empty, searches for the
     *                xml element with matching ID attribute are
     *                carried out from this XML node.
     */
    XML_Node* get_XML_Node(const string& src, XML_Node* root);

    /**
     * This routine will locate an XML node in either the input
     * XML tree or in another input file specified by the file
     * part of the file_ID string. Searches are based on the
     * XML element name and the ID attribute of the XML element.
     *
     * @param nameTarget This is the XML element name to look for.
     *                   
     * @param file_ID This is a concatenation of two strings seperated
     *                by the "#" character. The string before the
     *                pound character is the file name of an xml
     *                file to carry out the search. The string after
     *                the # character is the ID attribute 
     *                of the xml element to search for. 
     *                The string is interpreted as a file string if
     *                no # character is in the string.
     *
     * @param root    If the file string is empty, searches for the
     *                xml element with matching ID attribute are
     *                carried out from this XML node.
     */
    XML_Node* get_XML_NameID(const string& nameTarget,
			     const string& file_ID, XML_Node* root);

    bool installSpecies(int k, const XML_Node& s, thermo_t& p, 
		    			SpeciesThermo& spthermo, int rule); 

    bool importPhase(XML_Node& phase, ThermoPhase* th);
    bool importKinetics(const XML_Node& phase, vector<ThermoPhase*> th, 
        Kinetics* kin);
    bool installReactionArrays(const XML_Node& parent, Kinetics& kin, 
        string default_phase, bool check_for_duplicates = false);
    ThermoPhase* newPhase(XML_Node& phase);
    ThermoPhase* newPhase(string file, string id);
    bool buildSolutionFromXML(XML_Node& root, string id, string nm, 
        ThermoPhase* th, Kinetics* k);

}
 
#endif

