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

#include "ThermoPhase.h"
#include "Kinetics.h"

namespace Cantera {

  class Kinetics;
  class SpeciesThermoFactory;
  class XML_Node;

  //! Search for an XML_Node either wiithin an existing XML tree structure, or in another file,
  //! based on the file name or the XML id attribute.
  /*!
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
   *
   * @return
   *      This routine will process the XML file, creating an XML
   *      tree structure. It returns a pointer to the top of the tree.
   *
   *
   * For example,
   * @code
   *
   * XML_Node* xn = get_XML_Node("phase", "gri30.xml#gri30_mix", 0);
   *
   * @endcode
   *		  
   * will search in the file gri30.xml for an XML element of the following form, where
   * the XML element name, phase, is an optional hit:
   * @code
   *      <phase id="gri30_mix>
   *         . . .
   *      </phase>
   * @endcode
   *
   * It will return a pointer to an xml tree for the XML phase element.
   *
   * @ingroup inputfiles
   */
  XML_Node* get_XML_Node(const std::string& file_ID, XML_Node* root);
  
  //! Search for an XML node based on the XML element name, file name, or XML id attribute.
  /**
   * This routine will locate an XML node in either the input
   * XML tree or in another input file specified by the file
   * part of the file_ID string. Searches are based on the
   * XML element name and the ID attribute of the XML element.
   * An exact match of both is usually required. However, the
   * ID attribute may be set to "", in which case the first
   * xml element with the correct XML element name will be returned.
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
   *
   * @return
   *      This routine will process the XML file, possibly creating an XML
   *      tree structure. It returns a pointer to the XML node searched for.
   *
   * For example,
   * @code
   *
   * XML_Node* xn = get_XML_NameID("phase", "gri30.xml#gri30_mix", 0);
   *
   * @endcode
   *		  
   * will search in the file gri30.xml for an XML element of the following form:
   * @code
   *      <phase id="gri30_mix>
   *         . . .
   *      </phase>
   * @endcode
   *
   * It will return a pointer to an xml tree for the XML phase element.
   *
   * @ingroup inputfiles
   */
  XML_Node* get_XML_NameID(const std::string& nameTarget,
			   const std::string& file_ID, XML_Node* root);
  
  //! Install a species into a ThermoPhase object, which defines
  //! the phase thermodynamics and speciation.
  /*!
   *  This routine first gathers the information from the Species XML
   *  tree and calls addUniqueSpecies() to add it to the
   *  ThermoPhase object, p. This information consists of:
   *
   *         ecomp[] = element composition of species.
   *         chgr    = electric charge of species
   *         name    = string name of species
   *         sz      = size of the species 
   *                 (option double used a lot in thermo)
   *
   *  Then, the routine processes the "thermo" XML element and
   *  calls underlying utility routines to read the XML elements
   *  containing the thermodynamic information for the reference
   *  state of the species. Failures or lack of information trigger
   *  an "UnknownSpeciesThermoModel" exception being thrown.
   *
   *  Then the routine calls the factory routine for installing the
   *  species thermo into the species thermo manager. Note, this step
   *  may be customized (i.e., both the factory and the SpeciesThermo
   *  method may be inherited options, so that new standard states 
   *  may be input.
   *
   * @param k        Species index within the phase.
   * @param s        XML_Node containing the information about the species
   * @param p        ThermoPhase object
   * @param spthermo SpeciesThermo Manager
   * @param rule     If rule = 0, throw an error if the species
   *                 contains an undefined element. If rule ne 0,
   *                 just skip the species.
   * @param spfactory pointer to the SpeciesThermoFactory
   *
   * @return
   *      Returns true if the procedure succeeded. Returns false,
   *      if the species was not included in the ThermoPhase. Note
   *      this isn't necessarily a failure, because of the parameter, rule.
   *
   * @ingroup thermoprops
   */
  bool installSpecies(int k, const XML_Node& s, thermo_t& p, 
		      SpeciesThermo& spthermo, int rule,
		      SpeciesThermoFactory* spfactory);

  //! Import a phase information into an empty thermophase object
  /*!
   *   Here we read an XML description of the thermodynamic information
   *   for a phase. At the end of this routine, the phase should
   *   be ready to be used within applications. This routine contains
   *   some key routines that are used as pass back routines so that
   *   the phase (and the contents of the XML file) may contain
   *   variable paramerizations for the specification of the
   *   species standard states, the equation of state, and the
   *   specification of other nonidealities. Below, a description
   *   is presented of the main algorithm for bringing up a %ThermoPhase
   *   object, with care to present points where customizations 
   *   occur.
   *
   *   Before invoking this routine, either the ThermoPhase Factory routines
   *   are called or direct constructor routines are called that
   *   instantiate an inherited ThermoPhase object. This object is input
   *   to this routine, and therefore contains inherited routines that
   *   drive the custimation of the initialization process.
   *    
   *   At the start of the routine, we import descriptions of the elements 
   *   that make up the species in a phase.
   *
   *   We call setParametersFromXML(eos) to read parameters about
   *   the thermo phase before the species are read in.
   *
   *   We call addElementsFromXML() to add elements into the 
   *   description of the phase.
   *
   *   We create a new species thermo manager.  Function
   *   'newSpeciesThermoMgr' looks at the species in the database
   *   to see what thermodynamic property parameterizations are
   *   used, and selects a class that can handle the
   *   parameterizations found.
   *
   *   We import information about the species, including their
   *   reference state thermodynamic polynomials. We then freeze
   *   the state of the species in the element.
   *
   *   Finally, we call initThermoXML(),
   *   a member function of the ThermoPhase object, to "finish"
   *   the description. Now that the species are known, 
   *   additional information may be read in about the thermodynamics
   *   of the phase, (e.g.,  virial coefficients, which are 
   *   binary or ternary interaction parameters between species).
   *
   * @param phase This object must be the phase node of a
   *             complete XML tree
   *             description of the phase, including all of the
   *             species data. In other words while "phase" must
   *             point to an XML phase object, it must have
   *             sibling nodes "speciesData" that describe
   *             the species in the phase.
   * @param th   Pointer to the ThermoPhase object which will
   *             handle the thermodynamics for this phase.
   *             We initialize part of the Thermophase object
   *             here, especially for those objects which are
   *             part of the Cantera Kernel.
   *
   * @param spfactory species Thermo factory pointer, if
   *                  available. If not available, one will be
   *                  created.
   *
   * @ingroup thermoprops
   */
  bool importPhase(XML_Node& phase, ThermoPhase* th, 
		   SpeciesThermoFactory* spfactory = 0);
  
  //!This function returns a ratio if two reactions are duplicates of
  //!one another, and 0.0 otherwise. 
  /*!
   * The input arguments are two
   * maps from species number to stoichiometric coefficient, one for
   * each reaction. The reactions are considered duplicates if their
   * stoichiometric coefficients have the same ratio for all
   * species.
   *
   * @param r1 map 1
   * @param r2 map 2
   *
   * @return
   *    Returns 0.0 if the reactions are not the same.
   *    If the reactions are the same, it returns the ratio of the
   *    stoichiometric coefficients.
   *
   * @ingroup kineticsmgr
   */
  doublereal isDuplicateReaction(std::map<int, doublereal>& r1,
				 std::map<int, doublereal>& r2);
  
  
  //! This function will check a specific reaction to see if the elements balance.
  /*!
   *   @param kin Kinetics object
   *   @param rdata  Object containing the information about one reaction
   *   @param errorTolerance double containing the error tolerance.
   */
  void checkRxnElementBalance(Kinetics& kin, 
			      const ReactionData &rdata, 
			      doublereal errorTolerance = 1.0e-3);

  //! Read the rate coefficient data from the XML file. 
  /*!
   *  Extract the rate coefficient for a reaction from the xml node, kf.
   *  kf should point to a XML element named "rateCoeff".
   *  rdata is the partially filled ReactionData object for the reaction.
   *  This function will fill in more fields in the ReactionData object.
   *
   *  @param kf      XML_Node containing information about the rate coefficients.
   *  @param kin     kinetics manager
   *  @param rdata   ReactionData referece
   *  @param negA    Boolean indicating whether negative A's are ok.
   *
   *   Trigger anexception for negative A unless specifically authorized.
   */
  void getRateCoefficient(const XML_Node& kf, kinetics_t& kin, 
			  ReactionData& rdata, int negA);

  
  //! Create a new ThermoPhase object and initializes it according to the XML tree database.  
  /*!
   *  This routine first looks up the
   * identity of the model for the solution thermodynamics in the
   * model attribute of the thermo child of the xml phase
   * node. Then, it does a string lookup using Cantera's internal ThermoPhase Factory routines
   * on the model to figure out
   * what ThermoPhase derived class should be assigned. It creates a new
   * instance of that class, and then calls importPhase() to
   * populate that class with the correct parameters from the XML
   * tree.
   *
   * @param phase XML_Node reference pointing to the phase XML element.
   *
   * @return
   *    Returns a pointer to the completed and initialized ThermoPhase object. 
   *
   * @ingroup inputfiles
   */
  ThermoPhase* newPhase(XML_Node& phase);

  //! Create a new ThermoPhase object and initializes it according to a file name and phase id  
  /*!
   *  This routine will first find a matching XML phase description given a file name 
   *  and phase id string. Then, this routine will looks up the
   * identity of the model for the solution thermodynamics in the
   * model attribute of the thermo child of the xml phase
   * node. Then, it does a string lookup using Cantera's internal ThermoPhase Factory routines
   * on the model to figure out
   * what ThermoPhase derived class should be assigned. It creates a new
   * instance of that class, and then calls importPhase() to
   * populate that class with the correct parameters from the XML
   * tree.
   *
   * @param file String file name to look up the phase
   * @param id   XML id attribute of the phase.
   *
   * @return
   *    Returns a pointer to the completed and initialized ThermoPhase object. 
   *
   * @ingroup inputfiles
   */
  ThermoPhase* newPhase(std::string file, std::string id);

  //!  Install information about reactions into the kinetics object, kin.
  /*!
   *  At this point, parent usually refers to the phase xml element.
   *  One of the children of this element is reactionArray,
   *  the element which determines where in the xml file to
   *  look up the reaction rate data.
   *
   * This is a wrapper routine around the static function installReaction()
   *
   *  @param p             parent XML phase element
   *  @param kin           Kinetics object to install reactions into
   *  @param default_phase The default_phase is the default phase to assume when
   *                       looking up species.
   *  @param check_for_duplicates Check for reactions with exactly the same
   *                       reactants and products.
   *
   *  @return
   *    On return, if reaction instantiation goes correctly, return true.
   *    If there is a problem, return false.
   *
   * @ingroup kineticsmgr
   */
  bool installReactionArrays(const XML_Node& p, Kinetics& kin, 
			     std::string default_phase, 
			     bool check_for_duplicates = false);
  
  //! Import a reaction mechanism for a phase or an interface.
  /*!
   * This routine will import a reaction mechanism into a
   * kinetics object. The reaction
   * mechanism may either be homogeneous or heterogeneous,
   * involving multiple ThermoPhase objects. 
   * The hosting phase should be included as the first argument.
   * For example, if phase I is an interface phase between bulk
   * phases A and B. Then, the XML_Node for phase I should be
   * the first argument. 
   * The vector of %ThermoPhase objects should be consist of pointers
   * to phases I, A, and B. 
   *
   * @param phase This is an xml node containing a description
   *              of a phase. Within the phase is a XML element
   *              called reactionArray containing the location
   *              of the description of the reactions that make
   *              up the kinetics object. 
   *              Also within the phase is an XML element called
   *              phaseArray containing a listing of other phases
   *              that participate in the kinetics mechanism.
   *
   * @param th    This is a list of ThermoPhase pointers containing
   *              the phases that participate in the kinetics
   *              reactions. All of the phases must have already
   *              been initialized and formed within Cantera.
   *              However, their pointers should not have been
   *              added to the Kinetics object; this addition
   *              is carried out here.
   *
   * @param kin   This is a pointer to a bare kinetics manager class
   *              that will be initialized with the kinetics 
   *              mechanism.
   *
   * @ingroup kineticsmgr
   *
   */
  bool importKinetics(const XML_Node& phase, std::vector<ThermoPhase*> th, 
		      Kinetics* kin);
 
  //!Build a single-phase ThermoPhase object with associated kinetics mechanism.
  /*!
   *  In a single call, this routine initializes a ThermoPhase object and a 
   *  homogenous kinetics object for a phase.
   *
   * @param root pointer to the XML tree which will be searched to find the
   *             XML phase element.
   *
   * @param id   Name of the phase to be searched for.
   * @param nm   Name of the XML element. Should be "phase"
   * @param th   Pointer to a bare ThermoPhase object, which will be initialized
   *             by this operaton.
   * @param k    Pointer to a bare Kinetics object, which will be initialized
   *             by this operation to a homogeneous kinetics manager
   *
   * @return
   *    Returns true if all went well. If there are errors, it will return false.
   *
   * For Example
   *
   * @code
   *        ThermoPhase *th = new ThermoPhase();
   *        Kinetics    *k  = new Kinetics();
   *        XML_Node *root = get_XML_File("gri30.xml");
   *        ok =  buildSolutionFromXML(root, "gri30_mix", "phase", th, k)
   * @endcode
   *
   * @ingroup inputfiles
   * @see importKinetics()
   */
  bool buildSolutionFromXML(XML_Node& root, std::string id, std::string nm, 
			    ThermoPhase* th, Kinetics* k);
  
  //! Search an XML tree for species data.
  /*!
   *
   *   This utility routine will search the XML tree for the species
   *   named by the string, kname. It will return the XML_Node
   *   pointer.
   *   Failures of any kind return the null pointer.
   *
   * @param kname species Name
   * @param phaseSpeciesData Pointer to the phase XML node pertaining to the
   *                         species database for the phase to be found
   *
   * @return
   *    Returns a pointer to teh XML node containing the species data.
   *
   * @ingroup inputfiles
   */
  const XML_Node *speciesXML_Node(std::string kname,
				  const XML_Node *phaseSpeciesData);

}
 
#endif

