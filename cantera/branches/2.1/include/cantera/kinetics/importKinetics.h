/**
 *  @file importKinetics.h
 *   Definitions of global routines for the importing
 *   of data from XML files (see \ref inputfiles).
 *
 *     This file contains routines which are global routines, i.e.,
 *     not part of any object. These routine take as input, ctml
 *     pointers to data, and pointers to %Cantera objects. The purpose
 *     of these routines is to initialize the %Cantera objects with data
 *     from the ctml tree structures.
 */
// Copyright 2002  California Institute of Technology


#ifndef CT_IMPORTCTML_H
#define CT_IMPORTCTML_H

#include "cantera/thermo/ThermoPhase.h"
#include "Kinetics.h"

namespace Cantera
{

class Kinetics;
class SpeciesThermoFactory;
class XML_Node;

//! Rules for parsing and installing reactions
struct ReactionRules {
    ReactionRules();
    bool skipUndeclaredSpecies;
    bool skipUndeclaredThirdBodies;
    bool allowNegativeA;
};

//!This function returns a ratio if two reactions are duplicates of
//!one another, and 0.0 otherwise.
/*!
 * The input arguments are two maps from species number to stoichiometric
 * coefficient, one for each reaction. The reactions are considered duplicates
 * if their stoichiometric coefficients have the same ratio for all species.
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
 *
 * @ingroup kineticsmgr
 */
void checkRxnElementBalance(Kinetics& kin,
                            const ReactionData& rdata,
                            doublereal errorTolerance = 1.0e-3);

/**
 * Get the reactants or products of a reaction. The information is returned in
 * the spnum, stoich, and order vectors. The length of the vectors is the
 * number of different types of reactants or products found for the reaction.
 *
 *  @param[in] rxn xml node pointing to the reaction element in the xml tree.
 *  @param[in] kin Reference to the kinetics object to install the information
 *                 into.
 * @param[in] rp 1 -> Go get the reactants for a reaction; -1 -> Go get the
 *               products for a reaction
 * @param[in] default_phase Name for the default phase to loop up species in.
 * @param[out] spnum vector of species numbers found. Length is number of
 *                   reactants or products.
 * @param[out] stoich stoichiometric coefficient of the reactant or product.
 *                    Length is number of reactants or products.
 * @param[out] order Order of the reactant and product in the reaction rate
 *                   expression.
 * @param[in] rules If rules.skipUndeclaredSpecies is set and we fail to find
 *                  a species we simply return false, allowing the calling
 *                  routine to skip this reaction and continue. Otherwise, we
 *                  will throw an error.
 */
bool getReagents(const XML_Node& rxn, Kinetics& kin, int rp,
                 std::string default_phase,
                 std::vector<size_t>& spnum, vector_fp& stoich,
                 vector_fp& order, const ReactionRules& rules);

//! Read the rate coefficient data from the XML file.
/*!
 *  Extract the rate coefficient for a reaction from the xml node, kf.
 *  kf should point to a XML element named "rateCoeff".
 *  rdata is the partially filled ReactionData object for the reaction.
 *  This function will fill in more fields in the ReactionData object.
 *
 *  @param kf      XML_Node containing information about the rate coefficients.
 *  @param kin     kinetics manager
 *  @param rdata   ReactionData reference
 *  @param rules   Rules for parsing and installing reactions
 *
 *   Trigger an exception for negative A unless specifically authorized.
 *
 * @ingroup kineticsmgr
 */
void getRateCoefficient(const XML_Node& kf, Kinetics& kin, ReactionData& rdata,
                        const ReactionRules& rules);

//!  Install information about reactions into the kinetics object, kin.
/*!
 *  At this point, parent usually refers to the phase xml element.
 *  One of the children of this element is reactionArray,
 *  the element which determines where in the xml file to
 *  look up the reaction rate data.
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
 * This routine will import a reaction mechanism into a kinetics object. The
 * reaction mechanism may either be homogeneous or heterogeneous, involving
 * multiple ThermoPhase objects. The hosting phase should be included as the
 * first argument. For example, if phase I is an interface phase between bulk
 * phases A and B. Then, the XML_Node for phase I should be the first
 * argument. The vector of ThermoPhase objects should consist of pointers to
 * phases I, A, and B.
 *
 * @param phase This is an xml node containing a description of the owning
 *              phase for the kinetics object. Within the phase is a XML
 *              element called reactionArray containing the location of the
 *              description of the reactions that make up the kinetics object.
 *              Also within the phase is an XML element called phaseArray
 *              containing a listing of other phases that participate in the
 *              kinetics mechanism.
 *
 * @param th    This is a list of ThermoPhase pointers which must include all
 *              of the phases that participate in the kinetics operator. All
 *              of the phases must have already been initialized and formed
 *              within Cantera. However, their pointers should not have been
 *              added to the Kinetics object; this addition is carried out
 *              here. Additional phases may be include in the list; these have
 *              no effect.
 *
 * @param kin   This is a pointer to a kinetics manager class that will be
 *              initialized with the kinetics mechanism. Inherited Kinetics
 *              classes may be used here.
 *
 * @ingroup kineticsmgr
 *
 */
bool importKinetics(const XML_Node& phase, std::vector<ThermoPhase*> th,
                    Kinetics* kin);

//!Build a single-phase ThermoPhase object with associated kinetics mechanism.
/*!
 *  In a single call, this routine initializes a ThermoPhase object and a
 *  homogeneous kinetics object for a phase.
 *
 * @param root pointer to the XML tree which will be searched to find the
 *             XML phase element.
 *
 * @param id   Name of the phase to be searched for.
 * @param nm   Name of the XML element. Should be "phase"
 * @param th   Pointer to a bare ThermoPhase object, which will be initialized
 *             by this operation.
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
bool buildSolutionFromXML(XML_Node& root, const std::string& id,
                          const std::string& nm, ThermoPhase* th, Kinetics* k);

//! Search an XML tree for species data.
/*!
 *   This utility routine will search the XML tree for the species named by
 *   the string, kname. It will return the XML_Node pointer. Failures of any
 *   kind return the null pointer.
 *
 * @param kname species Name
 * @param phaseSpeciesData Pointer to the phase XML node pertaining to the
 *                         species database for the phase to be found
 *
 * @return
 *    Returns a pointer to the XML node containing the species data.
 *
 * @ingroup inputfiles
 */
//const XML_Node *speciesXML_Node(std::string kname,
//                                const XML_Node *phaseSpeciesData);

}

#endif
