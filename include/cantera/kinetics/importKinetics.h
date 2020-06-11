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

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IMPORTKINETICS_H
#define CT_IMPORTKINETICS_H

#include "Kinetics.h"

namespace Cantera
{

//!  Install information about reactions into the kinetics object, kin.
/*!
 *  At this point, parent usually refers to the phase XML element. One of the
 *  children of this element is reactionArray, the element which determines
 *  where in the XML file to look up the reaction rate data.
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
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
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
 * @param phase This is an XML node containing a description of the owning
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
 *              classes should be used here.
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
 *
 * @ingroup kineticsmgr
 */
bool importKinetics(const XML_Node& phase, std::vector<ThermoPhase*> th,
                    Kinetics* kin);

//! Build a single-phase ThermoPhase object with associated kinetics mechanism.
/*!
 *  In a single call, this routine initializes a ThermoPhase object and a
 *  homogeneous kinetics object for a phase. It returns the fully initialized
 *  ThermoPhase object pointer and kinetics pointer.
 *
 * @param root pointer to the XML tree which will be searched to find the
 *             XML phase element.
 * @param id   Name of the phase to be searched for.
 * @param nm   Name of the XML element. Should be "phase"
 * @param th   Pointer to a bare ThermoPhase object, which will be initialized
 *             by this operation.
 * @param kin  Pointer to a bare Kinetics object, which will be initialized
 *             by this operation to a homogeneous kinetics manager
 * @return
 *    Returns true if all went well. If there are errors, it will return false.
 *
 * For Example
 *
 * @code
 *        ThermoPhase *th   = new ThermoPhase();
 *        Kinetics    *kin  = new Kinetics();
 *        XML_Node *root = get_XML_File("gri30.xml");
 *        ok =  buildSolutionFromXML(root, "gri30_mix", "phase", th, kin)
 * @endcode
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
 *
 * @ingroup inputfiles
 * @see importKinetics()
 */
bool buildSolutionFromXML(XML_Node& root, const std::string& id,
                          const std::string& nm, ThermoPhase* th, Kinetics* kin);

//! Check to ensure that all electrochemical reactions are specified correctly
/*!
 *  This function ensures the user has correctly specified all electrochemical
 *  reactions. The routine counts the amount of charge (i.e. number of electron
 *  elements specified for each species in each phase) for both reactants and
 *  products. If net charge transfer phases during a reaction, the reaction is
 *  electrochemical. If not already specified as such, the function defines the
 *  reaction as electrochemical, corrects the reaction attributes, and sets
 *  beta = 0.5.
 *
 * @param p     This is an XML node containing a description of the owning
 *              phase for the kinetics object.
 * @param kin   This is a pointer to a kinetics manager class.
 * @param r     This is the reaction node that is being evaluated
 * @return      The function always returns true.
 *
 * @deprecated The XML input format is deprecated and will be removed in
 *     Cantera 3.0.
 */
bool checkElectrochemReaction(const XML_Node& p, Kinetics& kin, const XML_Node& r);


}

#endif
