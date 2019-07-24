/**
 *  @file Elements.h
 *  Contains the getElementWeight function and the definitions of element
 *  constraint types.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ELEMENTS_H
#define CT_ELEMENTS_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

/*!
 * @name Types of Element Constraint Equations
 *
 * There may be several different types of element constraints handled by the
 * equilibrium program and by Cantera in other contexts. These defines are used
 * to assign each constraint to one category.
 * @{
 */

//! An element constraint that is current turned off
#define CT_ELEM_TYPE_TURNEDOFF -1

//! Normal element constraint consisting of positive coefficients for the
//! formula matrix.
/*!
 * All species have positive coefficients within the formula matrix. With this
 * constraint, we may employ various strategies to handle small values of the
 * element number successfully.
 */
#define CT_ELEM_TYPE_ABSPOS 0

//! This refers to conservation of electrons
/*!
 * Electrons may have positive or negative values in the Formula matrix.
 */
#define CT_ELEM_TYPE_ELECTRONCHARGE 1

//! This refers to a charge neutrality of a single phase
/*!
 * Charge neutrality may have positive or negative values in the Formula matrix.
 */
#define CT_ELEM_TYPE_CHARGENEUTRALITY 2

//! Constraint associated with maintaining a fixed lattice stoichiometry in a solid
/*!
 * The constraint may have positive or negative values. The lattice 0 species
 * will have negative values while higher lattices will have positive values
 */
#define CT_ELEM_TYPE_LATTICERATIO 3

//! Constraint associated with maintaining frozen kinetic equilibria in
//! some functional groups within molecules
/*!
 * We seek here to say that some functional groups or ionic states should be
 * treated as if they are separate elements given the time scale of the problem.
 * This will be abs positive constraint. We have not implemented any examples
 * yet. A requirement will be that we must be able to add and subtract these
 * constraints.
 */
#define CT_ELEM_TYPE_KINETICFROZEN 4

//! Constraint associated with the maintenance of a surface phase
/*!
 * We don't have any examples of this yet either. However, surfaces only exist
 * because they are interfaces between bulk layers. If we want to treat surfaces
 * within thermodynamic systems we must come up with a way to constrain their
 * total number.
 */
#define CT_ELEM_TYPE_SURFACECONSTRAINT 5

//! Other constraint equations
/*!
 * currently there are none
 */
#define CT_ELEM_TYPE_OTHERCONSTRAINT 6
//@}

//! Number indicating we don't know the entropy of the element in its most
//! stable state at 298.15 K and 1 bar.
#define ENTROPY298_UNKNOWN -123456789.

//! Get the atomic weight of an element.
/*!
 * Get the atomic weight of an element defined in Cantera by its symbol
 * or by its name. This includes the named isotopes defined in Cantera.
 *
 * @param ename String, name or symbol of the element
 * @return The atomic weight of the element
 * @exception CanteraError if a match for ename is not found or the
 * element has no stable isotopes, and therefore no standard atomic weight
 */
double getElementWeight(const std::string& ename);

//! Get the atomic weight of an element.
/*!
 * Get the atomic weight of an element defined in Cantera by its atomic
 * number. The named isotopes cannot be accessed from this function,
 * since the atomic number of the isotopes is the same as the regular
 * element from which they are derived.
 *
 * @param atomicNumber Integer, atomic number of the element
 * @return The atomic weight of the element
 * @exception IndexError if the passed atomic number less than 1 or
 * larger than the number of elements defined
 * @exception CanteraError if the element has no stable isotopes, and
 * therefore no standard atomic weight
 */
double getElementWeight(int atomicNumber);

//! Get the symbol for an element
/*!
 * Get the symbol for an element defined in Cantera by its name. This
 * includes the named isotopes defined in Cantera.
 *
 * @param ename String, name of the element
 * @return The symbol of the element in a string
 * @exception CanteraError if a match for ename is not found
 */
std::string getElementSymbol(const std::string& ename);

//! Get the symbol for an element
/*!
 * Get the symbol for an element defined in Cantera by its atomic
 * number. The named isotopes cannot be accessed from this function,
 * since the atomic number of the isotopes is the same as the regular
 * element from which they are derived.
 *
 * @param atomicNumber Integer, atomic number of the element
 * @return The symbol of the element in a string
 * @exception IndexError if the passed atomic number less than 1 or
 * larger than the number of elements defined
 */
std::string getElementSymbol(int atomicNumber);

//! Get the name of an element
/*!
 * Get the name of an element defined in Cantera by its symbol. This
 * includes the named isotopes defined in Cantera.
 *
 * @param ename String, symbol for the element
 * @return The name of the element, in a string
 * @exception CanteraError if a match for ename is not found
 */
std::string getElementName(const std::string& ename);

//! Get the name of an element
/*!
 * Get the name of an element defined in Cantera by its atomic
 * number. The named isotopes cannot be accessed from this function,
 * since the atomic number of the isotopes is the same as the regular
 * element from which they are derived.
 *
 * @param atomicNumber Integer, atomic number of the element
 * @return The name of the element, in a string
 * @exception CanteraError IndexError if the passed atomic number less than 1 or
 * larger than the number of elements defined
 */
std::string getElementName(int atomicNumber);

//! Get the atomic number for an element
/*!
 * Get the atomic number of an element defined in Cantera by its symbol
 * or name. This includes the named isotopes included in Cantera.
 *
 *  @param ename String, name or symbol of the element
 *  @return The integer atomic number of the element
 *  @exception CanteraError if a match for ename is not found
 */
int getAtomicNumber(const std::string& ename);

//! Get the number of named elements defined in Cantera.
//! This array excludes named isotopes
int numElementsDefined();

//! Get the number of named isotopes defined in Cantera.
//! This array excludes the named elements
int numIsotopesDefined();

} // namespace

#endif
