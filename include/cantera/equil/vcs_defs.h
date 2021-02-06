/**
 * @file vcs_defs.h
 *    Defines and definitions within the vcs package
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef VCS_DEFS_H
#define VCS_DEFS_H

namespace Cantera
{
/*!
 *      ERROR CODES
 */
//@{
#define VCS_SUCCESS 0
#define VCS_NOMEMORY 1
#define VCS_FAILED_CONVERGENCE -1
#define VCS_SHOULDNT_BE_HERE -2
#define VCS_PUB_BAD -3
#define VCS_THERMO_OUTOFRANGE -4
#define VCS_FAILED_LOOKUP -5
#define VCS_MP_FAIL -6
//@}

/*!
 * @name  Sizes of Phases and Cutoff Mole Numbers
 *
 *      All size parameters are listed here
 * @{
 */

//! Cutoff relative mole fraction value, below which species are deleted from
//! the equilibrium problem.
#ifndef VCS_RELDELETE_SPECIES_CUTOFF
#define VCS_RELDELETE_SPECIES_CUTOFF 1.0e-64
#endif

//! Cutoff relative mole number value, below which species are deleted from the
//! equilibrium problem.
#ifndef VCS_DELETE_MINORSPECIES_CUTOFF
#define VCS_DELETE_MINORSPECIES_CUTOFF 1.0e-140
#endif

//! Relative value of multiphase species mole number for a multiphase species
//! which is small.
#ifndef VCS_SMALL_MULTIPHASE_SPECIES
#define VCS_SMALL_MULTIPHASE_SPECIES 1.0e-25
#endif

//! Cutoff relative moles below which a phase is deleted
//! from the equilibrium problem.
#ifndef VCS_DELETE_PHASE_CUTOFF
#define VCS_DELETE_PHASE_CUTOFF 1.0e-13
#endif

//! Relative mole number of species in a phase that is created We want this to
//! be comfortably larger than the VCS_DELETE_PHASE_CUTOFF value so that the
//! phase can have a chance to survive.
#ifndef VCS_POP_PHASE_MOLENUM
#define VCS_POP_PHASE_MOLENUM 1.0e-11
#endif


//! Cutoff moles below which a phase or species which comprises the bulk of an
//! element's total concentration is deleted.
#ifndef VCS_DELETE_ELEMENTABS_CUTOFF
#define VCS_DELETE_ELEMENTABS_CUTOFF 1.0e-280
#endif

//! Maximum steps in the inner loop
#ifndef VCS_MAXSTEPS
#define VCS_MAXSTEPS 50000
#endif

//@}

//! @name  Species Categories used during the iteration
/*!
 * These defines are valid values for spStatus()
 */
//@{

//! Species is a component which can never be nonzero because of a
//! stoichiometric constraint
/*!
 *  An example of this would be a species that contains Ni. But,
 *  the amount of Ni elements is exactly zero.
 */
#define VCS_SPECIES_COMPONENT_STOICHZERO 3

//! Species is a component which can be nonzero
#define VCS_SPECIES_COMPONENT 2

//! Species is a major species
/*!
 * A major species is either a species in a multicomponent phase with
 * significant concentration or it's a Stoich Phase
 */
#define VCS_SPECIES_MAJOR 1

//! Species is a major species
/*!
 * A major species is either a species in a multicomponent phase with
 * significant concentration or it's a Stoich Phase
 */
#define VCS_SPECIES_MINOR 0

//! Species lies in a multicomponent phase, with a small phase concentration
/*!
 * The species lies in a multicomponent phase that exists. It concentration is
 * currently very low, necessitating a different method of calculation.
 */
#define VCS_SPECIES_SMALLMS -1

//! Species lies in a multicomponent phase with concentration zero
/*!
 * The species lies in a multicomponent phase which currently doesn't exist.
 * It concentration is currently zero.
 */
#define VCS_SPECIES_ZEROEDMS -2

//! Species is a SS phase, that is currently zeroed out.
/*!
 * The species lies in a single-species phase which is currently zeroed out.
 */
#define VCS_SPECIES_ZEROEDSS -3

//! Species has such a small mole fraction it is deleted even though its
//! phase may possibly exist.
/*!
 * The species is believed to have such a small mole fraction that it best to
 * throw the calculation of it out. It will be added back in at the end of the
 * calculation.
 */
#define VCS_SPECIES_DELETED -4

//! Species refers to an electron in the metal.
/*!
 * The unknown is equal to the electric potential of the phase in which it
 * exists.
 */
#define VCS_SPECIES_INTERFACIALVOLTAGE -5

//! Species lies in a multicomponent phase that is zeroed atm
/*!
 * The species lies in a multicomponent phase that is currently deleted and will
 * stay deleted due to a choice from a higher level. These species will formally
 * always have zero mole numbers in the solution vector.
 */
#define VCS_SPECIES_ZEROEDPHASE -6

//! Species lies in a multicomponent phase that is active, but species
//! concentration is zero
/*!
 * The species lies in a multicomponent phase which currently does exist. It
 * concentration is currently identically zero, though the phase exists. Note,
 * this is a temporary condition that exists at the start of an equilibrium
 * problem. The species is soon "birthed" or "deleted".
 */
#define VCS_SPECIES_ACTIVEBUTZERO -7

//! Species lies in a multicomponent phase that is active,
//! but species concentration is zero due to stoich constraint
/*!
 * The species lies in a multicomponent phase which currently does exist.  Its
 * concentration is currently identically zero, though the phase exists. This is
 * a permanent condition due to stoich constraints.
 *
 * An example of this would be a species that contains Ni. But, the amount of Ni
 * elements in the current problem statement is exactly zero.
 */
#define VCS_SPECIES_STOICHZERO -8

//@}

//! @name  Phase Categories used during the iteration
/*!
 * These defines are valid values for the phase existence flag
 */
//@{
//! Always exists because it contains inerts which can't exist in any other phase
#define VCS_PHASE_EXIST_ALWAYS 3

//! Phase is a normal phase that currently exists
#define VCS_PHASE_EXIST_YES 2

//! Phase is a normal phase that exists in a small concentration
/*!
 * Concentration is so small that it must be calculated using an alternate
 * method
 */
#define VCS_PHASE_EXIST_MINORCONC 1

//! Phase doesn't currently exist in the mixture
#define VCS_PHASE_EXIST_NO 0

//! Phase currently is zeroed due to a programmatic issue
/*!
 * We zero phases because we want to follow phase stability boundaries.
 */
#define VCS_PHASE_EXIST_ZEROEDPHASE -6

//@}

/*!
 * @name Types of Element Constraint Equations
 *
 * There may be several different types of element constraints handled by the
 * equilibrium program.  These defines are used to assign each constraint to one
 * category.
 *   @{
 */


//! An element constraint that is current turned off
#define VCS_ELEM_TYPE_TURNEDOFF -1

//! Normal element constraint consisting of positive coefficients for the
//! formula matrix.
/*!
 * All species have positive coefficients within the formula matrix. With this
 * constraint, we may employ various strategies to handle small values of the
 * element number successfully.
 */
#define VCS_ELEM_TYPE_ABSPOS 0

//! This refers to conservation of electrons
/*!
 * Electrons may have positive or negative values in the Formula matrix.
 */
#define VCS_ELEM_TYPE_ELECTRONCHARGE 1

//! This refers to a charge neutrality of a single phase
/*!
 * Charge neutrality may have positive or negative values in the Formula matrix.
 */
#define VCS_ELEM_TYPE_CHARGENEUTRALITY 2

//! Constraint associated with maintaining a fixed lattice stoichiometry in the
//! solids
/*!
 * The constraint may have positive or negative values. The lattice 0 species
 * will have negative values while higher lattices will have positive values
 */
#define VCS_ELEM_TYPE_LATTICERATIO 3

//! Constraint associated with maintaining frozen kinetic equilibria in
//! some functional groups within molecules
/*!
 * We seek here to say that some functional groups or ionic states should be
 * treated as if they are separate elements given the time scale of the problem.
 * This will be abs positive constraint. We have not implemented any examples
 * yet. A requirement will be that we must be able to add and subtract these
 * constraints.
 */
#define VCS_ELEM_TYPE_KINETICFROZEN 4

//! Constraint associated with the maintenance of a surface phase
/*!
 * We don't have any examples of this yet either. However, surfaces only exist
 * because they are interfaces between bulk layers. If we want to treat surfaces
 * within thermodynamic systems we must come up with a way to constrain their
 * total number.
 */
#define VCS_ELEM_TYPE_SURFACECONSTRAINT 5
//! Other constraint equations
/*!
 * currently there are none
 */
#define VCS_ELEM_TYPE_OTHERCONSTRAINT 6
//@}

/*!
 * @name  Types of Species Unknowns in the problem
 * @{
 */
//! Unknown refers to mole number of a single species
#define VCS_SPECIES_TYPE_MOLNUM 0

//! Unknown refers to the voltage level of a phase
/*!
 * Typically, these species are electrons in metals. There is an infinite supply
 * of them. However, their electrical potential is sometimes allowed to vary,
 * for example if the open circuit voltage is sought after.
 */
#define VCS_SPECIES_TYPE_INTERFACIALVOLTAGE -5
//@}

/*!
 * @name  Types of State Calculations within VCS. These values determine where
 *        the results are stored within the VCS_SOLVE object.
 * @{
 */
//! State Calculation is currently in an unknown state
#define VCS_STATECALC_UNKNOWN -1
//! State Calculation based on the old or base mole numbers
#define VCS_STATECALC_OLD 0

//! State Calculation based on the new or tentative mole numbers
#define VCS_STATECALC_NEW 1

//! State Calculation based on tentative mole numbers for a phase which is
//! currently zeroed, but is being evaluated for whether it should pop back into
//! existence
#define VCS_STATECALC_PHASESTABILITY 2

//! State Calculation based on a temporary set of mole numbers
#define VCS_STATECALC_TMP 3
//@}

}

// namespace alias for backward compatibility
namespace VCSnonideal = Cantera;

#endif
