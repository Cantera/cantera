/**
 *  @file reaction_defs.h
 * This file defines some constants used to specify reaction types.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RXN_DEFS_H
#define CT_RXN_DEFS_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

const int INVALID_RXN = 0;
const int NONE = 0;

/// @name Reaction Types

//@{

//! A reaction with a rate coefficient that depends only on temperature and voltage
//! that also obeys mass-action kinetics.
/*!
 *  Here mass-action kinetics is defined as the reaction orders being equal to
 *  the reaction's stoichiometry.
 *
 * temperature. Example: O + OH <-> O2 + H
 */
const int ELEMENTARY_RXN = 1;

/**
 * A gas-phase reaction that requires a third-body collision partner. Example:
 * O2 + M <-> O + O + M
 */
const int THREE_BODY_RXN = 2;

/**
 * The general form for a gas-phase association or dissociation reaction, with a
 * pressure-dependent rate. Example: CH3 + H (+M) <-> CH4 (+M)
 */
const int FALLOFF_RXN = 4;

/**
 * A pressure-dependent rate expression consisting of several Arrhenius rate
 * expressions evaluated at different pressures. The final rate is calculated
 * by logarithmically interpolating between the two rates that bracket the
 * current pressure.
 */
const int PLOG_RXN = 5;

/**
 * A general gas-phase pressure-dependent reaction where k(T,P) is defined in
 * terms of a bivariate Chebyshev polynomial.
 */
const int CHEBYSHEV_RXN = 6;

/**
 * A chemical activation reaction. For these reactions, the rate falls
 * off as the pressure increases, due to collisional stabilization of
 * a reaction intermediate. Example: Si + SiH4 (+M) <-> Si2H2 + H2
 * (+M), which competes with Si + SiH4 (+M) <-> Si2H4 (+M).
 */
const int CHEMACT_RXN = 8;

/**
 * A reaction occurring on a surface.
 *  NOTE: This is a bit ambiguous, and will be taken out in the future
 *        The dimensionality of the interface is a separate concept from the type
 *        of the reaction.
 */
const int SURFACE_RXN = 20;

//! A reaction occurring on an interface, e.g a surface or edge.
const int INTERFACE_RXN = 20;

//!  This is a surface reaction that is formulated using the Butler-Volmer
//!  formulation and using concentrations instead of activity concentrations
//!  for its exchange current density formula.
//!  @deprecated To be removed after Cantera 2.5.
const int BUTLERVOLMER_NOACTIVITYCOEFFS_RXN = 25;

//!  This is a surface reaction that is formulated using the Butler-Volmer
//!  formulation. Note the B-V equations can be derived from the forward
//!  and reverse rate constants for a single step reaction. However, there
//!  are some advantages to using the formulation directly.
//!  @deprecated To be removed after Cantera 2.5.
const int BUTLERVOLMER_RXN = 26;

//!  This is a surface reaction that is formulated using the affinity
//!  representation, common in the geochemistry community.
//!  This is generally a global non-mass action reaction with an additional functional
//!  form dependence on delta G of reaction.
//!  @deprecated To be removed after Cantera 2.5.
const int SURFACEAFFINITY_RXN = 27;

/**
 * A reaction occurring at a one-dimensional interface between two surface phases.
 *  NOTE: This is a bit ambiguous, and will be taken out in the future
 *        The dimensionality of the interface is a separate concept from the type
 *        of the reaction.
 */
const int EDGE_RXN = 22;

/**
 * A global reaction. These may have non-mass action reaction orders,
 * and are not allowed to be reversible.
 * @deprecated To be removed after Cantera 2.5.
 */
const int GLOBAL_RXN = 30;

//@}

/** @name Rate Coefficient Types
 * These types define the supported rate coefficient types for
 * elementary reactions. Any of these may also be used as the high and
 * low-pressure limits of falloff and chemical activation reactions.
 *
 * Note that not all of these are currently implemented!
 * @todo Finish implementing reaction rate types.
 *
 * Magic numbers
 * @deprecated To be removed after Cantera 2.5.
 */
//@{

const int ARRHENIUS_REACTION_RATECOEFF_TYPE = 1;
const int LANDAUTELLER_REACTION_RATECOEFF_TYPE = 2;
const int TSTRATE_REACTION_RATECOEFF_TYPE = 3;
const int SURF_ARRHENIUS_REACTION_RATECOEFF_TYPE = 4;
const int ARRHENIUS_SUM_REACTION_RATECOEFF_TYPE = 5;
const int EXCHANGE_CURRENT_REACTION_RATECOEFF_TYPE = 6;
const int PLOG_REACTION_RATECOEFF_TYPE = 7;
const int CHEBYSHEV_REACTION_RATECOEFF_TYPE = 8;

//@}

/** @name Falloff Function Types
 *
 * Magic numbers
 * @deprecated To be removed after Cantera 2.5.
 */
//@{
const int SIMPLE_FALLOFF = 100;
const int TROE_FALLOFF = 110;
const int SRI_FALLOFF = 112;
//@}
}

#endif
