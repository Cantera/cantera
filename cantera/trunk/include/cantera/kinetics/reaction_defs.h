/**
 *  @file reaction_defs.h
 * This file defines some constants used to specify reaction types.
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_RXN_DEFS_H
#define CT_RXN_DEFS_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

const int NONE = 0;

/// @name Reaction Types

//@{

/**
 * A reaction with a rate coefficient that depends only on
 * temperature. Example: O + OH <-> O2 + H
 */
const int ELEMENTARY_RXN = 1;

/**
 * A reaction that requires a third-body collision partner. Example:
 * O2 + M <-> O + O + M
 */
const int THREE_BODY_RXN = 2;

/**
 * The general form for an association or dissociation reaction, with a
 * pressure-dependent rate. Example: CH3 + H (+M) <-> CH4 (+M)
 */
const int FALLOFF_RXN    = 4;

/**
 * A pressure-dependent rate expression consisting of several Arrhenius rate
 * expressions evaluated at different pressures. The final rate is calculated
 * by logarithmically interpolating between the two rates that bracket the
 * current pressure.
 */
const int PLOG_RXN = 5;

/**
 * A general pressure-dependent reaction where k(T,P) is defined in terms of
 * a bivariate Chebyshev polynomial.
 */
const int CHEBYSHEV_RXN = 6;

/**
 * A chemical activation reaction. For these reactions, the rate falls
 * off as the pressure increases, due to collisional stabilization of
 * a reaction intermediate. Example: Si + SiH4 (+M) <-> Si2H2 + H2
 * (+M), which competes with Si + SiH4 (+M) <-> Si2H4 (+M).
 */
const int CHEMACT_RXN    = 8;

/**
 * A reaction occurring on a surface.
 */
const int SURFACE_RXN    = 20;

/**
 * A reaction occurring at a one-dimensional interface between two
 * surface phases.
 */
const int EDGE_RXN  = 22;

/**
 * A global reaction. These may have non-integral reaction orders,
 * and are not allowed to be reversible.
 */
const int GLOBAL_RXN     = 30;

//@}

/** @name Rate Coefficient Types
 * These types define the supported rate coefficient types for
 * elementary reactions. Any of these may also be used as the high and
 * low-pressure limits of falloff and chemical activation reactions.
 *
 * Note that not all of these are currently implemented!
 * @todo Finish implementing reaction rate types.
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
 */
//@{
const int SIMPLE_FALLOFF = 100;
const int TROE3_FALLOFF = 110;
const int TROE4_FALLOFF = 111;
const int SRI3_FALLOFF  = 112;
const int SRI5_FALLOFF  = 113;
const int WF_FALLOFF = 114;
//@}
}

#endif
