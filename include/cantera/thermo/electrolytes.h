/**
 *  @file electrolytes.h
 *
 * Header file for a common definitions used in electrolytes
 * thermodynamics.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ELECTROLYTES_H
#define CT_ELECTROLYTES_H
namespace Cantera
{
/**
 * Electrolyte species type
 */
const int cEST_solvent = 0; // Solvent species (neutral)
const int cEST_chargedSpecies = 1; // Charged species (charged)
const int cEST_weakAcidAssociated = 2; // Species which can break
// apart into charged species.
// It may or may not be charged.
// These may or may not be
// be included in the
// species solution vector.
const int cEST_strongAcidAssociated = 3; // Species which always breaks
// apart into charged species.
// It may or may not be charged.
// Normally, these aren't included
// in the speciation vector.
const int cEST_polarNeutral = 4; // Polar neutral species
const int cEST_nonpolarNeutral = 5; // Nonpolar neutral species. These
// usually have activity coefficient
// corrections applied to them to
// account for salting-out effects

}
#endif
