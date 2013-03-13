/**
 *  @file electrolytes.h
 *
 * Header file for a common definitions used in electrolytes
 * thermodynamics.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_ELECTROLYTES_H
#define CT_ELECTROLYTES_H
namespace Cantera
{
/**
 * Electrolyte species type
 */
const int cEST_solvent  = 0;             // Solvent species (neutral)
const int cEST_chargedSpecies = 1;       // Charged species (charged)
const int cEST_weakAcidAssociated = 2;   // Species which can break
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
const int cEST_polarNeutral = 4;         // Polar neutral species
const int cEST_nonpolarNeutral = 5;      // Nonpolar neutral species. These
// usually have activity coefficnt
// corrections applied to them to
// account for salting-out effects

/**
 *  eosTypes returned for this ThermoPhase Object
 */
const int cHMWSoln0 = 45010;
const int cHMWSoln1 = 45011;
const int cHMWSoln2 = 45012;

/**
 *  eosTypes returned for this ThermoPhase Object
 */
const int cDebyeHuckel0 = 46010;
const int cDebyeHuckel1 = 46011;
const int cDebyeHuckel2 = 46012;
}
#endif
