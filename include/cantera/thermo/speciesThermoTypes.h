/**
 *  @file speciesThermoTypes.h Contains const definitions for types of species
 *       reference-state thermodynamics managers (see \ref spthermo)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef SPECIES_THERMO_TYPES_H
#define SPECIES_THERMO_TYPES_H

//! Constant Cp
#define CONSTANT_CP 1

//! Polynomial
#define POLYNOMIAL_4 2

//! Two regions of 7 coefficient NASA Polynomials
//! This is implemented in the class NasaPoly2 in NasaPoly2.h
#define NASA 4

//! Two regions of 7 coefficient NASA Polynomials
//! This is implemented in the class NasaPoly2 in NasaPoly2.h
#define NASA2 4

//! Two regions of Shomate Polynomials.
#define SHOMATE 8

//! Two regions of Shomate Polynomials.
#define SHOMATE2 8

//! Tiger Polynomials. Not implemented here.
#define TIGER 16

//! Constant Cp thermo.
//! This is implemented in ConstCpPoly in constCpPoly.h for one species.
//! If the whole phase is constcp, SimpleThermo in SimpleThermo.h
//! implements this for the whole phase.
#define SIMPLE 32

//! piecewise interpolation of mu0.
//! This is implemented in Mu0Poly in Mu0Poly.h
#define MU0_INTERP 64

//! one region of Shomate Polynomials used in NIST database
//! This is implemented in the NIST database.
//! This is implemented in ShomatePoly in ShomatePoly.h
#define SHOMATE1 128

//! 7 coefficient NASA Polynomials
//! This is implemented in the class NasaPoly1 in NasaPoly1.h
#define NASA1 256

//! 9 coefficient NASA Polynomials
//! This is implemented in the class Nasa9Poly1 in Nasa9Poly1.h
#define NASA9 512

//! 9 coefficient NASA Polynomials in multiple temperature regions
//! This is implemented in the class Nasa9PolyMultiTempRegion in Nasa9Poly1MultiTempRegion
#define NASA9MULTITEMP 513

#endif
