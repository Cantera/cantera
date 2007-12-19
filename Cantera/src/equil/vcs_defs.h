/**
 * @file vcs_defs.h
 *    Defines and definitions within the vcs package
 */
/*
 * $Id$
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef VCS_DEFS_H
#define VCS_DEFS_H

namespace VCSnonideal {

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *
 *      COMMON DEFINITIONS -> Protect them against redefinitions
 */

#ifndef TRUE
# define TRUE 1
#endif
   
#ifndef FALSE
# define FALSE 0
#endif

#ifndef BOOLEAN
# define BOOLEAN int
#endif

#ifndef MAX
# define MAX(x,y) (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
# define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif

#ifndef SWAP
# define SWAP(x1, x2, temp) ((temp) = (x1), (x1) = (x2), (x2) = (temp))
#endif

#ifndef SQUARE
# define SQUARE(x) ((x) * (x))
#endif

#ifndef DSIGN
# define DSIGN(x) (( (x) == (0.0) ) ? (0.0) : ( ((x) > 0.0) ? 1.0 : -1.0 ))
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *      ERROR CODES
 *
 */

#define VCS_SUCCESS             0
#define VCS_NOMEMORY            1
#define VCS_FAILED_CONVERGENCE -1
#define VCS_SHOULDNT_BE_HERE   -2
#define VCS_PUB_BAD            -3
#define VCS_THERMO_OUTOFRANGE  -4
#define VCS_FAILED_LOOKUP      -5
#define VCS_MP_FAIL            -6

#define VCS_ERR_NOVALUE       -1.1234E15

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *                   Problem Types
 */
#define VCS_PROBTYPE_TP      0
#define VCS_PROBTYPE_TV      1

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Maximum Length of any name in this package
 */
#define VCS_MAX_NAME_LEN 31

#define VCS_MAX_NAME_LEN_P1 32

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *               SIZES OF PHASES AND MOLE NUMBER CUTOFFS
 *
 *     VCS_DELETE_SPECIES_CUTOFF: Cutoff relative mole number value,
 *                                below which species are deleted 
 *                                from  the equilibrium problem.
 */
#ifndef VCS_DELETE_SPECIES_CUTOFF
#define VCS_DELETE_SPECIES_CUTOFF          1.0e-32
#endif
#ifndef VCS_DELETE_MINORSPECIES_CUTOFF
#define VCS_DELETE_MINORSPECIES_CUTOFF     1.0e-140
#endif

/*
* VCS_SMALL_MULTIPHASE_SPECIES:
*                            Relative value of multiphase
*                            species mole number for a multiphase
*                            species which is small.
*/
#define VCS_SMALL_MULTIPHASE_SPECIES  1.0e-25

/*
*     VCS_DELETE_PHASE_CUTOFF: Cutoff relative moles below
*                              which a phase is deleted 
*                              from  the equilibrium problem.
*/
#ifndef VCS_DELETE_PHASE_CUTOFF
#define VCS_DELETE_PHASE_CUTOFF     1.0e-11
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *      State of Dimensional Units for Gibbs free energies
 */
#define VCS_NONDIMENSIONAL_G       1
#define VCS_DIMENSIONAL_G          0

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *                  SPECIES CATEGORIES USED IN VCS_SOLVE_TP
 */
//! @name    Species Types during the iteration 
//! valid values for spStatus()
//@{
//! Species is a component
#define VCS_SPECIES_COMPONENT      2

//! Species is a major species
/*!
 * A major species is either a species in a multicomponent phase with
 * significant concentration or its a Stoich Phase
 */
#define VCS_SPECIES_MAJOR          1

//! Species is a major species
/*!
 * A major species is either a species in a multicomponent phase with
 * significant concentration or its a Stoich Phase
 */
#define VCS_SPECIES_MINOR          0

//! Species lies in a multicomponent phase that is zeroed atm
/*!
 * The species lies in a multicomponent phase that is currently
 * deleted. 
 */
#define VCS_SPECIES_ZEROEDPHASE   -1

//! Species lies in a multicomponent phase, with concentration zero
/*!
 *  The species lies in a multicomponent phase that exists.
 *  It concentration is currently zero, even though it may
 *  or may not actually have a low mole fraction in the phase
 *  this situation occurs when phases pop back into life.
 */
#define VCS_SPECIES_ZEROEDMS      -2

//! Species is a SS phase, that is currently zeroed out.
/*!
 *  The species lies in a single-species phase which
 *  is currently zereod out.
 */
#define VCS_SPECIES_ZEROEDSS      -3

//! Species has such a small mole fraction it is deleted.
/*!
 *  The species is believed to have such a small mole fraction
 *  that it best to throw the calculation of it out.
 *  It will be aded back in at the end of the calculation.
 */
#define VCS_SPECIES_DELETED       -4

//! Species refers to an electron in the metal
/*!
 *  The unknown is equal to the interfacial voltage
 *  drop across the interface on the SHE (standard
 *  hyrdogen electrode) scale (volts).
 */
#define VCS_SPECIES_INTERFACIALVOLTAGE  -5

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* 
 *      Units for the chemical potential data and pressure variables:
 *
 *                      Chem_Pot                 Pres      vol   moles
 *                        -------------------------------------------------
 *  VCS_UNITS_KCALMOL  = kcal/mol                  atm   cm**3   gmol
 *  VCS_UNITS_UNITLESS = MU / RT -> no units       atm   cm**3   gmol
 *  VCS_UNITS_KJMOL    = kJ / mol                  atm   cm**3   gmol
 *  VCS_UNITS_KELVIN   = KELVIN -> MU / R          atm   cm**3   gmol
 *  VCS_UNITS_MKS      = Joules / Kmol (Cantera)   Pa     m**3   kmol
 *
 *          Energy:
 *             VCS_UNITS_KCALMOL  = kcal/mol
 *             VCS_UNITS_UNITLESS = MU / RT -> no units
 *             VCS_UNITS_KJMOL    = kJ / mol
 *             VCS_UNITS_KELVIN   = KELVIN -> MU / R
 *             VCS_UNITS_MKS      = J / kmol
 *
 *          Pressure: (Pref and Pres)
 *             VCS_UNITS_KCALMOL  = atm
 *             VCS_UNITS_UNITLESS = no units
 *             VCS_UNITS_KJMOL    = atm
 *             VCS_UNITS_KELVIN   = atm
 *             VCS_UNITS_MKS      = Pa = kg / m s2
 */
#define VCS_UNITS_KCALMOL                   -1   
#define VCS_UNITS_UNITLESS                   0
#define VCS_UNITS_KJMOL                      1
#define VCS_UNITS_KELVIN                     2
#define VCS_UNITS_MKS                        3

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * Element type Defines
 */
#define VCS_ELEM_TYPE_ABSPOS           0                 
#define VCS_ELEM_TYPE_ELECTRONCHARGE   1                 
#define VCS_ELEM_TYPE_CHARGENEUTRALITY 2
#define VCS_ELEM_TYPE_OTHERCONSTRAINT  3
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * Species type Defines
 */
#define VCS_SPECIES_TYPE_MOLNUM              0
#define VCS_SPECIES_TYPE_INTERFACIALVOLTAGE -5

/****************************************************************************/
}

#endif
