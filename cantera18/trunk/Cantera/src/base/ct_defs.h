/**
 * @file ct_defs.h
 * This file contains definitions of terms that are used in internal 
 * routines and are unlikely to need modifying (text for module physConstants (see \ref physConstants) is found here). 
 * This file is included
 * in every file that is in the Cantera Namespace.
 *
 * All physical constants are stored here.
 * The module physConstants is defined here.
 */

/*
 * $Revision: 1.9 $
 * $Date: 2009/01/24 00:15:00 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_DEFS_H
#define CT_DEFS_H

#include <cmath>
#include "config.h"

// STL includes
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <numeric>
#include <string>
#include <algorithm>


//! creates a pointer to the start of the raw data for a ctvector
#define DATA_PTR(vec) &vec[0]

#ifdef WIN32
#define TYPENAME_KEYWORD
#pragma warning(disable:4267)
#else
//! create a define for the typename command
#define TYPENAME_KEYWORD typename
#endif

#undef CHEMKIN_COMPATIBILITY_MODE

/**
 * Namespace for the Cantera kernel.
 */
namespace Cantera {

 /*!
  * All physical constants are stored here.
  *
  * @defgroup physConstants Physical Constants
  * %Cantera uses the MKS system of units. The unit for moles 
  * is defined to be the kmol. All values of physical constants
  * are consistent with the 2006 CODATA recommendations.
  * @ingroup globalData
  * @{
  */
  //! Pi
    const doublereal Pi = 3.1415926;
  //! sqrt(Pi)
    const doublereal SqrtPi = std::sqrt(Pi);

  
  /*!
   *   @name Variations of the Gas Constant
   * %Cantera uses the MKS system of units. The unit for moles 
   * is defined to be the kmol.
   */
  //@{

    /// Avogadro's Number
    const doublereal Avogadro = 6.02214179e26; //   /kmol

    /// Universal Gas Constant. 2006 CODATA value.
    const doublereal GasConstant = 8314.47215; // J/kmol/K

    const doublereal logGasConstant = 9.025752908;

    /// One atmosphere
    const doublereal OneAtm = 1.01325e5;             // Pa

    /// Universal gas constant in cal/mol/K
    const doublereal GasConst_cal_mol_K = 1.987;

    /// Boltzmann's constant
    const doublereal Boltzmann = GasConstant / Avogadro;

    /// Planck's constant. Units of J-s
    const doublereal Planck = 6.62606896e-34;           // J-s
    const doublereal Planck_bar = 1.05457162853e-34;    // m2-kg/s

    /// log(k/h)
    const doublereal logBoltz_Planck = 23.7599032;  // ln(k_B/h)
    /// Stefan-Boltzmann constant
    const doublereal StefanBoltz = 5.6704004e-8;

    //@}
    /// @name Electron Properties
    //@{
    const doublereal ElectronCharge = 1.60217648740e-19; // C
    const doublereal ElectronMass = 9.1093821545e-31; // kg
    const doublereal Faraday = ElectronCharge * Avogadro;
    //@}

    /// @name Electromagnetism
    /// %Cantera uses the MKS unit system.
    //@{

    /// Permittivity of free space \f$ \epsilon_0 \f$ in F/m.
    const doublereal epsilon_0 = 8.85417817e-12;  // Farads/m = C^2/N/m^2

    /// Permeability of free space \f$ \mu_0 \f$ in N/A^2.
    const doublereal permeability_0 = 4.0e-7*Pi; // N/A^2

    /// Speed of Light (m/s).
    const doublereal lightSpeed = 1.0/std::sqrt(epsilon_0 * permeability_0);

    //@}
 //@}
 
  /*!
   * @name Thermodynamic Equilibrium Constraints
   * Integer numbers representing pairs of thermodynamic variables
   * which are held constant during equilibration.
   */
  //@{
    const int TV = 100, HP = 101, SP = 102, PV = 103, TP = 104, UV = 105, 
              ST = 106, SV = 107, UP = 108, VH = 109, TH = 110, SH = 111,
              PX = 112, TX = 113;
    const int VT = -100, PH = -101, PS = -102, VP = -103, PT = -104, 
              VU = -105, TS = -106, VS = -107, PU = -108, HV = -109, 
              HT = -110, HS = -111, XP = -112, XT = -113;
  //@}

  //! 1/3
    const doublereal OneThird = 1.0/3.0;
  //! 5/16
    const doublereal FiveSixteenths = 5.0/16.0;
  //! sqrt(10)
    const doublereal SqrtTen = std::sqrt(10.0);
  //! sqrt(8)
    const doublereal SqrtEight = std::sqrt(8.0);
  //! sqrt(2)
    const doublereal SqrtTwo = std::sqrt(2.0);

  //! smallest number to compare to zero.
    const doublereal SmallNumber = 1.e-300;
  //! largest number to compare to inf.
    const doublereal BigNumber = 1.e300;
  //! largest x such that exp(x) is valid
    const doublereal MaxExp = 690.775527898;         

  //! Fairly random number to be used to initialize variables against
  //! to see if they are subsequently defined.
    const doublereal Undef = -999.1234;
  //! Small number to compare differences of mole fractions against.
    const doublereal Tiny = 1.e-20;

  //! inline function to return the max value of two doubles.
  /*!
   * @param x  double value
   * @param y  second double value
   */
  inline doublereal fmaxx(doublereal x, doublereal y)
  { return (x > y) ? x : y; }

  //! inline function to return the min value of two doubles.
  /*!
   * @param x  double value
   * @param y  second double value
   */
  inline doublereal fminn(doublereal x, doublereal y)
  { return (x < y) ? x : y; }
  

  //! Map connecting a string name with a double.
  /*!
   * This is used mostly to assign concentrations and mole fractions
   * to species.
   */
    typedef std::map<std::string, doublereal> compositionMap;
//! Turn on the use of stl vectors for the basic array type within cantera
#define USE_STL_VECTOR
#ifdef USE_STL_VECTOR
  //! Vector of doubles.
    typedef std::vector<double>        array_fp;
  //! Vector of doubles.
    typedef std::vector<double>        vector_fp;
  //! Vector of ints
    typedef std::vector<int>           array_int;
  //! Vector of ints
    typedef std::vector<int>           vector_int;
#else
    typedef ct::ctvector_fp            array_fp;
    typedef ct::ctvector_fp            vector_fp;
    typedef ct::ctvector_int           array_int;
    typedef ct::ctvector_int           vector_int;
#endif
  //! typedef for a group of species.
  /*!
   * A group of species is a subset of the species in a phase.
   */
    typedef vector_int         group_t;
  //! typedef for a vector of groups of species.
  /*!
   * A grouplist of species is a vector of groups.
   */
  typedef std::vector<group_t>    grouplist_t;
 
  //! Typedef for a pointer to temporary work storage
    typedef doublereal* workPtr;
  //! typedef for a pointer to temporary work storage which is treated as constant
    typedef const doublereal* const_workPtr;


}  // namespace

#endif



