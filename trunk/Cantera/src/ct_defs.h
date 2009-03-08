/**
 * 
 * @file ct_defs.h
 *
 * This file contains definitions of terms that are used in internal 
 * routines and are unlikely to need modifying
 */

/* $Author: hkmoffa $
 * $Revision: 1.13 $
 * $Date: 2006/04/30 22:08:33 $
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

using namespace std;

//#include "ctvector.h"
//using namespace ct;

#define DATA_PTR(vec) &vec[0]

#ifdef WIN32
#define TYPENAME_KEYWORD
#pragma warning(disable:4267)
#else
#define TYPENAME_KEYWORD typename
#endif

#undef CHEMKIN_COMPATIBILITY_MODE

/**
 * Namespace for the Cantera kernel.
 */
namespace Cantera {

    const doublereal Pi = 3.1415926;
    const doublereal SqrtPi = sqrt(Pi);

    // use kg-moles, rather than g-moles.
    // Note: this constant is a relic of old versions,
    // and appears to be no longer used anywhere.
    const doublereal CtMoles_per_mole = 1.e-3;       // kmol

    /// @name Physical Constants
    //@{

    /// Avogadro's Number
    const doublereal Avogadro = 6.022136736e26;

    //const doublereal GasConstant = 8314.0;

    /// Universal Gas Constant. 1999 CODATA value.
    const doublereal GasConstant = 8314.47215; // J/kmol/K

    const doublereal logGasConstant = 9.025752908;

    /// One atmosphere
    const doublereal OneAtm = 1.01325e5;             // Pa

    /// Universal gas constant in cal/mol/K
    const doublereal GasConst_cal_mol_K = 1.987;

    /// Boltzmann's constant
    const doublereal Boltzmann = GasConstant / Avogadro;

    /// Planck's constant
    const doublereal Planck = 6.626068e-34;           // J-s
    const doublereal Planck_bar = 1.05457148e-34;    // m2-kg/s

    /// log(k/h)
    const doublereal logBoltz_Planck = 23.7599032;  // ln(k_B/h)
    /// Stefan-Boltzmann constant
    const doublereal StefanBoltz = 5.67e-8;

    /// @name Electron Properties
    //@{
    const doublereal ElectronCharge = 1.602e-19; // C
    const doublereal ElectronMass = 9.10938188e-31; // kg
    const doublereal Faraday = ElectronCharge * Avogadro;
    //@}

    /// @name Electromagnetism
    /// Cantera uses the MKS unit system.
    //@{

    /// Permittivity of free space \f$ \epsilon_0 \f$ in F/m.
    const doublereal epsilon_0 = 8.8542e-12;  // Farads/m = C^2/N/m^2

    /// Permeability of free space \f$ \mu_0 \f$ in N/A^2.
    const doublereal permeability_0 = 4.0e-7*Pi; // N/A^2

    //@}

    //@} 

    const doublereal OneThird = 1.0/3.0;
    const doublereal FiveSixteenths = 5.0/16.0;
    const doublereal SqrtTen = sqrt(10.0);
    const doublereal SqrtEight = sqrt(8.0);

    const doublereal SmallNumber = 1.e-300;
    const doublereal BigNumber = 1.e300;

    /// largest x such that exp(x) is valid
    const doublereal MaxExp = 690.775527898;         

    const int Undefined = -999;
    const doublereal Undef = -999.1234;
    const doublereal Cutoff = 1.e-12;
    const doublereal Tiny = 1.e-20;


    const int TV = 100, HP = 101, SP = 102, PV = 103, TP = 104, UV = 105, 
              ST = 106, SV = 107, UP = 108, VH = 109, TH = 110, SH = 111,
              PX = 112, TX = 113;
    const int VT = -100, PH = -101, PS = -102, VP = -103, PT = -104, 
              VU = -105, TS = -106, VS = -107, PU = -108, HV = -109, 
              HT = -110, HS = -111, XP = -112, XT = -113;

    inline doublereal fmaxx(doublereal x, doublereal y)
    { return (x > y) ? x : y; }
    inline doublereal fminn(doublereal x, doublereal y)
    { return (x < y) ? x : y; }

    const int GAS = 0;
    const int LIQUID = 1;
    const int SOLID = 2;
    const int PURE_FLUID = 3;

    //    enum Phase {GAS, LIQUID, SOLID, PURE_FLUID};

    const int Solid_Phase = 0,
             Liquid_Phase = 1,
              Vapor_Phase = 2,
                Gas_Phase = 2;

    const int None = 0;

    // typedefs
    typedef std::map<string, doublereal> compositionMap;
#define USE_STL_VECTOR
#ifdef USE_STL_VECTOR
    typedef std::vector<double>            array_fp;
    typedef std::vector<double>            vector_fp;
    typedef std::vector<int>           array_int;
    typedef std::vector<int>           vector_int;
#else
    typedef ct::ctvector_fp            array_fp;
    typedef ct::ctvector_fp            vector_fp;
    typedef ct::ctvector_int           array_int;
    typedef ct::ctvector_int           vector_int;
#endif
    typedef vector_int         group_t;
    typedef std::vector<group_t>    grouplist_t;
 
    typedef doublereal* workPtr;
    typedef const doublereal* const_workPtr;


    //    template<class A, class B>
    //inline doublereal operator*(const vector<A>& u, const vector<B>& v) {
    //    return inner_product(u.begin(), u.end(), v.begin(), 0.0);
    // }


}  // namespace

#endif



