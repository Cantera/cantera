/**
 * 
 * @file ct_defs.h
 *
 * This file contains definitions of terms that are used in internal 
 * routines and are unlikely to need modifying
 */

/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_DEFS_H
#define CT_DEFS_H

#include <math.h>
#include "config.h"

// STL includes
//#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <numeric>
#include <string>

using namespace std;

#include "ctvector.h"
using namespace ct;

#ifdef WIN32
#define TYPENAME_KEYWORD
#else
#define TYPENAME_KEYWORD typename
#endif

#undef CHEMKIN_COMPATIBILITY_MODE

/**
 * Namespace for the Cantera kernel.
 */
namespace Cantera {

#define CANTERA_VERSION 1.3

    // use kg-moles, rather than g-moles.
#define USE_KMOL

#ifdef USE_KMOL
    const doublereal CtMoles_per_mole = 1.e-3;       // kmol
    /// Avogadro's Number
    const doublereal Avogadro = 6.022136736e26;
    /// Universal Gas Constant
    const doublereal GasConstant = 8314.0;
    const doublereal logGasConstant = 9.02569612;
#else
    const doublereal CtMoles_per_mole = 1.0;         // mol
    const doublereal Avogadro = 6.022136736e23;
    const doublereal GasConstant = 8.3140;
    const doublereal logGasConstant = 2.117940841;
#endif

    /// One atmosphere
    const doublereal OneAtm = 1.01325e5;             // Pa
    /// Universal gas constant in cal/mol/K
    const doublereal GasConst_cal_mol_K = 1.987;
    /// Boltzmann's constant
    const doublereal Boltzmann = GasConstant / Avogadro;
    /// Planck's constant
    const doublereal Planck = 6.6262e-34;           // J-s
    /// log(k/h)
    const doublereal logBoltz_Planck = 23.7599032;  // ln(k_B/h)
    /// Stefan-Boltzmann constant
    const doublereal StefanBoltz = 5.67e-8;
    const doublereal ElectronCharge = 1.602e-19;
    const doublereal Faraday = ElectronCharge * Avogadro;
    const doublereal epsilon_0 = 8.85e-12;  // farads / m
    const doublereal Pi = 3.1415926;
    const doublereal SqrtPi = sqrt(Pi);

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
    typedef ct::ctvector_fp            array_fp;
    typedef ct::ctvector_fp            vector_fp;
    typedef ct::ctvector_int           array_int;
    typedef ct::ctvector_int           vector_int;
    typedef vector_int         group_t;
    typedef std::vector<group_t>    grouplist_t;
 
    typedef vector_fp::iterator workPtr;
    typedef vector_fp::const_iterator const_workPtr;


    template<class A, class B>
    inline doublereal operator*(const vector<A>& u, const vector<B>& v) {
        return inner_product(u.begin(), u.end(), v.begin(), 0.0);
    }

    template<class A>
    inline ostream& operator<<(ostream& s, const vector<A>& v) {
        int n = v.size();
        s << "<";
        for (int i = 0; i < n; i++) {
            s << v[i];
            if (i < n-1) s << ", ";
        }
        s << ">";
        return s;
    }

}  // namespace

#endif



