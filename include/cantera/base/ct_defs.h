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

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DEFS_H
#define CT_DEFS_H

#include "config.h"

#include <cmath>

// STL includes
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <memory>

/**
 * Namespace for the Cantera kernel.
 */
namespace Cantera
{

using std::shared_ptr;
using std::make_shared;
using std::unique_ptr;
using std::isnan; // workaround for bug in libstdc++ 4.8

/*!
 * All physical constants are stored here.
 *
 * @defgroup physConstants Physical Constants
 * %Cantera uses the MKS system of units. The unit for moles
 * is defined to be the kmol. All values of physical constants
 * are consistent with the 2018 CODATA recommendations.
 * @ingroup globalData
 * @{
 */

//! Pi
const double Pi = 3.14159265358979323846;

/*!
 * @name Defined Constants
 * These constants are defined by CODATA to have a particular value.
 * https://physics.nist.gov/cuu/Constants/index.html
 */
//@{

//! Avogadro's Number \f$ N_{\mathrm{A}} \f$ [number/kmol]
const double Avogadro = 6.02214076e26;

//! Boltzmann constant \f$ k \f$ [J/K]
const double Boltzmann = 1.380649e-23;

//! Planck constant \f$ h \f$ [J-s]
const double Planck = 6.62607015e-34;

//! Elementary charge \f$ e \f$ [C]
const double ElectronCharge = 1.602176634e-19;

/// Speed of Light in a vacuum \f$ c \f$ [m/s]
const double lightSpeed = 299792458.0;

//! One atmosphere [Pa]
const double OneAtm = 1.01325e5;

//! One bar [Pa]
const double OneBar = 1.0E5;

//@}

/*!
 * @name Measured Constants
 * These constants are measured and reported by CODATA
 */
//@{

//! Fine structure constant \f$ \alpha \f$ []
const double fineStructureConstant = 7.2973525693e-3;

//! Electron Mass \f$ m_e \f$ [kg]
const double ElectronMass = 9.1093837015e-31;

//@}

/*!
 * @name Derived Constants
 * These constants are found from the defined and measured constants
 */
//@{

//! Reduced Planck constant \f$ \hbar \f$ [m2-kg/s]
const double Planck_bar = Planck / (2 * Pi);

//! Universal Gas Constant \f$ R_u \f$ [J/kmol/K]
const double GasConstant = Avogadro * Boltzmann;

const double logGasConstant = std::log(GasConstant);

//! Universal gas constant in cal/mol/K
const double GasConst_cal_mol_K = GasConstant / 4184.0;

//! log(k_b/h)
const double logBoltz_Planck = std::log(Boltzmann / Planck);

//! Stefan-Boltzmann constant \f$ \sigma \f$ [W/m2/K4]
const double StefanBoltz = Pi * Pi * std::pow(Boltzmann, 4.0) / (60.0 * std::pow(Planck_bar, 3.0) * lightSpeed * lightSpeed); // 5.670374419e-8

//! Faraday constant \f$ F \f$ [C/kmol]
const double Faraday = ElectronCharge * Avogadro;

//! Permeability of free space \f$ \mu_0 \f$ [N/A2]
const double permeability_0 = 2 * fineStructureConstant * Planck / (ElectronCharge * ElectronCharge * lightSpeed);

//! Permittivity of free space \f$ \varepsilon_0 \f$ [F/m]
const double epsilon_0 = 1.0 / (lightSpeed * lightSpeed * permeability_0);

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

//! smallest number to compare to zero.
const double SmallNumber = 1.e-300;
//! largest number to compare to inf.
const double BigNumber = 1.e300;
//! largest x such that exp(x) is valid
const double MaxExp = 690.775527898;

//! Fairly random number to be used to initialize variables against
//! to see if they are subsequently defined.
const double Undef = -999.1234;

//! Small number to compare differences of mole fractions against.
/*!
 * This number is used for the interconversion of mole fraction and mass
 * fraction quantities when the molecular weight of a species is zero. It's also
 * used for the matrix inversion of transport properties when mole fractions
 * must be positive.
 */
const double Tiny = 1.e-20;

//! Map connecting a string name with a double.
/*!
 * This is used mostly to assign concentrations and mole fractions to species.
 */
typedef std::map<std::string, double> compositionMap;

//! Map from string names to doubles. Used for defining species mole/mass
//! fractions, elemental compositions, and reaction stoichiometries.
typedef std::map<std::string, double> Composition;

//! Turn on the use of stl vectors for the basic array type within cantera
//! Vector of doubles.
typedef std::vector<double> vector_fp;
//! Vector of ints
typedef std::vector<int> vector_int;

//! A grouplist is a vector of groups of species
typedef std::vector<std::vector<size_t> > grouplist_t;

//! index returned by functions to indicate "no position"
const size_t npos = static_cast<size_t>(-1);

} // namespace

#endif
