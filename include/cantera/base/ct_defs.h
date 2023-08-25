/**
 * @file ct_defs.h
 * This file contains definitions of constants, types and terms that are used
 * in internal routines and are unlikely to need modifying.
 *
 * All physical constants are stored here (see module @ref physConstants).
 *
 * This file is included in every file within the Cantera namespace.
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
#include <set>
#include <string>
#include <algorithm>
#include <memory>
#include <functional>

/**
 * Namespace for the Cantera kernel.
 */
namespace Cantera
{

using std::shared_ptr;
using std::make_shared;
using std::unique_ptr;
using std::make_unique;
using std::isnan; // workaround for bug in libstdc++ 4.8
using std::string;
using std::vector;
using std::map;
using std::set;
using std::function;
using std::pair;

/**
 * @defgroup physConstants Physical Constants
 * @brief Physical constants used by %Cantera.
 *
 * %Cantera uses the MKS system of units. The unit for moles
 * is defined to be the kmol. All values of physical constants
 * are consistent with the 2018 CODATA recommendations.
 * @ingroup globalData
 */

//! @addtogroup physConstants
//! @{

/**
 * @name Numerical Constants
 * @{
 */

//! Pi
const double Pi = 3.14159265358979323846;

//! Sqrt(2)
const double Sqrt2 = 1.41421356237309504880;

//! @}
//! @name Defined Constants
//!
//! These constants are defined by CODATA to have a particular value; see
//! [NIST Reference on Constants, Units, and Uncertainty](https://physics.nist.gov/cuu/Constants/index.html).
//! @{

//! Avogadro's Number @f$ N_{\mathrm{A}} @f$ [number/kmol]
const double Avogadro = 6.02214076e26;

//! Boltzmann constant @f$ k @f$ [J/K]
const double Boltzmann = 1.380649e-23;

//! Planck constant @f$ h @f$ [J-s]
const double Planck = 6.62607015e-34;

//! Elementary charge @f$ e @f$ [C]
const double ElectronCharge = 1.602176634e-19;

//! Speed of Light in a vacuum @f$ c @f$ [m/s]
const double lightSpeed = 299792458.0;

//! One atmosphere [Pa]
const double OneAtm = 1.01325e5;

//! One bar [Pa]
const double OneBar = 1.0E5;

//! @}
//! @name Measured Constants
//!
//! These constants are measured and reported by CODATA
//! @{

//! Fine structure constant @f$ \alpha @f$ []
const double fineStructureConstant = 7.2973525693e-3;

//! Electron Mass @f$ m_e @f$ [kg]
const double ElectronMass = 9.1093837015e-31;

//! @}
//! @name Derived Constants
//!
//! These constants are found from the defined and measured constants
//! @{

//! Universal Gas Constant @f$ R_u @f$ [J/kmol/K]
const double GasConstant = Avogadro * Boltzmann;

const double logGasConstant = std::log(GasConstant);

//! Universal gas constant in cal/mol/K
const double GasConst_cal_mol_K = GasConstant / 4184.0;

//! Stefan-Boltzmann constant @f$ \sigma @f$ [W/m2/K4]
const double StefanBoltz = 2.0 * std::pow(Pi, 5) * std::pow(Boltzmann, 4) / (15.0 * std::pow(Planck, 3) * lightSpeed * lightSpeed); // 5.670374419e-8

//! Faraday constant @f$ F @f$ [C/kmol]
const double Faraday = ElectronCharge * Avogadro;

//! Permeability of free space @f$ \mu_0 @f$ [N/A2]
const double permeability_0 = 2 * fineStructureConstant * Planck / (ElectronCharge * ElectronCharge * lightSpeed);

//! Permittivity of free space @f$ \varepsilon_0 @f$ [F/m]
const double epsilon_0 = 1.0 / (lightSpeed * lightSpeed * permeability_0);

//! @}

//! @}

//! @name Thermodynamic Equilibrium Constraints
//!
//! Integer numbers representing pairs of thermodynamic variables
//! which are held constant during equilibration.
//! @{

const int TV = 100, HP = 101, SP = 102, PV = 103, TP = 104, UV = 105,
          ST = 106, SV = 107, UP = 108, VH = 109, TH = 110, SH = 111,
          PX = 112, TX = 113;
const int VT = -100, PH = -101, PS = -102, VP = -103, PT = -104,
          VU = -105, TS = -106, VS = -107, PU = -108, HV = -109,
          HT = -110, HS = -111, XP = -112, XT = -113;
//! @}

//! smallest number to compare to zero.
const double SmallNumber = 1.e-300;
//! largest number to compare to inf.
const double BigNumber = 1.e300;

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

//! Map from string names to doubles. Used for defining species mole/mass
//! fractions, elemental compositions, and reaction stoichiometries.
typedef map<string, double> Composition;

//! index returned by functions to indicate "no position"
const size_t npos = static_cast<size_t>(-1);

} // namespace

#endif
