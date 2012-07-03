/**
 *  @file ckr_defs.h
 *
 */

// Copyright 2001  California Institute of Technology



#ifndef CKR_DEFS_H
#define CKR_DEFS_H

//#include "cantera/base/config.h"
#include <string>
#include <iostream>
#include <vector>

/// the namespace for the CKReader packaage
namespace ckr
{

typedef std::vector<double> vector_fp;
typedef std::vector<double> vector_int;

// exceptions
class CK_Exception
{
public:
    CK_Exception() {
        m_msg = "";
    }
    virtual ~CK_Exception() {}
    std::string errorMessage() {
        return m_msg;
    }
protected:
    std::string m_msg;
};

const double UNDEF = -9999.1234;

/**
 * @name Falloff Parameterizations
 * These constants are used to specify which falloff parameterization
 * to use for a pressure-dependent reaction.
 * @see Reaction.h
 */

//@{
const int Lindemann = 0;
const int Troe = 1;
const int SRI = 2;
//@}



/** @name Reaction Rate Types
 * These constant are used to specify the type of rate coefficient
 */
//@{
const int Arrhenius = 0, LandauTeller = 1, Jan = 2, Fit1 = 3;
//@}


/** @name Activation Energy Units
 *  These constants specify the supported units for the activation energy of
 *  a reaction
 *  The default is to assume Cal_per_Mole for unspecified units in the activation energy
 *  as this was the original default
 */
//@{
const int Cal_per_Mole = 1,
          Kcal_per_Mole = 2,
          Joules_per_Mole = 3,
          Kelvin= 4,
          Electron_Volts = 5,
          Kjoules_per_Mole = 6;
//@}

/**
 *  @name Quantity Units
 *  These constants define the supported units for number of molecules.
 */
//@{
const int Moles = 100;     ///< specify number of moles (6.023e23 molecules / mole)
const int Molecules = 101; ///< specify number of molecules
//@}

} // namespace


#endif




