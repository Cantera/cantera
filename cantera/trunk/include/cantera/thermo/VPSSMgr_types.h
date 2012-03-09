/**
 *  @file VPSSMgr_types.h
 *       Contains const definitions for types of calculation managers
 *       that are responsible for calculating the species standard
 *       state thermodynamic managers and
 *       reference-state thermodynamics managers
 *        (see
 *        class \link Cantera::VPSSMgr VPSSMgr\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

//! @deprecated remove include when UnknownVPSSMgr is removed
#include "cantera/base/stringUtils.h"

#ifndef VPSSMGR_TYPES_H
#define VPSSMGR_TYPES_H

//! Variable pressures SS calculator for ideal gas phases
#define VPSSMGR_IDEALGAS   1

//! Variable pressure SS calculate for phases consisting all
//!  species having a constant molar volume property
/*!
 *  This fits most solids
 */
#define VPSSMGR_CONSTVOL   2

//! Variable pressure SS calculate for phases consisting of real water
//! as the first species and species having a constant molar volume property
#define VPSSMGR_WATER_CONSTVOL   11

//! Variable pressure SS calculate for phases consisting of real water
//! as the first species and species obeying the HKFT standard state
#define VPSSMGR_WATER_HKFT   12

//! Variable pressure SS calculate for phases consisting of completing
//! general representations
#define VPSSMGR_GENERAL  22

namespace Cantera
{

//! Error for unknown thermo parameterization
class UnknownVPSSMgr : public CanteraError
{
public:
    //! Constructor
    /*!
     * @param func        String function id
     * @param thermotype  Integer specifying the thermo parameterization
     *
     * deprecated This class is unused
     */
    DEPRECATED(UnknownVPSSMgr(std::string func, int thermotype)) {
        CanteraError(func, std::string("\n ### ERROR ### \n") +
                     "Unknown species thermo parameterization ("
                     + int2str(thermotype) + ")\n\n");
    }
};


}

#endif


