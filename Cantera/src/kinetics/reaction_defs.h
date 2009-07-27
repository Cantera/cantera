/**
 *  @file reaction_defs.h
 * This file defines some constants used to specify reaction types.
 */
  
/*
 *  $Author: dggoodwin $
 *  $Date: 2007/05/04 14:27:24 $
 *  $Revision: 1.1 $
 *
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_RXN_DEFS_H
#define CT_RXN_DEFS_H

#include "ct_defs.h"

namespace Cantera {

    const int NONE = 0;

    /// @name Reaction Types

    //@{

    /**
     * A reaction with a rate coefficient that depends only on
     * temperature. Example: O + OH <-> O2 + H
     */
    const int ELEMENTARY_RXN = 1;

    /**
     * A reaction that requires a third-body collision partner. Example:
     * O2 + M <-> O + O + M
     */ 
    const int THREE_BODY_RXN = 2;
    
    /**
     * The general form for an association or dissociation reaction, with a
     * pressure-dependent rate. Example: CH3 + H (+M) <-> CH4 (+M)
     */
    const int FALLOFF_RXN    = 4;
    
    /**
     * A chemical activation reaction. For these reactions, the rate falls
     * off as the pressure increases, due to collisional stabilization of
     * a reaction intermediate. Example: Si + SiH4 (+M) <-> Si2H2 + H2
     * (+M), which competes with Si + SiH4 (+M) <-> Si2H4 (+M). 
     * @todo Implement chemical activation reactions.
     */
    const int CHEMACT_RXN    = 8;

    /**
     * A reaction occurring on a surface.
     */
    const int SURFACE_RXN    = 20;

    /**
     * A reaction occurring at a one-dimensional interface between two 
     * surface phases.
     */
    const int EDGE_RXN  = 22;

    /**
     * A global reaction. These may have non-integral reaction orders,
     * and are not allowed to be reversible.
     */
    const int GLOBAL_RXN     = 30;

    //@}
    
    /** @name Rate Coefficient Types
     * These types define the supported rate coefficient types for 
     * elementary reactions. Any of these may also be used as the high and
     * low-pressure limits of falloff and chemical activation reactions.
     *
     * Note that not all of these are currently implemented!
     * @todo Finish implementing reaction rate types.
     */
    //@{
    
    const int  ARRHENIUS = 1;
    const int  LANDAUTELLER = 2;
    const int  TSTRATE = 3;
    const int  SURF_ARRHENIUS = 4;
    const int  ARRHENIUS_SUM = 5;

    //@}

    /** @name Falloff Function Types
     */
    //@{    
    const int SIMPLE_FALLOFF = 100;
    const int TROE3_FALLOFF = 110;
    const int TROE4_FALLOFF = 111;
    const int SRI3_FALLOFF  = 112;
    const int SRI5_FALLOFF  = 113;
    const int WF_FALLOFF = 114;
    //@}

    // error flags
    const int NO_ERROR = 0;
    const int UNKNOWN_REACTION_TYPE = -100;
    const int UNKNOWN_RATE_COEFF_TYPE = -200;
    const int NOT_YET_IMPLEMENTED = -300;
    
}

#endif
