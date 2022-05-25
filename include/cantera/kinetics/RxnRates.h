/**
 *  @file RxnRates.h
 *  @deprecated To be removed after Cantera 2.6. See class Cantera::ReactionRate and
 *      derived classes for new reaction rate handlers.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RXNRATES_H
#define CT_RXNRATES_H

#include "cantera/kinetics/reaction_defs.h"
#include "cantera/kinetics/Arrhenius.h"
#include "cantera/kinetics/PlogRate.h"
#include "MultiRate.h"
#include "cantera/base/Array.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

namespace Cantera
{

typedef ArrheniusRate Arrhenius;

typedef PlogRate Plog;

}

#endif
