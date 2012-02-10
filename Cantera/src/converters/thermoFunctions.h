/**
 *  @file thermoFunctions.h
 *
 *  Thermodynamic properties. Note that these functions are used only
 *  for validation purposes by CKReader. They are not used by Cantera.
 */

// Copyright 2001  California Institute of Technology

#ifndef CKR_THERMOFUNCTIONS_H
#define CKR_THERMOFUNCTIONS_H

#include "Species.h"

namespace ckr
{

double enthalpy(double t, const Species& s);
double cp(double t, const Species& s);
double entropy(double t, const Species& s);
double gibbs(double t, const Species& s);

}

#endif

