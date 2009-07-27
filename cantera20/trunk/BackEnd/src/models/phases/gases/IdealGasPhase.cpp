/*
 *  IdealGasPhase.cpp
 *  Cantera
 *
 *  Created by David Goodwin on 4/21/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "IdealGasPhase.h"

namespace Cantera {

	IdealGasPhase::pressure( const State& s) {
	return s.N_ * GasConstant *  s.T_ / s.V_
	}
		
IdealGasPhase::	