/*
 *  IdealGasPhase.h
 *  Cantera
 *
 *  Created by David Goodwin on 4/21/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

namespace Cantera {
	
	class IdealGasPhase : public Phase {
	public: 
		Real pressure(const State& s);
		
		
	protected:
