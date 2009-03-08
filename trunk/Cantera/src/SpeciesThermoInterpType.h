/**
 *  @file SpeciesThermoInterpType.h
 *
 * $Author: dggoodwin $
 * $Revision: 1.4 $
 * $Date: 2006/11/06 15:20:34 $
 */

// Copyright 2001  California Institute of Technology

#include "speciesThermoTypes.h"

#ifndef CT_SPECIESTHERMOINTERPTYPE_H
#define CT_SPECIESTHERMOINTERPTYPE_H

namespace Cantera {

    /**
     * Base class.
     */    
    class SpeciesThermoInterpType {

    public:
    
        SpeciesThermoInterpType() {};
        virtual ~SpeciesThermoInterpType() {};

        virtual SpeciesThermoInterpType * 
          duplMyselfAsSpeciesThermoInterpType() const = 0;
       
        virtual doublereal minTemp() const = 0;
        virtual doublereal maxTemp() const = 0;
	virtual doublereal refPressure() const = 0;
        virtual int reportType() const = 0;

        virtual void updateProperties(const doublereal* tempPoly, 
				      doublereal* cp_R,
				      doublereal* h_RT,
				      doublereal* s_R) const = 0;

        virtual void updatePropertiesTemp(const doublereal temp, 
				      doublereal* cp_R,
				      doublereal* h_RT,
				      doublereal* s_R) const = 0;

	virtual void reportParameters(int &n, int &type,
				      doublereal &tlow, doublereal &thigh,
				      doublereal &pref,
				      doublereal* const coeffs) const = 0;

        virtual void modifyParameters(doublereal* coeffs) {}

    };

}

#endif

