/**
 *
 *  @file ThermoPhase.cpp
 */

/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"


namespace Cantera {


    void ThermoPhase::setState_HP(doublereal h, doublereal p, 
        doublereal tol) {
        doublereal dt;
        setPressure(p);
        for (int n = 0; n < 20; n++) {
            dt = (h - enthalpy_mass())/cp_mass();
            if (dt > 100.0) dt = 100.0;
            else if (dt < -100.0) dt = -100.0; 
            setState_TP(temperature() + dt, p);
            if (fabs(dt) < tol) {
                return;
            }
        }
        throw CanteraError("setState_HP","no convergence. dt = " + fp2str(dt));
    }

    void ThermoPhase::setState_UV(doublereal u, doublereal v, 
        doublereal tol) {
        doublereal dt;
        setDensity(1.0/v);
        for (int n = 0; n < 20; n++) {
            dt = (u - intEnergy_mass())/cv_mass();
            if (dt > 100.0) dt = 100.0;
            else if (dt < -100.0) dt = -100.0; 
            setTemperature(temperature() + dt);
            if (fabs(dt) < tol) {
                return;
            }
        }
        throw CanteraError("setState_UV","no convergence. dt = " + fp2str(dt));
    }

    void ThermoPhase::setState_SP(doublereal s, doublereal p, 
        doublereal tol) {
        doublereal dt;
        setPressure(p);
        for (int n = 0; n < 20; n++) {
            dt = (s - entropy_mass())*temperature()/cp_mass();
            if (dt > 100.0) dt = 100.0;
            else if (dt < -100.0) dt = -100.0; 
            setState_TP(temperature() + dt, p);
            if (fabs(dt) < tol) {
                return;
            }
        }
        throw CanteraError("setState_SP","no convergence. dt = " + fp2str(dt));
    }

    void ThermoPhase::setState_SV(doublereal s, doublereal v, 
        doublereal tol) {
        doublereal dt;
        setDensity(1.0/v);
        for (int n = 0; n < 20; n++) {
            dt = (s - entropy_mass())*temperature()/cv_mass();
            if (dt > 100.0) dt = 100.0;
            else if (dt < -100.0) dt = -100.0; 
            setTemperature(temperature() + dt);
            if (fabs(dt) < tol) {
                return;
            }
        }
        throw CanteraError("setState_SV","no convergence. dt = " + fp2str(dt));
    }

    doublereal ThermoPhase::err(string msg) const {
            throw CanteraError("ThermoPhase","Base class method "
                +msg+" called.");
            return 0;
    }

	/**
	 * Returns the units of the standard and general concentrations
	 * Note they have the same units, as their divisor is 
	 * defined to be equal to the activity of the kth species
	 * in the solution, which is unitless.
	 *
	 * This routine is used in print out applications where the
	 * units are needed. Usually, MKS units are assumed throughout
	 * the program and in the XML input files. 
	 *
	 * On return uA contains the powers of the units (MKS assumed)
	 * of the standard concentrations and generalized concentrations
	 * for the kth species.
	 *
	 *  uA[0] = kmol units - default  = 1
	 *  uA[1] = m    units - default  = -nDim(), the number of spatial
	 *                                dimensions in the Phase class.
	 *  uA[2] = kg   units - default  = 0;
	 *  uA[3] = Pa(pressure) units - default = 0;
	 *  uA[4] = Temperature units - default = 0;
	 *  uA[5] = time units - default = 0
	 */
    void ThermoPhase::getUnitsStandardConc(double *uA, int k, int sizeUA) {
	for (int i = 0; i < sizeUA; i++) {
	  if (i == 0) uA[0] = 1.0;
	  if (i == 1) uA[1] = -nDim();
	  if (i == 2) uA[2] = 0.0;
	  if (i == 3) uA[3] = 0.0;
	  if (i == 4) uA[4] = 0.0;
	  if (i == 5) uA[5] = 0.0;
	}
    }

}




