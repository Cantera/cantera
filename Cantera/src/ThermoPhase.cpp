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

}




