/**
 *  @file SurfPhase.cpp 
 *
 */

// Copyright 2002  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "SurfPhase.h"

//#include "ReactionData.h"
//#include "StoichManager.h"
//#include "RateCoeffMgr.h"

#include <iostream>
using namespace std;


    ///////////////////////////////////////////////////////////
    //
    //    class SurfPhase methods
    //
    ///////////////////////////////////////////////////////////

namespace Cantera {

    SurfPhase::
    SurfPhase(doublereal n0):
	ThermoPhase(),
	m_n0(n0),
	m_logn0(0.0),
	m_tmin(0.0),
	m_tmax(0.0),
	m_press(OneAtm),
	m_tlast(0.0) 
    {
	if (n0 > 0.0) m_logn0 = log(n0);
        setNDim(2);
    }

    doublereal SurfPhase::
    enthalpy_mole() const {
        if (m_n0 <= 0.0) return 0.0; 
        _updateThermo();
        return mean_X(m_h0.begin());
    }

    SurfPhase::
    ~SurfPhase() { }

    /**
     * For a surface phase, the pressure is not a relevant
     * thermodynamic variable, and so the enthalpy is equal to the
     * internal energy.
     */
    doublereal SurfPhase::
    intEnergy_mole() const { return enthalpy_mole(); }

    void SurfPhase::
    getStandardChemPotentials(doublereal* mu0) const {
        _updateThermo();
        copy(m_mu0.begin(), m_mu0.end(), mu0);
    }

    void SurfPhase::
    getActivityConcentrations(doublereal* c) const { 
        getConcentrations(c); 
    }

    doublereal SurfPhase::
    standardConcentration(int k) const {
        return m_n0/size(k); 
    }

    doublereal SurfPhase::
    logStandardConc(int k) const {
        return m_logn0 - m_logsize[k];
    }

    void SurfPhase::
    setParameters(int n, doublereal* c) {
        m_n0 = c[0];
	if (m_n0 <= 0.0) {
	  throw CanteraError("SurfPhase::setParameters",
			     "Bad value for parameter");
	}
        m_logn0 = log(m_n0);
    }

    void SurfPhase::
    initThermo() {
        m_h0.resize(m_kk);
        m_s0.resize(m_kk);
        m_cp0.resize(m_kk);
        m_mu0.resize(m_kk);
        m_work.resize(m_kk);
        m_pe.resize(m_kk, 0.0);
        vector_fp cov(m_kk, 0.0);
        cov[0] = 1.0;
        setCoverages(cov.begin());
        m_logsize.resize(m_kk);
        for (int k = 0; k < m_kk; k++) 
            m_logsize[k] = log(size(k));
    }

    void SurfPhase::
    setPotentialEnergy(int k, doublereal pe) {
        m_pe[k] = pe;
        _updateThermo(true);
    }

    void SurfPhase::
    setSiteDensity(doublereal n0) {
        doublereal x = n0;
        setParameters(1, &x);
    }


    //void SurfPhase::
    //setElectricPotential(doublereal V) {
    //    for (int k = 0; k < m_kk; k++) {
    //        m_pe[k] = charge(k)*Faraday*V;
    //    }
    //    _updateThermo(true);
    //}

    /**
     * Set the coverage fractions to a specified 
     * state. This routine converts to concentrations
     * in kmol/m2, using m_n0, the surface site density,
     * and size(k), which is defined to be the number of
     * surface sites occupied by the kth molecule.
     * It then calls State::setConcentrations to set the
     * internal concentration in the object.
     */
    void SurfPhase::
    setCoverages(const doublereal* theta) {
        for (int k = 0; k < m_kk; k++) {
            m_work[k] = m_n0*theta[k]/size(k);
        }
	/*
	 * Call the State:: class function
	 * setConcentrations.
	 */
        setConcentrations(m_work.begin());
    }

    void SurfPhase::
    getCoverages(doublereal* theta) const {
        getConcentrations(theta);
        for (int k = 0; k < m_kk; k++) {
            theta[k] *= size(k)/m_n0; 
        }
    }

    void SurfPhase::
    setCoveragesByName(string cov) {
        int kk = nSpecies();
        int k;
        compositionMap cc;
        for (k = 0; k < kk; k++) { 
            cc[speciesName(k)] = -1.0;
        }
        parseCompString(cov, cc);
        doublereal c;
        vector_fp cv(kk, 0.0);
        for (k = 0; k < kk; k++) { 
            c = cc[speciesName(k)];
            if (c > 0.0) cv[k] = c;
        }
        setCoverages(cv.begin());
    }


    void SurfPhase::
    _updateThermo(bool force) const {
        doublereal tnow = temperature();
        if (m_tlast != tnow || force) {
            m_spthermo->update(tnow, m_cp0.begin(), m_h0.begin(), 
                m_s0.begin());
            m_tlast = tnow;
            doublereal rt = GasConstant * tnow;
            int k;
            //doublereal deltaE;
            for (k = 0; k < m_kk; k++) {
                m_h0[k] *= rt;
                m_s0[k] *= GasConstant;
                m_cp0[k] *= GasConstant;
                //deltaE = m_pe[k];
                //m_h0[k] += deltaE;
                m_mu0[k] = m_h0[k] - tnow*m_s0[k];
            }
            m_tlast = tnow;
        }
    }

}
