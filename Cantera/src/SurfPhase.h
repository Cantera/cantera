/**
 *
 *  @file SurfPhase.h
 *
 */

/*  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */


#ifndef CT_SURFPHASE_H
#define CT_SURFPHASE_H

#include "mix_defs.h"
#include "ThermoPhase.h"

namespace Cantera {


    /**
     * A simple model for a surface phase. The surface consists of a
     * grid of equivalent sites. Surface species may be defined that
     * occupy one or more sites. The surface species are assumed to be
     * independent, and thus the species form an ideal solution.
     * The definitions of the member functions are located in 
     * InterfaceKinetics.cpp.
     */
    class SurfPhase : public ThermoPhase  {

    public:

        SurfPhase(doublereal n0 = 0.0);
        virtual ~SurfPhase();
        virtual int eosType() const { return cSurf; }
        virtual doublereal enthalpy_mole() const;
        virtual doublereal intEnergy_mole() const;
        virtual void getStandardChemPotentials(doublereal* mu0) const;
        virtual void getChemPotentials(doublereal* mu) const;
        virtual void getActivityConcentrations(doublereal* c) const;
        virtual doublereal standardConcentration(int k = 0) const;
        virtual doublereal logStandardConc(int k=0) const;
        virtual void setParameters(int n, doublereal* c);
        virtual void initThermo();
        doublereal siteDensity(){ return m_n0; }
        void setPotentialEnergy(int k, doublereal pe);
        doublereal potentialEnergy(int k) {return m_pe[k];}
        void setSiteDensity(doublereal n0);
        //void setElectricPotential(doublereal V);

        void getEnthalpy_RT(doublereal* hrt) const;
        void getEntropy_R(doublereal* sr) const;

	/**
         * Pressure. Units: Pa.
         */ 
        virtual doublereal pressure() const {
            return m_press;
        }
	
        /**
         * Set the pressure at constant temperature. Units: Pa.
         */
        virtual void setPressure(doublereal p) {
            m_press = p;
        }

        /**
         * Set the surface site fractions to a specified 
         * state. This routine converts to concentrations
         * in kmol/m2, using m_n0, the surface site density,
         * and size(k), which is defined to be the number of
         * surface sites occupied by the kth molecule.
         * It then calls State::setConcentrations to set the
         * internal concentration in the object.
	 *
	 * @param theta[k] This is the surface site fraction
	 *                 for the kth species in the surface phase.
	 *                 This is a dimensionless quantity.
         */
        void setCoverages(const doublereal* theta);
        void setCoveragesNoNorm(const doublereal* theta);
        void setCoveragesByName(string cov);
        void getCoverages(doublereal* theta) const;

    protected:

        doublereal m_n0;
        doublereal m_logn0;
        doublereal m_tmin, m_tmax;
	doublereal m_press;

        mutable doublereal    m_tlast;
        mutable array_fp      m_h0;
        mutable array_fp      m_s0;
        mutable array_fp      m_cp0;
        mutable array_fp      m_mu0;
        mutable array_fp      m_work;
        mutable array_fp      m_pe;
        mutable array_fp      m_logsize;

    private:

        void _updateThermo(bool force=false) const;

    };
}
        
#endif





