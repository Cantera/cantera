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

#include "ct_defs.h"
#include "mix_defs.h"
#include "ThermoPhase.h"
#include "ctvector.h"
#include <iostream>

namespace Cantera {


    /**
     * A simple model for a surface phase. The surface consists of a
     * grid of equivalent sites. Surface species may be defined that
     * occupy one or more sites. The surface species are assumed to be
     * independent, and thus the species form an ideal solution.
     */
    class SurfPhase : public ThermoPhase  {

    public:

        SurfPhase(doublereal n0 = 0.0);
        virtual ~SurfPhase() {}
        virtual int eosType() const { return cSurf; }
        virtual doublereal enthalpy_mole() const;
        virtual doublereal intEnergy_mole() const;
        virtual void getStandardChemPotentials(doublereal* mu0) const;
        virtual void getActivityConcentrations(doublereal* c) const;
        virtual doublereal standardConcentration(int k = 0) const;
        virtual doublereal logStandardConc(int k=0) const;
        virtual void setParameters(int n, doublereal* c);
        virtual void initThermo();
        doublereal siteDensity(){ return m_n0; }
        void setPotentialEnergy(int k, doublereal pe);
        doublereal potentialEnergy(int k) {return m_pe[k];}
        void setSiteDensity(doublereal n0);
        void setElectricPotential(doublereal V);
        void setCoverages(const doublereal* theta);
        void getCoverages(doublereal* theta) const;

    protected:

        doublereal m_n0, m_logn0;
        doublereal m_tmin, m_tmax;

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





