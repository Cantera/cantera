/**
 *
 *  @file IncompressibleThermo.h
 *
 *
 *  $Author: hkmoffa $
 *  $Date: 2003/08/21 22:39:16 $
 *  $Revision: 1.2 $
 *
 *  Copyright 2001 California Institute of Technology
 *
 */


#ifndef CT_INCOMPTHERMO_H
#define CT_INCOMPTHERMO_H

#include "ct_defs.h"
#include "mix_defs.h"
#include "Thermo.h"

namespace Cantera {

    /**
     * Overloads the virtual methods of class Thermo to implement the
     * ideal gas equation of state.
     */
    class IncompressibleThermo : public Thermo  {

    public:

        IncompressibleThermo(phase_t* phase=0) 
            : Thermo(phase), m_press(-1.0) {}
        virtual ~IncompressibleThermo() {}

        virtual bool eosType() const { return cIncompressible; }

        virtual doublereal enthalpy_mole() const {
            return GasConstant * m_s->temperature() * 
                m_s->mean_X(m_s->enthalpy_RT()) 
                + (pressure() - m_s->refPressure())/m_s->molarDensity();
        }

        virtual doublereal intEnergy_mole() const {
            return enthalpy_mole() - pressure()/m_s->molarDensity();
        }

        virtual doublereal entropy_mole() const {
            return GasConstant * (m_s->mean_X(m_s->entropy_R()) -
                m_s->sum_xlogx());
        }

        virtual doublereal gibbs_mole() const {
            return enthalpy_mole() - m_s->temperature() * entropy_mole();
        }

        virtual doublereal cp_mole() const {
            return GasConstant * m_s->mean_X(m_s->cp_R());
        }

        virtual doublereal cv_mole() const {
            return cp_mole();
        }

        virtual doublereal pressure() const {
            return m_press;
        }

        virtual void setPressure(doublereal p) {
            m_press = p;
        }

        virtual void getChemPotentials_RT(doublereal* mu) const {
            int kk = m_s->nSpecies();
            doublereal xx;
            const vector_fp& g_RT = m_s->gibbs_RT();
            doublereal vhat_dp = (pressure() - m_s->refPressure())/
                                 m_s->molarDensity();
            for (int k = 0; k < kk; k++) {
                xx = fmaxx(Tiny,m_s->moleFraction(k));
                mu[k] = g_RT[k] + log(xx);
            }
        }

    protected:
        doublereal m_press;
    private:

    };
}
        
#endif





