/**
 *  @file PureFluidPhase.h
 *
 * Declares class PureFluid
 */

/*  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2003 California Institute of Technology
 */

#ifndef CT_EOS_TPX_H
#define CT_EOS_TPX_H

#include "ThermoPhase.h"

#ifdef INCL_PURE_FLUIDS

#include "mix_defs.h"
#include "../../ext/tpx/Sub.h"
#include "../../ext/tpx/utils.h"

namespace Cantera {
    
    /// Class for single-component fluids
    class PureFluid  : public ThermoPhase {

    public:

        PureFluid() : ThermoPhase(), m_sub(0), m_subflag(0), 
                      m_mw(-1.0), m_verbose(false) {}

        virtual ~PureFluid() { delete m_sub; }
        
        
        virtual int eosType() const { return cPureFluid; }

        virtual doublereal enthalpy_mole() const;
        virtual doublereal intEnergy_mole() const;
        virtual doublereal entropy_mole() const;
        virtual doublereal gibbs_mole() const;
        virtual doublereal cp_mole() const;
        virtual doublereal cv_mole() const;
        virtual doublereal pressure() const;
        virtual void setPressure(doublereal p);

        virtual void getChemPotentials(doublereal* mu) const {
            mu[0] = gibbs_mole();
        }

        virtual doublereal isothermalCompressibility() {
            return m_sub->isothermalCompressibility();
        }

        virtual doublereal thermalExpansionCoeff() {
            return m_sub->thermalExpansionCoeff();
        }

        tpx::Substance& TPX_Substance() { return *m_sub; }

        /// critical temperature 
        virtual doublereal critTemperature() const { return m_sub->Tcrit(); }
        
        /// critical pressure
        virtual doublereal critPressure() const { return m_sub->Pcrit(); }
        
        /// critical density
        virtual doublereal critDensity() const { return 1.0/m_sub->Vcrit(); }
        
        
        /// saturation temperature
        virtual doublereal satTemperature(doublereal p) const { 
            try {
                doublereal ts = m_sub->Tsat(p);
                return ts;
            }
            catch(tpx::TPX_Error) {
                reportTPXError();
                return -1.0;
            }
        }
        
        virtual void setState_HP(doublereal h, doublereal p, 
            doublereal tol = 1.e-8) {
            Set(tpx::HP, h, p);
            setState_TR(m_sub->Temp(), 1.0/m_sub->v());
            check();
        }

        virtual void setState_UV(doublereal u, doublereal v, 
            doublereal tol = 1.e-8) {
            Set(tpx::UV, u, v);
            setState_TR(m_sub->Temp(), 1.0/m_sub->v());
            check();
        }

        virtual void setState_SV(doublereal s, doublereal v, 
            doublereal tol = 1.e-8) {
            Set(tpx::SV, s, v);
            setState_TR(m_sub->Temp(), 1.0/m_sub->v());
            check();
        }

        virtual void setState_SP(doublereal s, doublereal p, 
            doublereal tol = 1.e-8) {
            Set(tpx::SP, s, p);
            setState_TR(m_sub->Temp(), 1.0/m_sub->v());
            check();
        }

        /// saturation pressure
        virtual doublereal satPressure(doublereal t) const {
            //doublereal tsv = m_sub->Temp();
            doublereal vsv = m_sub->v();
            //if (t < 0.0)
            //    Set(tpx::TP, temperature(), 0.5*m_sub->Pcrit());
            //else
            //    Set(tpx::TP, t, 0.5*m_sub->Pcrit());
            try {
                Set(tpx::TV,t,vsv);
                doublereal ps = m_sub->Ps();
                //Set(tpx::TV,tsv,vsv);
                //check(ps);
                return ps;
            }
            catch(tpx::TPX_Error) {
                reportTPXError();
                return -1.0;
            }
        }
        
        virtual doublereal vaporFraction() const {
            setTPXState();
            doublereal x = m_sub->x();
            check(x);
            return x;
        }
        
        virtual void setState_Tsat(doublereal t, doublereal x) {
            setTemperature(t);
            setTPXState();
            Set(tpx::TX, t, x);
            setDensity(1.0/m_sub->v());
            check();
        }

        virtual void setState_Psat(doublereal p, doublereal x) {
            setTPXState();
            Set(tpx::PX, p, x);
            setTemperature(m_sub->Temp());
            setDensity(1.0/m_sub->v());
            check();
        }

        virtual void initThermo();
        virtual void setParametersFromXML(const XML_Node& eosdata);

protected:
        
        void Set(int n, double x, double y) const;
        void setTPXState() const;
        void check(doublereal v = 0.0) const;
        void reportTPXError() const;

private:
        mutable tpx::Substance* m_sub;
        int m_subflag;
        doublereal m_mw;
        bool m_verbose;
    };

}

#endif
#endif





