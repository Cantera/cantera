/**
 *  @file PureFluid.h
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
                      m_mw(-1.0), m_verbose(true) {}

        virtual ~PureFluid() { delete m_sub; }
        
        
        virtual int eosType() const { return cPureFluid; }

        /**
         * Mixture molar enthalpy. Units: J/mol. 
         */
        virtual doublereal enthalpy_mole() const {
            setTPXState();
            doublereal h = m_sub->h() * m_mw;
            check(h);
            return h;            
        }
        
        /**
         * Mixture molar internal energy. Units: J/mol. 
         */
        virtual doublereal intEnergy_mole() const {
            setTPXState();
            doublereal u = m_sub->u() * m_mw;
            check(u);
            return u;            
        }

        /**
         * Mixture molar entropy. Units: J/mol/K. 
         */
        virtual doublereal entropy_mole() const {
            setTPXState();
            doublereal s = m_sub->s() * m_mw;
            check(s);
            return s;            
        }

        /**
         * Mixture molar Gibbs function. Units: J/mol. 
         */
        virtual doublereal gibbs_mole() const {
            setTPXState();
            doublereal g = m_sub->g() * m_mw;
            check(g);
            return g;            
        }

        /**
         * Mixture molar heat capacity at constant pressure.
         * Units: J/mol/K. 
         */
        virtual doublereal cp_mole() const {
            setTPXState();
            doublereal cp = m_sub->cp() * m_mw;
            check(cp);
            return cp;            
        }

        /**
         * Mixture molar heat capacity at constant volume.
         * Units: J/mol/K. 
         */
        virtual doublereal cv_mole() const {
            setTPXState();
            doublereal cv = m_sub->cv() * m_mw;
            check(cv);
            return cv;
        }
        
        
        /**
         * Pressure. Units: Pa
         */
        virtual doublereal pressure() const {
            setTPXState();
            doublereal p = m_sub->P();
            check(p);
            return p;
        }
        
        /**
         * Set the pressure, holding temperature and composition
         * fixed.  
         */
        virtual void setPressure(doublereal p) {
            m_sub->Set(tpx::TP, temperature(), p);
            setDensity(1.0/m_sub->v());
            check();
        }

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
            doublereal ts = m_sub->Tsat(p);
            check(ts);
            return ts;
        }
        
        virtual void setState_HP(doublereal h, doublereal p, 
            doublereal tol = 1.e-8) {
            m_sub->Set(tpx::HP, h, p);
            setState_TR(m_sub->Temp(), 1.0/m_sub->v());
            check();
        }

        virtual void setState_UV(doublereal u, doublereal v, 
            doublereal tol = 1.e-8) {
            m_sub->Set(tpx::UV, u, v);
            setState_TR(m_sub->Temp(), 1.0/m_sub->v());
            check();
        }

        virtual void setState_SV(doublereal s, doublereal v, 
            doublereal tol = 1.e-8) {
            m_sub->Set(tpx::SV, s, v);
            setState_TR(m_sub->Temp(), 1.0/m_sub->v());
            check();
        }

        virtual void setState_SP(doublereal s, doublereal p, 
            doublereal tol = 1.e-8) {
            m_sub->Set(tpx::SP, s, p);
            setState_TR(m_sub->Temp(), 1.0/m_sub->v());
            check();
        }

        /// saturation pressure
        virtual doublereal satPressure(doublereal t) const {
            doublereal tsv = m_sub->Temp();
            doublereal vsv = m_sub->v();
            if (t < 0.0)
                m_sub->Set(tpx::TP, temperature(), 0.5*m_sub->Pcrit());
            else
                m_sub->Set(tpx::TP, t, 0.5*m_sub->Pcrit());

            doublereal ps = m_sub->Ps();
            m_sub->Set(tpx::TV,tsv,vsv);
            check(ps);
            return ps;
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
            m_sub->Set(tpx::TX, t, x);
            setDensity(1.0/m_sub->v());
            check();
        }

        virtual void setState_Psat(doublereal p, doublereal x) {
            setTPXState();
            m_sub->Set(tpx::PX, p, x);
            setTemperature(m_sub->Temp());
            setDensity(1.0/m_sub->v());
            check();
        }

        virtual void initThermo();
        virtual void setParametersFromXML(const XML_Node& eosdata);

protected:
        
        void setTPXState() const {
            m_sub->Set(tpx::TV, temperature(), 1.0/density());            
        }
        
        void check(doublereal v = 0.0) const {
            if (m_sub->Error() || v == tpx::Undef) {
                throw CanteraError("PureFluidPhase",string(tpx::errorMsg(
                                                               m_sub->Error())));
            }
        }
        
private:
        mutable tpx::Substance* m_sub;
        int m_subflag;
        doublereal m_mw;
        bool m_verbose;
    };

}

#endif
#endif





