/**
 *  @file PureFluid.h
 *
 * Declares class PureFluid
 */

// Copyright 2003  California Institute of Technology

#ifndef CT_EOS_TPX_H
#define CT_EOS_TPX_H

#include "ThermoPhase.h"

#ifdef INCL_PURE_FLUIDS

#include "../../ext/tpx/Sub.h"
#include "../../ext/tpx/utils.h"

namespace Cantera {
    
    /// Class for single-component fluids
    class PureFluid  : public ThermoPhase {

    public:

        PureFluid() : ThermoPhase(), m_sub(0) {}

        virtual ~PureFluid() { delete m_sub; }
        
        virtual void setParameters(int n, doublereal* c) {
            if (n == 1) {
                int subflag = int(c[0]);
                if (m_sub) delete m_sub;
                m_sub = tpx::GetSub(subflag);
                if (m_sub == 0) {
                    throw CanteraError("PureFluid::setParameters",
                        "could not create new substance object.");
                }
                m_subflag = subflag;
                m_mw = m_sub->MolWt();
                double cp0_R, h0_RT, s0_R, T0, p;
                T0 = 298.15;
                if (T0 < m_sub->Tcrit()) {
                    m_sub->Set(tpx::TX, T0, 1.0);
                    p = 0.01*m_sub->P();
                }
                else {
                    p = 0.001*m_sub->Pcrit();
                }
                m_sub->Set(tpx::TP, T0, p);
                
                m_spthermo->update_one(0, T0, &cp0_R, &h0_RT, &s0_R);
                double s_R = s0_R - log(p/refPressure());
                m_sub->setStdState(h0_RT*GasConstant*298.15/m_mw,
                    s_R*GasConstant/m_mw, T0, p);
            }
        }
        
        virtual int eosType() const { return cPureFluid; }

        /**
         * Mixture molar enthalpy. Units: J/mol. 
         */
        virtual doublereal enthalpy_mole() const {
            setTPXState();
            doublereal h = m_sub->h() * m_mw;
            check();
            return h;            
        }
        
        /**
         * Mixture molar internal energy. Units: J/mol. 
         */
        virtual doublereal intEnergy_mole() const {
            setTPXState();
            doublereal u = m_sub->u() * m_mw;
            check();
            return u;            
        }

        /**
         * Mixture molar entropy. Units: J/mol/K. 
         */
        virtual doublereal entropy_mole() const {
            setTPXState();
            doublereal s = m_sub->s() * m_mw;
            check();
            return s;            
        }

        /**
         * Mixture molar Gibbs function. Units: J/mol. 
         */
        virtual doublereal gibbs_mole() const {
            setTPXState();
            doublereal g = m_sub->g() * m_mw;
            check();
            return g;            
        }

        /**
         * Mixture molar heat capacity at constant pressure.
         * Units: J/mol/K. 
         */
        virtual doublereal cp_mole() const {
            setTPXState();
            doublereal cp = m_sub->cp() * m_mw;
            check();
            return cp;            
        }

        /**
         * Mixture molar heat capacity at constant volume.
         * Units: J/mol/K. 
         */
        virtual doublereal cv_mole() const {
            setTPXState();
            doublereal cv = m_sub->cv() * m_mw;
            check();
            return cv;
        }
        
        
        /**
         * Pressure. Units: Pa
         */
        virtual doublereal pressure() const {
            setTPXState();
            doublereal p = m_sub->P();
            check();
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
            check();
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
            m_sub->Set(tpx::TP, t, 0.5*m_sub->Pcrit());
            doublereal ps = m_sub->Ps();
            m_sub->Set(tpx::TV,tsv,vsv);
            check();
            return ps;
            }
        
        virtual doublereal vaporFraction() const {
            setTPXState();
            doublereal x = m_sub->x();
            check();
            return x;
        }
        
        virtual void setState_satLiquid() {
            setTPXState();
            m_sub->Set(tpx::TX, temperature(), 0.0);
            setDensity(1.0/m_sub->v());            
            check();            
        }
        
        virtual void setState_satVapor() {        
            setTPXState();
            m_sub->Set(tpx::TX, temperature(), 1.0);
            setDensity(1.0/m_sub->v());
            check();
        }        

protected:
        
        void setTPXState() const {
            m_sub->Set(tpx::TV, temperature(), 1.0/density());            
        }
        
        void check() const {
            if (m_sub->Error()) {
                throw CanteraError("PureFluidPhase",string(tpx::errorMsg(
                                                               m_sub->Error())));
            }
        }
        
private:
        mutable tpx::Substance* m_sub;
        int m_subflag;
        doublereal m_mw;
    };

}

#endif
#endif





