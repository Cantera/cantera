/**
 *  @file EOS.h
 *
 * Declares virtual base class EOS
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_EOS_TPX_H
#define CT_EOS_TPX_H

#include "EOS.h"
#include "../ext/tpx/Sub.h"
#include "../ext/tpx/utils.h"

namespace Cantera {

    class TPX_Error {
    public:
        TPX_Error(string proc, int err, int fatal = 1) {
            cerr << "Error in EOS_TPX::" << proc << ": "
                 << tpx::errorMsg(err) << endl;
            if (fatal) exit(-1);
        }
    };


    class EOS_TPX  : public EOS {

    public:

        EOS_TPX(int subflag, double h0 = 0.0, double s0 = 0.0) {
            m_sub = tpx::GetSub(subflag);
            m_mw = m_sub->MolWt();
            m_sub->setStdState(h0/m_mw, s0/m_mw);
        }

        virtual ~EOS_TPX() { delete m_sub; }


        /**
         * Mixture molar enthalpy. Units: J/mol. 
         */
        virtual doublereal enthalpy_mole(const  State& s, 
            const vector_fp& h0_RT) const {
            m_sub->Set(tpx::TV, s.temperature(), 1.0/s.density());
            return m_sub->h() * m_mw;
        }
        /**
         * Mixture molar internal energy. Units: J/mol. 
         */
        virtual doublereal intEnergy_mole(const State& s, 
            const vector_fp& h0_RT) const {
            m_sub->Set(tpx::TV, s.temperature(), 1.0/s.density());
            return m_sub->u() * m_mw;
        }

        /**
         * Mixture molar entropy. Units: J/mol/K. 
         */
        virtual doublereal entropy_mole(const State& s, 
            const vector_fp& s0_RT,
            doublereal log_Pp_bar = -999.0 ) const {
            m_sub->Set(tpx::TV, s.temperature(), 1.0/s.density());
            return m_sub->s() * m_mw;
        }

        /**
         * Mixture molar Gibbs function. Units: J/mol. 
         */
        virtual doublereal gibbs_mole(const State& s, 
            const vector_fp& g0_RT, 
            doublereal log_Pp_bar = -999.0 ) const {
            m_sub->Set(tpx::TV, s.temperature(), 1.0/s.density());
            return m_sub->g() * m_mw;
        }

        /**
         * Mixture molar heat capacity at constant pressure.
         * Units: J/mol/K. 
         */
        virtual doublereal cp_mole(const State& s, 
            const vector_fp& cp0_R ) const {
            m_sub->Set(tpx::TV, s.temperature(), 1.0/s.density());
            return m_sub->cp() * m_mw;
        }

        /**
         * Mixture molar heat capacity at constant volume.
         * Units: J/mol/K. 
         */
        virtual doublereal cv_mole(const State& s, 
            const vector_fp& cp0_R ) const {
            m_sub->Set(tpx::TV, s.temperature(), 1.0/s.density());
            return m_sub->cv() * m_mw;
        }

        /**
         * Mixture molar isothermal compressibility
         * \f$ -(1/V)(\partial V/\partial P)_T\f$.  Units: 1/Pa.
         */
        virtual doublereal compressibility_T(const State& s, 
            const vector_fp& cp0_R ) const {
            m_sub->Set(tpx::TV, s.temperature(), 1.0/s.density());
            return -999.0;
        }

        /**
         * Mixture molar volumetric thermal expansion coefficient
         * \f$ (1/V)(\partial V/\partial T)_P\f$.  Units: 1/K.
         */
        virtual doublereal thermalExpansionCoeff(const State& s, 
            const vector_fp& cp0_R ) const {
            m_sub->Set(tpx::TV, s.temperature(), 1.0/s.density());
            doublereal beta = m_sub->thermExpCoeff();
            if (m_sub->Error()) 
                throw TPX_Error("thermalExpansionCoeff", m_sub->Error());
            return beta;
        }


        /**
         * Pressure. Units: Pa
         */
        virtual doublereal pressure(const State& s) const {
            m_sub->Set(tpx::TV, s.temperature(), 1.0/s.density());
            doublereal pp = m_sub->P();
            if (m_sub->Error()) 
                throw TPX_Error("pressure", m_sub->Error());
            return pp;
        }


        /**
         * Set the pressure, holding temperature and composition
         * fixed.  
         *
         * @param s State instance defining the thermodynamic state.
         * The density attribute of s will be set to a value such that
         * the value of pressure() equals p.
         *  
         * @param p Pressure in Pa.
         */
        virtual void setPressure(State& s, doublereal p) const {
            m_sub->Set(tpx::TP, s.temperature(), p);
            s.setDensity(1.0/m_sub->v());
            if (m_sub->Error()) 
                throw TPX_Error("setPressure", m_sub->Error());
        }

        virtual void getChemPotentials_RT(const State& s, 
            const vector_fp& g0_RT, const vector_fp& x, 
            doublereal* mu) const {
            m_sub->Set(tpx::TV, s.temperature(), 1.0/s.density());
            mu[0] = gibbs_mole(s, g0_RT);
            if (m_sub->Error()) 
                throw TPX_Error("getChemPotentials_RT", m_sub->Error());
        }

        tpx::Substance& TPX_Substance() { return *m_sub; }
        doublereal Tmin() { return m_sub->Tmin(); }
        doublereal Tmax() { return m_sub->Tmax(); }

        /// critical state properties 
        virtual doublereal critTemperature() { return m_sub->Tcrit(); }
        virtual doublereal critPressure() { return m_sub->Pcrit(); }
        virtual doublereal critDensity() { return 1.0/m_sub->Vcrit(); }

        /// saturation properties
        virtual doublereal satTemperature(doublereal p) { 
            doublereal ts = m_sub->Tsat(p);
            if (ts == tpx::Undef) throw TPX_Error("satTemperature",m_sub->Error());
            return ts;
        }
        virtual doublereal satPressure(doublereal t) {
            doublereal tsv = m_sub->Temp();
            doublereal vsv = m_sub->v();
            m_sub->Set(tpx::TP, t, 0.5*m_sub->Pcrit());
            doublereal ps = m_sub->Ps();
            if (ps == tpx::Undef) throw TPX_Error("satPressure",m_sub->Error());
            m_sub->Set(tpx::TV,tsv,vsv);
            return ps;
            }

        virtual int phase() {
            doublereal xx = m_sub->x();
            if (xx > 0.99999) 
                return Vapor_Phase;
            else if (xx < 1.e-5)
                return Liquid_Phase;
            else
                return Liquid_Phase + Vapor_Phase;
        }


    protected:

        tpx::Substance* m_sub;
        doublereal m_mw;
    };

}
        
#endif





