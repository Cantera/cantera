/**
 *  @file PureFluidPhase.h
 *
 * Declares class PureFluid
 */

/*  $Author: hkmoffa $
 *  $Date: 2006/03/03 21:11:52 $
 *  $Revision: 1.13 $
 *
 *  Copyright 2003 California Institute of Technology
 */

#ifndef CT_EOS_TPX_H
#define CT_EOS_TPX_H

#include "ThermoPhase.h"

#ifdef WITH_PURE_FLUIDS

#include "mix_defs.h"

namespace tpx {
    class Substance;
}

namespace Cantera {


    /// Class for single-component fluids
    class PureFluidPhase  : public ThermoPhase {

    public:

        PureFluidPhase() : ThermoPhase(), m_sub(0), m_subflag(0), 
                      m_mw(-1.0), m_verbose(false) {}

        virtual ~PureFluidPhase();
        
        
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

        virtual doublereal isothermalCompressibility() const;
        virtual doublereal thermalExpansionCoeff() const;

        tpx::Substance& TPX_Substance();

        /// critical temperature 
        virtual doublereal critTemperature() const;
 
        /// critical pressure
        virtual doublereal critPressure() const;
        
        /// critical density
        virtual doublereal critDensity() const;
        
        /// saturation temperature
        virtual doublereal satTemperature(doublereal p) const;
        
        virtual void setState_HP(doublereal h, doublereal p, 
            doublereal tol = 1.e-8);

        virtual void setState_UV(doublereal u, doublereal v, 
            doublereal tol = 1.e-8);

        virtual void setState_SV(doublereal s, doublereal v, 
            doublereal tol = 1.e-8);

        virtual void setState_SP(doublereal s, doublereal p, 
            doublereal tol = 1.e-8);

        /// saturation pressure
        virtual doublereal satPressure(doublereal t) const;
        
        virtual doublereal vaporFraction() const;
        
        virtual void setState_Tsat(doublereal t, doublereal x);

        virtual void setState_Psat(doublereal p, doublereal x);

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





