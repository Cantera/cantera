//! @file RedlichKwong.h
#ifndef TPX_RK_H
#define TPX_RK_H

#include "cantera/tpx/Sub.h"
#include <math.h>

namespace tpx
{
const double GasConstant = 8314.3;

class RedlichKwong : public Substance
{

public:

    RedlichKwong() : Substance() {
        setParameters(1.0, 1.0, 1.0);
        m_name = "Redlich-Kwong";
        m_formula = "-";
    }

    void setParameters(double Tc, double Pc, double MolWt) {
        m_tcrit = Tc;
        m_pcrit = Pc;
        m_mw = MolWt;

        // compute the a and b parameters
        m_a = 0.42748*GasConstant*GasConstant*m_tcrit*m_tcrit*sqrt(m_tcrit)/m_pcrit;
        m_b = 0.08664*GasConstant*m_tcrit/m_pcrit;
    }

    double a() {
        return m_a;
    }
    double b() {
        return m_b;
    }

    double MolWt() {
        return m_mw;
    }
    double Tcrit() {
        return m_tcrit;
    }
    double Pcrit() {
        return m_pcrit;
    }
    double Vcrit() {
        return 0.3592725*GasConstant*T/(m_mw*m_pcrit);
    }
    double Tmin() {
        return 0.0;
    }
    double Tmax() {
        return 1.0e10;
    }
    char* name() {
        return (char*) m_name.c_str() ;
    }
    char* formula() {
        return (char*) m_formula.c_str() ;
    }

    double Pp();
    double up();
    double sp();
    double Psat();
    double dPsatdT();

    // compressibility
    double z();

    // enthalpy departure
    double hresid();

    // entropy departure
    double sresid();

    double ldens();

protected:

    double m_tcrit, m_pcrit, m_mw, m_a, m_b;
    //double m_tmin, m_tmax;
    //string m_name, m_formula;

private:

};
}

#endif // ! TPX_RK_H
