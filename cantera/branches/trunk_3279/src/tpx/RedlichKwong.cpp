//! @file RedlichKwong.cpp
#include "RedlichKwong.h"

namespace tpx
{
double RedlichKwong::up()
{
    return -Pp()/Rho + hresid() + m_energy_offset;
}

double RedlichKwong::hresid()
{
    double hh = m_b * (Rho/m_mw);
    double hresid_mol_RT =  z() - 1.0
                            - (1.5*m_a/(m_b*8314.3*T*sqrt(T)))*log(1.0 + hh);
    return 8314.3*T*hresid_mol_RT/m_mw;
}

double RedlichKwong::sresid()
{
    double hh = m_b * (Rho/m_mw);
    double sresid_mol_R =  log(z()*(1.0 - hh))
                           - (0.5*m_a/(m_b*8314.3*T*sqrt(T)))*log(1.0 + hh);
    return 8314.3*sresid_mol_R/m_mw;
}

double RedlichKwong::sp()
{
    const double Pref = 101325.0;
    double rgas = 8314.3/m_mw;
    double sr = sresid();
    double p = Pp();
    return rgas*(log(Pref/p)) + sr + m_entropy_offset;
}

double RedlichKwong::z()
{
    return Pp()*m_mw/(Rho*8314.3*T);
}

double RedlichKwong::Pp()
{
    double R = 8314.3;
    double V = m_mw/Rho;
    return R*T/(V - m_b) - m_a/(sqrt(T)*V*(V+m_b));
}

double RedlichKwong::Psat()
{
    double tt = m_tcrit/T;
    double lpr = -0.8734*tt*tt - 3.4522*tt + 4.2918;
    return m_pcrit*exp(lpr);
}

double RedlichKwong::ldens()
{
    double c;
    int i;
    double sqt = sqrt(T);
    double v = m_b, vnew;
    double pp = Psat();
    double Rhsave = Rho;
    for (i = 0; i < 50; i++) {
        c = m_b*m_b + m_b*GasConstant*T/pp - m_a/(pp*sqt);
        vnew = (1.0/c)*(v*v*v - GasConstant*T*v*v/pp - m_a*m_b/(pp*sqt));
        v = vnew;
    }
    Rho = Rhsave;
    return m_mw/vnew;
}

}
