// Lee-Kesler equation of state

#include "RedlichKwong.h"
#include <math.h>

namespace tpx {


    //--------------------------- member functions ------------------

    double RedlichKwong::up() {
	double u = -Pp()/Rho + hresid() + m_energy_offset;
        //cout << "up = " << u << endl;
        return u;
    }

    double RedlichKwong::hresid(){
        double hh = m_b * (Rho/m_mw);
	double hresid_mol_RT =  z() - 1.0 
                             - (1.5*m_a/(m_b*8314.3*pow(T,1.5)))*log(1.0 + hh);
        return 8314.3*T*hresid_mol_RT/m_mw;
    }

    double RedlichKwong::sresid(){
        double hh = m_b * (Rho/m_mw);
        //cout << "hh = " << hh << endl;
	double sresid_mol_R =  log(z()*(1.0 - hh)) 
                             - (0.5*m_a/(m_b*8314.3*pow(T,1.5)))*log(1.0 + hh);
        double sp = 8314.3*sresid_mol_R/m_mw;
        //cout << "sresid = " << sp << endl;
        return sp;
    }

    double RedlichKwong::sp() {
	const double Pref = 101325.0;
	double rgas = 8314.3/m_mw;
        //cout << "P = " << Rho*rgas*T << endl;
        double ss = rgas*(log(Pref/(Rho*rgas*T)));
        double sr = sresid(); 
        double p = Pp();
	double s = rgas*(log(Pref/p)) + sr + m_entropy_offset;
        //cout << "sp = " << s << " " << ss << " " << sr << " " << m_entropy_offset << endl;
        return s;
    }

    double RedlichKwong::z() {
        return Pp()*m_mw/(Rho*8314.3*T);
    }


    double RedlichKwong::Pp() {
	double R = 8314.3;
        double V = m_mw/Rho;
        double pp = R*T/(V - m_b) - m_a/(sqrt(T)*V*(V+m_b));
        //cout << "molar V, T, P = " << V << " " << T << " " << pp << endl;
        return pp;
        //cout << "Rho, T, Pp = " << pp << endl;
    }

    double RedlichKwong::Psat(){
	double tt = m_tcrit/T;
        double lpr = -0.8734*tt*tt - 3.4522*tt + 4.2918;
	return m_pcrit*exp(lpr);
    }

    double RedlichKwong::ldens(){
        double c;
        int i;
        double sqt = sqrt(T);
        double v = m_b, vnew;
        double pp = Psat();
        double Rhsave = Rho;
        for (i = 0; i < 50; i++) {
            //pp = Pp();
            c = m_b*m_b + m_b*GasConstant*T/pp - m_a/(pp*sqt);
            vnew = (1.0/c)*(v*v*v - GasConstant*T*v*v/pp - m_a*m_b/(pp*sqt));
            v = vnew;
            //Rho = m_mw/v;
            //cout << "ldens Rho = " << Rho << "  " << z() << "  " << Pp() << endl;
        }
        Rho = Rhsave;
        //cout << "ldens: " << m_mw/vnew << endl;
        return m_mw/vnew;
    }


}
