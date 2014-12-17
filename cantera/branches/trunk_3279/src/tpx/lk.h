//! @file lk.h Lee-Kesler equation of state
#ifndef TPX_LK_H
#define TPX_LK_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

class leekesler : public Substance
{

public:

    leekesler(double tc = 1.0, double pc = 1.0,
              double wt = 1.0, int itype = 0) {
        Tcr = tc;
        Pcr = pc;
        Mw = wt;
        Isr = itype;   // simple fluid or reference
        m_name = "Lee-Kesler";
        m_formula = "---";
    }

    double MolWt();
    double Tcrit();
    double Pcrit();
    double Vcrit();
    double Tmin();
    double Tmax();

    double Pp();
    double up();
    double sp();
    double Psat();
    double dPsatdT();

    // compressibility
    double z();

    // enthalpy departure
    double hdep();

    // entropy departure
    double sdep();

    double ldens();

protected:

    double Tcr, Pcr, Mw;
    int Isr;

private:

    double W(int n, double egrho, double gamma);
    double I();
    double J();
};
}

#endif // ! TPX_LK_H
