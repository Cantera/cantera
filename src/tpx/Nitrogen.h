
#ifndef TPX_NITROGEN_H
#define TPX_NITROGEN_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

class nitrogen : public Substance
{
public:
    nitrogen() {
        m_name = "nitrogen";
        m_formula = "N2";
    }
    virtual ~nitrogen() {}

    double MolWt();
    double Tcrit();
    double Pcrit();
    double Vcrit();
    double Tmin();
    double Tmax();
    char* name();
    char* formula();

    double Pp();
    double up();
    double sp();
    double Psat();

private:

    double ldens();
    double C(int i, double rt, double rt2);
    double Cprime(int i, double rt, double rt2, double rt3);
    double I(int i, double egrho);
    double H(int i, double egrho);
    double W(int i, double egrho);
};

}

#endif // ! TPX_NITROGEN_H

