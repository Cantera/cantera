#ifndef TPX_METHANE_H
#define TPX_METHANE_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

class methane : public Substance
{
public:
    methane() {
        m_name = "methane";
        m_formula = "CH4";
    }
    virtual ~methane() {}

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
#endif // ! METHANE_H
