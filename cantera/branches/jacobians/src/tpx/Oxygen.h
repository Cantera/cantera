#ifndef TPX_OXYGEN_H
#define TPX_OXYGEN_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

class oxygen : public Substance
{
public:
    oxygen() {
        m_name="oxygen";
        m_formula="O2";
    }
    virtual ~oxygen() {}

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
#endif // ! OXYGEN_H

