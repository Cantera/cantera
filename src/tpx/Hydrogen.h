#ifndef TPX_HYDROGEN_H
#define TPX_HYDROGEN_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

class hydrogen : public Substance
{
public:
    hydrogen() {
        m_name = "hydrogen";
        m_formula = "H2";
    }
    virtual ~hydrogen() {}

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
    double icv(int i, double x, double xlg);
};

}

#endif // ! HYDROGEN_H
