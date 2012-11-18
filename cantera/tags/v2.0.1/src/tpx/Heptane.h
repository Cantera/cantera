#ifndef TPX_HEPTANE_H
#define TPX_HEPTANE_H

#include "cantera/tpx/Sub.h"



/* FILE: Heptane.h
 * DESCRIPTION:
 *  representation of substance Heptane
 *  values and functions are from
 *  "Thermodynamic Properties in SI" bu W.C. Reynolds
 * AUTHOR: me@rebeccahhunt.com: GCEP, Stanford University
 * AUTHOR: jrh@stanford.edu: GCEP, Stanford University
 *
 */
namespace tpx
{

class Heptane : public Substance
{
public:
    Heptane() {
        m_name = "Heptane";
        m_formula = "C7H16";
    }
    virtual ~Heptane() {}

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
    double C(int jm, double, double, double, double);
    double Cprime(int i,  double, double, double);
    double I(int i,  double, double);
    double H(int i, double egrho);
};

}

#endif // ! TPX_HEPTANE_H


