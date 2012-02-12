#ifndef TPX_CARBONDIOXIDE_H
#define TPX_CARBONDIOXIDE_H

#include "cantera/tpx/Sub.h"



/* FILE: CarbonDioxide.h
 * DESCRIPTION:
 *  representation of substance Carbon Dioxide
 *  values and functions are from
 *  "Thermodynamic Properties in SI" bu W.C. Reynolds
 * AUTHOR: me@rebeccahhunt.com: GCEP, Stanford University
 *
 */
namespace tpx
{

class CarbonDioxide : public Substance
{
public:
    CarbonDioxide() :
        Substance() {
        m_name="CarbonDioxide";
        m_formula="CO2";
    }
    virtual ~CarbonDioxide() {}

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

#endif // ! TPX_CARBONDIOXIDE_H


