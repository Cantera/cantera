#ifndef CT_MMCOLLISIONINTCHARGED_H
#define CT_MMCOLLISIONINTCHARGED_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <iomanip>
#include "cantera/base/ct_defs.h"
#include "StringFunct.h"
#include "TransportReadExt.h"
#include "TransportCharged.h"


using namespace std;

namespace Cantera
{


//! Calculation of Collision integrals
class MMCollisionIntCharged
{
public:
    MMCollisionIntCharged();
    virtual ~MMCollisionIntCharged();

    doublereal omega11_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega22_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega12_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega13_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega14_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega15_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega23_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega24_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal bstar_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);

    vector<double> C_MMomega11;
    vector<double> C_MMomega22;
    vector<string> PairsMM;
    vector<double> C_MMbstar;

};
}
#endif
