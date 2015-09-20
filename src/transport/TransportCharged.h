#ifndef TRANSPORTCHARGED_H
#define TRANSPORTCHARGED_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include "StringFunct.h"
#include "cantera/numerics/DenseMatrix.h"


using namespace std;

namespace Cantera {

class TransportCharged
{
    public:
        TransportCharged();
	void load_data();
        double getDebyeLength(double T, double Te, double Xe, double P);
        double getNumberDensity(double T, double Te, double Xe, double P);
	std::vector<vector_fp> numberDensity_fit;
};

};

#endif // TRANSPORTCHARGED_H


