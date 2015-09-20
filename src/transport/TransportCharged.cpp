#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include "TransportCharged.h"
#include "StringFunct.h"
#include "cantera/numerics/polyfit.h"

using namespace std;

namespace Cantera {

TransportCharged::TransportCharged()
{	
}


void TransportCharged::load_data()
{
}



double TransportCharged::getDebyeLength(double T, double Te, double Xe, double P)
{

        const double Kb = 1.3806488*pow(10, -23);		// Boltzmann's constant (J/molecule-K)
	const double eps0 = 8.854187817*pow(10, -12);		// Vacuum permittivity (F/m)
	const double Qe = 1.602176565*pow(10, -19);		// Elementary positive charge (C)
	const double PI = 3.141592653;

	double debyeLength;
	double f;


	// Average closest impact parameters
        const double bfac = Qe * Qe / (8.0 * PI * eps0 * Kb);   //e-ion
        const double be   = bfac / Te;   // electron-electron
        const double bh   = bfac / T;   // ion-ion

        double number = 10000*(be+bh);

	double numberDensity;

	numberDensity = getNumberDensity(T, Te, Xe, P);

	f = (eps0*Kb*Te)/(2*numberDensity*Xe*Qe*Qe);

	debyeLength = pow(f,0.5);


	return min(debyeLength, number);

}


double TransportCharged::getNumberDensity(double T, double Te, double Xe, double P)
{
        const double Kb = 1.3806488*pow(10, -23);

	double numberDensity;

	numberDensity = P/(Kb* ( ((1-Xe)*T) + (Xe*Te)));

	return numberDensity;

}


}

