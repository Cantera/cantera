/**
 *  @file MMCollisionInt.cpp
 */
// Copyright 2001  California Institute of Technology */

#include "MMCollisionIntCharged.h"
#include "cantera/base/utilities.h"
#include "cantera/numerics/polyfit.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"
#include <math.h>
#include <vector>
#include "TransportReadExt.h"
#include "TransportCharged.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <iomanip>
#include <string>
#include <cstdio>

using namespace std;

namespace Cantera
{


MMCollisionIntCharged::MMCollisionIntCharged()
{

}

MMCollisionIntCharged::~MMCollisionIntCharged()
{
}


doublereal MMCollisionIntCharged::bstar_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P)
{

        const double Kb = 1.3806488*pow(10, -23);               // Boltzmann's constant (J/molecule-K)
        const double eps0 = 8.854187817*pow(10, -12);           // Vacuum permittivity (F/m)
        const double Qe = 1.602176565*pow(10, -19);             // Elementary positive charge (C)
        const double PI = 3.141592653;

	int charge = 0;
        double C_II[5] ={};      
        doublereal bstar;
        doublereal lambdaDeb;
        doublereal lambdaD;


	TransportCharged charged;
	lambdaDeb = charged.getDebyeLength(T, Te, Xe, P); 


        // Average closest impact parameters
        const double bfac = Qe * Qe / (8.0 * PI * eps0 * Kb);   //e-ion
        const double be   = bfac / Te;   // electron-electron
        const double bh   = bfac / T;   // ion-ion

        double number = 10000*(be+bh);
	charge = charge1*charge2;


        lambdaD = min(lambdaDeb, number);

        // Reduced temperatures
        const double Tste = max( (0.5 * lambdaD / be), 0.1);
        const double Tsth = max( (0.5 * lambdaD / bh), 0.1);

        // Factors used in the curve-fitting of the reduced collision integrals
        const double efac   = PI * lambdaD * lambdaD / (Tste * Tste);
        const double hfac   = PI * lambdaD * lambdaD / (Tsth * Tsth);
        const double lnTste = log(Tste);
        const double lnTsth = log(Tsth);

        //Repulsive interaction Ion-Ion   (REMARK: missing repulsive interaction ion-ion)
        C_II[0]=3.1204*pow(10, -6);
        C_II[1]=-0.0002626;
        C_II[2]=0.007332;
        C_II[3]=-0.0869;
        C_II[4]=1.4292;

        // repulsive interaction ion-ion
	if (charge > 0 )
                {
                        bstar = hfac * exp( C_II[0]*lnTsth*lnTsth*lnTsth + C_II[1]*lnTsth*lnTsth + C_II[2]*lnTsth + C_II[3]) ;
                }

	else
		{
			bstar = 1;
		}


                return bstar;

}



doublereal MMCollisionIntCharged::omega11_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P)
{

	const double Kb = 1.3806488*pow(10, -23);		// Boltzmann's constant (J/molecule-K)
	const double eps0 = 8.854187817*pow(10, -12);		// Vacuum permittivity (F/m)
	const double Qe = 1.602176565*pow(10, -19);		// Elementary positive charge (C)
	const double PI = 3.141592653;

	int charge = 0;
	double C_att_EI[5] ={};					// coefficients for actractive interaction electron-ion
	double C_rep_EI[5] ={};                                    // coefficients for repulsive interaction electron-ion
	double C_rep_II[4] ={};                                    // coefficients for repulsive interaction ion-ion


        doublereal omega11;
	doublereal lambdaDeb;
	doublereal lambdaD;
	charge = charge1*charge2;

	TransportCharged charged;
        lambdaDeb = charged.getDebyeLength(T, Te, Xe, P);

	// Average closest impact parameters
      	const double bfac = Qe * Qe / (8.0 * PI * eps0 * Kb);	//e-ion
      	const double be   = bfac / Te;   // electron-electron
      	const double bh   = bfac / T;   // ion-ion
	
	double number = 10000*(be+bh);

	lambdaD = min(lambdaDeb, number);

	// Reduced temperatures
     	const double Tste = max( (0.5 * lambdaD / be), 0.1);
     	const double Tsth = max( (0.5 * lambdaD / bh), 0.1);

	// Factors used in the curve-fitting of the reduced collision integrals
      	const double efac   = PI * lambdaD * lambdaD / (Tste * Tste);
      	const double hfac   = PI * lambdaD * lambdaD / (Tsth * Tsth);
      	const double lnTste = log(Tste);
      	const double lnTsth = log(Tsth);


	//Attractive interaction electron-ion	(REMARK: used also for attractive ion-ion)
        C_att_EI[0]=-1.67070199*pow(10, -4);
	C_att_EI[1]=4.87662165*pow(10, -3);
	C_att_EI[2]=-6.34929831*pow(10, -2);
        C_att_EI[3]=5.49816993*pow(10, -1);
	C_att_EI[4]=-7.91157696*pow(10, -1);

	//Repulsive interaction electron-ion	(REMARK: used also for repulsive electron-electron)
        C_rep_EI[0]=-3.25157374*pow(10, -4);
        C_rep_EI[1]=9.56207818*pow(10, -3);
        C_rep_EI[2]=-1.17684421*pow(10, -1);
        C_rep_EI[3]=8.41093450*pow(10, -1);
        C_rep_EI[4]=-1.39975085;


	//Repulsive interaction ion-ion
        C_rep_II[0]=3.41942*pow(10, -3);
        C_rep_II[1]=-8.01279*pow(10, -2);
        C_rep_II[2]=7.5847*pow(10, -1);
        C_rep_II[3]=-1.3513;


	// attractive interaction electron-ion
	if (    ( (species1 == "E") and ( charge < 0 ) ) or ( ( charge < 0 ) and (species2 == "E") ) )	
		{
			omega11 = efac * exp( C_att_EI[0]*lnTste*lnTste*lnTste*lnTste + C_att_EI[1]*lnTste*lnTste*lnTste + C_att_EI[2]*lnTste*lnTste + C_att_EI[3]*lnTste + C_att_EI[4]);
		}


	// attractive interaction ion-ion
	 if ( ( (species1 != "E") and ( charge < 0 ) ) or ( ( charge < 0 ) and (species2 != "E") ) )
                {
                        omega11 = hfac * exp( C_att_EI[0]*lnTsth*lnTsth*lnTsth*lnTsth + C_att_EI[1]*lnTsth*lnTsth*lnTsth + C_att_EI[2]*lnTsth*lnTsth + C_att_EI[3]*lnTsth + C_att_EI[4]);
                }

	// repulsive interaction electron-electron
	if (    ( (species1 == "E") and (species2 == "E") ) )
		
		{
			omega11 = efac * exp( C_rep_EI[0]*lnTste*lnTste*lnTste*lnTste + C_rep_EI[1]*lnTste*lnTste*lnTste + C_rep_EI[2]*lnTste*lnTste + C_rep_EI[3]*lnTste + C_rep_EI[4]);
		}

	// repulsive interaction electron-ion
	if (    ( (species1 == "E") and ( charge > 0 ) ) or ( ( charge > 0 ) and (species2 == "E") ) )	
		{
			omega11 = efac * exp( C_rep_EI[0]*lnTste*lnTste*lnTste*lnTste + C_rep_EI[1]*lnTste*lnTste*lnTste + C_rep_EI[2]*lnTste*lnTste + C_rep_EI[3]*lnTste + C_rep_EI[4]);
		}

	// repulsive interaction ion-ion
	if (    ( (species1 != "E") and ( charge > 0 ) ) or ( ( charge > 0 ) and (species2 != "E") ) )
		{
			omega11 = hfac * exp( C_rep_II[0]*lnTsth*lnTsth*lnTsth + C_rep_II[1]*lnTsth*lnTsth + C_rep_II[2]*lnTsth + C_rep_II[3]) ;
		}

		return omega11;
}


// add omega22 for actractive electron-ion actractive ion-ion repulsive electron-ion 
doublereal MMCollisionIntCharged::omega22_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P)
{

        const double Kb = 1.3806488*pow(10, -23);               // Boltzmann's constant (J/molecule-K)
        const double eps0 = 8.854187817*pow(10, -12);           // Vacuum permittivity (F/m)
        const double Qe = 1.602176565*pow(10, -19);             // Elementary positive charge (C)
        const double PI = 3.141592653;

	int charge = 0;
        double C_rep_EE[5] ={};                                    // coefficients for repulsive interaction electron-electron
        double C_rep_II[5] ={};                                    // coefficients for repulsive interaction ion-ion


        doublereal omega22;
        doublereal lambdaDeb;
        doublereal lambdaD;
	charge = charge1*charge2;


        TransportCharged charged;
        lambdaDeb = charged.getDebyeLength(T, Te, Xe, P);

        // Average closest impact parameters
        const double bfac = Qe * Qe / (8.0 * PI * eps0 * Kb);   //e-ion
        const double be   = bfac / Te;   // electron-electron
        const double bh   = bfac / T;   // ion-ion

        double number = 10000*(be+bh);


        lambdaD = min(lambdaDeb, number);


        // Reduced temperatures
        const double Tste = max( (0.5 * lambdaD / be), 0.1);
        const double Tsth = max( (0.5 * lambdaD / bh), 0.1);


        // Factors used in the curve-fitting of the reduced collision integrals
        const double efac   = PI * lambdaD * lambdaD / (Tste * Tste);
        const double hfac   = PI * lambdaD * lambdaD / (Tsth * Tsth);
        const double lnTste = log(Tste);
        const double lnTsth = log(Tsth);


        //Repulsive interaction electron-electron
        C_rep_EE[0]=-3.96530664*pow(10, -4);
        C_rep_EE[1]=1.08344320*pow(10, -2);
        C_rep_EE[2]=-1.22277658*pow(10, -1);
        C_rep_EE[3]=8.03918135*pow(10, -1);
        C_rep_EE[4]=-1.10681226;


        //Repulsive interaction ion-ion
        C_rep_II[0]=-3.96530664*pow(10, -4);
        C_rep_II[1]=1.08344320*pow(10, -2);
        C_rep_II[2]=-1.22277658*pow(10, -1);
        C_rep_II[3]=8.03918135*pow(10, -1);
        C_rep_II[4]=-1.10681226;


        // repulsive interaction electron-electron
        if (    ( (species1 == "E") and (species2 == "E") ) )

                {
                        omega22 = efac * exp( C_rep_EE[0]*lnTste*lnTste*lnTste*lnTste + C_rep_EE[1]*lnTste*lnTste*lnTste + C_rep_EE[2]*lnTste*lnTste + C_rep_EE[3]*lnTste + C_rep_EE[4]);
                }




        // repulsive interaction ion-ion
	if (    ( (species1 != "E") and ( charge > 0 ) ) or ( ( charge > 0 ) and (species2 != "E") ) )
                {
                        omega22 = hfac * exp( C_rep_II[0]*lnTsth*lnTsth*lnTsth*lnTsth + C_rep_II[1]*lnTsth*lnTsth*lnTsth + C_rep_II[2]*lnTsth*lnTsth + C_rep_II[3]*lnTsth + C_rep_II[4] ) ;
                }

                return omega22;
}



// add repulsive interaction electron-ion
doublereal MMCollisionIntCharged::omega12_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P)
{

        const double Kb = 1.3806488*pow(10, -23);               // Boltzmann's constant (J/molecule-K)
        const double eps0 = 8.854187817*pow(10, -12);           // Vacuum permittivity (F/m)
        const double Qe = 1.602176565*pow(10, -19);             // Elementary positive charge (C)
        const double PI = 3.141592653;
	const double oneThird = 1/3;
        double C_att_EI[5] ={};                                    // coefficients for attractive interaction electron-ion
	int charge = 0;


        doublereal omega12;
        doublereal lambdaDeb;
        doublereal lambdaD;
	doublereal omega11;
	doublereal num;
	omega11 = omega11_charged(species1, species2, charge1, charge2, T, Te, Xe, P);
	charge = charge1*charge2;

        TransportCharged charged;
        lambdaDeb = charged.getDebyeLength(T, Te, Xe, P);

        // Average closest impact parameters
        const double bfac = Qe * Qe / (8.0 * PI * eps0 * Kb);   //e-ion
        const double be   = bfac / Te;   // electron-electron
        const double bh   = bfac / T;   // ion-ion

        double number = 10000*(be+bh);


        lambdaD = min(lambdaDeb, number);

        // Reduced temperatures
        const double Tste = max( (0.5 * lambdaD / be), 0.1);

        // Factors used in the curve-fitting of the reduced collision integrals
        const double lnTste = log(Tste);

        //Attractive interaction electron-ion
        C_att_EI[0]=2.15877964*pow(10, -5);
        C_att_EI[1]=-6.06756574*pow(10, -4);
        C_att_EI[2]=7.22979820*pow(10, -3);
        C_att_EI[3]=-4.82251886*pow(10, -2);
        C_att_EI[4]=5.21304137*pow(10, -1);


	// attractive interaction electron-ion
	if (    ( (species1 == "E") and ( charge < 0 ) ) or ( ( charge < 0 ) and (species2 == "E") ) )
                {
                  	num = max( (C_att_EI[0]*lnTste*lnTste*lnTste*lnTste + C_att_EI[1]*lnTste*lnTste*lnTste + C_att_EI[2]*lnTste*lnTste + C_att_EI[3]*lnTste + C_att_EI[4]) , oneThird ) ;
			omega12 = omega11 * num;			
                }


	else		// omega12 set equal to omega11 for electron-neutral; missing data for repulsive interaction electron-ion;  
		{
			num = 1; // max( (C_att_EI[0]*lnTste*lnTste*lnTste*lnTste + C_att_EI[1]*lnTste*lnTste*lnTste + C_att_EI[2]*lnTste*lnTste + C_att_EI[3]*lnTste + C_att_EI[4]), oneThird ) ;
                        omega12 = 1;
		}

	return omega12;

}



doublereal MMCollisionIntCharged::omega13_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P)
{

        const double Kb = 1.3806488*pow(10, -23);               // Boltzmann's constant (J/molecule-K)
        const double eps0 = 8.854187817*pow(10, -12);           // Vacuum permittivity (F/m)
        const double Qe = 1.602176565*pow(10, -19);             // Elementary positive charge (C)
        const double PI = 3.141592653;
	int charge = 0;
        double C_att_EI[5] ={};                                    // coefficients for attractive interaction electron-ion


	vector<double> c;
	charge = charge1*charge2;
        doublereal omega13;
        doublereal lambdaDeb;
        doublereal lambdaD;
	doublereal omega11;
        doublereal omega12;
        doublereal num;
	omega11 = omega11_charged(species1, species2, charge1, charge2, T, Te, Xe, P);
        omega12 = omega12_charged(species1, species2, charge1, charge2, T, Te, Xe, P);

        TransportCharged charged;
        lambdaDeb = charged.getDebyeLength(T, Te, Xe, P);

        // Average closest impact parameters
        const double bfac = Qe * Qe / (8.0 * PI * eps0 * Kb);   //e-ion
        const double be   = bfac / Te;   // electron-electron
        const double bh   = bfac / T;   // ion-ion

        double number = 10000*(be+bh);


        lambdaD = min(lambdaDeb, number);

        // Reduced temperatures
        const double Tste = max( (0.5 * lambdaD / be), 0.1);


        // Factors used in the curve-fitting of the reduced collision integrals
        const double lnTste = log(Tste);

        //Attractive interaction electron-ion
        C_att_EI[0]=1.59680456*pow(10, -5);
        C_att_EI[1]=-5.13956309*pow(10, -4);
        C_att_EI[2]=7.74980901*pow(10, -3);
        C_att_EI[3]=-6.71047378*pow(10, -2);
        C_att_EI[4]=1.32373242;

        // attractive interaction electron-ion
	if (    ( (species1 == "E") and ( charge < 0 ) ) or ( ( charge < 0 ) and (species2 == "E") ) )
                {
                        num = max( C_att_EI[0]*lnTste*lnTste*lnTste*lnTste + C_att_EI[1]*lnTste*lnTste*lnTste + C_att_EI[2]*lnTste*lnTste + C_att_EI[3]*lnTste + C_att_EI[4], 1.0) ;
                        omega13 = 1.25*omega12 - 0.25*omega11*num;
                }

	else

		{
			num = 1;
                        omega13 = 1.25*omega12 - 0.25*omega11*num;
		}

        return omega13;

}




doublereal MMCollisionIntCharged::omega14_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P)
{

        const double Kb = 1.3806488*pow(10, -23);               // Boltzmann's constant (J/molecule-K)
        const double eps0 = 8.854187817*pow(10, -12);           // Vacuum permittivity (F/m)
        const double Qe = 1.602176565*pow(10, -19);             // Elementary positive charge (C)
        const double PI = 3.141592653;
	int charge = 0;

        double C_att_EI[5] ={};                                    // coefficients for attractive interaction electron-ion

	doublereal omega13;
	omega13 = omega13_charged(species1, species2, charge1, charge2, T, Te, Xe, P);

        doublereal omega14;
        doublereal lambdaDeb;
        doublereal lambdaD;

	charge = charge1*charge2;


        TransportCharged charged;
        lambdaDeb = charged.getDebyeLength(T, Te, Xe, P);

        // Average closest impact parameters
        const double bfac = Qe * Qe / (8.0 * PI * eps0 * Kb);   //e-ion
        const double be   = bfac / Te;   // electron-electron
        const double bh   = bfac / T;   // ion-ion

        double number = 10000*(be+bh);


        lambdaD = min(lambdaDeb, number);

        // Reduced temperatures
        const double Tste = max( (0.5 * lambdaD / be), 0.1);

        // Factors used in the curve-fitting of the reduced collision integrals
        const double efac   = PI * lambdaD * lambdaD / (Tste * Tste);
        const double lnTste = log(Tste);

        //Attractive interaction electron-ion
        C_att_EI[0]=-1.85869672*pow(10, -4);
        C_att_EI[1]=4.45192523*pow(10, -3);
        C_att_EI[2]=-4.77048273*pow(10, -2);
        C_att_EI[3]=3.89672330*pow(10, -1);
        C_att_EI[4]=-2.31134796;

        // attractive interaction electron-ion
	if (    ( (species1 == "E") and ( charge < 0 ) ) or ( ( charge < 0 ) and (species2 == "E") ) )
                {
                        omega14 = efac * exp( C_att_EI[0]*lnTste*lnTste*lnTste*lnTste + C_att_EI[1]*lnTste*lnTste*lnTste + C_att_EI[2]*lnTste*lnTste + C_att_EI[3]*lnTste + C_att_EI[4]);
                }


        else            // missing data for repulsive interaction electron-ion
                {
                        omega14 = omega13;
                }

        return omega14;

}


doublereal MMCollisionIntCharged::omega15_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P)
{

        const double Kb = 1.3806488*pow(10, -23);               // Boltzmann's constant (J/molecule-K)
        const double eps0 = 8.854187817*pow(10, -12);           // Vacuum permittivity (F/m)
        const double Qe = 1.602176565*pow(10, -19);             // Elementary positive charge (C)
        const double PI = 3.141592653;
	int charge = 0;

        double C_att_EI[5] ={};                                    // coefficients for attractive interaction electron-ion


        doublereal omega13;
        omega13 = omega13_charged(species1, species2, charge1, charge2, T, Te, Xe, P);

        doublereal omega15;
        doublereal lambdaDeb;
        doublereal lambdaD;
	charge = charge1*charge2;


        TransportCharged charged;
        lambdaDeb = charged.getDebyeLength(T, Te, Xe, P);


        // Average closest impact parameters
        const double bfac = Qe * Qe / (8.0 * PI * eps0 * Kb);   //e-ion
        const double be   = bfac / Te;   // electron-electron
        const double bh   = bfac / T;   // ion-ion

        double number = 10000*(be+bh);


        lambdaD = min(lambdaDeb, number);

        // Reduced temperatures
        const double Tste = max( (0.5 * lambdaD / be), 0.1);

        // Factors used in the curve-fitting of the reduced collision integrals
        const double efac   = PI * lambdaD * lambdaD / (Tste * Tste);
        const double lnTste = log(Tste);

        //Attractive interaction electron-ion
        C_att_EI[0]=-1.91341716*pow(10, -4);
        C_att_EI[1]=4.39527014*pow(10, -3);
        C_att_EI[2]=-4.53877276*pow(10, -2);
        C_att_EI[3]=3.69204436*pow(10, -1);
        C_att_EI[4]=-2.62434507;


        // attractive interaction electron-ion
	if (    ( (species1 == "E") and ( charge < 0 ) ) or ( ( charge < 0 ) and (species2 == "E") ) )
                {
                        omega15 = efac * exp( C_att_EI[0]*lnTste*lnTste*lnTste*lnTste + C_att_EI[1]*lnTste*lnTste*lnTste + C_att_EI[2]*lnTste*lnTste + C_att_EI[3]*lnTste + C_att_EI[4]);
                }


        else            // missing data for repulsive interaction electron-ion
                {
                        omega15 = omega13;
                }

        return omega15;

}



doublereal MMCollisionIntCharged::omega23_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P)
{

        const double Kb = 1.3806488*pow(10, -23);               // Boltzmann's constant (J/molecule-K)
        const double eps0 = 8.854187817*pow(10, -12);           // Vacuum permittivity (F/m)
        const double Qe = 1.602176565*pow(10, -19);             // Elementary positive charge (C)
        const double PI = 3.141592653;
	int charge = 0;
        double C_rep_EE[5] ={};                                    // coefficients for repulsive interaction electron-electron


	doublereal omega22;
        doublereal omega23;
        doublereal lambdaDeb;
        doublereal lambdaD;
	doublereal num;
	omega22 = omega22_charged(species1, species2, charge1, charge2, T, Te, Xe, P);
	charge = charge1*charge2;


        TransportCharged charged;
        lambdaDeb = charged.getDebyeLength(T, Te, Xe, P);


        // Average closest impact parameters
        const double bfac = Qe * Qe / (8.0 * PI * eps0 * Kb);   //e-ion
        const double be   = bfac / Te;   // electron-electron
        const double bh   = bfac / T;   // ion-ion

        double number = 10000*(be+bh);


        lambdaD = min(lambdaDeb, number);

        // Reduced temperatures
        const double Tste = max( (0.5 * lambdaD / be), 0.1);

        // Factors used in the curve-fitting of the reduced collision integrals
        const double lnTste = log(Tste);

        // Repulsive interaction electron-electron
        C_rep_EE[0]=2.14301076*pow(10, -5);
        C_rep_EE[1]=-7.35578975*pow(10, -4);
        C_rep_EE[2]=9.78964744*pow(10, -3);
        C_rep_EE[3]=-6.37201837*pow(10, -2);
        C_rep_EE[4]=7.01123442*pow(10, -1);



        // repulsive interaction electron-electron
        if (    ( (species1 == "E") and (species2 == "E") ) )

                {
			num = max( C_rep_EE[0]*lnTste*lnTste*lnTste*lnTste + C_rep_EE[1]*lnTste*lnTste*lnTste + C_rep_EE[2]*lnTste*lnTste + C_rep_EE[3]*lnTste + C_rep_EE[4], 0.5) ;
			omega23 = omega22*num;
                }

        else            // no interest on other species
                {
			num = max( C_rep_EE[0]*lnTste*lnTste*lnTste*lnTste + C_rep_EE[1]*lnTste*lnTste*lnTste + C_rep_EE[2]*lnTste*lnTste + C_rep_EE[3]*lnTste + C_rep_EE[4], 0.5) ;
                        omega23 = omega22*num;
                }

        return omega23;

}




doublereal MMCollisionIntCharged::omega24_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P)
{

        const double Kb = 1.3806488*pow(10, -23);               // Boltzmann's constant (J/molecule-K)
        const double eps0 = 8.854187817*pow(10, -12);           // Vacuum permittivity (F/m)
        const double Qe = 1.602176565*pow(10, -19);             // Elementary positive charge (C)
        const double PI = 3.141592653;
	int charge = 0;

        double C_rep_EE[5] ={};                                    // coefficients for repulsive interaction electron-electron


        doublereal omega24;
        doublereal lambdaDeb;
        doublereal lambdaD;
	charge = charge1*charge2;

        TransportCharged charged;
        lambdaDeb = charged.getDebyeLength(T, Te, Xe, P);

        // Average closest impact parameters
        const double bfac = Qe * Qe / (8.0 * PI * eps0 * Kb);   //e-ion
        const double be   = bfac / Te;   // electron-electron
        const double bh   = bfac / T;   // ion-ion

        double number = 10000*(be+bh);


        lambdaD = min(lambdaDeb, number);

        // Reduced temperatures
        const double Tste = max( (0.5 * lambdaD / be), 0.1);

        // Factors used in the curve-fitting of the reduced collision integrals
        const double efac   = PI * lambdaD * lambdaD / (Tste * Tste);
        const double lnTste = log(Tste);

        //Repulsive interaction electron-electron
        C_rep_EE[0]=-3.97192162*pow(10, -4);
        C_rep_EE[1]=1.00233388*pow(10, -2);
        C_rep_EE[2]=-1.03318360*pow(10, -1);
        C_rep_EE[3]=6.47397675*pow(10, -1);
        C_rep_EE[4]=-1.75549845;

        //Repulsive interaction electron-electron
        if (    ( (species1 == "E") and (species2 == "E") ) )

                {
                        omega24 = efac * exp( C_rep_EE[0]*lnTste*lnTste*lnTste*lnTste + C_rep_EE[1]*lnTste*lnTste*lnTste + C_rep_EE[2]*lnTste*lnTste + C_rep_EE[3]*lnTste + C_rep_EE[4]);
                }

        else            // no interest on other species
                {
                        omega24 = efac * exp( C_rep_EE[0]*lnTste*lnTste*lnTste*lnTste + C_rep_EE[1]*lnTste*lnTste*lnTste + C_rep_EE[2]*lnTste*lnTste + C_rep_EE[3]*lnTste + C_rep_EE[4]);
                }

        return omega24;

}


} // namespace
