#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include "TransportReadExt.h"
#include "StringFunct.h"

using namespace std;

namespace Cantera {


TransportReadExt::TransportReadExt()
{
	load_coeff();
}



bool TransportReadExt::isValidCollisionString(const std::string &str)
    {
      // Must be at least 5 characters long
      if (str.length() < 5)
        return false;

      // Must have periods bounding the string but nowhere else
      if (str[0] != '.' || str.find(".", 1) != str.length() - 1)
        return false;

      // Must have dash somewhere in the middle
      if (str.find("-", 2) > str.length() - 3)
        return false;

      // No spaces
      if (str.find(" ") != string::npos)
        return false;

      // Assume that above constraints are good enough for now...
      // Note we could use regex functions for a more accurate analysis
      // later when it is added to the standard library
      return true;
    }


void TransportReadExt::load_coeff()
{
ifstream file("heavy.dat", ios::in);

string str;
vector<string> CollPairs;

double coefficients_omega11[4];
double coefficients_omega22[4];

vector<double> coeffs_omega11;
vector<double> coeffs_omega22;

double coefficients_bstar[4];
vector<double> coeffs_bstar;      
file >> str;

int nS2 = 0;
int numS_ion2 = 0;

ifstream file2("list_neut.dat", ios::in);
file2 >> nS2;
ifstream file3("list_charge.dat", ios::in);
file3 >> numS_ion2;


const int nS = nS2;
const int numS_ion = numS_ion2;

 vector<string>  sp(nS);
 vector<string> sp_ion(numS_ion);


for (int i=0; i< nS; i++)
{
	file2 >> sp[i];
}

for (int i=0; i< numS_ion; i++)
{
        file3 >> sp_ion[i];
}


file2.close();
file3.close();



	while (str != "STOP"){
        // Check if we have landed on a collision identifier
        if (isValidCollisionString(str)) {

        CollisionPair CollPair(str);


	bool test = false;
	for (int i=0; i < nS; i++)
	{
		for (int i2=0; i2 < nS; i2++)
		{
        		for (int j=0; j< numS_ion; j++)
        		{
                		if (    ( (CollPair.speciesName1()==sp[i]) and (CollPair.speciesName2()==sp[i2])        )       or
                        		( (CollPair.speciesName1()==sp[i]) and (CollPair.speciesName2()==sp_ion[j])     )       or
                        		( (CollPair.speciesName1()==sp_ion[j]) and (CollPair.speciesName2()==sp[i])     )
                   		   )

                        	{
                                	test = true;

                        	}
        		}
		}
	}

	if (test)
	{

        CollPairs.push_back(CollPair.speciesName1());
        CollPairs.push_back(CollPair.speciesName2());

        for (int i=0; i<4; i++)
        {
        	file >> coefficients_omega11[i];
		coeffs_omega11.push_back(coefficients_omega11[i]);
		C_omega11.push_back(coefficients_omega11[i]);
        }



        for (int i=0; i<4; i++)
        {
        	file >> coefficients_omega22[i];
        	coeffs_omega22.push_back(coefficients_omega22[i]);
		C_omega22.push_back(coefficients_omega22[i]);

        }




        	for (int i=0; i<4; i++)
        	{
                	file >> coefficients_bstar[i];
                	coeffs_bstar.push_back(coefficients_bstar[i]);
                	C_bstar.push_back(coefficients_bstar[i]);

        	}

	}
      }

      file >> str;

      }

      file.close();

	Pairs = CollPairs;

}

}

