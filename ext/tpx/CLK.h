#ifndef CLK_H
#define CLK_H

#include "CFluid.h"

#include <fstream.h>

// custom Lee-Kesler fluid class

class CLK : public CFluid {

public:

	CLK(char* name, double tcrit, double pcrit, int itype, double MWt, int kmx, double *xmoles) :
		  	CFluid(name, MWt, kmx, xmoles) 
		{
			Tcr = tcrit;
			Pcr = pcrit;
			Isr = itype;
			}

	~CLK(){}

	double Tcrit();
	double Pcrit();
	double Vcrit();
	char * formula();

//protected:
	
	double Tcr;
	double Pcr;
	int Isr;

	double u_ni();
	double s_ni();
    double Pp();
	double z();

	double Psat();
	double ldens();

private:
	double W(int n, double egrho, double gamma);
	double I();
	double J();

};
#endif // ! CLK