#ifndef CFLUID_H
#define CFLUID_H

#include "sub.h"
#include "ck_gas.h"

#include <fstream.h>
#include <string.h>

// custom fluid base class

class CFluid : public Substance {

public:

	CFluid(char *name, double MWt, int kmx, double *xmoles) 
	{
		ig = new ck_gas(kmx, xmoles);
		T = Undef;
		Rho = Undef;
		Mw = MWt;
		CName = new char(strlen(name)+1);
		strcpy(CName,name);
	}


	~CFluid() { 
		delete ig; 
	}

	double R(){
		return 8314.3/Mw;
	}

    double MolWt() {
		return Mw;
	}

    double Tmin() {
		return ig->Tmin();
	}

    double Tmax() {
		return ig->Tmax();
	}

    char * name() {
		return CName;
	}

    char * formula() {
		return "doda";
	}

//protected:

	double Mw;
	ck_gas*  ig;
    char* CName;

    virtual double u_ni()=0;
	virtual double s_ni()=0;

	double up(){
		ig->Set(TP,T,101325.0);
		return ig->u() + u_ni();
	}

	double sp(){
		ig->Set(TP,T,101325.0);
		return ig->s() + s_ni();
	}
 
};
#endif // ! CFLUID