#ifndef LKW_H
#define LKW_H

#include "lk.h"

class lkw {
public:
	lkw(double tc = 1.0, double pc = 1.0, double wt = 8314.3, double omega = 0.0)
	{
		Tcr = tc;
		Pcr = pc;
		Mw = wt;
		Acent = omega;
		T = Undefined;
		Pressure = Undefined;
	};
	~lkw() {};
/*
	double MolWt();
	double Tcrit();
    double Pcrit();
	double Vcrit();
	double Tmin();
	double Tmax();
	char * name();
	char * formula();
*/
	double Ps();
//	double Pp();
//	double up();
//	double sp();
//	double Psat();
	double z(double tstar, double pstar);
//	double hdep();
//	double sdep();
//	double ldens();

protected:
	double Tcr, Pcr, Mw;
	double Acent;
	double Pressure;
	double T, Rho, Rhov, Rhof;
	lk f0; //(1.0, 1.0, 8314.3, 0);
	lk f1; //(1.0, 1.0, 8314.3, 1);

private:
	
};
#endif // ! LKW_H
