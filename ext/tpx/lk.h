#ifndef LK_H
#define LK_H

#include "sub.h"

class lk : public Substance{
public:
	lk(double tc = 1.0, double pc = 1.0, double wt = 1.0, int itype = 0)
	{
		Tcr = tc;
		Pcr = pc;
		Mw = wt;
		Isr = itype;   // simple fluid or reference 
	};
	~lk() {};

	double MolWt();
	double Tcrit();
        double Pcrit();
	double Vcrit();
	double Tmin();
	double Tmax();
	char * name();
	char * formula();
	
	double Pp();
	double up();
	double sp();
	double Psat();
	double z();
	double hdep();
	double sdep();
	double ldens();

protected:
	double Tcr, Pcr, Mw;
	int Isr;

private:
	double W(int n, double egrho, double gamma);
	double I();
	double J();
};
#endif // ! LK_H
