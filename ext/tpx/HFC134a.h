#ifndef HFC134_H
#define HFC134_H

#include "sub.h"

class HFC134a : public Substance{
public:
	HFC134a(){}
       ~HFC134a(){}

	double MolWt();
	double Tcrit();
    double Pcrit();
	double Vcrit();
	double Tmin();
	double Tmax();
	char * name();
	char * formula();

	double Pp();
	double fp();
	double up();
	double sp() 
	   { return (up() - fp())/T; }
	double Psat();
//    double dPsatdT();
private:
	double ldens();
    };
#endif // ! HFC134_H

