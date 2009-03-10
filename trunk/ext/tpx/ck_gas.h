#ifndef CK_GAS_H
#define CK_GAS_H

#include "Sub.h"

//#include <fstream.h>

extern "C" { 
	__declspec(dllimport) double __stdcall
    ckprop(int ijob, double temp, double dens, double *xmole);

	__declspec(dllimport) int __stdcall ckspecies();

	__declspec(dllimport) double __stdcall
    chem(double temp, double dens, int ks, double *xmole, int ijob);
}

namespace tpx {
class ck_gas : public Substance {
public:

	ck_gas() {
		kk = ckspecies();
		if (kk > 0) {
			T = 300.;
			Rho = 0.001;
			xm = new double[kk];
		}
	}


	ck_gas(int ki, double *xmoles)
	{
		kk = ckspecies();
		if (kk > 0 && kk == ki) {
			T = 300.;
			Rho = 1.;
			xm = new double[kk];
			double sum=0.0;
			for (int i = 0; i<kk; i++) { 
				xm[i] = xmoles[i];
				sum += xm[i];
			}
			for (i = 0; i<kk; i++) {
				xm[i] = xm[i]/sum;
			}
		}
		else {
			T = 300.;
		    Rho = 1.;
		    xm = new double[1];
			Err = CKError;
		}
	}


	~ck_gas() { 
		delete xm; 
	}

	double MolWt();
	double Tcrit() {return 1.0;}
    double Pcrit() {return 2.0;}
	double Vcrit() {return 1.0/300.0;}
	double Tmin()  {return 200.0;}
	double Tmax()  {return 5000.0;}
	char * name()  {return "Gas mixture";}
	char * formula() {return "Chemkin";}
	double Pp();
	double up();
	double sp();

	double Psat() { return 0.0; }
	double ldens() { return 1.e6; }

    double Cdot(int ks);
    double Ddot(int ks);
    double Wdot(int ks);
    double Tchem(int ks);

protected:
	int kk;  
	double *xm;

};
}
#endif // ! CK_GAS
