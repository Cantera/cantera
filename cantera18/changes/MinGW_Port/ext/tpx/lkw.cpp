// Lee-Kesler equation of state

#include "lkw.h"
#include <math.h>

const double omega_ref = 0.3978;

//--------------------------- member functions ------------------

double lkw::z(double Temp, double Pres) {
	T = Temp;
	f0.Set(TP,Temp,Pres);
	f1.Set(TP,Temp,Pres);
	if (Temp < Tcrit() && f0.x() != f1.x()) {  // need to find Psat
		if (Pres >= Ps()) {    //liquid
			f0.Set_meta(Liquid,Pres);
			f1.Set_meta(Liquid,Pres);
		}
		else {
			f0.Set_meta(Vapor,Pres);
			f1.Set_meta(Vapor,Pres);
		}
	}
	double z0 = f0.z();
	double zz = z0 + (f1.z() - z0)*Acent/omega_ref;
	return zz;
}


double lkw::Ps() {
	
// start with linear estimate of Psat
	double p0 = f0.Psat();
	double pse = p0 + (f1.Psat() - p0)/omega_ref;
	double lps = log(pse);
	double zf, zv;
// loop until g_f = g_v	
	for (int i = 0; i<20; i++) {
	
		f0.Set_meta(Liquid,pse);  // set both to liquid
		f1.Set_meta(Liquid,pse);
		double gf = f0.gp() + (f1.gp() - f0.gp())/omega_ref;
		zf = f0.z() + (f1.z() - f0.z())/omega_ref;

		f0.Set_meta(Vapor,pse);  // set both to vapor
		f1.Set_meta(Vapor,pse);
		double gv = f0.gp() + (f1.gp() - f0.gp())/omega_ref;
		zv = f0.z() + (f1.z() - f0.z())/omega_ref;

		double gfv = gv - gf;
		double dlp = gfv*Mw/(8314.3*T*(zv - zf));
		lps -= dlp;
		pse = exp(lps);
		if (fabs(gfv) < 0.001) break;
	}	
	if (i >= 20) {
		Pst = Undefined;
		Rhv = Undefined;
		Rhf = Undefined;
		Tslast = Undefined;
		set_Err(NoConverge);
	}
	else {
		Pst = pse;
		Tslast = T;
		Rhv = pse*Mw/(zv*8314.3*T);
		Rhf = pse*Mw/(zf*8314.3*T);
	}
	return Pst;
}
/*
double lk::Tcrit() {return Tcr;}
double lk::Pcrit() {return Pcr;}
double lk::Vcrit() {return 0.2901*8314.3*Tcr/(Pcr*Mw);}
double lk::Tmin() {return -100.0;}
double lk::Tmax() {return 10000.0;}
char * lk::name() {return "Lee-Kesler";}
char * lk::formula() {return "---";}
double lk::MolWt() {return Mw;}
*/
