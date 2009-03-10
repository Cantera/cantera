// Chemkin ideal gas mixture

#include "ck_gas.h"
#include <math.h>

double ck_gas::up(){
	return ckprop(2, T, Rho, xm);                   // + h_0(T) 
}

double ck_gas::sp() {return ckprop(3, T, Rho, xm);}


double ck_gas::Pp() {
	return ckprop(1, T, Rho, xm);
}

double ck_gas::MolWt() {return ckprop(0, T, Rho, xm);}

double ck_gas::Cdot(int ks) {return chem(T, Rho, ks, xm, 1);}

double ck_gas::Ddot(int ks) {return chem(T, Rho, ks, xm, 2);}

double ck_gas::Wdot(int ks) {return chem(T, Rho, ks, xm, 3);}

double ck_gas::Tchem(int ks) {return chem(T, Rho, ks, xm, 4);}

