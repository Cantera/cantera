// Lee-Kesler equation of state for use with custom fluids

#include "CLK.h"
#include <math.h>

const double b[2][4] = {{0.1181193, 0.265728, 0.154790, 0.030323},
                        {0.2026579, 0.331511, 0.027655, 0.203488}};
const double c[2][4] = {{0.0236744, 0.0186984, 0.0, 0.042724},
                        {0.0313385, 0.0503618, 0.016901, 0.041577}};
const double d[2][2] = {{1.55488e-5, 6.23689e-5},{4.8736e-5, 0.740336e-5}};
const double beta[2] = {0.65392, 1.226};
const double gamma[2] = {0.060167, 0.03754};

//--------------------------- member functions ------------------

double CLK::W(int n, double egrho, double Gamma) {
	return (n == 0 ? (1.0 - egrho)/(2.0*Gamma) : 
	   (n*W(n-1, egrho, Gamma) - 0.5*pow(Rho,2*n)*egrho)/Gamma);
}

double CLK::u_ni(){
	return -R()*T*T*I()/Tcr;                   // + u_0(T) 
}

double CLK::s_ni() {
	const double Pref = 101325.0;
	double rgas = R();
	return rgas*(log(Pref/(Rho*rgas*T)) - (T/Tcr)*I() - J());  // + s^0(T,p_0)
}

double CLK::I() {                    // \int_0^\rho_r (1/\rho_r)(dZ/dT_r) d\rho_r
	double Bp, Cp, Dp;
	double rtr = Tcr/T;
	double rtr2 = rtr*rtr;
	double rvr = 8314.3*Tcr*Rho/(Pcr*Mw);  //   1/v_r^\prime
	double rvr2 = rvr*rvr;
    double egrho;

	egrho = exp(-gamma[Isr]*rvr2);
	Bp = rtr2*b[Isr][1] + 2.0*rtr*rtr2*b[Isr][2] + 3.0*rtr2*rtr2*b[Isr][3];
	Cp = rtr2*c[Isr][1] - 3.0*c[Isr][2]*rtr2*rtr2;
	Dp = -d[Isr][1]*rtr2;
	double r = Bp*rvr + 0.5*rvr2*Cp + 0.2*pow(rvr,5)*Dp 
   		- 3.0*c[Isr][3]*rtr2*rtr2*(beta[Isr]*W(0,egrho,gamma[Isr]) 
		+ gamma[Isr]*W(1,egrho,gamma[Isr]));
	return r;
}

double CLK::J() {                    // \int_0^\rho_r (1/\rho_r)(Z - 1) d\rho_r
	double BB, CC, DD;
	double rtr = Tcr/T;
	double rtr2 = rtr*rtr;
	double rvr = 8314.3*Tcr*Rho/(Pcr*Mw);  //   1/v_r^\prime
	double rvr2 = rvr*rvr;
    double egrho;

	egrho = exp(-gamma[Isr]*rvr2);
	BB = b[Isr][0] - rtr*(b[Isr][1] 
		+ rtr*(b[Isr][2] + rtr*b[Isr][3]));
	CC = c[Isr][0] - rtr*(c[Isr][1] - c[Isr][2]*rtr*rtr);
	DD = d[Isr][0] + d[Isr][1]*rtr;
	double r = BB*rvr + 0.5*rvr2*CC + 0.2*pow(rvr,5)*DD 
		+ c[Isr][3]*rtr2*rtr*(beta[Isr]*W(0,egrho,gamma[Isr]) 
		+ gamma[Isr]*W(1,egrho,gamma[Isr]));
	return r;
}

double CLK::z() {
	double zz, rvr2, BB, CC, DD, EE;
	double rtr = Tcr/T;             //   1/T_r
	double rvr = Rho*8314.3*Tcr/(Pcr*Mw);
	rvr2 = rvr*rvr;
	BB = b[Isr][0] - rtr*(b[Isr][1] 
			+ rtr*(b[Isr][2] + rtr*b[Isr][3]));
	CC = c[Isr][0] - rtr*(c[Isr][1] - c[Isr][2]*rtr*rtr);
	DD = d[Isr][0] + d[Isr][1]*rtr;
	EE = exp(-gamma[Isr]*rvr2);

	zz = 1.0 + BB*rvr + CC*rvr2 + DD*pow(rvr,5)
			+ c[Isr][3]*pow(rtr,3)*rvr2*
			(beta[Isr] + gamma[Isr]*rvr2)*EE;
	return zz;
}


double CLK::Pp() {
	return 8314.3*z()*Rho*T/Mw;
}

double CLK::Psat(){
	double tr = 1.0 - Tcr/T;
	double lpr;

	if (Isr == 0)
		lpr = 5.395743797*tr + 0.05524287*tr*tr + 0.06853005*tr*tr*tr;
	else
		lpr = 7.259961465*tr - 0.549206092*tr*tr + 0.177581752*tr*tr*tr;
	return Pcr*exp(lpr);
}

double CLK::ldens(){
	double x = 1.0 - T/Tcr;

	double rho_r;
	if (Isr == 0)
		rho_r = 5.2307 + 15.16*x - 21.9778*x*x + 18.767*x*x*x;
	else {
		rho_r = 6.166930606 + 17.42866964*x - 18.62589833*x*x + 11.73957224*x*x*x;
		rho_r *= 1.0;
	}
	return Pcr*rho_r*Mw/(8314.3*Tcr);
}

double CLK::Tcrit() {return Tcr;}
double CLK::Pcrit() {return Pcr;}
double CLK::Vcrit() {return 0.2901*R()*Tcr/Pcr;}
char * CLK::formula() {return "---";}
