/**
 * @file WaterPropsIAPWSphi.cpp
 *  Definitions for Lowest level of the classes which support a real water model
 *  (see class #WaterPropsIAPWS and class  #WaterPropsIAPWSphi).
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: WaterPropsIAPWSphi.cpp,v 1.12 2008/12/29 22:53:36 hkmoffa Exp $
 */

#include "WaterPropsIAPWSphi.h"

#include <cstdio>
#include <cmath>

using std::printf;
using std::sqrt;
using std::log;
using std::exp;
using std::pow;
using std::fabs;

/*
 * Critical Point values in mks units: Note, these aren't used in this 
 * routine, except for internal checks. All calculations here are done
 * in dimensionless units.
 */
static const double T_c = 647.096;  // Kelvin
static const double P_c = 22.064E6; // Pascals
static const double Rho_c = 322.;    // kg m-3
static const double M_water = 18.015268; // kg kmol-1

/*
 * The added constants were calculated so that u = s = 0
 * for liquid at the triple point. These where determined
 * by the program testPress. I'm not quite satisfied with
 * the result, but will let it stand for the moment.
 * H didn't turn out to be .611872 J/kg, but .611782 J/kg.
 * There may be a slight error here somehow.
 */
static const double ni0[9] = {
  0.0,
  -8.32044648201 - 0.000000001739715,
  6.6832105268   + 0.000000000793232,
  3.00632,
  0.012436,
  0.97315,
  1.27950,
  0.96956,
  0.24873
};

static const double gammi0[9] = {
  0.0,
  0.0,
  0.0,
  0.0,
  1.28728967,
  3.53734222,
  7.74073708,
  9.24437796,
  27.5075105
};

static const int ciR[56] = {
  0,  //  0
  0,  //  1
  0,
  0,
  0,
  0,  //  5
  0,
  0,
  1,
  1,
  1,  // 10
  1,
  1,
  1,
  1,
  1,  // 15
  1,
  1,
  1,
  1,
  1,  // 20
  1,   
  1,
  2,
  2,
  2,  // 25
  2,
  2,
  2,
  2,
  2,  // 30
  2,
  2,
  2,
  2,
  2,  // 35
  2,
  2,
  2,
  2,
  2,  // 40
  2,
  2,
  3,
  3,
  3,  // 45
  3,
  4,
  6,
  6,
  6,  // 50
  6,
  0,
  0,
  0,
  0   // 55
};

static const int diR[55] = {
  0,  //  0
  1,  //  1
  1,
  1,
  2,
  2,  //  5
  3,
  4,
  1,
  1,
  1,  // 10
  2,
  2,
  3,
  4,
  4,  // 15
  5,
  7,
  9,
  10,
  11, // 20
  13,   
  15,
  1,
  2,
  2,  // 25
  2,
  3,
  4,
  4,
  4,  // 30
  5,
  6,
  6,
  7,
  9,  // 35
  9,
  9,
  9,
  9,
  10, // 40
  10,
  12,
  3,
  4,
  4,  // 45
  5,
  14,
  3,
  6,
  6,  // 50
  6,
  3,
  3,
  3   // 54
};

static const int tiR[55] = {
  0,  //  0
  0,  //  1
  0,
  0,
  0,
  0,  //  5
  0,
  0,
  4,  //  8
  6,
  12, // 10
  1,
  5,
  4,
  2,
  13, // 15
  9,
  3,
  4,
  11,
  4,  // 20
  13,   
  1,
  7,
  1,
  9,  // 25
  10,
  10,
  3,
  7,
  10,  // 30
  10,
  6,
  10,
  10,
  1,  // 35
  2,
  3,
  4,
  8,
  6, // 40
  9,
  8,
  16,
  22,
  23,  // 45
  23,
  10,
  50,
  44,
  46,  // 50
  50,
  0,
  1,
  4   // 54
};

static const double ni[57] = {
  +0.0,
  +0.12533547935523E-1, //  1
  +0.78957634722828E1,  //  2
  -0.87803203303561E1,  //  3
  +0.31802509345418E0,  //  4
  -0.26145533859358E0,  //  5
  -0.78199751687981E-2, //  6
  +0.88089493102134E-2, //  7
  -0.66856572307965E0,  //  8
  +0.20433810950965,    //  9
  -0.66212605039687E-4, // 10
  -0.19232721156002E0,  // 11
  -0.25709043003438E0,  // 12
  +0.16074868486251E0,  // 13
  -0.40092828925807E-1, // 14
  +0.39343422603254E-6, // 15
  -0.75941377088144E-5, // 16
  +0.56250979351888E-3, // 17
  -0.15608652257135E-4, // 18
  +0.11537996422951E-8, // 19
  +0.36582165144204E-6, // 20
  -0.13251180074668E-11,// 21
  -0.62639586912454E-9, // 22
  -0.10793600908932E0,  // 23
  +0.17611491008752E-1, // 24
  +0.22132295167546E0,  // 25
  -0.40247669763528E0,  // 26
  +0.58083399985759E0,  // 27
  +0.49969146990806E-2, // 28
  -0.31358700712549E-1, // 29
  -0.74315929710341E0,  // 30
  +0.47807329915480E0,  // 31
  +0.20527940895948E-1, // 32
  -0.13636435110343E0,  // 33
  +0.14180634400617E-1, // 34
  +0.83326504880713E-2, // 35
  -0.29052336009585E-1, // 36
  +0.38615085574206E-1, // 37
  -0.20393486513704E-1, // 38
  -0.16554050063734E-2, // 39
  +0.19955571979541E-2, // 40
  +0.15870308324157E-3, // 41
  -0.16388568342530E-4, // 42
  +0.43613615723811E-1, // 43
  +0.34994005463765E-1, // 44
  -0.76788197844621E-1, // 45
  +0.22446277332006E-1, // 46
  -0.62689710414685E-4, // 47
  -0.55711118565645E-9, // 48
  -0.19905718354408E0,  // 49
  +0.31777497330738E0,  // 50
  -0.11841182425981E0,  // 51
  -0.31306260323435E2,  // 52
  +0.31546140237781E2,  // 53
  -0.25213154341695E4,  // 54
  -0.14874640856724E0,  // 55
  +0.31806110878444E0   // 56
};


static const double alphai[3] = {
  +20.,
  +20.,
  +20.
};

static const double betai[3] = {
  +150.,
  +150.,
  +250.
};

static const double gammai[3] = {
  +1.21,
  +1.21,
  +1.25
};

static const double epsi[3] = {
  +1.0,
  +1.0,
  +1.0
};

static const double ai[2] = {
  +3.5,
  +3.5
};

static const double bi[2] = {
  +0.85,
  +0.95
};

static const double Bi[2] = {
  +0.2,
  +0.2
};

static const double Ci[2] = {
  +28.0,
  +32.0
};

static const double Di[2] = {
  +700.,
  +800.
};

static const double Ai[2] = {
  +0.32,
  +0.32
};

static const double Bbetai[2] = {
  +0.3,
  +0.3
};

/**
 * Constructor for the object. 
 */
WaterPropsIAPWSphi::WaterPropsIAPWSphi() :
  TAUsave(-1.0),
  TAUsqrt(-1.0),
  DELTAsave(-1.0)
{
  for (int i = 0; i < 52; i++) {
    TAUp[i] = 1.0;
  }
  for (int i = 0; i < 16; i++) {
    DELTAp[i] = 1.0;
  }
}

/*
 * intCheck() calculates all of the functions at a one point and
 * prints out the result. It's used for conducting the internal
 * check.
 */
void WaterPropsIAPWSphi::intCheck(double tau, double delta) {
  tdpolycalc(tau, delta);
  double nau    = phi0();
  double res    = phiR();
  double res_d  = phiR_d();
  double nau_d  = phi0_d();
  double res_dd = phiR_dd();
  double nau_dd = phi0_dd();
  double res_t  = phiR_t();
  double nau_t  = phi0_t();
  double res_tt = phiR_tt();
  double nau_tt = phi0_tt();
  double res_dt = phiR_dt();
  double nau_dt = phi0_dt();

  std::printf("nau    = %20.12e\t\tres    = %20.12e\n", nau,    res);
  std::printf("nau_d  = %20.12e\t\tres_d  = %20.12e\n", nau_d,  res_d);
  printf("nau_dd = %20.12e\t\tres_dd = %20.12e\n", nau_dd, res_dd);
  printf("nau_t  = %20.12e\t\tres_t  = %20.12e\n", nau_t,  res_t);
  printf("nau_tt = %20.12e\t\tres_tt = %20.12e\n", nau_tt, res_tt);
  printf("nau_dt = %20.12e\t\tres_dt = %20.12e\n", nau_dt, res_dt);
}

void WaterPropsIAPWSphi::check1() {
  double T = 500.;
  double rho = 838.025;
  double tau = T_c/T;
  double delta = rho / Rho_c;
  printf(" T = 500 K, rho = 838.025 kg m-3\n");
  intCheck(tau, delta);
}

void WaterPropsIAPWSphi::check2() {
  double T = 647;
  double rho = 358.0;
  double tau = T_c/T;
  double delta = rho / Rho_c;
  printf(" T = 647 K, rho = 358.0 kg m-3\n");
  intCheck(tau, delta);
}

/*
 * Calculate the polynomials in tau and delta, and store them in static 
 * storage.
 */
void WaterPropsIAPWSphi::tdpolycalc(double tau, double delta) {
  if ((tau != TAUsave) || 1) {
    TAUsave = tau;
    TAUsqrt = sqrt(tau);
    TAUp[0] = 1.0;
    for (int i = 1; i < 51; i++) {
      TAUp[i] = TAUp[i-1] * tau;
    }
  }
  if ((delta != DELTAsave) || 1) {
    DELTAsave = delta;
    DELTAp[0] = 1.0;
    for (int i = 1; i <= 15; i++) {
      DELTAp[i] = DELTAp[i-1] * delta;
    }
  }
}

/*
 * Calculate Eqn. 6.5 for phi0, the ideal gas part of the
 * dimensionless Helmholtz free energy.
 */
double WaterPropsIAPWSphi::phi0() const {
  double tau = TAUsave;
  double delta = DELTAsave;
  double retn = log(delta) + ni0[1] + ni0[2]*tau + ni0[3]*log(tau);

  retn += ni0[4] * log(1.0 - exp(-gammi0[4]*tau));
  retn += ni0[5] * log(1.0 - exp(-gammi0[5]*tau));
  retn += ni0[6] * log(1.0 - exp(-gammi0[6]*tau));
  retn += ni0[7] * log(1.0 - exp(-gammi0[7]*tau));
  retn += ni0[8] * log(1.0 - exp(-gammi0[8]*tau));
  return retn;
}

/*
 * Calculate Eqn. 6.6 for phiR, the residual part of the
 * dimensionless Helmholtz free energy.
 *
 *  tau = dimensionless temperature
 *  delta = dimensionless pressure
 */
double WaterPropsIAPWSphi::phiR() const {
  double tau = TAUsave;
  double delta = DELTAsave;
  int i, j;
    
  /*
   * Write out the first seven polynomials in the expression
   */
  double T375 = pow(tau, 0.375);
  double val = (ni[1] * delta / TAUsqrt +
		ni[2] * delta * TAUsqrt * T375 +
		ni[3] * delta * tau +
		ni[4] * DELTAp[2] * TAUsqrt +
		ni[5] * DELTAp[2] * T375 * T375 +
		ni[6] * DELTAp[3] * T375 +
		ni[7] * DELTAp[4] * tau);
  /*
   * Next, do polynomial contributions 8 to 51
   */
  for (i = 8; i <= 51; i++) {
    val += (ni[i] * DELTAp[diR[i]] * TAUp[tiR[i]] * exp(-DELTAp[ciR[i]]));
  }

  /*
   * Next do contributions 52 to 54
   */
  for (j = 0; j < 3; j++) {
    i = 52 + j;
    double dtmp = delta - epsi[j];
    double ttmp = tau - gammai[j];
    val += (ni[i] * DELTAp[diR[i]] * TAUp[tiR[i]] * 
	    exp(-alphai[j]*dtmp*dtmp - betai[j]*ttmp*ttmp));
  }

  /*
   * Next do contributions 55 and 56
   */
  for (j = 0; j < 2; j++) {
    i = 55 + j;
    double deltam1 = delta - 1.0;
    double dtmp2 = deltam1 * deltam1;
    double atmp = 0.5 / Bbetai[j];
    double theta = (1.0 - tau) + Ai[j] * pow(dtmp2, atmp);
    double triag = theta * theta + Bi[j] * pow(dtmp2, ai[j]);
    double ttmp = tau - 1.0;
   
    double triagtmp = pow(triag, bi[j]);
      
    double phi = exp(-Ci[j]*dtmp2 - Di[j]*ttmp*ttmp);
    val += (ni[i] * triagtmp * delta * phi);
  }

  return val;
}

/*
 * Calculate the Phi function, which is basically the helmholtz free energy
 * Eqn. (6.4)
 */
double WaterPropsIAPWSphi::phi(double tau, double delta) {
  tdpolycalc(tau, delta);
  double nau = phi0();
  double res = phiR();
  double retn = nau + res;
  return retn;
}


/*
 * Calculate d_phiR_d(delta), the first derivative of phiR
 * wrt delta
 *
 *  tau = dimensionless temperature
 *  delta = dimensionless pressure
 */
double WaterPropsIAPWSphi::phiR_d() const {
  double tau = TAUsave;
  double delta = DELTAsave;
  int i, j;
    
  /*
   * Write out the first seven polynomials in the expression
   */
  double T375 = pow(tau, 0.375);
  double val = (ni[1] / TAUsqrt +
		ni[2] * TAUsqrt * T375 +
		ni[3] * tau +
		ni[4] * 2.0 * delta * TAUsqrt +
		ni[5] * 2.0 * delta * T375 * T375 +
		ni[6] * 3.0 * DELTAp[2] * T375 +
		ni[7] * 4.0 * DELTAp[3] * tau);
  /*
   * Next, do polynomial contributions 8 to 51
   */
  for (i = 8; i <= 51; i++) {
    val += ( (ni[i] * exp(-DELTAp[ciR[i]]) *  DELTAp[diR[i] - 1] *
	      TAUp[tiR[i]]) * (diR[i] - ciR[i]* DELTAp[ciR[i]]));
  }

  /*
   * Next do contributions 52 to 54
   */
  for (j = 0; j < 3; j++) {
    i = 52 + j;
    double dtmp = delta - epsi[j];
    double ttmp = tau - gammai[j];
    double tmp = (ni[i] * DELTAp[diR[i]] * TAUp[tiR[i]] * 
		  exp(-alphai[j]*dtmp*dtmp - betai[j]*ttmp*ttmp));
    val += tmp * (diR[i]/delta - 2.0 * alphai[j] * dtmp);
  }

  /*
   * Next do contributions 55 and 56
   */
  for (j = 0; j < 2; j++) {
    i = 55 + j;
    double deltam1 = delta - 1.0;
    double dtmp2 = deltam1 * deltam1;
    double atmp = 0.5 / Bbetai[j];
    double theta = (1.0 - tau) + Ai[j] * pow(dtmp2, atmp);
    double triag = theta * theta + Bi[j] * pow(dtmp2, ai[j]);
    double ttmp = tau - 1.0;
   
    double triagtmp = pow(triag, bi[j]);
    double triagtmpm1 = pow(triag, bi[j]-1.0);
    double atmpM1 = atmp - 1.0;
    double ptmp = pow(dtmp2,atmpM1); 
    double p2tmp = pow(dtmp2, ai[j]-1.0);
    double dtriagddelta = 
      deltam1 *(Ai[j] * theta * 2.0 / Bbetai[j] * ptmp +
		2.0*Bi[j]*ai[j]*p2tmp);

    double phi = exp(-Ci[j]*dtmp2 - Di[j]*ttmp*ttmp);
    double dphiddelta = -2.0*Ci[j]*deltam1*phi;
    double dtriagtmpddelta = bi[j] * triagtmpm1 * dtriagddelta;

    double tmp = ni[i] * (triagtmp * (phi + delta*dphiddelta) +
			  dtriagtmpddelta * delta * phi);
    val += tmp;
  }

  return val;
}

/*
 * Calculate d_phi0_d(delta), the first derivative of phi0
 * wrt delta
 *
 *  tau = dimensionless temperature
 *  delta = dimensionless pressure
 */
double WaterPropsIAPWSphi::phi0_d() const {
  double delta = DELTAsave;
  return (1.0/delta);
}

/*
 * Calculate the dPhidDelta function, which is basically the derivative
 * of helmholtz free energy wrt delta
 * Eqn. (6.4)
 */
double WaterPropsIAPWSphi::phi_d(double tau, double delta) {
  tdpolycalc(tau, delta);
  double nau = phi0_d();
  double res = phiR_d();
  double retn = nau + res;
  return retn;
}

/*
 * Calculate the dimensionless pressure at tau and delta;
 *
 *       p/(rhoRT) = delta * phi_d()
 *
 * note: this is done so much, we have a seperate routine.
 */
double WaterPropsIAPWSphi::pressureM_rhoRT(double tau, double delta) {
  tdpolycalc(tau, delta);
  double res = phiR_d();
  double retn = 1.0 + delta * res; 
  return retn;
}

/*
 * Calculate d2_phiR_dd(delta), the second derivative of phiR
 * wrt delta
 *
 *  tau = dimensionless temperature
 *  delta = dimensionless pressure
 */
double WaterPropsIAPWSphi::phiR_dd() const {
  double tau = TAUsave;
  double delta = DELTAsave;
  int i, j;
  double atmp;
    
  /*
   * Write out the first seven polynomials in the expression
   */
  double T375 = pow(tau, 0.375);
  double val = (ni[4] * 2.0 * TAUsqrt +
		ni[5] * 2.0 * T375 * T375 +
		ni[6] * 6.0 * delta * T375 +
		ni[7] * 12.0 * DELTAp[2] * tau);
  /*
   * Next, do polynomial contributions 8 to 51
   */
  for (i = 8; i <= 51; i++) {
    double dtmp = DELTAp[ciR[i]];
    double tmp = ni[i] * exp(-dtmp) * TAUp[tiR[i]];
    if (diR[i] == 1) {
      atmp = 1.0/delta;
    } else {
      atmp = DELTAp[diR[i] - 2];
    }
    tmp *= atmp *((diR[i] - ciR[i]*dtmp)*(diR[i]-1.0-ciR[i]*dtmp) -
		  ciR[i]*ciR[i]*dtmp);
    val += tmp;
  }

  /*
   * Next do contributions 52 to 54
   */
  for (j = 0; j < 3; j++) {
    i = 52 + j;
    double dtmp = delta - epsi[j];
    double ttmp = tau - gammai[j];
    double tmp = (ni[i] * TAUp[tiR[i]] * 
		  exp(-alphai[j]*dtmp*dtmp - betai[j]*ttmp*ttmp));
    double deltmp =  DELTAp[diR[i]];
    double deltmpM1 = deltmp/delta;
    double deltmpM2 = deltmpM1 / delta;
    double d2tmp = dtmp * dtmp;
      
    val += tmp * (-2.0*alphai[j]*deltmp +
		  4.0 * alphai[j] * alphai[j] * deltmp * d2tmp -
		  4.0 * diR[i] * alphai[j] * deltmpM1 * dtmp +
		  diR[i] * (diR[i] - 1.0) * deltmpM2);
  }

  /*
   * Next do contributions 55 and 56
   */
  for (j = 0; j < 2; j++) {
    i = 55 + j;
    double deltam1 = delta - 1.0;
    double dtmp2 = deltam1 * deltam1;
    atmp = 0.5 / Bbetai[j];
    double theta = (1.0 - tau) + Ai[j] * pow(dtmp2, atmp);
    double triag = theta * theta + Bi[j] * pow(dtmp2, ai[j]);
    double ttmp = tau - 1.0;
   
    double triagtmp = pow(triag, bi[j]);
    double triagtmpm1 = pow(triag, bi[j]-1.0);
    double atmpM1 = atmp - 1.0;
    double ptmp = pow(dtmp2,atmpM1);
    double p2tmp = pow(dtmp2, ai[j]-1.0);
    double dtriagddelta = 
      deltam1 *(Ai[j] * theta * 2.0 / Bbetai[j] * ptmp +
		2.0*Bi[j]*ai[j]*p2tmp);

    double phi = exp(-Ci[j]*dtmp2 - Di[j]*ttmp*ttmp);
    double dphiddelta = -2.0*Ci[j]*deltam1*phi;
    double dtriagtmpddelta = bi[j] * triagtmpm1 * dtriagddelta;


    double d2phiddelta2 = 2.0 * Ci[j] * phi * (2.0*Ci[j]*dtmp2 - 1.0);

    double pptmp = ptmp / dtmp2;
    double d2triagddelta2 = dtriagddelta / deltam1;
    d2triagddelta2 += 
      dtmp2 *(4.0*Bi[j]*ai[j]*(ai[j]-1.0)*pow(dtmp2,ai[j]-2.0) +
	      2.0*Ai[j]*Ai[j]/(Bbetai[j]*Bbetai[j])*ptmp*ptmp +
	      Ai[j]*theta*4.0/Bbetai[j]*(atmp-1.0)*pptmp);

    double  d2triagtmpd2delta = 
      bi[j] * (triagtmpm1 * d2triagddelta2 +
	       (bi[j]-1.0)*triagtmpm1/triag*dtriagddelta*dtriagddelta);

    double ctmp = (triagtmp * (2.0*dphiddelta + delta*d2phiddelta2) +
		   2.0*dtriagtmpddelta*(phi + delta * dphiddelta) +
		   d2triagtmpd2delta * delta * phi);
 
    val += ni[i] * ctmp;
  }

  return val;
}

/*
 * Calculate d2_phi0_dd(delta), the second derivative of phi0
 * wrt delta
 *
 *  tau = dimensionless temperature
 *  delta = dimensionless pressure
 */
double WaterPropsIAPWSphi::phi0_dd() const {
  double delta = DELTAsave;
  return (-1.0/(delta*delta));
}

/*
 * Calculate the d2_PhidDelta2 function, which is the second derivative
 * of helmholtz free energy wrt delta
 * Eqn. (6.4)
 */
double WaterPropsIAPWSphi::phi_dd(double tau, double delta) {
  tdpolycalc(tau, delta);
  double nau = phi0_dd();
  double res = phiR_dd();
  double retn = nau + res;
  return retn;
}

double WaterPropsIAPWSphi::dimdpdrho(double tau, double delta) {
  tdpolycalc(tau, delta);
  double res1 = phiR_d();
  double res2 = phiR_dd();
  double retn = 1.0 + delta * (2.0*res1 + delta*res2);
  return retn;
}

double WaterPropsIAPWSphi::dimdpdT(double tau, double delta) {
  tdpolycalc(tau, delta);
  double res1 = phiR_d();
  double res2 = phiR_dt();
  double retn = (1.0 + delta * res1) - tau * delta * (res2);
  return retn;
}

/*
 * Calculate d_phi0/d(tau)
 */
double WaterPropsIAPWSphi::phi0_t() const {
  double tau = TAUsave;
  double retn = ni0[2] + ni0[3]/tau;;
  retn += (ni0[4] * gammi0[4] * (1.0/(1.0 - exp(-gammi0[4]*tau)) - 1.0));
  retn += (ni0[5] * gammi0[5] * (1.0/(1.0 - exp(-gammi0[5]*tau)) - 1.0));
  retn += (ni0[6] * gammi0[6] * (1.0/(1.0 - exp(-gammi0[6]*tau)) - 1.0));
  retn += (ni0[7] * gammi0[7] * (1.0/(1.0 - exp(-gammi0[7]*tau)) - 1.0));
  retn += (ni0[8] * gammi0[8] * (1.0/(1.0 - exp(-gammi0[8]*tau)) - 1.0));
  return retn;
}

/*
 * Calculate Eqn. 6.6 for dphiRdtau, the derivative residual part of the
 * dimensionless Helmholtz free energy wrt temperature
 *
 *  tau = dimensionless temperature
 *  delta = dimensionless pressure
 */
double WaterPropsIAPWSphi::phiR_t() const {
  double tau = TAUsave;
  double delta = DELTAsave;
  int i, j;
  double atmp, tmp;
    
  /*
   * Write out the first seven polynomials in the expression
   */
  double T375 = pow(tau, 0.375);
  double val = ((-0.5) *ni[1] * delta / TAUsqrt / tau +
		ni[2] * delta * 0.875 / TAUsqrt * T375 +
		ni[3] * delta +
		ni[4] * DELTAp[2] * 0.5 / TAUsqrt +
		ni[5] * DELTAp[2] * 0.75 * T375 * T375 / tau +
		ni[6] * DELTAp[3] * 0.375 * T375 / tau +
		ni[7] * DELTAp[4]);
  /*
   * Next, do polynomial contributions 8 to 51
   */
  for (i = 8; i <= 51; i++) {
    tmp = (ni[i] * DELTAp[diR[i]] * TAUp[tiR[i]-1] * exp(-DELTAp[ciR[i]]));
    val += tiR[i] * tmp;
  }

  /*
   * Next do contributions 52 to 54
   */
  for (j = 0; j < 3; j++) {
    i = 52 + j;
    double dtmp = delta - epsi[j];
    double ttmp = tau - gammai[j];
    tmp = (ni[i] * DELTAp[diR[i]] * TAUp[tiR[i]] * 
	   exp(-alphai[j]*dtmp*dtmp - betai[j]*ttmp*ttmp));
    val += tmp *(tiR[i]/tau - 2.0 * betai[j]*ttmp);
  }

  /*
   * Next do contributions 55 and 56
   */
  for (j = 0; j < 2; j++) {
    i = 55 + j;
    double deltam1 = delta - 1.0;
    double dtmp2 = deltam1 * deltam1;
    atmp = 0.5 / Bbetai[j];
    double theta = (1.0 - tau) + Ai[j] * pow(dtmp2, atmp);
    double triag = theta * theta + Bi[j] * pow(dtmp2, ai[j]);
    double ttmp = tau - 1.0;
   
    double triagtmp = pow(triag, bi[j]);
      
    double phi = exp(-Ci[j]*dtmp2 - Di[j]*ttmp*ttmp);
 

    double dtriagtmpdtau = -2.0*theta * bi[j] * triagtmp / triag;

    double dphidtau = - 2.0 * Di[j] * ttmp * phi;
      
    val += ni[i] * delta * (dtriagtmpdtau * phi + triagtmp * dphidtau);
  }

  return val;
}

/*
 * Calculate the dPhidtau function, which is basically the derivative
 * of helmholtz free energy wrt tau
 * Eqn. (6.4)
 */
double WaterPropsIAPWSphi::phi_t(double tau, double delta) {
  tdpolycalc(tau, delta);
  double nau = phi0_t();
  double res = phiR_t();
  double retn = nau + res;
  return retn;
}

/*
 * Calculate d2_phi0/dtau2
 */
double WaterPropsIAPWSphi::phi0_tt() const {
  double tau = TAUsave;
  double tmp, itmp;
  double retn = - ni0[3]/(tau * tau);
  for (int i = 4; i <= 8; i++) {
    tmp = exp(-gammi0[i]*tau);
    itmp = 1.0 - tmp;
    retn -= (ni0[i] * gammi0[i] * gammi0[i] * tmp / (itmp * itmp));
  }
  return retn;
}

/*
 * Calculate Eqn. 6.6 for dphiRdtau, the second derivative residual part of the
 * dimensionless Helmholtz free energy wrt temperature
 *
 *  tau = dimensionless temperature
 *  delta = dimensionless pressure
 */
double WaterPropsIAPWSphi::phiR_tt() const {
  double tau = TAUsave;
  double delta = DELTAsave;
  int i, j;
  double atmp, tmp;
    
  /*
   * Write out the first seven polynomials in the expression
   */
  double T375 = pow(tau, 0.375);
  double val = ((-0.5) * (-1.5) * ni[1] * delta / (TAUsqrt * tau * tau) +
		ni[2] * delta * 0.875 * (-0.125) * T375 / (TAUsqrt * tau) +
		ni[4] * DELTAp[2] * 0.5 * (-0.5)/ (TAUsqrt * tau) +
		ni[5] * DELTAp[2] * 0.75 *(-0.25) * T375 * T375 / (tau * tau) +
		ni[6] * DELTAp[3] * 0.375 *(-0.625) * T375 / (tau * tau));
  /*
   * Next, do polynomial contributions 8 to 51
   */
  for (i = 8; i <= 51; i++) {
    if (tiR[i] > 1) {
      tmp = (ni[i] * DELTAp[diR[i]] * TAUp[tiR[i]-2] * exp(-DELTAp[ciR[i]]));
      val += tiR[i] *  (tiR[i] - 1.0) * tmp;
    }
  }

  /*
   * Next do contributions 52 to 54
   */
  for (j = 0; j < 3; j++) {
    i = 52 + j;
    double dtmp = delta - epsi[j];
    double ttmp = tau - gammai[j];
    tmp = (ni[i] * DELTAp[diR[i]] * TAUp[tiR[i]] * 
	   exp(-alphai[j]*dtmp*dtmp - betai[j]*ttmp*ttmp));
    atmp = tiR[i]/tau - 2.0 * betai[j]*ttmp;
    val += tmp *(atmp * atmp - tiR[i]/(tau*tau) - 2.0*betai[j]);
  }

  /*
   * Next do contributions 55 and 56
   */
  for (j = 0; j < 2; j++) {
    i = 55 + j;
    double deltam1 = delta - 1.0;
    double dtmp2 = deltam1 * deltam1;
    atmp = 0.5 / Bbetai[j];
    double theta = (1.0 - tau) + Ai[j] * pow(dtmp2, atmp);
    double triag = theta * theta + Bi[j] * pow(dtmp2, ai[j]);
    double ttmp = tau - 1.0;
   
    double triagtmp = pow(triag, bi[j]);
    double triagtmpM1 = triagtmp / triag;
      
    double phi = exp(-Ci[j]*dtmp2 - Di[j]*ttmp*ttmp);
 

    double dtriagtmpdtau = -2.0*theta * bi[j] * triagtmp / triag;

    double dphidtau = - 2.0 * Di[j] * ttmp * phi;

    double d2triagtmpdtau2 = 
      (2 * bi[j] * triagtmpM1 +
       4 * theta * theta  * bi[j] * (bi[j]-1.0) * triagtmpM1 / triag);

    double d2phidtau2 = 2.0*Di[j]*phi *(2.0*Di[j]*ttmp*ttmp - 1.0);

    tmp =  (d2triagtmpdtau2 * phi +
	    2 * dtriagtmpdtau * dphidtau +
	    triagtmp * d2phidtau2);
    val += ni[i] * delta * tmp;
  }

  return val;
}

/*
 * Calculate the d2Phidtau2 function, which is basically the second derivative
 * of helmholtz free energy wrt tau
 * Eqn. (6.4)
 */
double WaterPropsIAPWSphi::phi_tt(double tau, double delta) {
  tdpolycalc(tau, delta);
  double nau = phi0_tt();
  double res = phiR_tt();
  double retn = nau + res;
  return retn;
}

/**
 * Calculate d2_phi0/dtauddelta
 */
double WaterPropsIAPWSphi::phi0_dt() const {
  return 0.0;
}

/*
 * Calculate d2_phiR_d(delta)d(tau), the mixed derivative of phi
 * wrt delta and tau.
 *
 *  tau = dimensionless temperature
 *  delta = dimensionless pressure
 */
double WaterPropsIAPWSphi::phiR_dt() const {
  double tau = TAUsave;
  double delta = DELTAsave;
  int i, j;
  double tmp;
  /*
   * Write out the first seven polynomials in the expression
   */
  double T375 = pow(tau, 0.375);
  double val = (ni[1] * (-0.5) / (TAUsqrt * tau) +
		ni[2] * (0.875) * T375 / TAUsqrt +
		ni[3] +
		ni[4] * 2.0 * delta * (0.5) / TAUsqrt +
		ni[5] * 2.0 * delta * (0.75) * T375 * T375 / tau +
		ni[6] * 3.0 * DELTAp[2] * 0.375 * T375 / tau +
		ni[7] * 4.0 * DELTAp[3]);
  /*
   * Next, do polynomial contributions 8 to 51
   */
  for (i = 8; i <= 51; i++) {
    tmp = (ni[i] * tiR[i] * exp(-DELTAp[ciR[i]]) *  DELTAp[diR[i] - 1] *
	   TAUp[tiR[i] - 1]);
    val += tmp * (diR[i] - ciR[i] * DELTAp[ciR[i]]);
  }

  /*
   * Next do contributions 52 to 54
   */
  for (j = 0; j < 3; j++) {
    i = 52 + j;
    double dtmp = delta - epsi[j];
    double ttmp = tau - gammai[j];
    tmp = (ni[i] * DELTAp[diR[i]] * TAUp[tiR[i]] * 
	   exp(-alphai[j]*dtmp*dtmp - betai[j]*ttmp*ttmp));
    val += tmp * ((diR[i]/delta - 2.0 * alphai[j] * dtmp) * 
		  (tiR[i]/tau   - 2.0 * betai[j]  * ttmp));
  }

  /*
   * Next do contributions 55 and 56
   */
  for (j = 0; j < 2; j++) {
    i = 55 + j;
    double deltam1 = delta - 1.0;
    double dtmp2 = deltam1 * deltam1;
    double atmp = 0.5 / Bbetai[j];
    double theta = (1.0 - tau) + Ai[j] * pow(dtmp2, atmp);
    double triag = theta * theta + Bi[j] * pow(dtmp2, ai[j]);
    double ttmp = tau - 1.0;
   
    double triagtmp = pow(triag, bi[j]);
    double triagtmpm1 = pow(triag, bi[j]-1.0);
    double atmpM1 = atmp - 1.0;
    double ptmp = pow(dtmp2,atmpM1); 
    double p2tmp = pow(dtmp2, ai[j]-1.0);
    double dtriagddelta = 
      deltam1 *(Ai[j] * theta * 2.0 / Bbetai[j] * ptmp +
		2.0*Bi[j]*ai[j]*p2tmp);

    double phi = exp(-Ci[j]*dtmp2 - Di[j]*ttmp*ttmp);
    double dphiddelta = -2.0*Ci[j]*deltam1*phi;
    double dtriagtmpddelta = bi[j] * triagtmpm1 * dtriagddelta;


    double dtriagtmpdtau = -2.0*theta * bi[j] * triagtmp / triag;

    double dphidtau = - 2.0 * Di[j] * ttmp * phi;

    double d2phiddeltadtau = 4.0 * Ci[j] * Di[j] * deltam1 * ttmp * phi;

    double  d2triagtmpddeltadtau =
      ( -Ai[j] * bi[j] * 2.0 / Bbetai[j] * triagtmpm1 * deltam1 * ptmp
	-2.0 * theta * bi[j] * (bi[j] - 1.0) * triagtmpm1 / triag * dtriagddelta);
	      

    double tmp = ni[i] * (triagtmp * (dphidtau + delta*d2phiddeltadtau) +
			  delta * dtriagtmpddelta * dphidtau +
			  dtriagtmpdtau * (phi + delta * dphiddelta) +
			  d2triagtmpddeltadtau * delta * phi);
    val += tmp;
  }

  return val;
}

/**
 * This program computes the reduced density, given the reduced pressure
 * and the reduced temperature, tau. It takes an initial guess, deltaGuess.
 * DeltaGuess is important as this is a multivalued function below the
 * critical point.
 * 
 */
double WaterPropsIAPWSphi::dfind(double p_red, double tau, double deltaGuess) {
  double dd = deltaGuess;
  bool conv = false;
  double deldd = dd;
  double pcheck = 1.0E-30 + 1.0E-8 * p_red;
  for (int n = 0; n < 200; n++) {
    /*
     * Calculate the internal polynomials, and then calculate the
     * phi deriv functions needed by this routine.
     */
    tdpolycalc(tau, dd);
    double q1 = phiR_d();
    double q2 = phiR_dd();

    /*
     * Calculate the predicted reduced pressure, pred0, based on the
     * current tau and dd. 
     */
    double pred0 = dd + dd * dd * q1;
    /*
     * Calculate the derivative of the predicted reduced pressure
     * wrt the reduced density, dd, This is dpddelta
     */
    double dpddelta = 1.0 + 2.0 * dd * q1 + dd * dd * q2; 
    /*
     * If dpddelta is negative, then we are in the middle of the
     * 2 phase region, beyond the stability curve. We need to adjust
     * the initial guess outwards and start a new iteration.
     */
    if (dpddelta <= 0.0) {
      if (deltaGuess > 1.0) dd = dd * 1.05;
      if (deltaGuess < 1.0) dd = dd * 0.95;
      continue;
    }
    /*
     * Check for convergence
     */
    if (fabs(pred0-p_red) < pcheck) {
      conv = true;
      break;
    }

    /*
     * Dampen and crop the update
     */
    double dpdx = dpddelta;
    if (n < 10) {
      dpdx = dpddelta * 1.1;
    }
    if (dpdx < 0.001) dpdx = 0.001;
      
    /*
     * Formulate the update to reduced density using 
     * Newton's method. Then, crop it to a max value
     * of 0.02
     */
    deldd = - (pred0 - p_red) / dpdx;
    if (fabs(deldd) > 0.05) {
      deldd = deldd * 0.05 / fabs(deldd);
    }
    /*
     * updated the reduced density value
     */
    dd = dd + deldd;
    if (fabs(deldd/dd) < 1.0E-14) {
      conv = true;
      break;
    }
    /*
     * Check for negative densities
     */
    if (dd <= 0.0) {
      dd = 1.0E-24;
    }
  }
  /*
   * Check for convergence, and return 0.0 if it wasn't achieved.
   */
  if (! conv) {
    dd = 0.0;
  }
  return dd;
}

/**
 * Calculate the dimensionless gibbs free energy g/RT.
 */
double WaterPropsIAPWSphi::gibbs_RT() const {
  double delta = DELTAsave;
  double rd = phiR_d();
  double g = 1.0 + phi0() + phiR() + delta * rd;
  return g;
}

/**
 * Calculate the dimensionless enthalpy h/RT.
 */
double WaterPropsIAPWSphi::enthalpy_RT() const {
  double delta = DELTAsave;
  double tau   = TAUsave;
  double rd = phiR_d();
  double nt = phi0_t();
  double rt = phiR_t();
  double hRT = 1.0 + tau * (nt + rt) + delta * rd;
  return hRT;
}

/*
 * Calculate the dimensionless entropy s/R.
 */
double WaterPropsIAPWSphi::entropy_R() const {
  double tau   = TAUsave;
  double nt = phi0_t();
  double rt = phiR_t();
  double p0 = phi0();
  double pR = phiR();
  double sR = tau * (nt + rt) - p0 - pR;
  return sR;
}

/*
 * Calculate the dimensionless internal energy, u/RT.
 */
double WaterPropsIAPWSphi::intEnergy_RT() const {
  double tau   = TAUsave;
  double nt = phi0_t();
  double rt = phiR_t();
  double uR = tau * (nt + rt);
  return uR;
}

/*
 * Calculate the dimensionless constant volume Heat Capacity, Cv/R
 */
double WaterPropsIAPWSphi::cv_R() const {
  double tau   = TAUsave;
  double ntt = phi0_tt();
  double rtt = phiR_tt();
  double cvR = - tau * tau * (ntt + rtt);
  return cvR;
}

/*
 * Calculate the dimensionless constant pressure Heat Capacity, Cp/R
 */
double WaterPropsIAPWSphi::cp_R() const {
  double tau   = TAUsave;
  double delta = DELTAsave;
  double cvR = cv_R();
  //double nd = phi0_d();
  double rd = phiR_d();
  double rdd = phiR_dd();
  double rdt = phiR_dt();
  double num = (1.0 + delta * rd - delta * tau * rdt);
  double cpR = cvR + (num * num / 
		      (1.0 + 2.0 * delta * rd + delta * delta * rdd));
  return cpR;
}

