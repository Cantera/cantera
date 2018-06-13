#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

/*
 *  Values of A_J/R : tabular form
 *       units  sqrt(kg/gmol),
 */
double A_JdR(double temp)
{
    double retn;
    if (temp == 323.15) {
        retn = 5.50274;
    } else if (temp == 473.15) {
        retn = 24.8263;
    } else {
        printf("A_JdR unknown temp value %g\n", temp);
        exit(-1);
    }
    return retn;
}

double Beta0(double temp, int ifunc)
{
    double q1 = 0.0765;
    double q2 = -777.03;
    double q3 = -4.4706;
    double q4 = 0.008946;
    double q5 = -3.3158E-6;
    double retn;
    double tref = 298.15;
    if (ifunc == 0) {
        retn = q1 + q2 * (1.0/temp - 1.0/tref)
               + q3 * (log(temp/tref)) + q4 * (temp - tref)
               + q5 * (temp * temp - tref * tref);
    } else if (ifunc == 1) {
        retn = (- q2 * 1.0/(temp* temp)
                + q3 / temp
                + q4
                + 2.0 * temp * q5);
    } else if (ifunc == 2) {
        retn = (2.0 * q2 * 1.0/(temp* temp*temp)
                - q3 / (temp*temp)
                + 2.0 * q5);
    } else {
        exit(-1);
    }
    return retn;
}

double Beta1(double temp, int ifunc)
{
    double q6 = 0.2664;
    double q9 = 6.1608E-5;
    double q10 = 1.0715E-6;
    double retn;
    double tref = 298.15;
    if (ifunc == 0) {
        retn = q6 + q9 * (temp - tref)
               + q10 * (temp * temp - tref * tref);
    } else if (ifunc == 1) {
        retn = q9 + 2.0 * q10 * temp;
    } else if (ifunc == 2) {
        retn = 2.0 * q10;
    } else {
        exit(-1);
    }
    return retn;
}

double Cphi(double temp, int ifunc)
{
    double q11 = 0.00127;
    double q12 = 33.317;
    double q13 = 0.09421;
    double q14 = -4.655E-5;
    double retn;
    double tref = 298.15;
    if (ifunc == 0) {
        retn = q11 + q12 * (1.0/temp - 1.0/tref)
               + q13 * (log(temp/tref)) + q14 * (temp - tref);
    } else if (ifunc == 1) {
        retn = - q12 / (temp * temp)
               + q13 / temp + q14;
    } else if (ifunc == 2) {
        retn = + 2.0 * q12 / (temp * temp * temp)
               - q13 / (temp * temp);
    } else {
        exit(-1);
    }
    return retn;
}

void calc(double temp, double Iionic)
{
    /*
     * Gas Constant in J gmol-1 K-1
     */
    double GasConst = 8.314472;

    /*
     * Calculate A_H in J gmol-1 sqrt(kg/gmol)
     */
    double A_J = A_JdR(temp);
    A_J *= GasConst;

    double beta0prime2 = Beta0(temp, 2);
    printf("   beta0prime2 = %g\n", beta0prime2);

    double beta1prime2 = Beta1(temp, 2);
    printf("   beta1prime2 = %g\n", beta1prime2);

    double cphiprime2 = Cphi(temp, 2);
    printf("   Cphiprime2 = %g\n", cphiprime2);

    double beta0prime = Beta0(temp, 1);
    printf("   beta0prime2 = %g\n", beta0prime);

    double beta1prime = Beta1(temp, 1);
    printf("   beta1prime = %g\n", beta1prime);

    double cphiprime = Cphi(temp, 1);
    printf("   Cphiprime = %g\n", cphiprime);

    double vm = 1.0;
    double vx = 1.0;
    double v = vm + vx;
    double m = Iionic;
    double zm = 1.;
    double zx = 1.0;

    double sqrtI = sqrt(Iionic);

    double alpha = 2.0;
    double a2 = alpha * alpha;
    double b = 1.2;

    double Bpmx = beta0prime + 2.0 * beta1prime / (a2* Iionic) *
                  (1.0 - (1.0 + alpha * sqrtI) * exp(-alpha*sqrtI));

    double Bppmx = beta0prime2 + 2.0 * beta1prime2 / (a2* Iionic) *
                   (1.0 - (1.0 + alpha * sqrtI) * exp(-alpha*sqrtI));

    double Cpmx = 0.5 * sqrt(vm * vx) * cphiprime;

    double Cppmx = 0.5 * sqrt(vm * vx) * cphiprime2;

    double Bmx = Bppmx + 2.0 / temp * Bpmx;
    double Cmx = Cppmx + 2.0 / temp * Cpmx;

    double phiJ = v * zm * zx * (A_J/(2.*b)) * log(1 + 1.2 * sqrtI) -
                  2 * vm * vx * GasConst * temp * temp * (m * Bmx + m * m * Cmx);
    phiJ *= 1.0E-3;
    printf("   phiJ = %15.8g kJ/gmolSalt\n", phiJ);

    double molecWeight = 18.01528;

    double xo = 1.0 / (molecWeight/1000. * 2 * m + 1.0);
    printf(" no = %g\n", xo);
}

int main()
{
    printf("Standalone test of the apparent relative molal enthalpy, phiL:\n");
    printf(" (Check against simple formula in Silvester&Pitzer, J. Phys. Chem. 81, 1822 (1977)\n");

    printf("T = 50C\n");
    double Iionic = 6.146;
    printf("Ionic Strength = %g\n", Iionic);

    calc(273.15 + 50., Iionic);
    printf("T = 200C\n");
    printf("Ionic Strength = %g\n", Iionic);

    calc(273.15 + 200., Iionic);
    return 0;
}
