#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

/*
 *  Values of A_V : tabular form
 *       units  sqrt(kg/gmol)cm3/gmol
 */
double A_V(double temp)
{
    double retn;
    if (temp == 323.15) {
        retn = 2.37356;
    } else if (temp == 473.15) {
        retn = 15.0766;
    } else {
        printf("A_V unknown temp value %g\n", temp);
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
    } else if (ifunc == 3) {
        retn = 0.0;
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
        retn = q6  + q9 * (temp - tref)
               + q10 * (temp * temp - tref * tref);
    } else if (ifunc == 1) {
        retn = q9 + 2.0 * q10 * temp;
    } else if (ifunc == 2) {
        retn = 2.0 * q10;
    } else if (ifunc == 3) {
        retn = 0.0;
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
    } else if (ifunc == 3) {
        retn = 0.0;
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
     * Calculate A_V in sqrt(kg/gmol)cm3/gmol
     */
    double Av = A_V(temp);

    double beta0prime3 = Beta0(temp, 3);
    printf("   beta0prime3 = %g\n", beta0prime3);

    double beta1prime3 = Beta1(temp, 3);
    printf("   beta1prime3 = %g\n", beta1prime3);

    double cphiprime3= Cphi(temp, 3);
    printf("   Cphiprime = %g\n", cphiprime3);

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

    double BVmx = beta0prime3 + 2.0 * beta1prime3 / (a2* Iionic) *
                  (1.0 - (1.0 + alpha * sqrtI) * exp(-alpha*sqrtI));

    double CVmx = 0.5 * sqrt(vm * vx) * cphiprime3;

    double phiV = v * zm * zx * (Av/(2.*b)) * log(1 + 1.2 * sqrtI) +
                  2 * vm * vx * GasConst * temp *1.0E3 * (m * BVmx + m * m * CVmx);

    printf("   phiV = %15.8g cm3/gmolSalt\n", phiV);

    double molecWeight = 18.01528;

    double xo = 1.0 / (molecWeight/1000. * 2 * m + 1.0);
    printf(" no = %g\n", xo);
}

int main()
{
    printf("Standalone test of the apparent relative molal excess volume, phiV:\n");
    printf(" (Check against simple formula in \n");
    printf("    Activity Coefficients in Eletrolyte Solutions, 2nd Ed K. S. Pitzer, "
           "CRC Press, Boca Raton, 1991   \n");
    printf("T = 50C\n");
    double Iionic = 6.146;
    printf("Ionic Strength = %g\n", Iionic);

    calc(273.15 + 50., Iionic);
    printf("T = 200C\n");
    printf("Ionic Strength = %g\n", Iionic);

    calc(273.15 + 200., Iionic);
    return 0;
}
