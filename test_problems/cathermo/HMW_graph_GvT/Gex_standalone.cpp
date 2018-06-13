#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

double Beta0(double temp)
{
    double q1 = 0.0765;
    double q2 = -777.03;
    double q3 = -4.4706;
    double q4 = 0.008946;
    double q5 = -3.3158E-6;
    double tref = 298.15;
    return q1 + q2 * (1.0/temp - 1.0/tref)
           + q3 * (log(temp/tref)) + q4 * (temp - tref)
           + q5 * (temp * temp - tref * tref);
}

double Beta1(double temp)
{
    double q6 = 0.2664;
    double q9 = 6.1608E-5;
    double q10 = 1.0715E-6;
    double tref = 298.15;
    return q6 + q9 * (temp - tref)
           + q10 * (temp * temp - tref * tref);
}

double Cphi(double temp)
{
    double q11 = 0.00127;
    double q12 = 33.317;
    double q13 = 0.09421;
    double q14 = -4.655E-5;
    double tref = 298.15;
    return q11 + q12 * (1.0/temp - 1.0/tref)
           + q13 * (log(temp/tref)) + q14 * (temp - tref);
}

void calc(double temp, double Iionic)
{
    double Aphi = 0.0;
    if (temp == 323.15) {
        Aphi = 0.4102995331359;
    } else if (temp == 473.15) {
        Aphi = 0.622777;
    } else {
        printf("ERROR: unknown temp\n");
        exit(-1);
    }

    printf("   Aphi = %g\n", Aphi);

    double beta0 = Beta0(temp);
    printf("   beta0 = %g\n", beta0);

    double beta1 = Beta1(temp);
    printf("   beta1 = %g\n", beta1);

    double cphi = Cphi(temp);
    printf("   Cphi = %g\n", cphi);

    double vm = 1.0;
    double vx = 1.0;
    double v = vm + vx;
    double m = Iionic;
    double zm = 1.;
    double zx = 1.0;

    double sqrtI = sqrt(Iionic);

    double alpha = 2.0;
    double b = 1.2;

    double osm1 = - zm * zx * Aphi * sqrtI / (1.0 + b * sqrtI)
                  + m * 2.0 * vm * vx / v * (beta0 + beta1 * exp(- alpha * sqrtI))
                  + m * m * 2 * pow((vm*vx), 1.5) / v * cphi;

    double os = osm1 + 1.0;
    double a2 = alpha * alpha;

    printf("osmotic coeff = %20.13g\n", os);

    double lnmeanAct = - zm * zx * Aphi *
                       (sqrtI / (1.0 + b * sqrtI) + 2 / b * log(1.0 + b * sqrtI)) +
                       m * 2 * vm * vx / v * (2.0 * beta0 +
                               2.0 * beta1/ (a2 * Iionic)
                               * (1.0- (1.0 + alpha* sqrtI - a2*Iionic/2.0)* exp(- alpha * sqrtI))
                                             )
                       + 3.0 * m * m / 2.0 * (2 * sqrt(vm * vx) * vm * vx / v) * cphi;

    printf("ln(meanac) = %20.13g\n", lnmeanAct);

    double actCoeff = exp(lnmeanAct);

    printf("actCoeff = %20.13g\n", actCoeff);

    /*
     * Gas constant in J gmol-1 K-1
     */
    double GasConst = 8.314472;
    double gex = v * m * GasConst * 1.0E-3 * temp * (-osm1 + lnmeanAct);
    printf("   Gex = %20.13g kJ/kg_water\n", gex);

    double RT = GasConst * temp * 1.0E-3;
    double IdealMixing = 2.0 * RT * m * (log(m) - 1.0);
    printf("   IdealMixing = %20.13g kJ/kg_water\n", IdealMixing);

    double DelG = gex + IdealMixing;
    printf("   G - G0 = %20.13g kJ/kg_water\n", DelG);

    double mu0[6], mu[6];
    mu0[0] = -307.76256;
    double molecWeight = 18.01528;
    double diff = - RT * molecWeight / 1000. * 2. * m * os;
    mu[0] = mu0[0] + diff;
    printf("mus_kJ/gmol - H2O(L) - %20.13g   %20.13g\n", mu0[0], mu[0]);
    printf(" diff = %20.14g\n", diff);
    double xo = 1.0 / (molecWeight/1000. * 2 * m + 1.0);
    printf(" no = %g\n", xo);
    double tmp = diff / RT;
    double actCoefWater = exp(tmp) / xo;
    printf("actCoefWater = %g\n", actCoefWater);
}

int main()
{
    printf("standalone test of Gibbs excess free energy:\n");
    printf("T = 50C\n");
    double Iionic = 6.146;
    printf("Ionic Strength = %g\n", Iionic);

    calc(273.15 + 50., Iionic);
    printf("T = 200C\n");

    calc(273.15 + 200., Iionic);
    return 0;
}
