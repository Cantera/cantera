
#include "cantera/thermo/WaterSSTP.h"
#include "cantera/transport/WaterTransport.h"

#include <cstdio>

using namespace std;
using namespace Cantera;

double tvalue(double val, double atol = 1.0E-9)
{
    double rval = val;
    if (fabs(val) < atol) {
        rval = 0.0;
    }
    return rval;
}

int main()
{

    try {
        double lambda;
        WaterSSTP* w = new WaterSSTP("waterTPphase.xml", "");

        WaterTransport* wtTran = new WaterTransport(w, 3);
        printf("------------------------------------------------------------------------------------\n");
        printf("   T(C)    MPa     Phase         Visc     Visc(paper)   lambda     lambda(paper)\n");
        printf("                                   10-6 kg/m/s          10-3 W/m/s \n");
        printf("------------------------------------------------------------------------------------\n");

        double T = 273.15 + 25.0;
        double pres = 1.0E5;
        w->setState_TP(T, pres);
        double visc = wtTran->viscosity();
        lambda = wtTran->thermalConductivity();
        printf("%8g %10.3g  L  %13.4g      890.5    %13.4g     607.2\n",
               T - 273.15, pres * 1.0E-6, visc * 1.0E6, lambda * 1.0E3);


        T = 273.15 + 100.0;
        pres = 1.0E5;
        w->setState_TP(T, pres);
        visc = wtTran->viscosity();
        lambda = wtTran->thermalConductivity();
        printf("%8g %10.3g  L  %13.4g      281.9    %13.4g     679.1\n",
               T - 273.15, pres * 1.0E-6, visc * 1.0E6, lambda * 1.0E3);



        T = 273.15 + 100.0;
        pres = 1.0E7;
        w->setState_TP(T, pres);
        visc = wtTran->viscosity();
        lambda = wtTran->thermalConductivity();
        printf("%8g %10.3g  L  %13.6g      284.5   %13.6g      684.5\n",
               T - 273.15, pres * 1.0E-6, visc * 1.0E6, lambda * 1.0E3);

        T = 273.15 + 250.0;
        pres = 5.0E6;
        w->setState_TP(T, pres);
        visc = wtTran->viscosity();
        lambda = wtTran->thermalConductivity();
        printf("%8g %10.3g  L  %13.6g      106.4    %13.6g     622.7\n",
               T - 273.15, pres * 1.0E-6, visc * 1.0E6, lambda * 1.0E3);

        T = 273.15 + 250.0;
        pres = 5.0E7;
        w->setState_TP(T, pres);
        visc = wtTran->viscosity();
        lambda = wtTran->thermalConductivity();
        printf("%8g %10.3g  L  %13.6g      117.5    %13.6g     672.1\n",
               T - 273.15, pres * 1.0E-6, visc * 1.0E6, lambda * 1.0E3);


        T = 273.15 + 350.0;
        pres = 1.75E7;
        w->setState_TP(T, pres);
        visc = wtTran->viscosity();
        lambda = wtTran->thermalConductivity();
        printf("%8g %10.3g  L  %13.6g      67.0     %13.6g     452.3\n",
               T - 273.15, pres * 1.0E-6, visc * 1.0E6, lambda * 1.0E3);


        T = 273.15 + 400.0;
        pres = 1.50E7;
        w->setState_TP(T, pres);
        visc = wtTran->viscosity();
        lambda = wtTran->thermalConductivity();
        printf("%8g %10.3g  SC %13.6g      24.93     %13.6g     80.68\n",
               T - 273.15, pres * 1.0E-6, visc * 1.0E6, lambda * 1.0E3);

        printf("---------------------------------------------------------------------------------\n");

        delete w;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        Cantera::appdelete();
        return -1;
    }


    return 0;
}
