#include "cantera/thermo/PureFluidPhase.h"
#include "cantera/thermo.h"
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
#if defined(_MSC_VER) && _MSC_VER < 1900
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    double pres;
    try {

        unique_ptr<ThermoPhase> w(newPhase("liquidvapor.yaml", "water"));

        /*
         * Print out the triple point conditions
         */
        double temp = 273.16;
        pres = w->satPressure(temp);
        writelog("psat({:g}) = {:.4g}\n", temp, pres);
        double presLow = 1.0E-2;
        temp = 298.15;
        double oneBar = 1.0E5;

        writelog("Comparisons to NIST: (see http://webbook.nist.gov):\n\n");

        w->setDensity(1.0E-8);
        w->setState_TP(temp, presLow);
        double h = w->enthalpy_mole();
        writelog("H0(298.15) = {:g} J/kmol\n", h);
        double h298 = h;

        double s = w->entropy_mole();
        s -= GasConstant * log(oneBar/presLow);
        writelog("S0(298.15) = {:g} J/kmolK\n", s);


        double T[20];
        T[0] = 298.15;
        T[1] = 500.;
        T[2] = 600.;
        T[3] = 1000.;

        double Cp0, delh0, delg0, g;

        writelog("\nIdeal Gas Standard State:\n");
        writelog("        T      Cp0           S0     "
               " -(G0-H298)/T       H0-H298\n");
        writelog("       (K)   (J/molK)     (J/molK)  "
               "   (J/molK)        (kJ/mol)\n");
        for (int i = 0; i < 4; i++) {
            temp = T[i];
            w->setState_TP(temp, presLow);
            h = w->enthalpy_mole();
            delh0 = tvalue(h - h298, 1.0E-6);
            g = w->gibbs_mole();
            delg0 = (g - h298)/temp + GasConstant * log(oneBar/presLow);
            Cp0 = w->cp_mole();
            s = w->entropy_mole();
            s -= GasConstant * log(oneBar/presLow);
            writelog("{:10g} {:10g} {:13.4g} {:13.4g} {:13.4g}\n", temp, Cp0*1.0E-3, s*1.0E-3,
                   -delg0*1.0E-3, delh0*1.0E-6);
        }
        writelog("\n\n");

        temp = 298.15;
        w->setDensity(1000.);
        w->setState_TP(temp, oneBar);
        h = w->enthalpy_mole();
        writelog("H_liq(298.15, onebar) = {:g} J/kmol\n", h);
        double h298l = h;
        s = w->entropy_mole();
        writelog("S_liq(298.15, onebar) = {:g} J/kmolK\n", s);


        T[0] = 273.19;
        T[1] = 298.15;
        T[2] = 300.;
        T[3] = 373.15;
        T[4] = 400.;
        T[5] = 500.;
        writelog("\nLiquid 1bar or psat Standard State\n");
        writelog("       T     press         psat            Cp0            S0   "
               "  -(G0-H298)/T       H0-H298\n");
        writelog("      (K)     (bar)        (bar)        (J/molK)       (J/molK)"
               "     (J/molK)        (kJ/mol)\n");

        for (int i = 0; i < 6; i++) {
            temp = T[i];
            double psat = w->satPressure(temp);
            double press = oneBar;
            if (psat > press) {
                press = psat*1.002;
            }
            w->setState_TP(temp, press);
            h = w->enthalpy_mole();
            delh0 = tvalue(h - h298l, 1.0E-6);
            g = w->gibbs_mole();
            delg0 = (g - h298l)/temp;
            Cp0 = w->cp_mole();
            s = w->entropy_mole();
            writelog("{:10g} {:10g} {:12g} {:13.4g} {:13.4g} {:13.4g} {:13.4g}\n", temp, press*1.0E-5,
                   psat*1.0E-5,
                   Cp0*1.0E-3, s*1.0E-3,
                   -delg0*1.0E-3, delh0*1.0E-6);
        }

        writelog("\nLiquid Densities:\n");
        writelog("       T     press         psat        Density          molarVol   "
               "\n");
        writelog("      (K)     (bar)        (bar)      (kg/m3)          (m3/kmol)"
               "\n");
        for (int i = 0; i < 6; i++) {
            temp = T[i];
            double psat = w->satPressure(temp);
            double press = oneBar;
            if (psat > press) {
                press = psat*1.002;
            }
            w->setState_TP(temp, press);
            double d = w->density();
            double mw = w->molecularWeight(0);
            double vbar = mw/d;
            writelog("{:10g} {:10g} {:12g} {:13.4g} {:13.4g}\n", temp, press*1.0E-5,
                   psat*1.0E-5, d, vbar);

        }

        writelog("\n\nTable of increasing Enthalpy at 1 atm\n\n");
        double dens;
        writelog("  Enthalpy,   Temperature,     x_Vapor,    Density, Entropy_mass, Gibbs_mass\n");
        w->setState_TP(298., OneAtm);
        double Hset = w->enthalpy_mass();
        double vapFrac = w->vaporFraction();
        double Tcalc = w->temperature();
        double Scalc = w->entropy_mass();
        double Gcalc = w->gibbs_mass();
        dens = w->density();
        writelog(" {:10g}, {:10g}, {:10g}, {:11.5g}, {:11.5g}, {:11.5g}\n", Hset , Tcalc, vapFrac, dens, Scalc, Gcalc);
        w->setState_HP(Hset, OneAtm);
        vapFrac = w->vaporFraction();
        Tcalc = w->temperature();
        dens = w->density();
        Scalc = w->entropy_mass();
        Gcalc = w->gibbs_mass();
        writelog(" {:10g}, {:10g}, {:10g}, {:11.5g}, {:11.5g}, {:11.5g}\n", Hset , Tcalc, vapFrac, dens, Scalc, Gcalc);

        double deltaH = 100000.;
        for (int i = 0; i < 40; i++) {
            Hset += deltaH;
            w->setState_HP(Hset, OneAtm);
            vapFrac = w->vaporFraction();
            Tcalc = w->temperature();
            dens = w->density();
            Scalc = w->entropy_mass();
            Gcalc = w->gibbs_mass();
            writelog(" {:10g}, {:10g}, {:10g}, {:11.5g}, {:11.5g}, {:11.5g}\n", Hset , Tcalc, vapFrac, dens, Scalc, Gcalc);
        }


        writelog("Critical Temp     = {:10.3g} K\n", w->critTemperature());
        writelog("Critical Pressure = {:10.3g} atm\n", w->critPressure()/OneAtm);
        writelog("Critical Dens     = {:10.3g} kg/m3\n", w->critDensity());

    } catch (CanteraError& err) {
        writelog(err.what());
        Cantera::appdelete();
        return -1;
    }


    return 0;
}
