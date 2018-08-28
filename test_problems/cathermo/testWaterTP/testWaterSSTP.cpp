#include "cantera/thermo/ThermoFactory.h"
#include <iostream>

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
        ThermoPhase* w = newPhase("liquid-water.xml");

        /*
         * Print out the triple point conditions
         */
        double temp = 273.16;
        pres = w->satPressure(temp);
        printf("psat(%g) = %g\n", temp, pres);
        double presLow = 1.0E-2;
        temp = 298.15;
        double oneBar = 1.0E5;
        double vol;

        printf("Comparisons to NIST: (see http://webbook.nist.gov):\n\n");

        w->setDensity(1.0E-8);
        w->setState_TP(temp, presLow);
        double h = w->enthalpy_mole();
        printf("H0(298.15) = %g J/kmol\n", h);
        double h298 = h;

        double s = w->entropy_mole();
        s -= GasConstant * log(oneBar/presLow);
        printf("S0(298.15) = %g J/kmolK\n", s);


        double T[20];
        T[0] = 298.15;
        T[1] = 500.;
        T[2] = 600.;
        T[3] = 1000.;

        double Cp0, delh0, delg0, g;
        double Cp0_ss;
        printf("\nIdeal Gas Standard State:\n");
        printf("        T      Cp0           S0     "
               " -(G0-H298)/T       H0-H298\n");
        printf("       (K)   (J/molK)     (J/molK)  "
               "   (J/molK)        (kJ/mol)\n");
        for (int i = 0; i < 4; i++) {
            temp = T[i];
            w->setState_TP(temp, presLow);
            h = w->enthalpy_mole();
            delh0 = tvalue(h - h298, 1.0E-6);
            g = w->gibbs_mole();
            delg0 = (g - h298)/temp + GasConstant * log(oneBar/presLow);
            Cp0 = w->cp_mole();
            {
                w->getCp_R(&Cp0_ss);
                Cp0_ss *= GasConstant;
                if (fabs(Cp0_ss - Cp0) > 1.0E-5) {
                    printf("Inconsistency!\n");
                    exit(-1);
                }
            }
            s = w->entropy_mole();
            s -= GasConstant * log(oneBar/presLow);
            printf("%10g %10g %13.4f %13.4f %13.4f\n", temp, Cp0*1.0E-3, s*1.0E-3,
                   -delg0*1.0E-3, delh0*1.0E-6);
        }
        printf("\n\n");

        temp = 298.15;
        w->setDensity(1000.);
        w->setState_TP(temp, oneBar);
        h = w->enthalpy_mole();
        printf("H_liq(298.15, onebar) = %g J/kmol\n", h);
        double h298l = h;
        s = w->entropy_mole();
        printf("S_liq(298.15, onebar) = %g J/kmolK\n", s);


        T[0] = 273.19;
        T[1] = 298.15;
        T[2] = 300.;
        T[3] = 373.15;
        T[4] = 400.;
        T[5] = 500.;
        printf("\nLiquid 1bar or psat Standard State\n");
        printf("       T     press         psat            Cp0            S0   "
               "  -(G0-H298)/T       H0-H298\n");
        printf("      (K)     (bar)        (bar)        (J/molK)       (J/molK)"
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
            printf("%10g %10g %12g %13.4f %13.4f %13.4f %13.4f\n", temp, press*1.0E-5,
                   psat*1.0E-5,
                   Cp0*1.0E-3, s*1.0E-3,
                   -delg0*1.0E-3, delh0*1.0E-6);
        }

        printf("\nLiquid Densities:\n");
        printf("       T     press         psat        Density          molarVol   "
               "\n");
        printf("      (K)     (bar)        (bar)      (kg/m3)          (m3/kmol)"
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
            printf("%10g %10g %12g %13.4f %13.4f\n", temp, press*1.0E-5,
                   psat*1.0E-5, d, vbar);

        }


        printf("\nLiquid 1bar or psat State: Partial Molar Quantities\n");
        printf("       T     press         psat            Cpbar          Sbar "
               "  -(G0-H298)/T       H0-H298       Volume\n");
        printf("      (K)     (bar)        (bar)        (J/molK)       (J/molK)"
               "     (J/molK)        (kJ/mol)     m3/kmol\n");

        for (int i = 0; i < 6; i++) {
            temp = T[i];
            double psat = w->satPressure(temp);
            double press = oneBar;
            if (psat > press) {
                press = psat*1.002;
            }
            w->setState_TP(temp, press);
            w->getPartialMolarEnthalpies(&h);
            delh0 = tvalue(h - h298l, 1.0E-6);
            w->getChemPotentials(&g);
            delg0 = (g - h298l)/temp;
            w->getPartialMolarCp(&Cp0);
            w->getPartialMolarEntropies(&s);
            w->getPartialMolarVolumes(&vol);
            printf("%10g %10g %12g %13.4f %13.4f %13.4f %13.4f %13.4f\n", temp, press*1.0E-5,
                   psat*1.0E-5,
                   Cp0*1.0E-3, s*1.0E-3,
                   -delg0*1.0E-3, delh0*1.0E-6, vol);
        }

        printf("\nLiquid 1bar or psat State: Standard State Quantities\n");
        printf("       T     press         psat            Cpbar          Sbar "
               "  -(G0-H298)/T       H0-H298       Volume\n");
        printf("      (K)     (bar)        (bar)        (J/molK)       (J/molK)"
               "     (J/molK)        (kJ/mol)     m3/kmol\n");

        for (int i = 0; i < 6; i++) {
            temp = T[i];
            double psat = w->satPressure(temp);
            double press = oneBar;
            if (psat > press) {
                press = psat*1.002;
            }
            w->setState_TP(temp, press);
            w->getEnthalpy_RT(&h);
            h *= temp * GasConstant;
            delh0 = tvalue(h - h298l, 1.0E-6);
            w->getStandardChemPotentials(&g);
            delg0 = (g - h298l)/temp;
            w->getCp_R(&Cp0);
            Cp0 *= GasConstant;
            w->getEntropy_R(&s);
            s *= GasConstant;
            w->getStandardVolumes(&vol);
            printf("%10g %10g %12g %13.4f %13.4f %13.4f %13.4f %13.4f\n", temp, press*1.0E-5,
                   psat*1.0E-5,
                   Cp0*1.0E-3, s*1.0E-3,
                   -delg0*1.0E-3, delh0*1.0E-6, vol);
        }

        printf("\nLiquid 1bar or psat State: Reference State Quantities (Always 1 atm no matter what system pressure is)\n");
        printf("       T     press         psat            Cpbar          Sbar "
               "  -(G0-H298)/T       H0-H298       Volume\n");
        printf("      (K)     (bar)        (bar)        (J/molK)       (J/molK)"
               "     (J/molK)        (kJ/mol)     m3/kmol\n");

        for (int i = 0; i < 6; i++) {
            temp = T[i];
            double psat = w->satPressure(temp);
            double press = oneBar;
            if (psat > press) {
                press = psat*1.002;
            }
            w->setState_TP(temp, press);
            w->getEnthalpy_RT_ref(&h);
            h *= temp * GasConstant;
            delh0 = tvalue(h - h298l, 1.0E-6);
            w->getGibbs_ref(&g);
            delg0 = (g - h298l)/temp;
            w->getCp_R_ref(&Cp0);
            Cp0 *= GasConstant;
            w->getEntropy_R_ref(&s);
            s *= GasConstant;
            w->getStandardVolumes_ref(&vol);
            printf("%10g %10g %12g %13.4f %13.4f %13.4f %13.4f %13.4f\n", temp, press*1.0E-5,
                   psat*1.0E-5,
                   Cp0*1.0E-3, s*1.0E-3,
                   -delg0*1.0E-3, delh0*1.0E-6, vol);
        }

        printf("\nLiquid 1 atm: Standard State Quantities - Should agree with table above\n");
        printf("       T     press         psat            Cpbar          Sbar "
               "  -(G0-H298)/T       H0-H298       Volume\n");
        printf("      (K)     (bar)        (bar)        (J/molK)       (J/molK)"
               "     (J/molK)        (kJ/mol)     m3/kmol\n");

        for (int i = 0; i < 6; i++) {
            temp = T[i];
            double psat = w->satPressure(temp);
            double press = OneAtm;
            w->setState_TP(temp, press);
            w->getEnthalpy_RT(&h);
            h *= temp * GasConstant;
            delh0 = tvalue(h - h298l, 1.0E-6);
            w->getStandardChemPotentials(&g);
            delg0 = (g - h298l)/temp;
            w->getCp_R(&Cp0);
            Cp0 *= GasConstant;
            w->getEntropy_R(&s);
            s *= GasConstant;
            w->getStandardVolumes(&vol);
            printf("%10g %10g %12g %13.4f %13.4f %13.4f %13.4f %13.4f\n", temp, press*1.0E-5,
                   psat*1.0E-5,
                   Cp0*1.0E-3, s*1.0E-3,
                   -delg0*1.0E-3, delh0*1.0E-6, vol);
        }
        printf("\n\nTable of increasing Enthalpy at 1 atm\n\n");
        double dens;
        printf("  Enthalpy,   Temperature,     x_Vapor,    Density, Entropy_mass, Gibbs_mass\n");
        w->setState_TP(298., OneAtm);
        double Hset = w->enthalpy_mass();
        double vapFrac = w->vaporFraction();
        double Tcalc = w->temperature();
        double Scalc = w->entropy_mass();
        double Gcalc = w->gibbs_mass();
        dens = w->density();
        printf(" %10g, %10g, %10g, %11.5g, %11.5g, %11.5g\n", Hset , Tcalc, vapFrac, dens, Scalc, Gcalc);
        w->setState_HP(Hset, OneAtm);
        vapFrac = w->vaporFraction();
        Tcalc = w->temperature();
        dens = w->density();
        Scalc = w->entropy_mass();
        Gcalc = w->gibbs_mass();
        printf(" %10g, %10g, %10g, %11.5g, %11.5g, %11.5g\n", Hset , Tcalc, vapFrac, dens, Scalc, Gcalc);

        double deltaH = 100000.;
        for (int i = 0; i < 40; i++) {
            Hset += deltaH;
            try {
                w->setState_HP(Hset, OneAtm);
            } catch (CanteraError&) {
                printf(" %10g, -> Failed to converge, beyond the spinodal probably \n\n", Hset);
                break;
            }
            vapFrac = w->vaporFraction();
            Tcalc = w->temperature();
            dens = w->density();
            Scalc = w->entropy_mass();
            Gcalc = w->gibbs_mass();
            printf(" %10g, %10g, %10g, %11.5g, %11.5g, %11.5g\n", Hset , Tcalc, vapFrac, dens, Scalc, Gcalc);
        }



        printf("Critical Temp     = %10.3g K\n", w->critTemperature());
        printf("Critical Pressure = %10.3g atm\n", w->critPressure()/OneAtm);
        printf("Critical Dens     = %10.3g kg/m3\n", w->critDensity());

        delete w;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        Cantera::appdelete();
        return -1;
    }


    return 0;
}
