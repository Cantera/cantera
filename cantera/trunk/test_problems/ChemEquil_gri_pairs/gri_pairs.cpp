/*
 *  Copyright 2002 California Institute of Technology
 */

#include "cantera/IdealGasMix.h"
#include "cantera/equilibrium.h"

#include <cstdio>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    try {
        IdealGasMix g("gri30.xml", "gri30_mix");
        double pres = 1.0E5;
        g.setState_TPX(1500.0, pres, "CH4:0.3, O2:0.3, N2:0.4");
        equilibrate(g, "TP", -1);
        //cout << g;
        double enth = g.enthalpy_mass();
        printf(" enth = %g\n", enth);
        enth -= 1.0E5;
        printf("attempted equil at (H,P) = %10.5g, %10.5g\n", enth, pres);
        g.setState_HP(enth, pres);
        equilibrate(g, "HP", -1);
        //cout << g;

        double entrop = g.entropy_mass();
        printf(" entropy = %g\n", entrop);
        entrop += 1.0E2;
        printf("attempted equil at (S,P) = %10.5g, %10.5g\n", entrop, pres);
        g.setState_SP(entrop, pres);
        equilibrate(g, "SP", -1);
        //cout << g;

        double dens = g.density();
        printf(" dens = %g\n", dens);
        dens *= 0.9;
        double vol = 1.0/dens;
        printf("attempted equil at (S,V,dens) = %10.5g, %10.5g, %10.5g\n", entrop, vol, 1.0/vol);
        g.setState_SV(entrop, vol);
        equilibrate(g, "SV", -1);
        //cout << g;

        double temp = 1800.;
        printf("attempted equil at (T,V, dens) = %10.5g, %10.5g, %10.5g\n", temp, vol, 1.0/vol);
        g.setTemperature(temp);
        equilibrate(g, "TV", -1);
        //cout << g;

        double inte = g.intEnergy_mass();
        printf(" inte = %g\n", inte);
        inte -= 1.0E5;
        printf("attempted equil at (U,V,dens) = %10.5g, %10.5g, %10.5g\n", inte, vol, 1.0/vol);
        g.setState_UV(inte, vol);
        equilibrate(g, "UV", -1);
        //cout << g;



        return 0;
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        cerr << "program terminating." << endl;
        return -1;
    }
}
