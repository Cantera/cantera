// An open Rankine cycle

#include "cantera/thermo/PureFluidPhase.h"

using namespace Cantera;

std::map<std::string,double> h, s, T, P, x;
std::vector<std::string> states;

template<class F>
void saveState(F& fluid, std::string name)
{
    h[name] = fluid.enthalpy_mass();
    s[name] = fluid.entropy_mass();
    T[name] = fluid.temperature();
    P[name] = fluid.pressure();
    x[name] = fluid.vaporFraction();
    states.push_back(name);
}

void printStates()
{
    int nStates = states.size();
    for (int n = 0; n < nStates; n++) {
        std::string name = states[n];
        writelog(" {:5s} {:10.6g} {:10.6g} {:12.6g} {:12.6g} {:5.2g}\n",
                 name, T[name], P[name], h[name], s[name], x[name]);
    }
}

int openRankine()
{
    double etap = 0.6; // pump isentropic efficiency
    double etat = 0.8; // turbine isentropic efficiency
    double phigh = 8.0e5; // high pressure

    PureFluidPhase w;
    w.initThermoFile("liquidvapor.yaml", "water");

    // begin with water at 300 K, 1 atm
    w.setState_TP(300.0, OneAtm);
    saveState(w,"1");

    // pump water to 0.8 MPa
    w.setState_SP(s["1"], phigh);
    saveState(w,"2s");
    double h2 = (h["2s"] - h["1"])/etap + h["1"];
    w.setState_HP(h2, phigh);
    saveState(w,"2");

    // heat to saturated vapor
    w.setState_Psat(phigh, 1.0);
    saveState(w,"3");

    // expand to 1 atm
    w.setState_SP(s["3"], OneAtm);
    saveState(w,"4s");
    double work_s = h["3"] - h["4s"];
    double work = etat*work_s;
    w.setState_HP(h["3"] - work, OneAtm);
    saveState(w,"4");

    printStates();

    double heat_in = h["3"] - h["2"];
    double efficiency = work/heat_in;

    writelog("efficiency = {:8.6g}\n", efficiency);
    return 0;
}

int main()
{
    try {
        return openRankine();
    } catch (CanteraError& err) {
        writelog(err.what());
        return -1;
    }
}
