/**
 * @file PureFluidPhase.cpp Definitions for a ThermoPhase object for a pure
 *     fluid phase consisting of gas, liquid, mixed-gas-liquid and supercritical
 *     fluid (see \ref thermoprops and class \link Cantera::PureFluidPhase
 *     PureFluidPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/xml.h"
#include "cantera/thermo/PureFluidPhase.h"

#include "cantera/tpx/Sub.h"
#include "cantera/tpx/utils.h"
#include "cantera/base/stringUtils.h"

#include <cstdio>

using std::string;

namespace Cantera
{

PureFluidPhase::PureFluidPhase() :
    m_subflag(-1),
    m_mw(-1.0),
    m_verbose(false)
{
}

void PureFluidPhase::initThermo()
{
    if (m_input.hasKey("pure-fluid-name")) {
        setSubstance(m_input["pure-fluid-name"].asString());
    }

    if (m_tpx_name != "") {
        m_sub.reset(tpx::newSubstance(m_tpx_name));
    } else {
        m_sub.reset(tpx::GetSub(m_subflag));
    }

    m_mw = m_sub->MolWt();
    setMolecularWeight(0,m_mw);

    double cp0_R, h0_RT, s0_R, p;
    double T0 = 298.15;
    if (T0 < m_sub->Tcrit()) {
        m_sub->Set(tpx::PropertyPair::TX, T0, 1.0);
        p = 0.01*m_sub->P();
    } else {
        p = 0.001*m_sub->Pcrit();
    }
    p = 0.001 * p;
    m_sub->Set(tpx::PropertyPair::TP, T0, p);

    m_spthermo.update_single(0, T0, &cp0_R, &h0_RT, &s0_R);
    double s_R = s0_R - log(p/refPressure());
    m_sub->setStdState(h0_RT*GasConstant*298.15/m_mw,
                       s_R*GasConstant/m_mw, T0, p);
    debuglog("PureFluidPhase::initThermo: initialized phase "
             +name()+"\n", m_verbose);
}

void PureFluidPhase::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","PureFluid");
    m_subflag = atoi(eosdata["fluid_type"].c_str());
    if (m_subflag < 0) {
        throw CanteraError("PureFluidPhase::setParametersFromXML",
                           "missing or negative substance flag");
    }
}

std::vector<std::string> PureFluidPhase::fullStates() const
{
    return {"TD", "UV", "DP", "HP", "SP", "SV",
            "ST", "TV", "PV", "UP", "VH", "TH", "SH", "TPQ"};
}

std::vector<std::string> PureFluidPhase::partialStates() const
{
    return {"TP", "TQ", "PQ"};
}

std::string PureFluidPhase::phaseOfMatter() const
{
    if (temperature() >= critTemperature() || pressure() >= critPressure()) {
        return "supercritical";
    } else if (m_sub->TwoPhase() == 1) {
        return "liquid-gas-mix";
    } else if (pressure() < m_sub->Ps()) {
        return "gas";
    } else {
        return "liquid";
    }
}

double PureFluidPhase::minTemp(size_t k) const
{
    return m_sub->Tmin();
}

double PureFluidPhase::maxTemp(size_t k) const
{
    return m_sub->Tmax();
}

doublereal PureFluidPhase::enthalpy_mole() const
{
    return m_sub->h() * m_mw;
}

doublereal PureFluidPhase::intEnergy_mole() const
{
    return m_sub->u() * m_mw;
}

doublereal PureFluidPhase::entropy_mole() const
{
    return m_sub->s() * m_mw;
}

doublereal PureFluidPhase::gibbs_mole() const
{
    return m_sub->g() * m_mw;
}

doublereal PureFluidPhase::cp_mole() const
{
    return m_sub->cp() * m_mw;
}

doublereal PureFluidPhase::cv_mole() const
{
    return m_sub->cv() * m_mw;
}

doublereal PureFluidPhase::pressure() const
{
    return m_sub->P();
}

void PureFluidPhase::setPressure(doublereal p)
{
    Set(tpx::PropertyPair::TP, temperature(), p);
    ThermoPhase::setDensity(1.0/m_sub->v());
}

void PureFluidPhase::setTemperature(double T)
{
    ThermoPhase::setTemperature(T);
    Set(tpx::PropertyPair::TV, T, m_sub->v());
}

void PureFluidPhase::setDensity(double rho)
{
    ThermoPhase::setDensity(rho);
    Set(tpx::PropertyPair::TV, m_sub->Temp(), 1.0/rho);
}

void PureFluidPhase::Set(tpx::PropertyPair::type n, double x, double y) const
{
    m_sub->Set(n, x, y);
}

doublereal PureFluidPhase::isothermalCompressibility() const
{
    return m_sub->isothermalCompressibility();
}

doublereal PureFluidPhase::thermalExpansionCoeff() const
{
    return m_sub->thermalExpansionCoeff();
}

tpx::Substance& PureFluidPhase::TPX_Substance()
{
    return *m_sub;
}

void PureFluidPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    hbar[0] = enthalpy_mole();
}

void PureFluidPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    sbar[0] = entropy_mole();
}

void PureFluidPhase::getPartialMolarIntEnergies(doublereal* ubar) const
{
    ubar[0] = intEnergy_mole();
}

void PureFluidPhase::getPartialMolarCp(doublereal* cpbar) const
{
    cpbar[0] = cp_mole();
}

void PureFluidPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    vbar[0] = 1.0 / molarDensity();
}

Units PureFluidPhase::standardConcentrationUnits() const
{
    return Units(1.0);
}

void PureFluidPhase::getActivityConcentrations(doublereal* c) const
{
    c[0] = 1.0;
}

doublereal PureFluidPhase::standardConcentration(size_t k) const
{
    return 1.0;
}

void PureFluidPhase::getActivities(doublereal* a) const
{
    a[0] = 1.0;
}

void PureFluidPhase::getStandardChemPotentials(doublereal* mu) const
{
    mu[0] = gibbs_mole();
}

void PureFluidPhase::getEnthalpy_RT(doublereal* hrt) const
{
    hrt[0] = enthalpy_mole() / RT();
}

void PureFluidPhase::getEntropy_R(doublereal* sr) const
{
    sr[0] = entropy_mole() / GasConstant;
}

void PureFluidPhase::getGibbs_RT(doublereal* grt) const
{
    grt[0] = gibbs_mole() / RT();
}

void PureFluidPhase::getEnthalpy_RT_ref(doublereal* hrt) const
{
    double psave = pressure();
    double t = temperature();
    double plow = 1.0E-8;
    Set(tpx::PropertyPair::TP, t, plow);
    getEnthalpy_RT(hrt);
    Set(tpx::PropertyPair::TP, t, psave);

}

void PureFluidPhase::getGibbs_RT_ref(doublereal* grt) const
{
    double psave = pressure();
    double t = temperature();
    double pref = refPressure();
    double plow = 1.0E-8;
    Set(tpx::PropertyPair::TP, t, plow);
    getGibbs_RT(grt);
    grt[0] += log(pref/plow);
    Set(tpx::PropertyPair::TP, t, psave);
}

void PureFluidPhase::getGibbs_ref(doublereal* g) const
{
    getGibbs_RT_ref(g);
    g[0] *= RT();
}

void PureFluidPhase::getEntropy_R_ref(doublereal* er) const
{
    double psave = pressure();
    double t = temperature();
    double pref = refPressure();
    double plow = 1.0E-8;
    Set(tpx::PropertyPair::TP, t, plow);
    getEntropy_R(er);
    er[0] -= log(pref/plow);
    Set(tpx::PropertyPair::TP, t, psave);
}

doublereal PureFluidPhase::critTemperature() const
{
    return m_sub->Tcrit();
}

doublereal PureFluidPhase::critPressure() const
{
    return m_sub->Pcrit();
}

doublereal PureFluidPhase::critDensity() const
{
    return 1.0/m_sub->Vcrit();
}

doublereal PureFluidPhase::satTemperature(doublereal p) const
{
    return m_sub->Tsat(p);
}

/* The next several functions set the state. They run the Substance::Set
 * function, followed by setting the state of the ThermoPhase object
 * to the newly computed temperature and density of the Substance.
 */

void PureFluidPhase::setState_HP(double h, double p, double tol)
{
    Set(tpx::PropertyPair::HP, h, p);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_UV(double u, double v, double tol)
{
    Set(tpx::PropertyPair::UV, u, v);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_SV(double s, double v, double tol)
{
    Set(tpx::PropertyPair::SV, s, v);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_SP(double s, double p, double tol)
{
    Set(tpx::PropertyPair::SP, s, p);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_ST(double s, double t, double tol)
{
    Set(tpx::PropertyPair::ST, s, t);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_TV(double t, double v, double tol)
{
    Set(tpx::PropertyPair::TV, t, v);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_PV(double p, double v, double tol)
{
    Set(tpx::PropertyPair::PV, p, v);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_UP(double u, double p, double tol)
{
    Set(tpx::PropertyPair::UP, u, p);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_VH(double v, double h, double tol)
{
    Set(tpx::PropertyPair::VH, v, h);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_TH(double t, double h, double tol)
{
    Set(tpx::PropertyPair::TH, t, h);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_SH(double s, double h, double tol)
{
    Set(tpx::PropertyPair::SH, s, h);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

doublereal PureFluidPhase::satPressure(doublereal t)
{
    Set(tpx::PropertyPair::TV, t, m_sub->v());
    return m_sub->Ps();
}

doublereal PureFluidPhase::vaporFraction() const
{
    return m_sub->x();
}

void PureFluidPhase::setState_Tsat(doublereal t, doublereal x)
{
    Set(tpx::PropertyPair::TX, t, x);
    ThermoPhase::setTemperature(t);
    ThermoPhase::setDensity(1.0/m_sub->v());
}

void PureFluidPhase::setState_Psat(doublereal p, doublereal x)
{
    Set(tpx::PropertyPair::PX, p, x);
    ThermoPhase::setTemperature(m_sub->Temp());
    ThermoPhase::setDensity(1.0/m_sub->v());
}

std::string PureFluidPhase::report(bool show_thermo, doublereal threshold) const
{
    fmt::memory_buffer b;
    // This is the width of the first column of names in the report.
    int name_width = 18;

    string blank_leader = fmt::format("{:{}}", "", name_width);

    string one_property = "{:>{}}   {:<.5g} {}\n";

    string two_prop_header = "{}   {:^15}   {:^15}\n";
    string kg_kmol_header = fmt::format(
        two_prop_header, blank_leader, "1 kg", "1 kmol"
    );
    string Y_X_header = fmt::format(
        two_prop_header, blank_leader, "mass frac. Y", "mole frac. X"
    );
    string two_prop_sep = fmt::format(
        "{}   {:-^15}   {:-^15}\n", blank_leader, "", ""
    );
    string two_property = "{:>{}}   {:15.5g}   {:15.5g}  {}\n";

    string three_prop_header = fmt::format(
        "{}   {:^15}   {:^15}   {:^15}\n", blank_leader, "mass frac. Y",
        "mole frac. X", "chem. pot. / RT"
    );
    string three_prop_sep = fmt::format(
        "{}   {:-^15}   {:-^15}   {:-^15}\n", blank_leader, "", "", ""
    );
    string three_property = "{:>{}}   {:15.5g}   {:15.5g}   {:15.5g}\n";

    if (name() != "") {
        format_to(b, "\n  {}:\n", name());
    }
    format_to(b, "\n");
    format_to(b, one_property, "temperature", name_width, temperature(), "K");
    format_to(b, one_property, "pressure", name_width, pressure(), "Pa");
    format_to(b, one_property, "density", name_width, density(), "kg/m^3");
    format_to(b, one_property, "mean mol. weight", name_width, meanMolecularWeight(), "kg/kmol");
    format_to(b, "{:>{}}   {:<.5g}\n", "vapor fraction", name_width, vaporFraction());
    format_to(b, "{:>{}}   {}\n", "phase of matter", name_width, phaseOfMatter());

    if (show_thermo) {
        format_to(b, "\n");
        format_to(b, kg_kmol_header);
        format_to(b, two_prop_sep);
        format_to(b, two_property, "enthalpy", name_width, enthalpy_mass(), enthalpy_mole(), "J");
        format_to(b, two_property, "internal energy", name_width, intEnergy_mass(), intEnergy_mole(), "J");
        format_to(b, two_property, "entropy", name_width, entropy_mass(), entropy_mole(), "J/K");
        format_to(b, two_property, "Gibbs function", name_width, gibbs_mass(), gibbs_mole(), "J");
        format_to(b, two_property, "heat capacity c_p", name_width, cp_mass(), cp_mole(), "J/K");
        format_to(b, two_property, "heat capacity c_v", name_width, cv_mass(), cv_mole(), "J/K");
    }

    return to_string(b);
}

}
