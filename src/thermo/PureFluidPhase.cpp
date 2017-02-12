/**
 * @file PureFluidPhase.cpp Definitions for a ThermoPhase object for a pure
 *     fluid phase consisting of gas, liquid, mixed-gas-liquid and supercritical
 *     fluid (see \ref thermoprops and class \link Cantera::PureFluidPhase
 *     PureFluidPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

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
    m_subflag(0),
    m_mw(-1.0),
    m_verbose(false)
{
}

void PureFluidPhase::initThermo()
{
    if (m_tpx_name != "") {
        m_sub.reset(tpx::newSubstance(m_tpx_name));
    } else {
        m_sub.reset(tpx::GetSub(m_subflag));
    }
    if (!m_sub) {
        throw CanteraError("PureFluidPhase::initThermo",
                           "could not create new substance object.");
    }
    m_mw = m_sub->MolWt();
    setMolecularWeight(0,m_mw);
    double one = 1.0;
    setMoleFractions(&one);
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

    m_spthermo->update_single(0, T0, &cp0_R, &h0_RT, &s0_R);
    double s_R = s0_R - log(p/refPressure());
    m_sub->setStdState(h0_RT*GasConstant*298.15/m_mw,
                       s_R*GasConstant/m_mw, T0, p);
    debuglog("PureFluidPhase::initThermo: initialized phase "
             +id()+"\n", m_verbose);
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

doublereal PureFluidPhase::enthalpy_mole() const
{
    setTPXState();
    return m_sub->h() * m_mw;
}

doublereal PureFluidPhase::intEnergy_mole() const
{
    setTPXState();
    return m_sub->u() * m_mw;
}

doublereal PureFluidPhase::entropy_mole() const
{
    setTPXState();
    return m_sub->s() * m_mw;
}

doublereal PureFluidPhase::gibbs_mole() const
{
    setTPXState();
    return m_sub->g() * m_mw;
}

doublereal PureFluidPhase::cp_mole() const
{
    setTPXState();
    return m_sub->cp() * m_mw;
}

doublereal PureFluidPhase::cv_mole() const
{
    setTPXState();
    return m_sub->cv() * m_mw;
}

doublereal PureFluidPhase::pressure() const
{
    setTPXState();
    return m_sub->P();
}

void PureFluidPhase::setPressure(doublereal p)
{
    Set(tpx::PropertyPair::TP, temperature(), p);
    setDensity(1.0/m_sub->v());
}

void PureFluidPhase::Set(tpx::PropertyPair::type n, double x, double y) const
{
    m_sub->Set(n, x, y);
}

void PureFluidPhase::setTPXState() const
{
    Set(tpx::PropertyPair::TV, temperature(), 1.0/density());
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
    double pref = m_spthermo->refPressure();
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
    double pref = m_spthermo->refPressure();
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
    setTPXState();
    return m_sub->x();
}

void PureFluidPhase::setState_Tsat(doublereal t, doublereal x)
{
    setTemperature(t);
    setTPXState();
    Set(tpx::PropertyPair::TX, t, x);
    setDensity(1.0/m_sub->v());
}

void PureFluidPhase::setState_Psat(doublereal p, doublereal x)
{
    setTPXState();
    Set(tpx::PropertyPair::PX, p, x);
    setTemperature(m_sub->Temp());
    setDensity(1.0/m_sub->v());
}

std::string PureFluidPhase::report(bool show_thermo, doublereal threshold) const
{
    fmt::MemoryWriter b;
    if (name() != "") {
        b.write("\n  {}:\n", name());
    }
    b.write("\n");
    b.write("       temperature    {:12.6g}  K\n", temperature());
    b.write("          pressure    {:12.6g}  Pa\n", pressure());
    b.write("           density    {:12.6g}  kg/m^3\n", density());
    b.write("  mean mol. weight    {:12.6g}  amu\n", meanMolecularWeight());
    b.write("    vapor fraction    {:12.6g}\n", vaporFraction());

    doublereal phi = electricPotential();
    if (phi != 0.0) {
        b.write("         potential    {:12.6g}  V\n", phi);
    }
    if (show_thermo) {
        b.write("\n");
        b.write("                          1 kg            1 kmol\n");
        b.write("                       -----------      ------------\n");
        b.write("          enthalpy    {:12.6g}     {:12.4g}     J\n",
                enthalpy_mass(), enthalpy_mole());
        b.write("   internal energy    {:12.6g}     {:12.4g}     J\n",
                intEnergy_mass(), intEnergy_mole());
        b.write("           entropy    {:12.6g}     {:12.4g}     J/K\n",
                entropy_mass(), entropy_mole());
        b.write("    Gibbs function    {:12.6g}     {:12.4g}     J\n",
                gibbs_mass(), gibbs_mole());
        b.write(" heat capacity c_p    {:12.6g}     {:12.4g}     J/K\n",
                cp_mass(), cp_mole());
        try {
            b.write(" heat capacity c_v    {:12.6g}     {:12.4g}     J/K\n",
                    cv_mass(), cv_mole());
        } catch (NotImplementedError&) {
            b.write(" heat capacity c_v    <not implemented>\n");
        }
    }
    return b.str();
}

}
