/**
 *   @file PureFluidPhase.cpp
 *   Definitions for a ThermoPhase object for a pure fluid phase consisting
 *   of gas, liquid, mixed-gas-liquid
 *   and supercritical fluid (see \ref thermoprops
 *   and class \link Cantera::PureFluidPhase PureFluidPhase\endlink).
 */
#include "cantera/base/xml.h"
#include "cantera/thermo/PureFluidPhase.h"

#include "cantera/tpx/Sub.h"
#include "cantera/tpx/utils.h"

#include <cstdio>

using std::string;
using std::endl;

namespace Cantera
{

PureFluidPhase::PureFluidPhase() :
    ThermoPhase(),
    m_sub(0),
    m_subflag(0),
    m_mw(-1.0),
    m_verbose(false)
{
}

PureFluidPhase::PureFluidPhase(const PureFluidPhase& right) :
    ThermoPhase(),
    m_sub(0),
    m_subflag(0),
    m_mw(-1.0),
    m_verbose(false)
{
    *this = right;
}

PureFluidPhase& PureFluidPhase::operator=(const PureFluidPhase& right)
{
    if (&right != this) {
        ThermoPhase::operator=(right);
        delete m_sub;
        m_subflag    = right.m_subflag;
        m_sub        = tpx::GetSub(m_subflag);
        m_mw         = right.m_mw;
        m_verbose    = right.m_verbose;
    }
    return *this;
}

ThermoPhase* PureFluidPhase::duplMyselfAsThermoPhase() const
{
    return new PureFluidPhase(*this);
}

PureFluidPhase::~PureFluidPhase()
{
    delete m_sub;
}

void PureFluidPhase::
initThermo()
{
    delete m_sub;
    m_sub = tpx::GetSub(m_subflag);
    if (m_sub == 0) {
        throw CanteraError("PureFluidPhase::initThermo",
                           "could not create new substance object.");
    }
    m_mw = m_sub->MolWt();
    setMolecularWeight(0,m_mw);
    double one = 1.0;
    setMoleFractions(&one);
    double cp0_R, h0_RT, s0_R, T0, p;
    T0 = 298.15;
    if (T0 < m_sub->Tcrit()) {
        m_sub->Set(tpx::PropertyPair::TX, T0, 1.0);
        p = 0.01*m_sub->P();
    } else {
        p = 0.001*m_sub->Pcrit();
    }
    p = 0.001 * p;
    m_sub->Set(tpx::PropertyPair::TP, T0, p);

    m_spthermo->update_one(0, T0, &cp0_R, &h0_RT, &s0_R);
    double s_R = s0_R - log(p/refPressure());
    m_sub->setStdState(h0_RT*GasConstant*298.15/m_mw,
                       s_R*GasConstant/m_mw, T0, p);
    writelog("PureFluidPhase::initThermo: initialized phase "
             +id()+"\n", m_verbose);
}

void PureFluidPhase::
setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","PureFluid");
    m_subflag = atoi(eosdata["fluid_type"].c_str());
    if (m_subflag < 0)
        throw CanteraError("PureFluidPhase::setParametersFromXML",
                           "missing or negative substance flag");
}

doublereal PureFluidPhase::
enthalpy_mole() const
{
    setTPXState();
    return m_sub->h() * m_mw;
}

doublereal PureFluidPhase::
intEnergy_mole() const
{
    setTPXState();
    return m_sub->u() * m_mw;
}

doublereal PureFluidPhase::
entropy_mole() const
{
    setTPXState();
    return m_sub->s() * m_mw;
}

doublereal PureFluidPhase::
gibbs_mole() const
{
    setTPXState();
    return m_sub->g() * m_mw;
}

doublereal PureFluidPhase::
cp_mole() const
{
    setTPXState();
    return m_sub->cp() * m_mw;
}

doublereal PureFluidPhase::
cv_mole() const
{
    setTPXState();
    return m_sub->cv() * m_mw;
}

doublereal PureFluidPhase::
pressure() const
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

void  PureFluidPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    hbar[0] = enthalpy_mole();
}

void  PureFluidPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    sbar[0] = entropy_mole();
}

void  PureFluidPhase::getPartialMolarIntEnergies(doublereal* ubar) const
{
    ubar[0] = intEnergy_mole();
}

void  PureFluidPhase::getPartialMolarCp(doublereal* cpbar) const
{
    cpbar[0] = cp_mole();
}

void  PureFluidPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    vbar[0] = 1.0 / molarDensity();
}

int PureFluidPhase::standardStateConvention() const
{
    return cSS_CONVENTION_TEMPERATURE;
}

void  PureFluidPhase::getActivityConcentrations(doublereal* c) const
{
    c[0] = 1.0;
}

doublereal PureFluidPhase::standardConcentration(size_t k) const
{
    return 1.0;
}

void  PureFluidPhase::getActivities(doublereal* a) const
{
    a[0] = 1.0;
}

void  PureFluidPhase::getStandardChemPotentials(doublereal* mu) const
{
    mu[0] = gibbs_mole();
}

void PureFluidPhase::getEnthalpy_RT(doublereal* hrt) const
{
    doublereal rt = _RT();
    doublereal h = enthalpy_mole();
    hrt[0] = h / rt;
}

void PureFluidPhase::getEntropy_R(doublereal* sr) const
{
    doublereal s = entropy_mole();
    sr[0] = s / GasConstant;
}

void  PureFluidPhase::getGibbs_RT(doublereal* grt) const
{
    doublereal rt = _RT();
    doublereal g = gibbs_mole();
    grt[0] = g / rt;
}

void PureFluidPhase::getEnthalpy_RT_ref(doublereal* hrt) const
{
    double psave = pressure();
    double t = temperature();
    //double pref = m_spthermo->refPressure();
    double plow = 1.0E-8;
    Set(tpx::PropertyPair::TP, t, plow);
    getEnthalpy_RT(hrt);
    Set(tpx::PropertyPair::TP, t, psave);

}

void  PureFluidPhase::getGibbs_RT_ref(doublereal* grt) const
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
    g[0] *= (GasConstant * temperature());
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

void PureFluidPhase::setState_HP(doublereal h, doublereal p,
                                 doublereal tol)
{
    Set(tpx::PropertyPair::HP, h, p);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_UV(doublereal u, doublereal v,
                                 doublereal tol)
{
    Set(tpx::PropertyPair::UV, u, v);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_SV(doublereal s, doublereal v,
                                 doublereal tol)
{
    Set(tpx::PropertyPair::SV, s, v);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_SP(doublereal s, doublereal p,
                                 doublereal tol)
{
    Set(tpx::PropertyPair::SP, s, p);
    setState_TR(m_sub->Temp(), 1.0/m_sub->v());
}

doublereal PureFluidPhase::satPressure(doublereal t) 
{
    doublereal vsv = m_sub->v();
    Set(tpx::PropertyPair::TV,t,vsv);
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

std::string PureFluidPhase::report(bool show_thermo) const
{
    char p[800];
    string s = "";
    if (name() != "") {
        sprintf(p, " \n  %s:\n", name().c_str());
        s += p;
    }
    sprintf(p, " \n       temperature    %12.6g  K\n", temperature());
    s += p;
    sprintf(p, "          pressure    %12.6g  Pa\n", pressure());
    s += p;
    sprintf(p, "           density    %12.6g  kg/m^3\n", density());
    s += p;
    sprintf(p, "  mean mol. weight    %12.6g  amu\n", meanMolecularWeight());
    s += p;
    sprintf(p, "    vapor fraction    %12.6g \n", vaporFraction());
    s += p;

    doublereal phi = electricPotential();
    if (phi != 0.0) {
        sprintf(p, "         potential    %12.6g  V\n", phi);
        s += p;
    }
    if (show_thermo) {
        sprintf(p, " \n");
        s += p;
        sprintf(p, "                          1 kg            1 kmol\n");
        s += p;
        sprintf(p, "                       -----------      ------------\n");
        s += p;
        sprintf(p, "          enthalpy    %12.6g     %12.4g     J\n",
                enthalpy_mass(), enthalpy_mole());
        s += p;
        sprintf(p, "   internal energy    %12.6g     %12.4g     J\n",
                intEnergy_mass(), intEnergy_mole());
        s += p;
        sprintf(p, "           entropy    %12.6g     %12.4g     J/K\n",
                entropy_mass(), entropy_mole());
        s += p;
        sprintf(p, "    Gibbs function    %12.6g     %12.4g     J\n",
                gibbs_mass(), gibbs_mole());
        s += p;
        sprintf(p, " heat capacity c_p    %12.6g     %12.4g     J/K\n",
                cp_mass(), cp_mole());
        s += p;
        try {
            sprintf(p, " heat capacity c_v    %12.6g     %12.4g     J/K\n",
                    cv_mass(), cv_mole());
            s += p;
        } catch (CanteraError& e) {
            e.save();
            sprintf(p, " heat capacity c_v    <not implemented>       \n");
            s += p;
        }
    }
    return s;
}

}
