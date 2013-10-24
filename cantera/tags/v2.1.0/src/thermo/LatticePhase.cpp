/**
 *
 *  @file LatticePhase.cpp
 *  Definitions for a simple thermodynamics model of a bulk phase
 *  derived from ThermoPhase,
 *  assuming a lattice of solid atoms
 *  (see \ref thermoprops and class \link Cantera::LatticePhase LatticePhase\endlink).
 *
 */
#include "cantera/base/config.h"
#include "cantera/base/ct_defs.h"
#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/LatticePhase.h"
#include "cantera/thermo/SpeciesThermo.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{

LatticePhase::LatticePhase() :
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm),
    m_tlast(0.0),
    m_speciesMolarVolume(0),
    m_site_density(0.0)
{
}

LatticePhase::LatticePhase(const LatticePhase& right) :
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm),
    m_tlast(0.0),
    m_speciesMolarVolume(0),
    m_site_density(0.0)
{
    *this = operator=(right);
}

LatticePhase& LatticePhase::operator=(const LatticePhase& right)
{
    if (&right != this) {
        ThermoPhase::operator=(right);
        m_Pref       = right.m_Pref;
        m_Pcurrent     = right.m_Pcurrent;
        m_tlast      = right.m_tlast;
        m_h0_RT      = right.m_h0_RT;
        m_cp0_R      = right.m_cp0_R;
        m_g0_RT      = right.m_g0_RT;
        m_s0_R       = right.m_s0_R;
        m_vacancy    = right.m_vacancy;
        m_speciesMolarVolume = right.m_speciesMolarVolume;
        m_site_density = right.m_site_density;
    }
    return *this;
}

LatticePhase::LatticePhase(const std::string& inputFile, const std::string& id_)
{
    initThermoFile(inputFile, id_);
}

LatticePhase::LatticePhase(XML_Node& phaseRef, const std::string& id_)
{
    importPhase(*findXMLPhase(&phaseRef, id_), this);
}

ThermoPhase* LatticePhase::duplMyselfAsThermoPhase() const
{
    return new LatticePhase(*this);
}

doublereal LatticePhase::enthalpy_mole() const
{
    doublereal p0 = m_spthermo->refPressure();
    return GasConstant * temperature() *
           mean_X(&enthalpy_RT_ref()[0])
           + (pressure() - p0)/molarDensity();
}

doublereal LatticePhase::intEnergy_mole() const
{
    doublereal p0 = m_spthermo->refPressure();
    return GasConstant * temperature() *
           mean_X(&enthalpy_RT_ref()[0])
           - p0/molarDensity();
}

doublereal LatticePhase::entropy_mole() const
{
    return GasConstant * (mean_X(&entropy_R_ref()[0]) -
                          sum_xlogx());
}

doublereal LatticePhase::gibbs_mole() const
{
    return enthalpy_mole() - temperature() * entropy_mole();
}

doublereal LatticePhase::cp_mole() const
{
    return GasConstant * mean_X(&cp_R_ref()[0]);
}

doublereal LatticePhase::cv_mole() const
{
    return cp_mole();
}

doublereal LatticePhase::calcDensity()
{
    setMolarDensity(m_site_density);
    doublereal mw = meanMolecularWeight();
    doublereal dens = mw * m_site_density;
    /*
     * Calculate the molarVolume of the solution (m**3 kmol-1)
     */
    // const doublereal * const dtmp = moleFractdivMMW();
    // doublereal invDens = dot(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), dtmp);
    /*
     * Set the density in the parent State object directly,
     * by calling the Phase::setDensity() function.
     */
    // doublereal dens = 1.0/invDens;
    //  Phase::setDensity(dens);
    return dens;
}

void LatticePhase::setPressure(doublereal p)
{
    m_Pcurrent = p;
    calcDensity();
}

void LatticePhase::setMoleFractions(const doublereal* const x)
{
    Phase::setMoleFractions(x);
    calcDensity();
}

void LatticePhase::setMoleFractions_NoNorm(const doublereal* const x)
{
    Phase::setMoleFractions(x);
    calcDensity();
}

void LatticePhase::setMassFractions(const doublereal* const y)
{
    Phase::setMassFractions(y);
    calcDensity();
}

void LatticePhase::setMassFractions_NoNorm(const doublereal* const y)
{
    Phase::setMassFractions_NoNorm(y);
    calcDensity();
}

void LatticePhase::setConcentrations(const doublereal* const c)
{
    Phase::setConcentrations(c);
    calcDensity();
}

void LatticePhase::getActivityConcentrations(doublereal* c) const
{
    getMoleFractions(c);
}

void LatticePhase::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

doublereal LatticePhase::standardConcentration(size_t k) const
{
    return 1.0;
}

doublereal LatticePhase::logStandardConc(size_t k) const
{
    return 0.0;
}

void LatticePhase::getChemPotentials(doublereal* mu) const
{
    doublereal delta_p = m_Pcurrent - m_Pref;
    doublereal xx;
    doublereal RT = temperature() * GasConstant;
    const vector_fp& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(SmallNumber, moleFraction(k));
        mu[k] = RT * (g_RT[k] + log(xx))
                + delta_p * m_speciesMolarVolume[k];
    }

}

void LatticePhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal rt = GasConstant * temperature();
    scale(_h.begin(), _h.end(), hbar, rt);
}

void LatticePhase::getPartialMolarEntropies(doublereal* sbar) const
{
    const vector_fp& _s = entropy_R_ref();
    doublereal r = GasConstant;
    doublereal xx;
    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] = r * (_s[k] - log(xx));
    }
}

void LatticePhase::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

void LatticePhase::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
}

void LatticePhase::getStandardChemPotentials(doublereal* mu0) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), mu0, _RT());
}

void LatticePhase::getPureGibbs(doublereal* gpure) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    doublereal delta_p = (m_Pcurrent - m_Pref);
    double RT = GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        gpure[k] = RT * gibbsrt[k] + delta_p * m_speciesMolarVolume[k];
    }
}

void LatticePhase::getEnthalpy_RT(doublereal* hrt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal delta_prt = ((m_Pcurrent - m_Pref) / (GasConstant * temperature()));
    for (size_t k = 0; k < m_kk; k++) {
        hrt[k] = _h[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

void LatticePhase::getEntropy_R(doublereal* sr) const
{
    const vector_fp& _s = entropy_R_ref();
    std::copy(_s.begin(), _s.end(), sr);
}

void LatticePhase::getGibbs_RT(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    doublereal RT = _RT();
    doublereal delta_prt = (m_Pcurrent - m_Pref)/ RT;
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] = gibbsrt[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

void LatticePhase::getGibbs_ref(doublereal* g) const
{
    getGibbs_RT_ref(g);
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= GasConstant * temperature();
    }
}

void LatticePhase::getCp_R(doublereal* cpr) const
{
    const vector_fp& _cpr = cp_R_ref();
    std::copy(_cpr.begin(), _cpr.end(), cpr);
}

void LatticePhase::getStandardVolumes(doublereal* vbar) const
{
    copy(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), vbar);
}

const vector_fp& LatticePhase::enthalpy_RT_ref() const
{
    _updateThermo();
    return m_h0_RT;
}

const vector_fp& LatticePhase::gibbs_RT_ref() const
{
    _updateThermo();
    return m_g0_RT;
}

void LatticePhase::getGibbs_RT_ref(doublereal* grt) const
{
    _updateThermo();
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] = m_g0_RT[k];
    }
}

const vector_fp& LatticePhase::entropy_R_ref() const
{
    _updateThermo();
    return m_s0_R;
}

const vector_fp& LatticePhase::cp_R_ref() const
{
    _updateThermo();
    return m_cp0_R;
}

void LatticePhase::initThermo()
{
    m_Pref = refPressure();
    size_t leng = m_kk;
    m_h0_RT.resize(leng);
    m_g0_RT.resize(leng);
    m_cp0_R.resize(leng);
    m_s0_R.resize(leng);
    m_speciesMolarVolume.resize(leng, 0.0);

    ThermoPhase::initThermo();
}

void LatticePhase::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    std::string idattrib = phaseNode.id();
    if (!id_.empty() && id_ != idattrib) {
        throw CanteraError("LatticePhase::initThermoXML",
                           "ids don't match");
    }

    std::string subname = "LatticePhase::initThermoXML";
    /*
     * Check on the thermo field. Must have:
     * <thermo model="Lattice" />
     */
    if (phaseNode.hasChild("thermo")) {
        XML_Node& thNode = phaseNode.child("thermo");
        std::string mStringa = thNode.attrib("model");
        std::string mString = lowercase(mStringa);
        if (mString != "lattice") {
            throw CanteraError(subname.c_str(),
                               "Unknown thermo model: " + mStringa);
        }
    } else {
        throw CanteraError(subname.c_str(),
                           "Unspecified thermo model");
    }
    /*
     * Now go get the molar volumes. use the default if not found
     */
    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"], &phaseNode.root());
    const std::vector<std::string> &sss = speciesNames();

    for (size_t k = 0; k < m_kk; k++) {
        m_speciesMolarVolume[k] = m_site_density;
        XML_Node* s =  speciesDB->findByAttr("name", sss[k]);
        if (!s) {
            throw CanteraError(" LatticePhase::initThermoXML", "database problems");
        }
        XML_Node* ss = s->findByName("standardState");
        if (ss) {
            if (ss->findByName("molarVolume")) {
                m_speciesMolarVolume[k] = ctml::getFloat(*ss, "molarVolume", "toSI");
            }
        }
    }

    /*
     * Call the base initThermo, which handles setting the initial
     * state.
     */
    ThermoPhase::initThermoXML(phaseNode, id_);
}

void LatticePhase::_updateThermo() const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow) {
        m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0], &m_s0_R[0]);
        m_tlast = tnow;
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_tlast = tnow;
    }
}

void LatticePhase::setParameters(int n, doublereal* const c)
{
    warn_deprecated("LatticePhase::setParameters");
    m_site_density = c[0];
    setMolarDensity(m_site_density);
}

void LatticePhase::getParameters(int& n, doublereal* const c) const
{
    warn_deprecated("LatticePhase::getParameters");
    double d = molarDensity();
    c[0] = d;
    n = 1;
}

void LatticePhase::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model", "Lattice");
    m_site_density = ctml::getFloat(eosdata, "site_density", "toSI");
    m_vacancy = ctml::getChildValue(eosdata, "vacancy_species");
}

}
