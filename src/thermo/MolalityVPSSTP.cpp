/**
 *  @file MolalityVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ molality based activity coefficient formulations
 *  (see \ref thermoprops
 * and class \link Cantera::MolalityVPSSTP MolalityVPSSTP\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MolalityVPSSTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/base/utilities.h"

#include <fstream>

using namespace std;

namespace Cantera
{

MolalityVPSSTP::MolalityVPSSTP() :
    m_pHScalingType(PHSCALE_PITZER),
    m_indexCLM(npos),
    m_weightSolvent(18.01528),
    m_xmolSolventMIN(0.01),
    m_Mnaught(18.01528E-3)
{
    // Change the default to be that charge neutrality in the phase is necessary
    // condition for the proper specification of thermodynamic functions within
    // the phase
    m_chargeNeutralityNecessary = true;
}

// -------------- Utilities -------------------------------

void MolalityVPSSTP::setpHScale(const int pHscaleType)
{
    m_pHScalingType = pHscaleType;
    if (pHscaleType != PHSCALE_PITZER && pHscaleType != PHSCALE_NBS) {
        throw CanteraError("MolalityVPSSTP::setpHScale",
                           "Unknown scale type: {}", pHscaleType);
    }
}

int MolalityVPSSTP::pHScale() const
{
    return m_pHScalingType;
}

void MolalityVPSSTP::setMoleFSolventMin(doublereal xmolSolventMIN)
{
    if (xmolSolventMIN <= 0.0) {
        throw CanteraError("MolalityVPSSTP::setMoleFSolventMin ", "trouble");
    } else if (xmolSolventMIN > 0.9) {
        throw CanteraError("MolalityVPSSTP::setMoleFSolventMin ", "trouble");
    }
    m_xmolSolventMIN = xmolSolventMIN;
}

doublereal MolalityVPSSTP::moleFSolventMin() const
{
    return m_xmolSolventMIN;
}

void MolalityVPSSTP::calcMolalities() const
{
    getMoleFractions(m_molalities.data());
    double xmolSolvent = std::max(m_molalities[0], m_xmolSolventMIN);
    double denomInv = 1.0/ (m_Mnaught * xmolSolvent);
    for (size_t k = 0; k < m_kk; k++) {
        m_molalities[k] *= denomInv;
    }
}

void MolalityVPSSTP::getMolalities(doublereal* const molal) const
{
    calcMolalities();
    for (size_t k = 0; k < m_kk; k++) {
        molal[k] = m_molalities[k];
    }
}

void MolalityVPSSTP::setMolalities(const doublereal* const molal)
{
    double Lsum = 1.0 / m_Mnaught;
    for (size_t k = 1; k < m_kk; k++) {
        m_molalities[k] = molal[k];
        Lsum += molal[k];
    }
    double tmp = 1.0 / Lsum;
    m_molalities[0] = tmp / m_Mnaught;
    double sum = m_molalities[0];
    for (size_t k = 1; k < m_kk; k++) {
        m_molalities[k] = tmp * molal[k];
        sum += m_molalities[k];
    }
    if (sum != 1.0) {
        tmp = 1.0 / sum;
        for (size_t k = 0; k < m_kk; k++) {
            m_molalities[k] *= tmp;
        }
    }
    setMoleFractions(m_molalities.data());

    // Essentially we don't trust the input: We calculate the molalities from
    // the mole fractions that we just obtained.
    calcMolalities();
}

void MolalityVPSSTP::setMolalitiesByName(const compositionMap& mMap)
{
    // HKM -> Might need to be more complicated here, setting neutrals so that
    //        the existing mole fractions are preserved.

    // Get a vector of mole fractions
    vector_fp mf(m_kk, 0.0);
    getMoleFractions(mf.data());
    double xmolSmin = std::max(mf[0], m_xmolSolventMIN);
    for (size_t k = 0; k < m_kk; k++) {
        double mol_k = getValue(mMap, speciesName(k), 0.0);
        if (mol_k > 0) {
            mf[k] = mol_k * m_Mnaught * xmolSmin;
        }
    }

    // check charge neutrality
    size_t largePos = npos;
    double cPos = 0.0;
    size_t largeNeg = npos;
    double cNeg = 0.0;
    double sum = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        double ch = charge(k);
        if (mf[k] > 0.0) {
            if (ch > 0.0 && ch * mf[k] > cPos) {
                largePos = k;
                cPos = ch * mf[k];
            }
            if (ch < 0.0 && fabs(ch) * mf[k] > cNeg) {
                largeNeg = k;
                cNeg = fabs(ch) * mf[k];
            }
        }
        sum += mf[k] * ch;
    }
    if (sum != 0.0) {
        if (sum > 0.0) {
            if (cPos > sum) {
                mf[largePos] -= sum / charge(largePos);
            } else {
                throw CanteraError("MolalityVPSSTP:setMolalitiesbyName",
                                   "unbalanced charges");
            }
        } else {
            if (cNeg > (-sum)) {
                mf[largeNeg] -= (-sum) / fabs(charge(largeNeg));
            } else {
                throw CanteraError("MolalityVPSSTP:setMolalitiesbyName",
                                   "unbalanced charges");
            }
        }
    }
    sum = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum += mf[k];
    }
    sum = 1.0/sum;
    for (size_t k = 0; k < m_kk; k++) {
        mf[k] *= sum;
    }
    setMoleFractions(mf.data());

    // After we formally set the mole fractions, we calculate the molalities
    // again and store it in this object.
    calcMolalities();
}

void MolalityVPSSTP::setMolalitiesByName(const std::string& x)
{
    compositionMap xx = parseCompString(x, speciesNames());
    setMolalitiesByName(xx);
}

// - Activities, Standard States, Activity Concentrations -----------

int MolalityVPSSTP::activityConvention() const
{
    return cAC_CONVENTION_MOLALITY;
}

void MolalityVPSSTP::getActivityConcentrations(doublereal* c) const
{
    throw NotImplementedError("MolalityVPSSTP::getActivityConcentrations");
}

doublereal MolalityVPSSTP::standardConcentration(size_t k) const
{
    throw NotImplementedError("MolalityVPSSTP::standardConcentration");
}

void MolalityVPSSTP::getActivities(doublereal* ac) const
{
    throw NotImplementedError("MolalityVPSSTP::getActivities");
}

void MolalityVPSSTP::getActivityCoefficients(doublereal* ac) const
{
    getMolalityActivityCoefficients(ac);
    double xmolSolvent = std::max(moleFraction(0), m_xmolSolventMIN);
    for (size_t k = 1; k < m_kk; k++) {
        ac[k] /= xmolSolvent;
    }
}

void MolalityVPSSTP::getMolalityActivityCoefficients(doublereal* acMolality) const
{
    getUnscaledMolalityActivityCoefficients(acMolality);
    applyphScale(acMolality);
}

doublereal MolalityVPSSTP::osmoticCoefficient() const
{
    // First, we calculate the activities all over again
    vector_fp act(m_kk);
    getActivities(act.data());

    // Then, we calculate the sum of the solvent molalities
    double sum = 0;
    for (size_t k = 1; k < m_kk; k++) {
        sum += std::max(m_molalities[k], 0.0);
    }
    double oc = 1.0;
    if (sum > 1.0E-200) {
        oc = - log(act[0]) / (m_Mnaught * sum);
    }
    return oc;
}

void MolalityVPSSTP::setStateFromXML(const XML_Node& state)
{
    VPStandardStateTP::setStateFromXML(state);
    string comp = getChildValue(state,"soluteMolalities");
    if (comp != "") {
        setMolalitiesByName(comp);
    }
    if (state.hasChild("pressure")) {
        double p = getFloat(state, "pressure", "pressure");
        setPressure(p);
    }
}

void MolalityVPSSTP::setState_TPM(doublereal t, doublereal p,
                                  const doublereal* const molalities)
{
    setMolalities(molalities);
    setState_TP(t, p);
}

void MolalityVPSSTP::setState_TPM(doublereal t, doublereal p, const compositionMap& m)
{
    setMolalitiesByName(m);
    setState_TP(t, p);
}

void MolalityVPSSTP::setState_TPM(doublereal t, doublereal p, const std::string& m)
{
    setMolalitiesByName(m);
    setState_TP(t, p);
}

void MolalityVPSSTP::setState(const AnyMap& state) {
    AnyValue molalities;
    if (state.hasKey("molalities")) {
        molalities = state["molalities"];
    } else if (state.hasKey("M")) {
        molalities = state["M"];
    }

    if (molalities.is<string>()) {
        setMolalitiesByName(molalities.asString());
    } else if (molalities.is<AnyMap>()) {
        setMolalitiesByName(molalities.asMap<double>());
    }

    VPStandardStateTP::setState(state);
}

void MolalityVPSSTP::initThermo()
{
    VPStandardStateTP::initThermo();

    // Find the Cl- species
    m_indexCLM = findCLMIndex();
}

void MolalityVPSSTP::getUnscaledMolalityActivityCoefficients(doublereal* acMolality) const
{
    throw NotImplementedError("MolalityVPSSTP::getUnscaledMolalityActivityCoefficients");
}

void MolalityVPSSTP::applyphScale(doublereal* acMolality) const
{
    throw NotImplementedError("MolalityVPSSTP::applyphScale");
}

size_t MolalityVPSSTP::findCLMIndex() const
{
    size_t indexCLM = npos;
    size_t eCl = npos;
    size_t eE = npos;
    size_t ne = nElements();
    for (size_t e = 0; e < ne; e++) {
        string sn = elementName(e);
        if (sn == "Cl" || sn == "CL") {
            eCl = e;
            break;
        }
    }
    // We have failed if we can't find the Cl element index
    if (eCl == npos) {
        return npos;
    }
    for (size_t e = 0; e < ne; e++) {
        string sn = elementName(e);
        if (sn == "E" || sn == "e") {
            eE = e;
            break;
        }
    }
    // We have failed if we can't find the E element index
    if (eE == npos) {
        return npos;
    }
    for (size_t k = 1; k < m_kk; k++) {
        doublereal nCl = nAtoms(k, eCl);
        if (nCl != 1.0) {
            continue;
        }
        doublereal nE = nAtoms(k, eE);
        if (nE != 1.0) {
            continue;
        }
        for (size_t e = 0; e < ne; e++) {
            if (e != eE && e != eCl) {
                doublereal nA = nAtoms(k, e);
                if (nA != 0.0) {
                    continue;
                }
            }
        }
        string sn = speciesName(k);
        if (sn != "Cl-" && sn != "CL-") {
            continue;
        }

        indexCLM = k;
        break;
    }
    return indexCLM;
}

bool MolalityVPSSTP::addSpecies(shared_ptr<Species> spec)
{
    bool added = VPStandardStateTP::addSpecies(spec);
    if (added) {
        if (m_kk == 1) {
            // The solvent defaults to species 0
            m_weightSolvent = molecularWeight(0);
            m_Mnaught = m_weightSolvent / 1000.;
        }
        m_molalities.push_back(0.0);
    }
    return added;
}

std::string MolalityVPSSTP::report(bool show_thermo, doublereal threshold) const
{
    fmt::memory_buffer b;
    try {
        if (name() != "") {
            format_to(b, "\n  {}:\n", name());
        }
        format_to(b, "\n");
        format_to(b, "       temperature    {:12.6g}  K\n", temperature());
        format_to(b, "          pressure    {:12.6g}  Pa\n", pressure());
        format_to(b, "           density    {:12.6g}  kg/m^3\n", density());
        format_to(b, "  mean mol. weight    {:12.6g}  amu\n", meanMolecularWeight());

        doublereal phi = electricPotential();
        format_to(b, "         potential    {:12.6g}  V\n", phi);

        vector_fp x(m_kk);
        vector_fp molal(m_kk);
        vector_fp mu(m_kk);
        vector_fp muss(m_kk);
        vector_fp acMolal(m_kk);
        vector_fp actMolal(m_kk);
        getMoleFractions(&x[0]);
        getMolalities(&molal[0]);
        getChemPotentials(&mu[0]);
        getStandardChemPotentials(&muss[0]);
        getMolalityActivityCoefficients(&acMolal[0]);
        getActivities(&actMolal[0]);

        size_t iHp = speciesIndex("H+");
        if (iHp != npos) {
            double pH = -log(actMolal[iHp]) / log(10.0);
            format_to(b, "                pH    {:12.4g}\n", pH);
        }

        if (show_thermo) {
            format_to(b, "\n");
            format_to(b, "                          1 kg            1 kmol\n");
            format_to(b, "                       -----------      ------------\n");
            format_to(b, "          enthalpy    {:12.6g}     {:12.4g}     J\n",
                    enthalpy_mass(), enthalpy_mole());
            format_to(b, "   internal energy    {:12.6g}     {:12.4g}     J\n",
                    intEnergy_mass(), intEnergy_mole());
            format_to(b, "           entropy    {:12.6g}     {:12.4g}     J/K\n",
                    entropy_mass(), entropy_mole());
            format_to(b, "    Gibbs function    {:12.6g}     {:12.4g}     J\n",
                    gibbs_mass(), gibbs_mole());
            format_to(b, " heat capacity c_p    {:12.6g}     {:12.4g}     J/K\n",
                    cp_mass(), cp_mole());
            try {
                format_to(b, " heat capacity c_v    {:12.6g}     {:12.4g}     J/K\n",
                        cv_mass(), cv_mole());
            } catch (NotImplementedError&) {
                format_to(b, " heat capacity c_v    <not implemented>\n");
            }
        }

        format_to(b, "\n");
        int nMinor = 0;
        doublereal xMinor = 0.0;
        if (show_thermo) {
            format_to(b, "                           X        "
                    "   Molalities         Chem.Pot.    ChemPotSS    ActCoeffMolal\n");
            format_to(b, "                                    "
                    "                      (J/kmol)      (J/kmol)\n");
            format_to(b, "                     -------------  "
                    "  ------------     ------------  ------------    ------------\n");
            for (size_t k = 0; k < m_kk; k++) {
                if (x[k] > threshold) {
                    if (x[k] > SmallNumber) {
                        format_to(b, "{:>18s}  {:12.6g}     {:12.6g}     {:12.6g}   {:12.6g}   {:12.6g}\n",
                                speciesName(k), x[k], molal[k], mu[k], muss[k], acMolal[k]);
                    } else {
                        format_to(b, "{:>18s}  {:12.6g}     {:12.6g}          N/A      {:12.6g}   {:12.6g}\n",
                                speciesName(k), x[k], molal[k], muss[k], acMolal[k]);
                    }
                } else {
                    nMinor++;
                    xMinor += x[k];
                }
            }
        } else {
            format_to(b, "                           X"
                    "Molalities\n");
            format_to(b, "                     -------------"
                    "     ------------\n");
            for (size_t k = 0; k < m_kk; k++) {
                if (x[k] > threshold) {
                    format_to(b, "{:>18s}   {:12.6g}     {:12.6g}\n",
                            speciesName(k), x[k], molal[k]);
                } else {
                    nMinor++;
                    xMinor += x[k];
                }
            }
        }
        if (nMinor) {
            format_to(b, "     [{:+5d} minor] {:12.6g}\n", nMinor, xMinor);
        }
    } catch (CanteraError& err) {
        return to_string(b) + err.what();
    }
    return to_string(b);
}

void MolalityVPSSTP::getCsvReportData(std::vector<std::string>& names,
                                      std::vector<vector_fp>& data) const
{
    names.clear();
    data.assign(10, vector_fp(nSpecies()));

    names.push_back("X");
    getMoleFractions(&data[0][0]);

    names.push_back("Molal");
    getMolalities(&data[1][0]);

    names.push_back("Chem. Pot. (J/kmol)");
    getChemPotentials(&data[2][0]);

    names.push_back("Chem. Pot. SS (J/kmol)");
    getStandardChemPotentials(&data[3][0]);

    names.push_back("Molal Act. Coeff.");
    getMolalityActivityCoefficients(&data[4][0]);

    names.push_back("Molal Activity");
    getActivities(&data[5][0]);

    names.push_back("Part. Mol Enthalpy (J/kmol)");
    getPartialMolarEnthalpies(&data[5][0]);

    names.push_back("Part. Mol. Entropy (J/K/kmol)");
    getPartialMolarEntropies(&data[6][0]);

    names.push_back("Part. Mol. Energy (J/kmol)");
    getPartialMolarIntEnergies(&data[7][0]);

    names.push_back("Part. Mol. Cp (J/K/kmol");
    getPartialMolarCp(&data[8][0]);

    names.push_back("Part. Mol. Cv (J/K/kmol)");
    getPartialMolarVolumes(&data[9][0]);
}

}
