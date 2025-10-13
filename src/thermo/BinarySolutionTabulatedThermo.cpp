/**
 *  @file BinarySolutionTabulatedThermo.cpp Implementation file for an binary
 *      solution model with tabulated standard state thermodynamic data (see
 *       @ref thermoprops and class
 *       @link Cantera::BinarySolutionTabulatedThermo BinarySolutionTabulatedThermo@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/BinarySolutionTabulatedThermo.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"

namespace Cantera
{

BinarySolutionTabulatedThermo::BinarySolutionTabulatedThermo(const string& inputFile,
                                                             const string& id_)
{
    initThermoFile(inputFile, id_);
}

void BinarySolutionTabulatedThermo::compositionChanged()
{
    IdealSolidSolnPhase::compositionChanged();
    _updateThermo();
}

void BinarySolutionTabulatedThermo::_updateThermo() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    bool x_changed = !cached.validate(stateMFNumber());

    if (x_changed) {
        double x_tab = moleFraction(m_kk_tab);
        double x_other = moleFraction(1 - m_kk_tab);
        m_h0_tab = interpolate(x_tab, m_enthalpy_tab);
        m_s0_tab = interpolate(x_tab, m_entropy_tab);
        if (x_tab == 0) {
            m_s0_tab = -BigNumber;
        } else if (x_other == 0) {
            m_s0_tab = BigNumber;
        } else {
            m_s0_tab += GasConstant*std::log(x_tab / x_other) +
                    GasConstant/Faraday*std::log(standardConcentration(1-m_kk_tab)
                    /standardConcentration(m_kk_tab));
        }
    }

    double tnow = temperature();
    if (x_changed || m_tlast != tnow) {
        // Update the thermodynamic functions of the reference state.
        m_spthermo.update(tnow, m_cp0_R.data(), m_h0_RT.data(), m_s0_R.data());
        double rrt = 1.0 / RT();
        m_h0_RT[m_kk_tab] += m_h0_tab * rrt;
        m_s0_R[m_kk_tab] += m_s0_tab / GasConstant;
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_tlast = tnow;
    }
}

bool BinarySolutionTabulatedThermo::addSpecies(shared_ptr<Species> spec)
{
    if  (m_kk == 2) {
        throw CanteraError("BinarySolutionTabulatedThermo::addSpecies",
                           "No. of species should be equal to 2");
    }
    bool added = IdealSolidSolnPhase::addSpecies(spec);
    return added;
}

void BinarySolutionTabulatedThermo::initThermo()
{
    if (m_input.hasKey("tabulated-thermo")) {
        m_kk_tab = speciesIndex(m_input["tabulated-species"].asString(), true);
        if (nSpecies() != 2) {
            throw InputFileError("BinarySolutionTabulatedThermo::initThermo",
                m_input["species"],
                "No. of species should be equal to 2 in phase '{}'!",name());
        }
        if (m_kk_tab == npos) {
            throw InputFileError("BinarySolutionTabulatedThermo::initThermo",
                m_input["tabulated-species"],
                "Species '{}' is not in phase '{}'",
                m_input["tabulated-species"].asString(), name());
        }
        const AnyMap& table = m_input["tabulated-thermo"].as<AnyMap>();
        vector<double> x = table["mole-fractions"].asVector<double>();
        size_t N = x.size();
        vector<double> h = table.convertVector("enthalpy", "J/kmol", N);
        vector<double> s = table.convertVector("entropy", "J/kmol/K", N);
        vector<double> vmol(N);
        // Check for molar-volume key in tabulatedThermo table,
        // otherwise calculate molar volume from pure species molar volumes
        if (table.hasKey("molar-volume")) {
            vmol = table.convertVector("molar-volume", "m^3/kmol", N);
        } else {
            for(size_t i = 0; i < N; i++) {
                vmol[i] = x[i] * m_speciesMolarVolume[m_kk_tab] + (1-x[i])
                        * m_speciesMolarVolume[1-m_kk_tab];
            }
        }

        // Sort the x, h, s, vmol data in the order of increasing x
        vector<pair<double,double>> x_h(N), x_s(N), x_vmol(N);
        for(size_t i = 0; i < N; i++) {
            x_h[i] = {x[i], h[i]};
            x_s[i] = {x[i], s[i]};
            x_vmol[i] = {x[i], vmol[i]};
        }
        std::sort(x_h.begin(), x_h.end());
        std::sort(x_s.begin(), x_s.end());
        std::sort(x_vmol.begin(), x_vmol.end());

        // Store the sorted values in different arrays
        m_molefrac_tab.resize(N);
        m_enthalpy_tab.resize(N);
        m_entropy_tab.resize(N);
        m_molar_volume_tab.resize(N);
        m_derived_molar_volume_tab.resize(N);

        for (size_t i = 0; i < N; i++) {
            m_molefrac_tab[i] = x_h[i].first;
            m_enthalpy_tab[i] = x_h[i].second;
            m_entropy_tab[i] = x_s[i].second;
            m_molar_volume_tab[i] = x_vmol[i].second;
        }

        diff(m_molar_volume_tab, m_derived_molar_volume_tab);
    }
    IdealSolidSolnPhase::initThermo();
}

bool BinarySolutionTabulatedThermo::ready() const
{
    return !m_molefrac_tab.empty();
}

void BinarySolutionTabulatedThermo::getParameters(AnyMap& phaseNode) const
{
    IdealSolidSolnPhase::getParameters(phaseNode);
    phaseNode["tabulated-species"] = speciesName(m_kk_tab);
    AnyMap tabThermo;
    tabThermo["mole-fractions"] = m_molefrac_tab;
    tabThermo["enthalpy"].setQuantity(m_enthalpy_tab, "J/kmol");
    tabThermo["entropy"].setQuantity(m_entropy_tab, "J/kmol/K");
    tabThermo["molar-volume"].setQuantity(m_molar_volume_tab, "m^3/kmol");
    phaseNode["tabulated-thermo"] = std::move(tabThermo);
}

double BinarySolutionTabulatedThermo::interpolate(const double x,
                                                  const vector<double>& inputData) const
{
    double c;
    // Check if x is out of bound
    if (x > m_molefrac_tab.back()) {
        c = inputData.back();
        return c;
    }
    if (x < m_molefrac_tab.front()) {
        c = inputData.front();
        return c;
    }
    size_t i = std::distance(m_molefrac_tab.begin(),
            std::lower_bound(m_molefrac_tab.begin(), m_molefrac_tab.end(), x));
    c = inputData[i-1] + (inputData[i] - inputData[i-1])
            * (x - m_molefrac_tab[i-1]) / (m_molefrac_tab[i] - m_molefrac_tab[i-1]);
    return c;
}

void BinarySolutionTabulatedThermo::diff(const vector<double>& inputData,
                                         vector<double>& derivedData) const
{
    if (inputData.size() > 1) {
        derivedData[0] = (inputData[1] - inputData[0]) /
                (m_molefrac_tab[1] - m_molefrac_tab[0]);
        derivedData.back() = (inputData.back() - inputData[inputData.size()-2]) /
                (m_molefrac_tab.back() - m_molefrac_tab[m_molefrac_tab.size()-2]);

        if (inputData.size() > 2) {
            for (size_t i = 1; i < inputData.size()-1; i++) {
                derivedData[i] = (inputData[i+1] - inputData[i-1]) /
                        (m_molefrac_tab[i+1] - m_molefrac_tab[i-1]);
            }
        }
    } else {
        derivedData.front() = 0;
    }
}

void BinarySolutionTabulatedThermo::getPartialMolarVolumes(double* vbar) const
{
    std::copy(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), vbar);
}

void BinarySolutionTabulatedThermo::calcDensity()
{
    double Xtab = moleFraction(m_kk_tab);
    double Vm = interpolate(Xtab, m_molar_volume_tab);
    double dVdX_tab = interpolate(Xtab, m_derived_molar_volume_tab);
    m_speciesMolarVolume[m_kk_tab] = Vm + (1 - Xtab) * dVdX_tab;
    m_speciesMolarVolume[1-m_kk_tab] = Vm - Xtab * dVdX_tab;

    double dens = meanMolecularWeight() / Vm;

    // Set the density in the parent State object directly, by calling the
    // Phase::assignDensity() function.
    Phase::assignDensity(dens);
}
}
