/**
 *  @file BinarySolutionTabulatedThermo.cpp Implementation file for an binary
 *      solution model with tabulated standard state thermodynamic data (see
 *       \ref thermoprops and class
 *       \link Cantera::BinarySolutionTabulatedThermo BinarySolutionTabulatedThermo\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/BinarySolutionTabulatedThermo.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"

namespace Cantera
{

BinarySolutionTabulatedThermo::BinarySolutionTabulatedThermo()
    : m_kk_tab(npos)
    , m_xlast(-1)
{
}

BinarySolutionTabulatedThermo::BinarySolutionTabulatedThermo(const std::string& inputFile,
                                                             const std::string& id_)
    : m_kk_tab(npos)
    , m_xlast(-1)

{
    initThermoFile(inputFile, id_);
}

BinarySolutionTabulatedThermo::BinarySolutionTabulatedThermo(XML_Node& root,
                                                             const std::string& id_)
    : m_kk_tab(npos)
    , m_xlast(-1)
{
    importPhase(root, this);
}

void BinarySolutionTabulatedThermo::compositionChanged()
{
    IdealSolidSolnPhase::compositionChanged();
    _updateThermo();
}

void BinarySolutionTabulatedThermo::_updateThermo() const
{
    double xnow = moleFraction(m_kk_tab);
    bool x_changed = (m_xlast != xnow);

    if (x_changed) {
        std::tie(m_h0_tab, m_s0_tab) = interpolate(xnow);
        if (xnow == 0) {
            m_s0_tab = -BigNumber;
        } else if (xnow == 1) {
            m_s0_tab = BigNumber;
        } else {
            m_s0_tab += GasConstant*std::log(xnow/(1.0-xnow)) +
              GasConstant/Faraday*std::log(standardConcentration(1-m_kk_tab)
              /standardConcentration(m_kk_tab));
        }
        m_xlast = xnow;
    }

    double tnow = temperature();
    if (x_changed || m_tlast != tnow) {
        // Update the thermodynamic functions of the reference state.
        m_spthermo.update(tnow, m_cp0_R.data(), m_h0_RT.data(), m_s0_R.data());
        double rrt = 1.0 / RT();
        m_h0_RT[m_kk_tab] += m_h0_tab * rrt;
        m_s0_R[m_kk_tab] += m_s0_tab / GasConstant;
        for (size_t k = 0; k < m_kk; k++) {
            double deltaE = rrt * m_pe[k];
            m_h0_RT[k] += deltaE;
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_tlast = tnow;
    }
}

void BinarySolutionTabulatedThermo::initThermo()
{
    if (m_input.hasKey("tabulated-thermo")) {
        m_kk_tab = speciesIndex(m_input["tabulated-species"].asString());
        if (m_kk_tab == npos) {
            throw InputFileError("BinarySolutionTabulatedThermo::initThermo",
                m_input["tabulated-species"],
                "Species '{}' is not in phase '{}'",
                m_input["tabulated-species"].asString(), name());
        }
        const AnyMap& table = m_input["tabulated-thermo"].as<AnyMap>();
        vector_fp x = table["mole-fractions"].asVector<double>();
        size_t N = x.size();
        vector_fp h = table.convertVector("enthalpy", "J/kmol", N);
        vector_fp s = table.convertVector("entropy", "J/kmol/K", N);

        // Sort the x, h, s data in the order of increasing x
        std::vector<std::pair<double,double>> x_h(N), x_s(N);
        for(size_t i = 0; i < N; i++){
            x_h[i] = {x[i], h[i]};
            x_s[i] = {x[i], s[i]};
        }
        std::sort(x_h.begin(), x_h.end());
        std::sort(x_s.begin(), x_s.end());

        // Store the sorted values in different arrays
        m_molefrac_tab.resize(N);
        m_enthalpy_tab.resize(N);
        m_entropy_tab.resize(N);
        for (size_t i = 0; i < N; i++) {
            m_molefrac_tab[i] = x_h[i].first;
            m_enthalpy_tab[i] = x_h[i].second;
            m_entropy_tab[i] = x_s[i].second;
        }
    }
    IdealSolidSolnPhase::initThermo();
}

void BinarySolutionTabulatedThermo::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    vector_fp x, h, s;
    std::vector<std::pair<double,double>> x_h_temp, x_s_temp;

    if (id_.size() > 0) {
        if (phaseNode.id() != id_) {
            throw CanteraError("BinarySolutionTabulatedThermo::initThermoXML",
                    "phasenode and Id are incompatible");
        }
    }
    if (nSpecies()!=2) {
        throw CanteraError("BinarySolutionTabulatedThermo::initThermoXML",
                            "No. of species should be equal to 2!");
    }
    if (phaseNode.hasChild("thermo")) {
        XML_Node& thermoNode = phaseNode.child("thermo");
        std::string mString = thermoNode["model"];
        if (!caseInsensitiveEquals(mString, "binarysolutiontabulatedthermo")) {
            throw CanteraError("BinarySolutionTabulatedThermo::initThermoXML",
                    "Unknown thermo model: " + mString);
        }
        if (thermoNode.hasChild("tabulatedSpecies")) {
            XML_Node& speciesNode = thermoNode.child("tabulatedSpecies");
            std::string tabulated_species_name = speciesNode["name"];
            m_kk_tab = speciesIndex(tabulated_species_name);
            if (m_kk_tab == npos) {
                throw CanteraError("BinarySolutionTabulatedThermo::initThermoXML",
                        "Species " + tabulated_species_name + " not found.");
            }
        }
        if (thermoNode.hasChild("tabulatedThermo")) {
            XML_Node& dataNode = thermoNode.child("tabulatedThermo");
            getFloatArray(dataNode, x, true, "", "moleFraction");
            getFloatArray(dataNode, h, true, "J/kmol", "enthalpy");
            getFloatArray(dataNode, s, true, "J/kmol/K", "entropy");

            // Check for data length consistency
            if ((x.size() != h.size()) || (x.size() != s.size()) || (h.size() != s.size())) {
                throw CanteraError("BinarySolutionTabulatedThermo::initThermoXML",
                        "Species tabulated thermo data has different lengths.");
            }
            // Sort the x, h, s data in the order of increasing x
            for(size_t i = 0; i < x.size(); i++){
                x_h_temp.push_back(std::make_pair(x[i],h[i]));
                x_s_temp.push_back(std::make_pair(x[i],s[i]));
            }
            std::sort(x_h_temp.begin(), x_h_temp.end());
            std::sort(x_s_temp.begin(), x_s_temp.end());

            // Store the sorted values in different arrays
            m_molefrac_tab.resize(x_h_temp.size());
            m_enthalpy_tab.resize(x_h_temp.size());
            m_entropy_tab.resize(x_h_temp.size());
            for (size_t i = 0; i < x_h_temp.size(); i++) {
                m_molefrac_tab[i] = x_h_temp[i].first;
                m_enthalpy_tab[i] = x_h_temp[i].second;
                m_entropy_tab[i] = x_s_temp[i].second;
            }
        } else {
            throw CanteraError("BinarySolutionTabulatedThermo::initThermoXML",
                    "Unspecified tabulated species or thermo");
        }
    } else {
        throw CanteraError("BinarySolutionTabulatedThermo::initThermoXML",
                "Unspecified thermo model");
    }

    /*
     * Form of the standard concentrations. Must have one of:
     *
     *     <standardConc model="unity" />
     *     <standardConc model="molar_volume" />
     *     <standardConc model="solvent_volume" />
     */
    if (phaseNode.hasChild("standardConc")) {
        XML_Node& scNode = phaseNode.child("standardConc");
        setStandardConcentrationModel(scNode.attrib("model"));
    } else {
        throw CanteraError("BinarySolutionTabulatedThermo::initThermoXML",
                "Unspecified standardConc model");
    }

    // Call the base initThermo, which handles setting the initial state
    ThermoPhase::initThermoXML(phaseNode, id_);
}

std::pair<double, double> BinarySolutionTabulatedThermo::interpolate(double x) const
{
    std::pair<double, double> c;
    // Check if x is out of bound
    if (x > m_molefrac_tab.back()) {
        c.first = m_enthalpy_tab.back();
        c.second = m_entropy_tab.back();
        return c;
    }
    if (x < m_molefrac_tab[0]) {
        c.first = m_enthalpy_tab[0];
        c.second = m_entropy_tab[0];
        return c;
    }
    size_t i = std::distance(m_molefrac_tab.begin(),
      std::lower_bound(m_molefrac_tab.begin(), m_molefrac_tab.end(), x));
    c.first = m_enthalpy_tab[i-1] + (m_enthalpy_tab[i] - m_enthalpy_tab[i-1])
      * (x - m_molefrac_tab[i-1])/(m_molefrac_tab[i]- m_molefrac_tab[i-1]);
    c.second = m_entropy_tab[i-1] + (m_entropy_tab[i] - m_entropy_tab[i-1])
      * (x - m_molefrac_tab[i-1])/(m_molefrac_tab[i]- m_molefrac_tab[i-1]);
    return c;
}

}
