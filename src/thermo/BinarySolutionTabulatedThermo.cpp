/**
 *  @file BinarySolutionTabulatedThermo.cpp Implementation file for an binary solution model
 *      with tabulated standard state thermodynamic data (see \ref thermoprops and
 *      class \link Cantera::BinarySolutionTabulatedThermo BinarySolutionTabulatedThermo\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/BinarySolutionTabulatedThermo.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"

namespace Cantera
{

BinarySolutionTabulatedThermo::BinarySolutionTabulatedThermo(int formGC) :
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError(" BinarySolutionTabulatedThermo Constructor",
                           " Illegal value of formGC");
    }
}

BinarySolutionTabulatedThermo::BinarySolutionTabulatedThermo(const std::string& inputFile,
        const std::string& id_, int formGC) :
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError(" BinarySolutionTabulatedThermo Constructor",
                           " Illegal value of formGC");
    }
    initThermoFile(inputFile, id_);
}

BinarySolutionTabulatedThermo::BinarySolutionTabulatedThermo(XML_Node& root, const std::string& id_,
        int formGC) :
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError(" BinarySolutionTabulatedThermo Constructor",
                           " Illegal value of formGC");
    }
    importPhase(root, this);
}

void BinarySolutionTabulatedThermo::compositionChanged()
{
    IdealSolidSolnPhase::compositionChanged();
    _updateThermo();
}

void BinarySolutionTabulatedThermo::_updateThermo()
{
    double tnow = temperature();
    double xnow = moleFraction(m_kk_tab);
    double c[4];
    double *d;
    double dS_corr = 0.0;
    double tlow = 0.0, thigh = 0.0;
    int type = 0;
    if (m_tlast != tnow || m_xlast != xnow) {
        c[0] = tnow;
        d = interpolate(xnow);
        c[1] = d[0] * 1e3; // 1e3 for conversion J/mol -> J/kmol
        if (xnow == 0)
        {
            dS_corr = -BigNumber;
        } else if (xnow == 1)
        {
            dS_corr = BigNumber;
        } else
        {
            dS_corr = GasConstant*std::log(xnow/(1.0-xnow)) + GasConstant/Faraday*std::log(this->standardConcentration(1-m_kk_tab)/this->standardConcentration(m_kk_tab));
        }
        c[2] = d[1] * 1e3 + dS_corr; // 1e3 for conversion J/K/mol -> J/K/kmol
        c[3] = 0.0;
        type = m_spthermo.reportType(m_kk_tab);
        tlow = m_spthermo.minTemp(m_kk_tab);
        thigh = m_spthermo.maxTemp(m_kk_tab);
        shared_ptr<SpeciesThermoInterpType> stit(
                newSpeciesThermoInterpType(type, tlow, thigh, OneAtm, c));
        m_spthermo.modifySpecies(m_kk_tab, stit);
        // Update the thermodynamic functions of the reference state.
        m_spthermo.update(tnow, m_cp0_R.data(), m_h0_RT.data(), m_s0_R.data());
        doublereal rrt = 1.0 / RT();
        for (size_t k = 0; k < m_kk; k++) {
            double deltaE = rrt * m_pe[k];
            m_h0_RT[k] += deltaE;
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_xlast = xnow;
        m_tlast = tnow;
    }
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
            m_xlast = moleFraction(m_kk_tab);
        }
        if (thermoNode.hasChild("tabulatedThermo")) {
            XML_Node& dataNode = thermoNode.child("tabulatedThermo");
            getFloatArray(dataNode, x, false, "", "moleFraction");
            getFloatArray(dataNode, h, false, "", "enthalpy");
            getFloatArray(dataNode, s, false, "", "entropy");

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
        std::string formString = scNode.attrib("model");
        if (caseInsensitiveEquals(formString, "unity")) {
            m_formGC = 0;
        } else if (caseInsensitiveEquals(formString, "molar_volume")) {
            m_formGC = 1;
        } else if (caseInsensitiveEquals(formString, "solvent_volume")) {
            m_formGC = 2;
        } else {
            throw CanteraError("BinarySolutionTabulatedThermo::initThermoXML",
                    "Unknown standardConc model: " + formString);
        }
    } else {
        throw CanteraError("BinarySolutionTabulatedThermo::initThermoXML",
                "Unspecified standardConc model");
    }

    // Call the base initThermo, which handles setting the initial state
    ThermoPhase::initThermoXML(phaseNode, id_);
}

double* BinarySolutionTabulatedThermo::interpolate(double x) const
{
    static double c[2];
    // Check if x is out of bound
    if (x > m_molefrac_tab.back()) {
        c[0] = m_enthalpy_tab.back();
        c[1] = m_entropy_tab.back();
        return c;
    }
    if (x < m_molefrac_tab[0]) {
        c[0] = m_enthalpy_tab[0];
        c[1] = m_entropy_tab[0];
        return c;
    }
    size_t i = std::distance(m_molefrac_tab.begin(), std::lower_bound(m_molefrac_tab.begin(), m_molefrac_tab.end(), x));
    c[0] = m_enthalpy_tab[i-1] + (m_enthalpy_tab[i] - m_enthalpy_tab[i-1]) * (x - m_molefrac_tab[i-1])/(m_molefrac_tab[i]- m_molefrac_tab[i-1]);
    c[1] = m_entropy_tab[i-1] + (m_entropy_tab[i] - m_entropy_tab[i-1]) * (x - m_molefrac_tab[i-1])/(m_molefrac_tab[i]- m_molefrac_tab[i-1]);
    return c;
}

}
