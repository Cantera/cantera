//! @file Radiation1D.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.


#include "cantera/oneD/Radiation1D.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/global.h"


namespace Cantera
{


Radiation1D::Radiation1D(ThermoPhase* thermo, double pressure, size_t points,
                         std::function<double(const double*, size_t)> temperatureFunction,
                         std::function<double(const double*, size_t, size_t)> moleFractionFunction)
    : m_thermo(thermo), m_press(pressure), m_points(points),
    m_T(temperatureFunction), m_X(moleFractionFunction)
{
    parseRadiationData();
}

void Radiation1D::setBoundaryEmissivities(double e_left, double e_right)
{
    if (e_left < 0 || e_left > 1) {
        throw CanteraError("Radiation1D::setBoundaryEmissivities",
            "The left boundary emissivity must be between 0.0 and 1.0!");
    } else if (e_right < 0 || e_right > 1) {
        throw CanteraError("Radiation1D::setBoundaryEmissivities",
            "The right boundary emissivity must be between 0.0 and 1.0!");
    } else {
        m_epsilon_left = e_left;
        m_epsilon_right = e_right;
    }
}

void Radiation1D::parseRadiationData() {
    AnyMap radiationPropertiesDB;
    // Search 'crit-properties.yaml' to find Tc and Pc. Load data if needed.
    if (radiationPropertiesDB.empty()) {
        radiationPropertiesDB = AnyMap::fromYamlFile("radiation-properties.yaml");
    }

    if( radiationPropertiesDB.hasKey("PMAC")) {
        auto& data = radiationPropertiesDB["PMAC"].as<AnyMap>();

        // Needs to loop over only the species that are in the input yaml data
        for (const auto& name : m_thermo->speciesNames()) {
            if (data.hasKey("radiation")) {
                std::cout << "Radiation data found for species " << name << std::endl;
                m_absorptionSpecies.insert({name, m_thermo->speciesIndex(name)});
                if (data["radiation"].hasKey("fit-type")) {
                    m_PMAC[name]["fit-type"] = data["radiation"]["fit-type"].asString();
                } else {
                    throw InputFileError("Flow1D::Flow1D", data,
                    "No 'fit-type' entry found for species '{}'", name);
                }

                // This is the direct tabulation of the optical path length
                if (data["radiation"]["fit-type"] == "table") {
                    if (data["radiation"].hasKey("temperatures")) {
                        std::cout << "Storing temperatures for species " << name << std::endl;
                        // Each species may have a specific set of temperatures that are used
                        m_PMAC[name]["temperatures"] = data["radiation"]["temperatures"].asVector<double>();
                    } else {
                        throw InputFileError("Flow1D::Flow1D", data,
                        "No 'temperatures' entry found for species '{}'", name);
                    }

                    if (data["radiation"].hasKey("data")) {
                        std::cout << "Storing data for species " << name << std::endl;
                        // This data is the Plank mean absorption coefficient
                        m_PMAC[name]["coefficients"] = data["radiation"]["data"].asVector<double>();
                    } else {
                        throw InputFileError("Flow1D::Flow1D", data,
                        "No 'data' entry found for species '{}'", name);
                    }
                } else if (data["radiation"]["fit-type"] == "polynomial") {
                    std::cout << "Polynomial fit found for species " << name << std::endl;
                    if (data["radiation"].hasKey("data")) {
                        std::cout << "Storing data for species " << name << std::endl;
                        m_PMAC[name]["coefficients"] = data["radiation"]["data"].asVector<double>();
                    } else {
                        throw InputFileError("Flow1D::Flow1D", data,
                        "No 'data' entry found for species '{}'", name);
                    }
                } else {
                    throw InputFileError("Flow1D::Flow1D", data,
                    "Invalid 'fit-type' entry found for species '{}'", name);
                }
            }
        }

        // Polynomial coefficients for CO2 and H2O (backwards compatibility)
        // Check if "CO2" is already in the map, if not, add the polynomial fit data
        if (!m_PMAC.hasKey("CO2")) {
            const std::vector<double> c_CO2 = {18.741, -121.310, 273.500, -194.050, 56.310,
                                                -5.8169};
            m_PMAC["CO2"]["fit-type"] = "polynomial";
            m_PMAC["CO2"]["coefficients"] = c_CO2;
        }

        // Check if "H2O" is already in the map, if not, add the polynomial fit data
        if (!m_PMAC.hasKey("H2O")) {
            const std::vector<double> c_H2O = {-0.23093, -1.12390, 9.41530, -2.99880,
                                                0.51382, -1.86840e-5};
            m_PMAC["H2O"]["fit-type"] = "polynomial";
            m_PMAC["H2O"]["coefficients"] = c_H2O;
        }
    }
}

void Radiation1D::computeRadiation(double* x, size_t jmin, size_t jmax, std::vector<double>& qdotRadiation) {
    const double k_P_ref = 1.0 * OneAtm;
    const double StefanBoltz = 5.67e-8;

    double boundary_Rad_left = m_epsilon_left * StefanBoltz * std::pow(m_T(x, 0), 4);
    double boundary_Rad_right = m_epsilon_right * StefanBoltz * std::pow(m_T(x, m_points - 1), 4);

    for (size_t j = jmin; j < jmax; j++) {
        double k_P = 0;
        for (const auto& [sp_name, sp_idx] : m_absorptionSpecies) {
            const auto& fit_type = m_PMAC[sp_name]["fit-type"].asString();
            if (fit_type == "table") {
                k_P += m_press * m_X(x, sp_idx, j) *
                       interpolateTable(m_PMAC[sp_name]["temperatures"].asVector<double>(), m_PMAC[sp_name]["coefficients"].asVector<double>(), m_T(x, j));
            } else if (fit_type == "polynomial") {
                k_P += m_press * m_X(x, sp_idx, j) *
                       calculatePolynomial(m_PMAC[sp_name]["coefficients"].asVector<double>(), m_T(x, j));
            }
        }
        qdotRadiation[j] = 2 * k_P * (2 * StefanBoltz * std::pow(m_T(x, j), 4) - boundary_Rad_left - boundary_Rad_right);
    }
}

double Radiation1D::calculatePolynomial(const std::vector<double>& coefficients, double temperature) {
    double result = 0.0;
    for (size_t n = 0; n < coefficients.size(); ++n) {
        result += coefficients[n] * std::pow(1000 / temperature, static_cast<double>(n));
    }
    return result / (1.0 * OneAtm);
}

double Radiation1D::interpolateTable(const std::vector<double>& temperatures, const std::vector<double>& data, double temperature) {
    size_t index = 0;
    for (size_t i = 1; i < temperatures.size(); ++i) {
        if (temperature < temperatures[i]) {
            index = i - 1;
            break;
        }
    }
    double t1 = temperatures[index], t2 = temperatures[index + 1];
    double d1 = data[index], d2 = data[index + 1];
    return d1 + (d2 - d1) * (temperature - t1) / (t2 - t1);
}


} // namespace Cantera