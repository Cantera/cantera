//! @file Radiation1D.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.


#include "cantera/oneD/Radiation1D.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/global.h"


namespace Cantera
{

TabularPlanckMean::TabularPlanckMean(ThermoPhase* thermo)
    : m_thermo(thermo)
{
    parseRadiationData();
}

void TabularPlanckMean::parseRadiationData()
{
    AnyMap radiationPropertiesDB;

    try {
        radiationPropertiesDB = AnyMap::fromYamlFile("radiation-parameters.yaml");
    } catch (CanteraError& err) {
        warn_user("TabularPlanckMean::parseRadiationData",
            "Failed to load 'radiation-parameters.yaml':\n{}"
            "\nFalling back to default polynomial data for CO2, H2O.", err.what());
    }

    if(!radiationPropertiesDB.empty() && radiationPropertiesDB.hasKey("PMAC")) {
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
    }

    // Polynomial coefficients for CO2 and H2O (backwards compatibility)
    // Check if "CO2" is already in the map, if not, add the polynomial fit data
    if (!m_PMAC.hasKey("CO2")) {
        const std::vector<double> c_CO2 = {18.741, -121.310, 273.500, -194.050, 56.310,
                                            -5.8169};
        m_PMAC["CO2"]["fit-type"] = "polynomial";
        m_PMAC["CO2"]["coefficients"] = c_CO2;
    }
    if (m_absorptionSpecies.find("CO2") == m_absorptionSpecies.end()) {
         size_t k = m_thermo->speciesIndex("CO2");
         if (k != npos) m_absorptionSpecies.emplace("CO2", k);
     }

    // Check if "H2O" is already in the map, if not, add the polynomial fit data
    if (!m_PMAC.hasKey("H2O")) {
        const std::vector<double> c_H2O = {-0.23093, -1.12390, 9.41530, -2.99880,
                                            0.51382, -1.86840e-5};
        m_PMAC["H2O"]["fit-type"] = "polynomial";
        m_PMAC["H2O"]["coefficients"] = c_H2O;
    }
    if (m_absorptionSpecies.find("H2O") == m_absorptionSpecies.end()) {
         size_t k = m_thermo->speciesIndex("H2O");
         if (k != npos) m_absorptionSpecies.emplace("H2O", k);
     }
}

double TabularPlanckMean::calculatePolynomial(const std::vector<double>& coefficients,
                                              double temperature)
{
    double result = 0.0;
    for (size_t n = 0; n < coefficients.size(); ++n) {
        result += coefficients[n] * std::pow(1000 / temperature, static_cast<double>(n));
    }
    return result / (1.0 * OneAtm);
}

double TabularPlanckMean::interpolateTable(const std::vector<double>& temperatures,
                                           const std::vector<double>& data,
                                           double temperature)
{
    // Handle edge cases first
    if (temperature <= temperatures.front()) {
        // alpha = 1.0 / data[0]
        return 1.0 / data.front();
    } else if (temperature >= temperatures.back()) {
        // alpha = 1.0 / data[last]
        return 1.0 / data.back();
    }

    // Find the interval [t1, t2] where t1 <= T < t2
    // so that temperatures[i-1] <= T < temperatures[i]
    size_t idx = 1;
    for (; idx < temperatures.size(); ++idx) {
        if (temperature < temperatures[idx]) {
            break;
        }
    }

    // Perform linear interpolation
    double t1 = temperatures[idx - 1];
    double t2 = temperatures[idx];
    double v1 = data[idx - 1];
    double v2 = data[idx];

    // ln(alpha) = ln(1/v1) + ( ln(1/v2) - ln(1/v1) ) * (T - t1)/(t2 - t1)
    double frac = (temperature - t1) / (t2 - t1);
    double lnAlpha = log(1.0 / v1) + (log(1.0 / v2) - log(1.0 / v1)) * frac;
    return exp(lnAlpha) / (1.0 * OneAtm);
}

void TabularPlanckMean::getBandProperties(std::vector<double>& kabs,
                                    std::vector<double>& awts,
                                    const RadComposition& comp)
{
    double k_P = 0;

    // Loop over absorbing species
    for (const auto& [sp_name, sp_idx] : m_absorptionSpecies) {
        const auto& fit_type = m_PMAC[sp_name]["fit-type"].asString();

        // Get the species mole fraction from the Composition
        // If the species doesn't exist in comp.X, error out
        double x_sp = 0;
        if (comp.X.find(sp_name) == comp.X.end()) {
            throw CanteraError("TabularPlanckMean::getBandProperties",
                "Species '{}' not found in composition data", sp_name);
        } else {
            x_sp = comp.X.at(sp_name);
        }

        if (fit_type == "table") {
            double kVal = interpolateTable(
                m_PMAC[sp_name]["temperatures"].asVector<double>(),
                m_PMAC[sp_name]["coefficients"].asVector<double>(),
                comp.T);
            k_P += comp.P * x_sp * kVal;
        } else if (fit_type == "polynomial") {
            double kVal = calculatePolynomial(
                m_PMAC[sp_name]["coefficients"].asVector<double>(),
                comp.T);
            k_P += comp.P * x_sp * kVal;
        }
    }

    // Store the single-band result
    kabs.resize(1);
    awts.resize(1);
    kabs[0] = k_P;
    awts[0] = 1.0; // single “band” weighting
}


Radiation1D::Radiation1D(ThermoPhase* thermo, double pressure, size_t points,
                         std::function<double(const double*, size_t)> temperatureFunction,
                         std::function<double(const double*, size_t, size_t)> moleFractionFunction,
                         std::unique_ptr<RadiationPropertyCalculator> props,
                         std::unique_ptr<RadiationSolver> solver)
    : m_thermo(thermo), m_press(pressure), m_points(points),
    m_T(temperatureFunction), m_X(moleFractionFunction),
    m_props(std::move(props)), m_solver(std::move(solver))
{
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


void Radiation1D::computeRadiation(double* x, size_t jmin, size_t jmax,
                                   std::vector<double>& qdotRadiation) {
    double boundary_Rad_left = m_epsilon_left * StefanBoltz * std::pow(m_T(x, 0), 4);
    double boundary_Rad_right = m_epsilon_right * StefanBoltz * std::pow(m_T(x, m_points - 1), 4);

    for (size_t j = jmin; j < jmax; j++) {
        RadComposition comp;
        comp.T = m_T(x, j);
        comp.P = m_press;
        // Build local mole-fraction map at j
        comp.X.clear();
        const auto& names = m_thermo->speciesNames();
        for (const auto& nm : names) {
            size_t k = m_thermo->speciesIndex(nm);
            comp.X[nm] = m_X(x, k, j);
        }

        // Get the band absorption coefficients and weighting factors
        std::vector<double> kabs, awts;
        m_props->getBandProperties(kabs, awts, comp);

        // Solve for radiative heat loss
        qdotRadiation[j] = m_solver->computeHeatLoss(kabs, awts, comp.T,
                                                     boundary_Rad_left,
                                                     boundary_Rad_right);
    }
}




} // namespace Cantera