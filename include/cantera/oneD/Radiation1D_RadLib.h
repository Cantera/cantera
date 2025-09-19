//
// Optional RadLib-backed radiation property calculators for Cantera 1-D domains.
// Enabled when CT_ENABLE_RADLIB is defined in Cantera's config.h.
//
// This file is intended to be compiled only when RadLib is available.
// RadLib: https://github.com/BYUignite/radlib  (MIT license)
//
// The adapter classes here implement Cantera's RadiationPropertyCalculator
// interface and forward getBandProperties() calls to RadLib model objects.
// They return band-wise absorption coefficients k_i [1/m] and weights a_i [-],
// which are then consumed by a RadiationSolver (e.g. OpticallyThinSolver).
//
// Notes:
//  - Units expected by RadLib are SI: T [K], P [Pa], fvsoot [dimensionless volume fraction].
//  - Current RadLib models use only {H2O, CO2, CO, CH4} for gas absorption,
//    plus soot via fvsoot. Mole fractions for other species are ignored.
//  - For RCSLW, RadLib requires the number of gray gases (nGG) and a reference
//    state at construction; users should choose a representative T_ref for the domain.
//
// ----------------------------------------------------------------------------

#pragma once

#include "Radiation1D.h"

#ifdef CT_ENABLE_RADLIB

#include <rad.h>              // base class
#include <rad_planck_mean.h>  // rad_planck_mean
#include <rad_wsgg.h>         // rad_wsgg
#include <rad_rcslw.h>        // rad_rcslw

namespace Cantera {

// Helper to fetch mole fraction from RadComposition::X by key; returns 0 if missing
inline double _getX(const std::map<std::string, double>& X, const char* name) {
    auto it = X.find(name);
    return it == X.end() ? 0.0 : it->second;
}

// -------------------------- Base adapter --------------------------

class RadLibAdapterBase : public RadiationPropertyCalculator {
public:
    explicit RadLibAdapterBase(ThermoPhase* thermo)
        : m_thermo(thermo) {}

    ~RadLibAdapterBase() override = default;

protected:
    ThermoPhase* m_thermo;
    double m_fvsoot {0.0}; // uniform soot volume fraction, default 0
};

// -------------------------- Planck Mean ---------------------------

class RadLibPlanckMean : public RadLibAdapterBase {
public:
    explicit RadLibPlanckMean(ThermoPhase* thermo, double fvsoot = 0.0)
        : RadLibAdapterBase(thermo)
        , m_rad(new rad_planck_mean())
    {
        m_fvsoot = fvsoot;
    }

    void getBandProperties(std::vector<double>& kabs,
                           std::vector<double>& awts,
                           const RadComposition& comp) override
    {
        // Extract required mole fractions from comp.X (others are ignored by RadLib)
        const double xH2O = _getX(comp.X, "H2O");
        const double xCO2 = _getX(comp.X, "CO2");
        const double xCO  = _getX(comp.X, "CO");
        const double xCH4 = _getX(comp.X, "CH4");

        // Ask RadLib for k and a (Planck mean -> one band)
        kabs.clear(); awts.clear();
        m_rad->get_k_a(kabs, awts, comp.T, comp.P, m_fvsoot,
                       xH2O, xCO2, xCO, xCH4);
    }

    std::vector<std::string> requiredSpecies() const override {
        return {"H2O", "CO2", "CO", "CH4"};
    }

private:
    std::unique_ptr<rad_planck_mean> m_rad;
};

// -------------------------- WSGG ----------------------------------

class RadLibWSGG : public RadLibAdapterBase {
public:
    explicit RadLibWSGG(ThermoPhase* thermo, double fvsoot = 0.0)
        : RadLibAdapterBase(thermo)
        , m_rad(new rad_wsgg())
    {
        m_fvsoot = fvsoot;
    }

    void getBandProperties(std::vector<double>& kabs,
                           std::vector<double>& awts,
                           const RadComposition& comp) override
    {
        const double xH2O = _getX(comp.X, "H2O");
        const double xCO2 = _getX(comp.X, "CO2");
        const double xCO  = _getX(comp.X, "CO");
        const double xCH4 = _getX(comp.X, "CH4");

        kabs.clear(); awts.clear();
        m_rad->get_k_a(kabs, awts, comp.T, comp.P, m_fvsoot,
                       xH2O, xCO2, xCO, xCH4);
    }

    std::vector<std::string> requiredSpecies() const override {
        return {"H2O", "CO2", "CO", "CH4"};
    }

private:
    std::unique_ptr<rad_wsgg> m_rad;
};

// -------------------------- RCSLW ---------------------------------

class RadLibRCSLW : public RadLibAdapterBase {
public:
    // nGrayGases: number of grey gases (not including the clear gas). Typical 15–30.
    // Tref: recommended representative temperature in K (e.g., 1500)
    RadLibRCSLW(ThermoPhase* thermo,
                int nGrayGases = 25,
                double Tref = 1500.0,
                double Pref = OneAtm,
                double fvsoot = 0.0,
                double xH2O_ref = 0.1,
                double xCO2_ref = 0.1,
                double xCO_ref  = 0.0)
        : RadLibAdapterBase(thermo)
        , m_rad(new rad_rcslw(nGrayGases, Tref, Pref, fvsoot,
                              xH2O_ref, xCO2_ref, xCO_ref))
    {
        m_fvsoot = fvsoot;
    }

    void getBandProperties(std::vector<double>& kabs,
                           std::vector<double>& awts,
                           const RadComposition& comp) override
    {
        const double xH2O = _getX(comp.X, "H2O");
        const double xCO2 = _getX(comp.X, "CO2");
        const double xCO  = _getX(comp.X, "CO");
        const double xCH4 = _getX(comp.X, "CH4");

        kabs.clear(); awts.clear();
        m_rad->get_k_a(kabs, awts, comp.T, comp.P, m_fvsoot,
                       xH2O, xCO2, xCO, xCH4);
    }

    std::vector<std::string> requiredSpecies() const override {
        return {"H2O", "CO2", "CO", "CH4"};
    }

private:
    std::unique_ptr<rad_rcslw> m_rad;
};

// ------------------- Factory utility (optional) -------------------
// You can call these helpers from createRadiation1D(...)

inline std::unique_ptr<RadiationPropertyCalculator>
makeRadLibProps(const std::string& model,
                ThermoPhase* thermo,
                double fvsoot = 0.0,
                int nGray = 25,
                double Tref = 1500.0,
                double Pref = OneAtm)
{
    if (model == "RadLib.PlanckMean" || model == "radlib-pm") {
        return std::unique_ptr<RadiationPropertyCalculator>(
            new RadLibPlanckMean(thermo, fvsoot));
    } else if (model == "RadLib.WSGG" || model == "radlib-wsgg") {
        return std::unique_ptr<RadiationPropertyCalculator>(
            new RadLibWSGG(thermo, fvsoot));
    } else if (model == "RadLib.RCSLW" || model == "radlib-rcslw") {
        return std::unique_ptr<RadiationPropertyCalculator>(
            new RadLibRCSLW(thermo, nGray, Tref, Pref, fvsoot));
    } else {
        throw CanteraError("makeRadLibProps",
            "Unknown RadLib model '{}'", model);
    }
}

} // namespace Cantera

#else  // CT_ENABLE_RADLIB not defined

// If a build without RadLib tries to select a RadLib model, throw a clear error.
namespace Cantera {
inline std::unique_ptr<RadiationPropertyCalculator>
makeRadLibProps(const std::string&, ThermoPhase*, double fvsoot = 0.0,
                int nGray = 25, double Tref = 1500.0, double Pref = OneAtm)
{
    throw CanteraError("makeRadLibProps",
        "RadLib support is not enabled in this Cantera build. "
        "Rebuild Cantera with RadLib enabled (set 'system_radlib=y' for a system install, "
        "or 'system_radlib=n' to use the bundled RadLib submodule).");
}
} // namespace Cantera

#endif // CT_ENABLE_RADLIB
