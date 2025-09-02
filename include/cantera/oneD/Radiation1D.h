//! @file Radiation1D.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef RADIATION1D_H
#define RADIATION1D_H

#include "Domain1D.h"
#include "cantera/base/Array.h"
#include "cantera/thermo/ThermoPhase.h"

#include <functional>

namespace Cantera
{

/** Stores the temperature, pressure, an optional soot fraction (fvSoot),
 * and a map of species mole fractions. Allows the property calculator
 * classes to retrieve the local state information they need .
 *
 * The temperature is given in Kelvin [K], the pressure in Pascals [Pa],
 * and fvSoot is a dimensionless volume fraction for soot. The map 'X'
 * holds species names as keys and their mole fractions (unitless) as values.
 */
struct RadComposition {
    double T = 0.0;      //! Temperature (K)
    double P = 0.0;      //! Pressure (Pa)
    double fvSoot = 0.0; //! Soot volume fraction

    Composition X; //! Map of name->mole fraction
};

/** Base class for radiation property calculators.
 *
 * Responsible for computing the spectral absorption coefficients (kabs) and
 * weighting factors (awts) for a given thermodynamic state. Different models
 * e.g. polynomial fits, tabular data, or external libraries such as RadLib
 * are implemented by deriving from this class.
 *
 * The data produced by getBandProperties() are used by a RadiationSolver to compute
 * the net radiative heat loss at each point in a 1D domain.
 */
class RadiationPropertyCalculator {
public:
    virtual ~RadiationPropertyCalculator() = default;

    /** Calculate absorption coefficients and weighting factors for each spectral band.
     *
     * The size of `kabs` and `awts` determine how many "gray gases" or bands are
     * used. For a simple Planck mean approach, there may be only one band. For
     * multi-band or weighted-sum-of-grey-gases (WSGG) models, there could be multiple.
     *
     * @param kabs  A vector to be filled with absorption coefficients (k_i).
     * @param awts  A vector to be filled with weighting factors (a_i).
     * @param comp  A RadComposition containing T, P, composition, etc.
     */
    virtual void getBandProperties(std::vector<double>& kabs,
                                   std::vector<double>& awts,
                                   const RadComposition& comp) = 0;

    // Optional: list of species names required by this calculator.
    // If empty, the caller may provide all species; otherwise, the caller can
    // optimize by providing only these species in RadComposition::X.
    // @since New in Cantera (radiation plugin)
    virtual std::vector<std::string> requiredSpecies() const { return {}; }
};


/* Commented out for now until RadLib dependency is resolved
class RadLibPlanckMean : public RadiationPropertyCalculator {
public:
    RadLibPlanckMean() {
        m_rad = new rad_planck_mean();
    }
    ~RadLibPlanckMean() {
        delete m_rad;
    }

    void getBandProperties(std::vector<double>& kabs,
                           std::vector<double>& awts,
                           double T, double P, const RadComposition& comp) override
    {
        m_rad->get_k_a(kabs, awts, T, P, comp.fvSoot, comp.xH2O, comp.xCO2, comp.xCO, comp.xCH4);
    }

private:
    rad* m_rad; // pointer to rad_planck_mean
};
*/

/* Reads species-specific Planck-mean absorption coefficient data from a YAML file
 * (radiation-properties.yaml) if available, falling back to polynomial approximations
 * for CO2 and H2O otherwise. A `fit-type` of "table" or "polynomial" can be specified
 * in the YAML data.
 *
 * The table-based data uses interpolation for a discrete set of temperatures,
 * whereas the polynomial data uses functional fits.
 *
 * The `fit-type` of `polynomial` is uses the model described below:
 *
 * Polynomial lines calculate the species Planck coefficients for H2O and CO2. The
 * data for the lines are taken from the RADCAL program @cite RADCAL.
 * The coefficients for the polynomials are taken from
 * [TNF Workshop](https://tnfworkshop.org/radiation/) material.
 *
 * The `fit-type` of `table` is uses the model described below.
 *
 * Spectra for molecules are downloaded with HAPI library from // https://hitran.org/hapi/
 * [R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski,
 * HITRAN Application Programming Interface (HAPI): A comprehensive approach
 * to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177,
 * 15-30 (2016), https://doi.org/10.1016/j.jqsrt.2016.03.005].
 *
 * Planck mean optical path lengths are what are read in from a YAML input file.
*/
class TabularPlanckMean : public RadiationPropertyCalculator {
public:
    /**
     * The constructor will attempt to parse radiation data from
     * "radiation-properties.yaml". If that file doesn't exist, a warning is
     * issued and polynomial defaults for CO2 and H2O are used.
     *
     * @param thermo Pointer to a ThermoPhase object which provides species names
     *        and other properties. This is needed to match species found in the
     *        YAML database to the actual species in the simulation.
     */
    TabularPlanckMean(ThermoPhase* thermo);

    /** Calculate absorption coefficients and weighting factors for each band.
     * This method sums absorption contributions from all absorbing species for
     * which the table or polynomial data is defined. The final result is stored
     * as a single-band coefficient (kabs.size()==1), with awts.size()==1=1.0,
     * representing a gray approximation.
     *
     * @param kabs  A vector to store absorption coefficients (k_i).
     * @param awts  A vector to store weighting factors (a_i).
     * @param comp  The RadComposition struct with T, P, and species mole fractions.
     */
    void getBandProperties(std::vector<double>& kabs, std::vector<double>& awts,
                           const RadComposition& comp) override;

    std::vector<std::string> requiredSpecies() const override;

private:
    /** Parse optional YAML data from "radiation-properties.yaml".
     * If the file is not found, a warning is issued and default polynomial data
     * for H2O and CO2 is used. If it is found, then species listed in the file
     * are read into 'm_PMAC' along with their "fit-type" and associated
     * coefficients. This might be polynomial or tabulated data.
     *
     * The method also ensures that H2O and CO2 have some default data even if
     * the file does not provide them.
     */
    void parseRadiationData();

    /** Compute polynomial-based absorption coefficient.
     *
     * This evaluates a polynomial that has a form given as:
     *
     * kabs = c0 + c1(1000/T) + c2(1000/T)^2 + c3(1000/T)^3 + c4(1000/T)^4 + c5(1000/T)^5
     *
     * The value computed is the Plank mean absorption coefficient. This is just one way
     * to represent the variation of the absorption coefficient with temperature, and
     * it used often for species such as H2O and CO2. The units of the coefficients are
     * (m-1 atm-1) and the temperature is in Kelvin( see https://tnfworkshop.org/radiation/).
     * This function converts the units to m-1 Pa-1.
     *
     * @param coefficients  The polynomial coefficients for the species.
     * @param temperature   The local temperature in K.
     * @return The Planck-mean absorption coefficient for that species
     *         at the given temperature in units of m-1 Pa-1.
     */
    double calculatePolynomial(const vector<double>& coefficients, double temperature);

    /** Compute an absorption coefficient using a log-linear interpolation.
     *
     * Use log-linear interpolation to compute an absorption  coefficient from a
     * table of Planck mean optical-path-length values ('data') at discrete
     * 'temperatures'. Units of the optical path length are meters and the units of
     * temperature are Kelvin.
     *
     * The tables hold values of the optical path length (OPL) for a gas at
     * different temperatures. The absorption coefficient is the inverse of the
     * OPL i.e. kabs = 1.0 / OPL
     *
     * The method uses the following algorithm:
     *  alpha = 1.0 / data[i]
     *  ln(alpha) is interpolated linearly vs. T in the bracket
     *  [temperatures[i-1], temperatures[i]] using the following formula:
     *
     * ln(alpha) = ln(1/v1) + ( ln(1/v2) - ln(1/v1) ) * (T - t1)/(t2 - t1)
     *
     * If T is below the lowest or above the highest table entry, the boundary value
     * without interpolation is used.
     *
     * @param temperatures  Sorted vector of tabulated temperatures
     * @param data          Corresponding tabulated OPL data
     *                      (optical path lengths)
     * @param temperature   Query temperature
     * @returns The absorption coefficient at 'temperature'
     *
     * @param data          A vector of corresponding absorption data.
     * @return The interpolated absorption coefficient at the given temperature
     *         in units of m-1 Pa-1.
     */
    double interpolateTable(const vector<double>& temperatures,
                            const vector<double>& data, double temperature);

private:
    ThermoPhase* m_thermo; //!< Pointer to the ThermoPhase object

    map<string, size_t> m_absorptionSpecies; //!< Absorbing species mapping names to indices
    AnyMap m_PMAC; //!< Absorption coefficient data for each species
};

/** Base class for radiation solvers.
 *
 * Computes the net radiative heat loss (or gain) from
 * absorption coefficients, weighting factors, boundary emissivities,
 * and local temperature. Different solver implementations should derive from
 * this class.
 */
class RadiationSolver {
public:
    virtual ~RadiationSolver() = default;

    // compute the radiative heat loss given boundary conditions, geometry, etc.
    virtual double computeHeatLoss(const std::vector<double>& kabs,
                                   const std::vector<double>& awts,
                                   double T,
                                   double boundaryRadLeft,
                                   double boundaryRadRight) = 0;
};

/*
 * The simple radiation model used was established by Liu and Rogg
 * @cite liu1991. This model considers the radiation of CO2 and H2O.
 *
 * This model uses the *optically thin limit* and the *gray-gas approximation*,
 * calculating a volumetric heat loss (qdot) using Planck-mean absorption
 * coefficients, boundary emissivities, and local temperature.
 *
 * Typically, qdot for each spectral band i is given by:
 *
 * \f[
 *   \dot{q}_i = 2 k_i \left( 2 \sigma T^4 - E_{\text{left}} - E_{\text{right}} \right)
 * \f]
 *
 * Summed over all bands i, weighted by a_i. The 2 factor comes from the assumption
 * that radiation can escape in both directions in a 1D domain, ignoring scattering.
 */
class OpticallyThinSolver : public RadiationSolver {
public:

    //! Calculate optically thin radiative heat loss for each band, summing the
    //! contributions.
    double computeHeatLoss(const std::vector<double>& kabs,
                           const std::vector<double>& awts,
                           double T,
                           double boundaryRadLeft,
                           double boundaryRadRight) override
    {
        // Sum over each band.
        double sum = 0.0;
        double sigma = StefanBoltz;
        for (size_t i=0; i<kabs.size(); ++i) {
            sum += awts[i] * 2.0 * kabs[i] *
                   (2.0 * sigma * std::pow(T,4) - boundaryRadLeft - boundaryRadRight);
        }
        return sum;
    }
};


/**
 * The Radiation1D class ties together a RadiationPropertyCalculator and a
 * RadiationSolver to compute volumetric radiative heat losses at each grid
 * point (j). It fetches the local temperature and species mole fractions
 * from the solution vector (using the provided lambda functions) and
 * optionally applies boundary emissivities to represent radiative flux
 * at the domain edges.
 *
*/
class Radiation1D {
public:

    /**
     * Constructor for a Radiation1D object.
     *
     * @param thermo          A pointer to the ThermoPhase object for species info
     * @param pressure        The operating pressure [Pa] (assumed uniform)
     * @param points          Number of grid points in the domain
     * @param temperatureFunction  A lambda to retrieve T(x, j) from the solution
     * @param moleFractionFunction A lambda to retrieve X_k(x, j) from the solution
     * @param props           A unique pointer to a RadiationPropertyCalculator
     * @param solver          A unique pointer to a RadiationSolver
     */
    Radiation1D(ThermoPhase* thermo, double pressure, size_t points,
                std::function<double(const double*, size_t)> temperatureFunction,
                std::function<double(const double*, size_t, size_t)> moleFractionFunction,
                std::unique_ptr<RadiationPropertyCalculator> props,
                std::unique_ptr<RadiationSolver> solver);

    /** Compute radiative heat loss from jmin to jmax and fill the qdotRadiation array.
     *
     *  This method extracts T and mole fractions from the solution at each grid point.
     *  It then uses the RadiationPropertyCalculator to get the band properties, and then
     *  uses the RadiationSolver to get the heat loss and stores it in qdotRadiation.
     *
     * @param x             Pointer to the solution vector for the 1D domain.
     * @param jmin          The first grid index to compute.
     * @param jmax          One past the last grid index to compute.
     * @param qdotRadiation A vector of size at least jmax, which will be filled
     *                      with the volumetric radiative heat loss at each grid point.
     */
    void computeRadiation(double* x, size_t jmin, size_t jmax,
                          vector<double>& qdotRadiation);

    /**
     * Sets the emissivities for the left and right boundary values in the
     * radiative term.
     */
    void setBoundaryEmissivities(double e_left, double e_right);

    //! Return emissivity at left boundary
    double leftEmissivity() const {
        return m_epsilon_left;
    }

    //! Return emissivity at right boundary
    double rightEmissivity() const {
        return m_epsilon_right;
    }

private:
    ThermoPhase* m_thermo; //!< Pointer to the ThermoPhase object
    double m_press; //!< Pressure in Pa
    size_t m_points; //!< Number of grid points

    //! Property calculator for absorption coefficients
    std::unique_ptr<RadiationPropertyCalculator> m_props;

    //! Solver for radiative heat loss
    std::unique_ptr<RadiationSolver> m_solver;

    //! Emissivity of the surface to the left and right of the domain. Used for calculating
    //! radiative heat loss.
    double m_epsilon_left = 0.0;
    double m_epsilon_right = 0.0;

    //! Lambda function to get temperature at a given point
    std::function<double(const double*, size_t)> m_T;

    //! Lambda function to get mole fraction at a given point
    std::function<double(const double*, size_t, size_t)> m_X;
};

} // namespace Cantera

// Bring in RadLib-backed property adapters and factory helper. This header is
// safe to include unconditionally: if CT_ENABLE_RADLIB is not defined,
// it provides a stub that throws a clear error when RadLib models are selected.
#include "Radiation1D_RadLib.h"

namespace Cantera
{
/**
 * Create a Radiation1D instance based on the selected property model and solver.
 *
 * @param propertyModel  String specifying which property calculator to use
 * @param solverModel    String specifying which solver to use
 * @param thermo         Pointer to ThermoPhase
 * @param pressure       Pressure in Pa
 * @param points         Number of grid points
 * @param Tfunc          Lambda for temperature
 * @param Xfunc          Lambda for mole fractions
 * @param e_left         Emissivity at left boundary
 * @param e_right        Emissivity at right boundary
 *
 * @return Radiation1D object
 */
inline std::unique_ptr<Radiation1D> createRadiation1D(
    const std::string& propertyModel,
    const std::string& solverModel,
    ThermoPhase* thermo,
    double pressure,
    size_t points,
    std::function<double(const double*, size_t)> Tfunc,
    std::function<double(const double*, size_t, size_t)> Xfunc,
    double e_left,
    double e_right
)
{
    // Create the RadiationPropertyCalculator
    std::unique_ptr<RadiationPropertyCalculator> props;
    if (propertyModel == "TabularPlanckMean") {
        props = std::make_unique<TabularPlanckMean>(thermo);
    }
    else if (propertyModel == "RadLib.PlanckMean" || propertyModel == "radlib-pm") {
       props = makeRadLibProps("RadLib.PlanckMean", thermo, 0.0);
    } else if( propertyModel == "RadLib.WSGG" || propertyModel == "radlib-wsgg") {
        props = makeRadLibProps("RadLib.WSGG", thermo, 0.0);
    } else if( propertyModel == "RadLib.RCSLW" || propertyModel == "radlib-rcslw") {
        // Default to 25 gray gases and Tref=1500 K for RCSLW if not specified elsewhere
        props = makeRadLibProps("RadLib.RCSLW", thermo, 0.0, 25, 1500.0, pressure);
    }
    else {
        throw CanteraError("createRadiation1D",
            "Unknown property model: " + propertyModel);
    }

    // Create the RadiationSolver
    std::unique_ptr<RadiationSolver> solver;
    if (solverModel == "OpticallyThin") {
        solver = std::make_unique<OpticallyThinSolver>();
    }
    else {
        throw CanteraError("createRadiation1D",
            "Unknown solver model: " + solverModel);
    }

    // Build the Radiation1D object
    auto rad = std::make_unique<Radiation1D>(
        thermo, pressure, points,
        std::move(Tfunc), std::move(Xfunc),
        std::move(props), std::move(solver)
    );

    rad->setBoundaryEmissivities(e_left, e_right);

    return rad;
};


} // namespace Cantera

#endif // RADIATION1D_H
