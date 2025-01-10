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

/**
 * Computes the radiative heat loss vector over points jmin to jmax and stores
 * the data in the qdotRadiation variable.
 *
 * The `fit-type` of `polynomial` is uses the model described below.
 *
 * The simple radiation model used was established by Liu and Rogg
 * @cite liu1991. This model considers the radiation of CO2 and H2O.
 *
 * This model uses the optically thin limit and the gray-gas approximation to
 * simply calculate a volume specified heat flux out of the Planck absorption
 * coefficients, the boundary emissivities and the temperature. Polynomial lines
 * calculate the species Planck coefficients for H2O and CO2. The data for the
 * lines are taken from the RADCAL program @cite RADCAL.
 * The coefficients for the polynomials are taken from
 * [TNF Workshop](https://tnfworkshop.org/radiation/) material.
 *
 *
 * The `fit-type` of `table` is uses the model described below.
 *
 * Spectra for molecules are downloaded with HAPI library from // https://hitran.org/hapi/
 * [R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski,
 * HITRAN Application Programming Interface (HAPI): A comprehensive approach
 *  to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177,
 *   15-30 (2016), https://doi.org/10.1016/j.jqsrt.2016.03.005].
 *
 * Planck mean optical path lengths are what are read in from a YAML input file.
*/
class Radiation1D {
public:
    Radiation1D(ThermoPhase* thermo, double pressure, size_t points,
                std::function<double(const double*, size_t)> temperatureFunction,
                std::function<double(const double*, size_t, size_t)> moleFractionFunction);

    // Parse radiation data from YAML input
    void parseRadiationData();

    // Compute radiative heat loss
    void computeRadiation(double* x, size_t jmin, size_t jmax, vector<double>& qdotRadiation);

     //! Set the emissivities for the boundary values
    /*!
     * Reads the emissivities for the left and right boundary values in the
     * radiative term and writes them into the variables, which are used for the
     * calculation.
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
    ThermoPhase* m_thermo;
    double m_press;
    size_t m_points;

    map<string, size_t> m_absorptionSpecies; //!< Absorbing species
    AnyMap m_PMAC; //!< Absorption coefficient data for each species

    //! Emissivity of the surface to the left of the domain. Used for calculating
    //! radiative heat loss.
    double m_epsilon_left = 0.0;

    //! Emissivity of the surface to the right of the domain. Used for calculating
    //! radiative heat loss.
    double m_epsilon_right = 0.0;

    // Helper functions
    double calculatePolynomial(const vector<double>& coefficients, double temperature);
    double interpolateTable(const vector<double>& temperatures, const vector<double>& data, double temperature);

    //! Lambda function to get temperature at a given point
    std::function<double(const double*, size_t)> m_T;

    //! Lambda function to get mole fraction at a given point
    std::function<double(const double*, size_t, size_t)> m_X;
};

} // namespace Cantera

#endif // RADIATION1D_H
