/**
 * @file Coverage.h
 * Header for reaction rates that occur at interfaces.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_COVERAGE_H
#define CT_COVERAGE_H

#include "cantera/kinetics/Arrhenius.h"
#include "MultiRate.h"

namespace Cantera
{

class AnyMap;

/**
 *  @defgroup surfaceGroup  Coverage-dependent rate parameterizations
 *
 *  This section describes the parameterizations used to describe rate
 *  parameterization that involve interfaces.
 *
 *  @ingroup chemkinetics
 */

//! Base class for reaction rate parameterizations that involve interfaces
class Coverage
{
public:
    Coverage();

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param dependencies  Coverage dependencies
     *  @param units  Unit system
     */
    void setCoverageDependencies(const AnyMap& dependencies,
                                 const UnitSystem& units=UnitSystem());

    void getCoverageDependencies(AnyMap& dependencies, bool asVector=false) const;

    //! Add a coverage dependency for species *sp*, with exponential dependence
    //! *a*, power-law exponent *m*, and activation energy dependence *e*,
    //! where *e* is in Kelvin, i.e. energy divided by the molar gas constant.
    void addCoverageDependence(std::string sp, double a, double m, double e);

    //! Set species indices within coverage array
    void setSpecies(const Kinetics& kin);
    void setSpecies(const std::vector<std::string>& species);

    //! Update reaction rate parameters
    //! @param shared_data  data shared by all reactions of a given type
    void updateFromStruct(const CoverageData& shared_data) {
        if (shared_data.ready) {
            m_siteDensity = shared_data.siteDensity;
        }
        if (m_indices.size() != m_cov.size()) {
            // object is not set up correctly (setSpecies needs to be run)
            m_acov = NAN;
            m_ecov = NAN;
            m_mcov = NAN;
            return;
        }
        m_acov = 0.0;
        m_ecov = 0.0;
        m_mcov = 0.0;
        for (auto& item : m_indices) {
            m_acov += m_ac[item.first] * shared_data.coverages[item.second];
            m_ecov += m_ec[item.first] * shared_data.coverages[item.second];
            m_mcov += m_mc[item.first] * shared_data.logCoverages[item.second];
        }
    }

    //! Return site density [kmol/m^2]
    double siteDensity() const {
        return m_siteDensity;
    }

    //! Set site density [kmol/m^2]
    //! @internal  used for testing purposes only
    //! Note that this quantity is not an independent variable and will be
    //! overwritten during an update of the thermodynamic state.
    void setSiteDensity(double siteDensity) {
        m_siteDensity = siteDensity;
    }

protected:
    double m_siteDensity; //!< Site density [kmol/m^2]
    double m_acov; //!< Coverage contribution to pre-exponential factor
    double m_ecov; //!< Coverage contribution to activation energy
    double m_mcov; //!< Coverage term in reaction rate
    std::map<size_t, size_t> m_indices; //!< Map holding indices of coverage species
    std::vector<std::string> m_cov; //!< Vector holding names of coverage species
    vector_fp m_ac; //!< Vector holding coverage-specific exponential dependence
    vector_fp m_ec; //!< Vector holding coverage-specific activation energy dependence
    vector_fp m_mc; //!< Vector holding coverage-specific power-law exponents
};


}
#endif
