
/**
 * @file PlasmaPhase.h
 * Header file for class PlasmaPhase.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAPHASE_H
#define CT_PLASMAPHASE_H

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/numerics/eigen_sparse.h"

namespace Cantera
{

/**
 * Base class for a phase with plasma properties. This class manages the
 * plasma properties such as electron energy distribution function (EEDF).
 * For isotropic velocity space such as Maxweillian and Druyvesteyn, a
 * generalized electron energy distribution can be expressed as [1,2],
 *   \f[
 *          f(\epsilon) = c_1 * \frac{\sqrt{epsilon}}{\epsilon_m^{3/2}}
 *          exp(-c_2 * (\frac{\epsilon}{\epsilon_m})^x),
 *   \f]
 * where \f$ x = 1 \f$ and \f$ x = 2 \f$ correspond to the Maxwellian and
 * Druyvesteyn (default) electron energy distribution, respectively.
 * Refereces:
 * [1] J. T. Gudmundsson. On the effect of the electron energy distribution on the
 * plasma parameters of an argon discharge: a global (volume-averaged) model study.
 * Plasma Sources Science and Technology, 10.1 (2001): 76.
 * doi: https://doi.org/10.1088/0963-0252/10/1/310
 * [2] H. Khalilpour and G. Foroutan. The effects of electron energy distribution
 * function on the plasma sheath structure in the presence of charged nanoparticles
 * Journal of Plasma Physics 86.2 (2020).
 * doi: https://doi.org/10.1017/S0022377820000161
 * @todo Implement electron Boltzmann equation solver to solve EEDF.
 * @ingroup phase
 */
class PlasmaPhase: public IdealGasPhase
{
public:
    //! Construct and initialize an PlasmaPhase ThermoPhase object
    //! directly from an ASCII input file. The constructor initializes the electron
    //! energy distribution to be Druyvesteyn distribution (m_x = 2.0). The initial
    //! electron energy grid is set to a linear space which starts at 0.01 eV and ends
    //! at 1 eV with 1000 points.
    /*!
     * @param inputFile Name of the input file containing the phase definition
     *                  to set up the object. If blank, an empty phase will be
     *                  created.
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     */
    explicit PlasmaPhase(const std::string& inputFile="",
                         const std::string& id="");

    virtual std::string type() const {
        return "plasma";
    }

    //! Set electron energy grid.
    //!     @param grid The vector of electron energy grid. [eV]
    void setElectronEnergyGrid(const vector_fp& grid);

    //! Get electron energy grid.
    //!     @param grid The vector of electron energy grid. [eV]
    void getElectronEnergyGrid(vector_fp& grid) const;

    //! Update electron energy distribution.
    void updateIsotropicElectronEnergyDistrb();

    //! Set the internally stored electron temperature of the phase (K).
    //!     @param Te Electron temperature in Kelvin
    virtual void setElectronTemperature(const double Te);

protected:
    // Method to get electron energy distribution function
    std::string m_electronEnergyDistrbMethod;

    // Electron energy order in the exponential term
    double m_x;

    //! Number of points of electron energy grid
    size_t m_nPoints;

    //! electron energy distribution vector
    Eigen::VectorXd m_EE_grid;

    //! Normalized electron energy distribution vector
    Eigen::VectorXd m_EE_distrb_vec;

    //! Mean electron energy
    double m_EE_mean;
};

}

#endif
