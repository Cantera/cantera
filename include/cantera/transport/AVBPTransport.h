/**
 *
 *  @file AVBPTransport.h
 *   Header file defining class AVBPTransport which implements
 *   the simplified transport model used in the solver AVBP
 */

/* $Author: B. Franzelli (v. 1.7) $
 * $Revision: A. Felden (v 2.1-2.3) $
 * $Date: 01/2018 $
 */


#ifndef CT_AVBPTRAN_H
#define CT_AVBPTRAN_H

#include "GasTransport.h"
#include "cantera/numerics/DenseMatrix.h"

// STL includes
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>

namespace Cantera
{

// class GasTransportParams;

/**
 * Class AVBPTransport
 * Constant Sch for each species, Pr, and simplified viscosity
 */
class AVBPTransport : public GasTransport
{

public:
    //! Default constructor.
    AVBPTransport() = default;

    //! Return the model id for transport
    /*!
     * @return cAVBPAverage
     */
    // virtual int model() const {
    //     warn_deprecated("AVBPTransport::model",
    //                     "To be removed after Cantera 2.3.");
    //     return cAVBPTransport;
    // }
    string transportModel() const override {
        return (m_mode == CK_Mode) ? "AVBP-CK" : "AVBP";
    }
    // virtual std::string transportType() const {
    //     return "AVBP";
    // }

    //! Return the thermal diffusion coefficients
    //virtual void getThermalDiffCoeffs(double* const dt);

    //! Returns the mixture thermal conductivity
    double thermalConductivity() override;

    //! Get the Electrical mobilities (m^2/V/s).
    //virtual void getMobilities(double* const mobil);

    //! Update the internal parameters whenever the temperature has changed
    void update_T() override;

    //! Update the internal parameters whenever the concentrations have changed
    void update_C() override;

    //virtual void getSpeciesFluxes(size_t ndim,
    //                              const double* const grad_T,
    //                              size_t ldx,
    //                              const double* const grad_X,
    //                              size_t ldf, double* const fluxes);

    //! Initialize the transport object
    //virtual bool initGas(GasTransportParams& tr);
    void init(ThermoPhase* thermo, int mode=0, int log_level=0) override;

    //! Viscosity of the mixture
    virtual double viscosity();

    virtual void getSpeciesViscosities(double* const visc) {
        update_T();
        updateViscosity_T();
        std::copy(m_visc.begin(), m_visc.end(), visc);
    }

    //! Mixture diffusion coefficients [m^2/s].
    virtual void getMixDiffCoeffs(double* const d);

    virtual void read_mixture(std::string s);

    size_t avbp_ipea;
    vector<double> avbp_pea_coeffs;

protected:

    //! Calculate the pressure from the ideal gas law
    double pressure_ig();

    double m_lambda;
    //bool m_debug;

    // AVBP variables
    vector<double> avbp_Sch;
    vector<double> avbp_Le;
    double avbp_Prandtl;
    double avbp_mu0;
    double avbp_T0;
    double avbp_beta;
    std::string avbp_fuel;

};
}
#endif
