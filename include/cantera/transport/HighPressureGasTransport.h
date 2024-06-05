/**
 *  @file HighPressureGasTransport.h
 *  Interface for class HighPressureGasTransport
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_HIGHPRESSUREGASTRAN_H
#define CT_HIGHPRESSUREGASTRAN_H

// Cantera includes
#include "GasTransport.h"
#include "cantera/transport/MixTransport.h"

namespace Cantera
{

//! Class MultiTransport implements transport properties for
//! high pressure gas mixtures.
/*!
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to %Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of  %Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * The implementation employs a method of corresponding states, using the Takahashi
 * @cite takahashi1975 approach for binary diffusion coefficients (using multicomponent
 * averaging rules for the mixture properties), and the Lucas method for the viscosity
 * of a high-pressure gas mixture. All methods are described in Poling et al.
 * @cite poling2001 (viscosity in Ch. 9, thermal conductivity in Ch. 10, and diffusion
 * coefficients in Ch. 11).
 *
 * @ingroup tranprops
 */
class HighPressureGasTransport : public MixTransport
{
protected:
    //! default constructor
    HighPressureGasTransport() = default;

public:
    string transportModel() const override {
        return "high-pressure";
    }

    double thermalConductivity() override;
    double viscosity() override;

    /**
     * Returns the matrix of binary diffusion coefficients
     *
     *      d[ld*j +  i] = rp*m_bdiff(i,j)*(DP)_R;
     *
     * @param ld    offset of rows in the storage
     * @param d     output vector of diffusion coefficients.  Units of m**2 / s
     */
    void getBinaryDiffCoeffs(const size_t ld, double* const d) override;

    friend class TransportFactory;

protected:
    double Tcrit_i(size_t i);
    double Pcrit_i(size_t i);
    double Vcrit_i(size_t i);
    double Zcrit_i(size_t i);

    double low_pressure_nondimensional_viscosity(double Tr, double FP, double FQ);
    double high_pressure_nondimensional_viscosity(double Tr, double Pr, double FP_low, double FQ_low, double P_vap, double P_crit);

    /**
     * @brief  Returns the quantum correction term for a species based on Tr
     * and MW, used in viscosity calculation.
     *
     * @param Q
     * @param Tr // Reduced temperature
     * @param MW // Molecular weight
     * @return double
     */
    double FQ_i(double Q, double Tr, double MW);

    /**
     * @brief  Returns the polarity correction term for a species based on Tr
     * and MW, used in viscosity calculation.
     *
     * @param Q
     * @param Tr // Reduced temperature
     * @param MW // Molecular weight
     * @return double
     */
    double FP_i(double mu_r, double Tr, double Z_c);

    /**
     * @brief Returns interpolated value of (DP)_R obtained from the data
     * in Table 2 of the Takahashi 1975 paper, given a value of the reduced
     * pressure (Pr) and reduced temperature (Tr).
     *
     * @param Pr  Reduced pressure
     * @param Tr  Reduced temperature
\    */
    double compute_correction_factor(double Pr, double Tr);
};



//These are the parameters that are needed to calculate the viscosity using the Chung method.
struct ChungMixtureParameters
{
    // Mixture critical properties used by the Chung viscosity model.
    double Vc_mix = 0;
    double Tc_mix = 0;

    // Values associated with the calculation of sigma and the molecular weight used in the Chung viscosity model.
    double sigma_mix = 0;
    double epsilon_over_k_mix = 0;
    double MW_mix = 0;

    // Values associated with the calculation of the Fc factor in the Chung viscosity model.
    double mu_mix = 0;
    double mu_r_mix = 0;
    double acentric_factor_mix = 0;
    double kappa_mix = 0;
};

//! Class ChungHighPressureTransport implements transport properties for
//! high pressure gas mixtures.
/*!
 * The implementation employs a method of corresponding states, using the Takahashi
 * @cite takahashi1975 approach for binary diffusion coefficients (using mixture
 * averaging rules for the mixture properties), and the Chung method for the viscosity
 *  and thermal conductivity of a high-pressure gas mixture. All methods are described in Poling et al.
 * @cite poling2001 (viscosity in Ch. 9, thermal conductivity in Ch. 10, and diffusion
 * coefficients in Ch. 11).
 *
 *
 *
 * @ingroup tranprops
 */
class ChungHighPressureGasTransport : public MixTransport
{
protected:
    //! default constructor
    ChungHighPressureGasTransport() = default;

public:
    string transportModel() const override {
        return "high-pressure-chung";
    }

    double viscosity() override;
    double thermalConductivity() override;

    /**
     * Returns the matrix of binary diffusion coefficients
     *
     *      d[ld*j +  i] = rp*m_bdiff(i,j)*(DP)_R;
     *
     * @param ld    offset of rows in the storage
     * @param d     output vector of diffusion coefficients.  Units of m**2 / s
     */
    void getBinaryDiffCoeffs(const size_t ld, double* const d) override;

    friend class TransportFactory;

protected:
    double Tcrit_i(size_t i);
    double Pcrit_i(size_t i);
    double Vcrit_i(size_t i);
    double Zcrit_i(size_t i);

    // Uses the low-pressure Chung viscosity model to calculate the viscosity
    // Defined by equation 9-4.10 in Poling et al.
    // Gives viscosity in units of micropoise
    // T must be units of K
    // MW must be units of kg/kmol
    // sigma must be units of Angstroms
    double low_pressure_viscosity(double T, double T_star, double MW, double acentric_factor, double mu_r, double sigma, double kappa);

    // Computes the high-pressure viscosity using the Chung method (Equation 9-6.18).
    // Gives viscosity in units of micropoise
    double high_pressure_viscosity(double T_star, double MW, double rho, double Vc, double Tc, double acentric_factor, double mu_r, double kappa);

    // Computes and store composition-dependent values of the parameters needed for the Chung viscosity model.
    void compute_mixture_parameters(ChungMixtureParameters& params);



    /**
     * @brief Returns interpolated value of (DP)_R obtained from the data
     * in Table 2 of the Takahashi 1975 paper, given a value of the reduced
     * pressure (Pr) and reduced temperature (Tr).
     *
     * @param Pr  Reduced pressure
     * @param Tr  Reduced temperature
\    */
    double compute_correction_factor(double Pr, double Tr);
};






}
#endif
