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
 * @note: All equations that are cited in this implementation are from the 5th edition
 * of the book "The Properties of Gases and Liquids" by Poling, Prausnitz, and O'Connell.
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

    /**
    * Calculates the high-pressure mixture thermal conductivity using the Chung method.
    *
    * The Chung method is described in on page 10.23.
    *
    * @f[
    *    \lambda = \frac{31.2 \eta^0 \Psi}{M'} \left( G_1^{-1} + B_6 y \right) + q B_7 y^2 T_r^{1/2} G_2
    * @f]
    ** where:
    *
    * @f[
    *   \begin{align*}
    *       \lambda &= \text{thermal conductivity, W/(m·K)} \\
    *       \eta^0 &= \text{low-pressure gas viscosity, N·s/m}^2 \\
    *       M' &= \text{molecular weight, kg/mol} \\
    *       \Psi &= f(C_v, \omega, T_r) \text{ [as defined under Eq. (10-3.14)]} \\
    *       q &= 3.586 \times 10^{-3} \left( \frac{T_c}{M'} \right)^{1/2} V_c^{2/3} \\
    *       T &= \text{temperature, K} \\
    *       T_c &= \text{critical temperature, K} \\
    *       T_r &= \text{reduced temperature, } \frac{T}{T_c} \\
    *       V_c &= \text{critical volume, cm}^3/\text{mol} \\
    *       \gamma &= \frac{V_c}{6V}
    *   \end{align*}
    * @f]
    *
    * The details about the parameters and constants used in the expression are
    * found in the Poling book.
    *
    * The mixture values of the pseudo-critical temperature and other model parameters
    * are calculated using the Chung mixing rules defined on page 9.25.
    *
    * The mixture value of the specific heat is computed using equation 10-6.6, which
    * is the mole fraction weighted sum of the pure species specific heats.
    *
    * @f[
    *   C_{v,m} = \sum_i y_i C_{v,i}
    * #f]
    *
    */
    double thermalConductivity() override;

    /**
     * Returns the matrix of binary diffusion coefficients, augmented by the
     * Takahashi correction factor.
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

    /**
    * Returns the low-pressure mixture viscosity using the Chung method in micropoise.
    *
    * Defined by equation 9-4.10.
    *
    * @f[
    *   \eta = 26.69 F_c \frac{(M*T)^(\frac{1}{2})}{\sigma^2 \Omega}
    * @f]
    *
    * T must be in units of K, MW must be units of kg/kmol, and sigma must be units of Angstroms
    *
    * This function is structured such that it can be used for pure species or mixtures, with the
    * only difference being the values that are passed to the function (pure values versus mixture values).
    */
    double low_pressure_viscosity(double T, double T_star, double MW, double acentric_factor,
                                 double mu_r, double sigma, double kappa);

    // Computes the high-pressure viscosity using the Chung method (Equation 9-6.18).
    // Gives viscosity in units of micropoise
    double high_pressure_viscosity(double T_star, double MW, double rho, double Vc, double Tc,
                                   double acentric_factor, double mu_r, double kappa);

    /**
    * Returns the composition-dependent values of the parameters needed for the Chung viscosity model.
    *
    * The equations for the mixing rules defined on page 9.25 for the Chung method's
    * composition dependent parameters.
    */
    void compute_mixture_parameters(ChungMixtureParameters& params);

    /**
    * Computes the high-pressure thermal conductivity using the Chung method (Equation 10-5.5).
    *
    * Gives thermal conductivity in units of W/m/K.
    *
    * This function is structured such that it can be used for pure species or mixtures, with the
    * only difference being the values that are passed to the function (pure values versus mixture values).
    *
    * This method utilizes the low-pressure Chung viscosity as that is a required parameter in the model, and
    * thus makes a call to the low pressure viscosity implementation. This is why it requires parameters
    * typically associated with the viscosity calculation.
    *
    * M_prime (M' in the model) has units of kg/mol, and is just the molecular weight (kg/kmol) divided by 1000.
    */
    double high_pressure_thermal_conductivity(double T, double T_star, double MW, double rho, double Cv, double Vc, double Tc, double sigma, double acentric_factor, double mu_r, double kappa);

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
