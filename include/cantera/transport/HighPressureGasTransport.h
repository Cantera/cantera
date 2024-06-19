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

/**
 * @brief Returns interpolated value of (DP)_R obtained from the data
 * in Table 2 of the Takahashi 1975 paper, given a value of the reduced
 * pressure (Pr) and reduced temperature (Tr).
 *
 * @param Pr  Reduced pressure
 * @param Tr  Reduced temperature
 */
double takahashi_correction_factor(double Pr, double Tr);


/**
 * @brief Returns the value of the Neufeld collision integral for a given
 * dimensionless temperature. Implementation of equation 9-4.3.
 * Applicable over the range of 0.3 <= T_star <= 100.
 *
 * @param T_star  Dimensionless temperature (Defined in Equation 9-4.1)
 */
double neufeld_collision_integral(double T_star);


//These are the parameters that are needed to calculate the viscosity using the Lucas method.
struct LucasMixtureParameters
{
    double FQ_mix_o;
    double FP_mix_o;
    double Tr_mix;
    double Pr_mix;
    double Pc_mix;
    double Tc_mix;
    double MW_mix;
    double P_vap_mix;
};

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
 * @cite takahashi1975 approach for binary diffusion coefficients (using mixture
 * averaging rules for the mixture properties), the Lucas method for the viscosity, and
 * a method from Ely and Hanley. All methods are described in Poling et al.
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

    /**
     * Returns the mixture high-pressure thermal conductivity in W/m/K
     * a method by Ely and Hanley.
     *
     */
    double thermalConductivity() override;

     /**
     * Returns the mixture high-pressure viscosity in Pa*s using the Lucas method.
     *
     * This uses the approach described in chapter 9-7.
     *
     * The mixture pseudo-critical temperature and pressure are calculated using
     * Equations 9-5.18 and 9-5.19. The mixture molecular weight is computed using
     * Equation 9-5.20. The mixture values of the low-pressure polarity and quantum
     * correction factors are computed using Equations 9-5.21 and 9-5.22.
     *
     */
    double viscosity() override;

    /**
     * Returns the matrix of binary diffusion coefficients using the Takahashi
     * correction factor.
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
     * Returns the composition-dependent values of the parameters needed for
     * the Lucas viscosity model.
     *
     * The equations for the mixing rules defined on page 9.23 for the Lucas method's
     * composition dependent parameters. The primary mixing rules are defined below,
     * and the reduced properties are just the properties divided by the pseudo-critical
     * mixture properties defined below.
     *
     * @f[
     *  T_{\t{c,m}} = \sum_i y_i T_{\t{c,i}}
     * @f]
     *
     * @f[
     *  P_{\t{c,m}} = R T_{\t{c,m}} \frac{\sum_i y_i Z_{\t{c,i}}}{\sum_i y_i V_{\t{c,i}}}
     * @f]
     *
     * @f[
     *  M_m = \sum y_i M_i
     * @f]
     *
     * @f[
     *  F_{P,m}^{\t{o}} = \sum y_i F_{P,i}^{\t{o}}
     * @f]
     *
     * @f[
     *   F_{Q,m}^{\t{o}} = \left ( \sum y_i F_{Q,i}^{\t{o}} \right ) A
     * @f]
     *
     * @f[
     *   A = 1 - 0.01 \left ( \frac{M_H}{M_L} \right )^{0.87}
     * @f]
     *
     * For $\frac{M_H}{M_L} > 9$ and $ 0.05 < y_H < 0.7$, otherwise A = 1. In the
     * above equation, $M_H$ and $M_L$ are the molecular weights of the heaviest
     * and lightest components in the mixture, and $y_H$ is the mole fraction of
     * the heaviest component.
     *
     * While it isn't returned as a parameter, the species-specific reduced dipole
     * moment is used to compute the mixture polarity correction factor. It is
     * defined as:
     *
     * @f[
     *   \mu_r = 52.46 \frac{\mu^2 P_{\t{c,i}}}{T_{\t{c,i}}
     * @f]
     *
     */
    void compute_mixture_parameters(LucasMixtureParameters& params);

    /**
     * Returns the non-dimensional low-pressure mixture viscosity in using the Lucas method.
     *
     * Defined by equation 9-4.16.
     *
     * @f[
     *   \eta \xi = [0.807 T_r^{0.618} - 0.357 e^{-0.449 T_r} + 0.340e^{-4.058 T_r} + 0.018] F_P^{\t{o}} F_Q^{\t{o}}
     * @f]
     *
     * This function is structured such that it can be used for pure species or mixtures, with the
     * only difference being the values that are passed to the function (pure values versus mixture values).
     *
     * @param Tr Reduced temperature [unitless]
     * @param FP Polarity correction factor [unitless]
     * @param FQ Quantum correction factor [unitless]
     * @return double
     */
    double low_pressure_nondimensional_viscosity(double Tr, double FP, double FQ);

    /**
     * Returns the non-dimensional high-pressure mixture viscosity in using the Lucas method.
     *
     * Defined by equation 9-6.12.
     *
     * @f[
     *   \eta \xi = Z_2 F_P F_Q
     * @f]
     *
     * This returns the value of η*ξ (by multiplying both sides of 9-6.12 by ξ and
     * simply returning the right-side of the resulting equation).
     *
     * This function is structured such that it can be used for pure species or mixtures, with the
     * only difference being the values that are passed to the function (pure values versus mixture values).
     *
     * @param Tr Reduced temperature [unitless]
     * @param Pr Reduced pressure [unitless]
     * @param FP_low Low-pressure polarity correction factor [unitless]
     * @param FQ_low Low-pressure quantum correction factor [unitless]
     * @param P_vap Vapor pressure [Pa]
     * @param P_crit Critical pressure [Pa]
     * @return double
     */
    double high_pressure_nondimensional_viscosity(double Tr, double Pr, double FP_low,
                                                 double FQ_low, double P_vap, double P_crit);

    /**
     * @brief  Returns the quantum correction term for a species based on Tr
     * and MW, used in viscosity calculation.
     *
     * Calculates quantum correction term of the Lucas method for a species based
     * on the reduced temperature(Tr) and molecular weight(MW), used in viscosity
     * calculation from equation 9-4.19.
     *
     * @f[
     *    F_{Q}^{\text{o}} = 1.22 Q^{0.15} \left( 1 + 0.00385 \left ( \left ( T_r - 12 \right ) ^2 \right ) ^{\frac{1}{MW}} sign(T_r - 12 \right ) \right )
     * @f]
     *
     * @param Q  Species-specific constant
     * @param Tr  Reduced temperature
     * @param MW  Molecular weight
     * @return double
     */
    double FQ_i(double Q, double Tr, double MW);

    /**
     * @brief  Returns the polarity correction term for a species based on reduced
     * temperature, reduced dipole moment, and critical compressibility. Used in
     * the viscosity calculation.
     *
     * Calculates polarity correction term of the Lucas method for a species based
     * on the reduced temperature(Tr) and molecular weight(MW). Equation 9.4.18.
     *
     * @f[
     *  \begin{equation}
     *   F_P^0 =
     *   \begin{cases}
     *       1 & 0 \leq \mu_r < 0.022 \\
     *      1 + 30.55(0.292 - Z_c)^{1.72} & 0.022 \leq \mu_r < 0.075 \\
     *      1 + 30.55(0.292 - Z_c)^{1.72} \times 0.96 + 0.1(T_r - 0.7) & 0.075 \leq \mu_r
     *   \end{cases}
     *  \end{equation}
     *

     * @note The original description in Poling(2001) neglects to mention what happens
     * when the quantity raised to the 1.72 power goes negative. That is an undefined
     * operation that generates real+imaginary numbers. For now, we
     * take the absolute value of the argument.
     *
     * @param mu_r  Species Reduced dipole moment
     * @param Tr  Reduced temperature
     * @param Z_c  Species Critical compressibility
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
 * and thermal conductivity of a high-pressure gas mixture. All methods are described
 * in Poling et al. @cite poling2001 (viscosity in Ch. 9, thermal conductivity in
 * Ch. 10, and diffusion coefficients in Ch. 11).
 *
 * @note: All equations that are cited in this implementation are from the 5th edition
 * of the book "The Properties of Gases and Liquids" by Poling, Prausnitz, and O'Connell.
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

    /**
     * @brief  Returns the estimate of the critical temperature that is given
     * from the thermo object for species i.
     *
     * This method sets the species composition vector to unity for
     * species i and zero for all other species, and then queries the
     * thermo object for the critical temperature. It then resets the
     * composition vector to the original state.
     *
     * @param i  Species index
     * @return double
     */
    double Tcrit_i(size_t i);

    /**
     * @brief  Returns the estimate of the critical pressure that is given
     * from the thermo object for species i.
     *
     * This method sets the species composition vector to unity for
     * species i and zero for all other species, and then queries the
     * thermo object for the critical temperature. It then resets the
     * composition vector to the original state.
     *
     * @param i  Species index
     * @return double
     */
    double Pcrit_i(size_t i);

    /**
     * @brief  Returns the estimate of the critical volume that is given
     * from the thermo object for species i.
     *
     * This method sets the species composition vector to unity for
     * species i and zero for all other species, and then queries the
     * thermo object for the critical temperature. It then resets the
     * composition vector to the original state.
     *
     * @param i  Species index
     * @return double
     */
    double Vcrit_i(size_t i);

    /**
     * @brief  Returns the estimate of the critical compressibility that is given
     * from the thermo object for species i.
     *
     * This method sets the species composition vector to unity for
     * species i and zero for all other species, and then queries the
     * thermo object for the critical temperature. It then resets the
     * composition vector to the original state.
     *
     * @param i  Species index
     * @return double
     */
    double Zcrit_i(size_t i);

    /**
     * Returns the composition-dependent values of the parameters needed for
     * the Chung viscosity model.
     *
     * The equations for the mixing rules defined on page 9.25 for the Chung method's
     * composition dependent parameters. The primary mixing rules are defined below.
     *
     * @f[
     *  \sigma_m^3 = \sum_{i} \sum_{j} y_i y_j \sigma_{ij}^3
     * @f]
     *
     * @f[
     *  T_m^* = \frac{T}{\left( \frac{\epsilon}{k} \right )_m}
     * @f]
     *
     * @f[
     *  \left ( \frac{\epsilon}{k} \right )_m =  \frac{\sum_{i} \sum_{j} y_i y_j \left ( \frac{\epsilon_{ij}}{k} \right ) \sigma_{ij}^3}{\sigma_m^3}
     * @f]
     *
     * @f[
     *  MW_m = \left [ \frac{\sum_{i} \sum_{j} y_i y_j \left ( \frac{\epsilon_{ij}}{k} \right ) \sigma_{ij}^2 MW_{ij}^{\frac{1}{2}}}{\left ( \frac{\epsilon}{k} \right )_m \sigma_m^2} \right ]^2
     * @f]
     *
     * @f[
     *   \omega_m = \frac{\sum_{i} \sum_{j} y_i y_j \omega_{ij} \sigma+{ij}^3}{\sigma_m^3}
     * @f]
     *
     * @f[
     *   \mu_m^4 = \sigma_m^3 \sum_{i} \sum_{j} \left( \frac{y_i y_j \mu_i^2 \mu_j^2}{\sigma_{ij}^3} \right)
     * @f]
     *
     * @f[
     *  \kappa_m = \sum_{i} \sum_{j} y_i y_j \kappa_{ij}
     * @f]
     *
     * The combining rules are defined as:
     *
     * @f[
     *   \sigma_{i} = 0.809 V_{c,i}^{1/3}
     * @f]
     *
     * @f[
     *  \sigma_{ij} =  \xi_{ij} \left( \sigma_{i} \sigma_{j} \right)^{1/2}
     * @f]
     *
     * @f[
     *  \left( \frac{\epsilon_i}{k} \right) = \frac{T_{c,i}}{1.2593}
     * @f]
     *
     * @f[
     *  \left( \frac{\epsilon_{ij}}{k} \right) = \zeta_{ij} \left( \right) ^{\frac{1}{2}}
     * @f]
     *
     * @f[
     *  \omega_{ij} = \frac{\omega_i + \omega_j}{2}
     * @f]
     *
     * @f[
     *  \kappa_{ij} = \left( \kappa_i \kappa_j \right)^{1/2}
     * @f]
     *
     * @f[
     *  MW_{ij} = \frac{2 MW_i MW_j}{MW_i + MW_j}
     * @f]
     *
     * $\xi and \zeta$ are the binary interaction parameters, and are assumed to be unity
     * in this implementation, in keeping with the Chung method.
     *
     * The Chung viscosity correction factor is defined as:
     *
     * @f[
     *  F_{c,m} = 1 - 0.275 \omega_m + 0.059035 \mu_{r,m}^4 + \kappa_m
     * @f]
     *
     * The reduced dipole moment computed using:
     *
     * @f[
     * \mu_{r,m} = \frac{131.3 \mu_m}{\left( V_{c,m} T_{c,m}\right)^{\frac{1}{2}}}
     * @f]
     *
     * @f[
     *  V_{c,m} = \left( \frac{\sigma_m}{0.809} \right)
     * @f]
     *
     * @f[
     * T_{c,m} = 1.2593 \left( \frac{\epsilon}{k} \right)_m
     * @f]
     *
     * In the equations, $T_c$ must be in units of K, $V_c$ must be in units of cm^3/mol,
     * and $\mu$ must be in units of Debye.
     *
     */
    void compute_mixture_parameters(ChungMixtureParameters& params);

    /**
     * Returns the low-pressure mixture viscosity in micropoise using the Chung method.
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
     *
     * @param T  Temperature [K]
     * @param T_star  Reduced temperature [unitless]
     * @param MW  Molecular weight [kg/kmol]
     * @param acentric_factor  Acentric factor [unitless]
     * @param mu_r  Dipole moment [Debye]
     * @param sigma  Lennard-Jones collision diameter [Angstroms]
     * @param kappa  Polar correction factor [unitless]
     * @return double
     */
    double low_pressure_viscosity(double T, double T_star, double MW, double acentric_factor,
                                 double mu_r, double sigma, double kappa);

    /**
     * Returns the high-pressure mixture viscosity in micropoise using the Chung method.
     *
     * Defined by equation 9-6.18.
     *
     * @f[
     *   \eta = \eta^* \frac{36.344 (M*T_c)^(\frac{1}{2})}{V^{\frac{2}{3}}}
     * @f]
     *
     * $T_c$ must be in units of K, MW must be units of kg/kmol, and sigma must be units of Angstroms
     *
     * This function is structured such that it can be used for pure species or mixtures, with the
     * only difference being the values that are passed to the function (pure values versus mixture values).
     *
     * @param T_star  Reduced temperature [unitless]
     * @param MW  Molecular weight [kg/kmol]
     * @param rho  Density [mol/cm^3]
     * @param Vc  Critical volume [cm^3/mol]
     * @param Tc  Critical temperature [K]
     * @param acentric_factor  Acentric factor [unitless]
     * @param mu_r  Dipole moment [Debye]
     * @param kappa  Polar correction factor [unitless]
     * @return double
     */
    double high_pressure_viscosity(double T_star, double MW, double rho, double Vc, double Tc,
                                   double acentric_factor, double mu_r, double kappa);

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
     *
     * @param T  Temperature [K]
     * @param T_star  Reduced temperature [unitless]
     * @param MW  Molecular weight [kg/kmol]
     * @param rho  Density [mol/cm^3]
     * @param Vc  Critical volume [cm^3/mol]
     * @param Tc  Critical temperature [K]
     * @param acentric_factor  Acentric factor [unitless]
     * @param mu_r  Dipole moment [Debye]
     * @param kappa  Polar correction factor [unitless]
     * @return double
     */
    double high_pressure_thermal_conductivity(double T, double T_star, double MW,
                                              double rho, double Cv, double Vc,
                                              double Tc, double sigma,
                                              double acentric_factor, double mu_r,
                                              double kappa);

    /**
     * @brief Returns interpolated value of (DP)_R obtained from the data
     * in Table 2 of the Takahashi 1975 paper, given a value of the reduced
     * pressure (Pr) and reduced temperature (Tr).
     *
     * @param Pr  Reduced pressure
     * @param Tr  Reduced temperature
     */
};

}
#endif
