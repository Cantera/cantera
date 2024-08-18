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
 * The implementation employs a method of corresponding states, using the Takahashi
 * @cite takahashi1975 approach for binary diffusion coefficients (using mixture
 * averaging rules for the mixture properties), the Lucas method for the viscosity, and
 * a method from Ely and Hanley for the thermal conductivity. All methods are described
 * in Poling et al. @cite poling2001 (viscosity in Ch. 9, thermal conductivity in
 * Ch. 10, and diffusion coefficients in Ch. 11).
 *
 * @ingroup tranprops
 */
class HighPressureGasTransport : public MixTransport
{
protected:
    //! default constructor
    HighPressureGasTransport()=default;

public:
    string transportModel() const override {
        return "high-pressure";
    }

    void init(ThermoPhase* thermo, int mode=0, int log_level=0) override;

    /**
     * Returns the mixture high-pressure thermal conductivity in W/m/K
     * using a method of corresponding states by Ely and Hanley.
     *
     * This method for computing the thermal conductivity is described in
     * @cite ely-hanley1983 . There are references to earlier work in
     * @cite ely-hanley1981 in the description of the model for thermal conductivity.
     * The equations referenced in this description primarily are from
     * @cite ely-hanley1983 , and equations that come from @cite ely-hanley1981 are
     * marked as such. The equations referenced here are the same as the ones from
     * the papers from Ely and Hanley.
     *
     * The thermal conductivity is assumed to have two components: a
     * translational/collisional part and an internal part that is related to the
     * transfer of energy to internal degrees of freedom.
     *
     * @f[
     *   \lambda_{mix}(\rho, T) = \lambda^{''}_{mix}(T) + \lambda^{'}_{mix}(\rho, T)
     * @f]
     *
     * The first term on the right-hand side is the internal part and the second term
     * is the translational part. The internal part is only a function of temperature,
     * and the translational part is a function of both temperature and density.
     *
     * For the internal component of the thermal conductivity the following equations
     * are used.
     *
     * @f[
     *  \lambda^{''}_{mix}(T) = \sum_i \sum_j X_i X_j \lambda^{''}_{ij}(T)
     *
     *  \quad \text{(Equation 2)}
     * @f]
     *
     * The mixing rule for combining pure-species internal thermal conductivity
     * components is given by:
     *
     * @f[
     *  (\lambda^{''}_{ij})^{-1} = 2[(\lambda^{''}_{i})^{-1}
     *                              + (\lambda^{''}_{j})^{-1}]
     *
     *  \quad \text{(Equation 3)}
     * @f]
     *
     * The pure species internal thermal conductivity is given by:
     *
     * @f[
     *  \lambda^{''}_{i} = \frac{\eta_i^*}{MW_i}f_{int} \left(C_{p,i}^0 - 2.5R \right)
     *
     *  \quad \text{(Equation 1)}
     * @f]
     *
     * Where @f$ \eta_i^* @f$ is the referred to as the dilute gas viscosity and is the
     * component that is temperature dependent described in @cite ely-hanley1981.
     * @f$ MW_i @f$ is the molecular weight of species i, @f$ f_int @f$ is a constant
     * and is 1.32, @f$ C_{p,i}^0 @f$ is the ideal gas heat capacity of species i,
     * and R is the gas constant (J/kmol/K).
     *
     * For the translational component of the thermal conductivity the following
     * equations are used.
     *
     * @f[
     *   \lambda^{'}_{mix}(\rho, T) = \lambda^{'}_{0}(\rho_0, T_0) F_{\lambda}
     *
     *   \quad \text{(Equation 5)}
     * @f]
     *
     * Where @f$ \lambda^{'}_{0}(\rho_0, T_0) @f$ is the internal component of the
     * thermal conductivity of a reference fluid that is at a reference temperature
     * and density, @f$ T_0 @f$ and @f$ \rho_0 @f$. These reference conditions are not
     * the same as the conditions that the mixture is at (@f$ T @f$ and @f$ \rho @f$).
     * The subscript 0 refers to the reference fluid.
     *
     * @f[
     *  F_{\lambda} = \left( \frac{MW_0}{MW_m} \right)^{1/2} f_{m,0}^{1/2}
     *                h_{m,0}^{-2/3}
     *
     *  \quad \text{(Equation 6)}
     * @f]
     *
     * Where @f$ MW_0 @f$ is the molecular weight of the reference fluid, @f$ MW_m @f$
     * is the molecular weight of the mixture, @f$ f_{m,0} @f$ and @f$ h_{m,0} @f$ are
     * called reducing ratios. The @f$ h_{m,0} @f$ quantity is defined as:
     *
     * @f[
     *   h_{m,0} = \sum_i \sum_j X_i X_j h_{ij}
     *
     *   \quad \text{(Equation 9)}
     * @f]
     *
     * with,
     *
     * @f[
     *   h_{ij} = \frac{1}{8} (h_i^{1/3} + h_j^{1/3})^3 (1 - l_{ij})
     *
     *   \quad \text{(Equation 11)}
     * @f]
     *
     * The variable @f$ l_{ij} @f$ is a binary interaction parameter and is assumed to
     * be zero as is done by Ely and Hanley.
     *
     * @f[
     *   h_i = \frac{V_c^i}{V_c^0} \phi_i(T_R^i, V_R^i, \omega_i)
     *
     *   \quad \text{(Equation 13)}
     * @f]
     *
     * Where @f$ \phi_i @f$ is a shape factor of Leach and Leland, @f$ T_R^i @f$ is
     * the reduced temperature of species i, @f$ V_R^i @f$ is the reduced volume of
     * species i, and @f$ \omega_i @f$ is the acentric factor of species i The shape
     * factor @f$ \phi_i @f$ is defined in @cite ely-hanley1981 as follows:
     *
     * @f[
     *   \phi_i(T_R^i, V_R^i, \omega_i) = \frac{Z_0^c}{Z_i^c} [1 + (\omega_i -
     *                                    \omega_0)G(T_R^i, V_R^i)]
     *
     *   \quad \text{(Equation 12 of @cite ely-hanley1981)}
     * @f]
     *
     * Where @f$ Z_0^c @f$ is the critical compressibility of the reference fluid,
     * @f$ Z_i^c @f$ is the critical compressibility of species i, @f$ \omega_0 @f$
     * is the acentric factor of the reference fluid, and @f$ G(T_R^i, V_R^i) @f$ is a
     * function of the reduced temperature and reduced volume of species i. The
     * function @f$ G(T_R^i, V_R^i) @f$ is defined as:
     *
     * @f[
     *   G(T_R^i, V_R^i) = a_2(V_i^+ + b_2) + c_2(V_i^+ + d_2)ln(T_R^i)
     *
     *   \quad \text{(Equation 14 of @cite ely-hanley1981)}
     * @f]
     *
     * Where @f$ T_i^+ @f$ and @f$ V_i^+ @f$ are limited values of the reduced
     * temperature and pressure and are defined as:
     *
     * @f[
     *   T_i^+ = min(2, max(T_R^i, 0.5))
     *
     *   \quad \text{(Equation 15 of @cite ely-hanley1981)}
     * @f]
     *
     * and
     *
     * @f[
     *  V_i^+ = min(2, max(V_R^i, 0.5))
     *
     *  \quad \text{(Equation 16 of @cite ely-hanley1981)}
     * @f]
     *
     * The values of @f$ a_2 @f$, @f$ b_2 @f$, @f$ c_2 @f$, and @f$ d_2 @f$ are
     * given in Table II of @cite ely-hanley1981. They are:
     *
     * @f[
     *   a_2 = 0.394901, b_2 = -1.023545, c_2 = -0.932813, d_2 = -0.754639
     * @f]
     *
     * With the definitions above, the value of @f$ h_{m,0} @f$ can be computed. This
     * leaves the other reducing ratio, @f$ f_{m,0} @f$, to be defined. The value of
     * that reducing ratio is defined by the following equation.
     *
     * @f[
     *   f_{m,0} = \frac{1}{h_{m,0}} \sum_i \sum_j X_i X_j f_{ij} h_{ij}
     *
     *   \quad \text{(Equation 8)}
     * @f]
     *
     *
     * The value of @f$ h_{ij} @f$ is the same as the one defined earlier. The
     * combining rule for @f$ f_{ij} @f$ is given by:
     *
     * @f[
     *   f_{ij} = (f_i f_j)^{1/2} (1-\kappa_{ij})
     *
     *   \quad \text{(Equation 10)}
     * @f]
     *
     * @f$ \kappa_{ij} @f$ is binary interaction parameters and is assumed to be zero
     * as was done in the work of Ely and Hanley. The species-specific value of
     * @f$ f_i @f$ is defined as:
     *
     * @f[
     *   f_i = \frac{T_c^i}{T_c^0} \theta_i(T_R^i, V_R^i, \omega_i)
     *
     *   \quad \text{(Equation 12 @cite ely-hanley1981)}
     * @f]
     *
     * Where @f$ \theta_i @f$ is a shape factor of Leach and Leland, @f$ T_c^i @f$ is
     * the critical temperature of species i, and @f$ T_c^0 @f$ is the critical
     * temperature of the reference fluid. The shape factor @f$ \theta_i @f$ is defined
     * in @cite ely-hanley1981 as:
     *
     * @f[
     *   \theta_i(T_R^i, V_R^i, \omega_i) = 1 + (\omega_i - \omega_0)F(T_R^i, V_R^i)]
     *
     *   \quad \text{(Equation 11 of @cite ely-hanley1981)}
     * @f]
     *
     * Where @f$ Z_0^c @f$ is the critical compressibility of the reference fluid,
     * @f$ Z_i^c @f$ is the critical compressibility of species i, @f$ \omega_0 @f$ is
     * the acentric factor of the reference fluid, and @f$ G(T_R^i, V_R^i) @f$ is a
     * function of the reduced temperature and reduced volume of species i. The
     * function @f$ G(T_R^i, V_R^i) @f$ is defined as:
     *
     * @f[
     *   F(T_R^i, V_R^i) = a_1 + b_1 ln(T_R^i) + (c_1 + d_1/T_R^i) (V_R^i - 0.5)
     *
     *   \quad \text{(Equation 13 of @cite ely-hanley1981)}
     * @f]
     *
     * Where @f$ T_i^+ @f$ and @f$ V_i^+ @f$ are limited values of the reduced
     * temperature and pressure and have the same definition as was shown earlier.
     *
     * The values of @f$ a_1 @f$, @f$ b_1 @f$, @f$ c_1 @f$, and @f$ d_1 @f$ are
     * given in Table II of @cite ely-hanley1981. They are:
     *
     * @f[
     *   a_1 = 0.090569, b_1 = -0.862762, c_1 = 0.316636, d_1 = -0.465684
     * @f]
     *
     * With the definitions above, the value of @f$ h_{m,0} @f$ can be computed from
     * the equation that defined the value of @f$ f_{m,0} h_{m,0} @f$. The value of
     * @f$ h_{m,0} @f$ must be computed first and then the value of @f$ f_{m,0} @f$
     * can be computed.
     *
     * The value of @f$ MW_m @f$ is the molecular weight of the mixture and is computed
     * using the following equation. This expression is somewhat complex, but the
     * method to obtain the value of @f$ MW_m @f$ is to compute the right-hand-side
     * first with the definitions of @f$ MW_{ij} @f$, @f$ f_{ij} @f$, and
     * @f$ h_{ij} @f$ . Then once the right-hand-side is known, along with the
     * previously computed values of @f$ f_{m,0} @f$ and @f$ h_{m,0} @f$,we can use
     * simple algebra to isolate @f$ MW_m @f$.
     *
     * @f[
     *   MW_m^{1/2} f_{m,0}^{1/2} h_{m,0}^{-4/3} = \sum_i \sum_j X_i X_j MW_{ij}^{-1/2}
     *                                             f_{ij}^{1/2} h_{ij}^{-4/3}
     *
     *   \quad \text{(Equation 14)}
     * @f]
     *
     * where,
     *
     * @f[
     *   MW_{ij} = \frac{1}{2} (MW_i^{-1} + MW_j^{-1})
     *
     *   \quad \text{(Equation 15)}
     * @f]
     *
     * For clarity, if we assume the right-hand-side of the molecular weight mixture
     * equation is A, then the mixture molecular weight is given by:
     *
     * @f[
     *   MW_m = A^{-2} f_{m,0} h_{m,0}^{-8/3}
     * @f]
     *
     * All of the above descriptions allow for the calculation of @f$ F_{\lambda} @f$
     * in the expression for the mixture value of the translational component of the
     * thermal conductivity. The reference fluid value of the translational component
     * of the thermal conductivity is still needed. The first thing that needs to be
     * done is the obtain the temperature and density of the reference fluid. The
     * following relations are used to compute the reference fluid temperature and
     * density.
     *
     * @f[
     *   \rho_0 = \rho h_{m,0}
     *
     *   \quad \text{(Equation 7)}
     * @f]
     *
     * and
     *
     * @f[
     *   T_0 = \frac{T}{f_{m,0}}
     *   \quad \text{(Equation 7)}
     * @f]
     *
     * The reference fluid translational component of thermal conductivity can now be
     * computed using the following equation.
     *
     * @f[
     *   \lambda^{'}_{0}(\rho_0, T_0) = \frac{15R}{4MW_0} \eta_0^* +
     *                                  \lambda_0^{(1)}\rho_0 +
     *                                  \Delta\lambda_0(\rho_0, T_0)
     *
     *   \quad \text{(Equation 18)}
     * @f]
     *
     * Where @f$ \eta_0^* @f$ is the dilute gas viscosity, @f$ \lambda_0^{(1)} @f$ is
     * the first density correction to the thermal conductivity, and
     * @f$ \Delta\lambda_0 @f$ is the high-density contribution. Ely and Hanley provide
     * a correlation equation to solve for this quantity for the methane reference
     * fluid. The details of the correlation and the parameters are shown in Table I
     * of @cite ely-hanley1983.
     *
     * For completeness, the relations are reproduced here.
     *
     * @f[
     *   \lambda_0^*(T_0) = \frac{15R}{4MW_0} \eta_0^*
     * @f]
     *
     * @f[
     *   \eta_0^* = \sum_{n=1}^{9} C_n T_0^{(n-4)/3}
     * @f]
     *
     * @f[
     *   \lambda_0^{(1)} = b_1 + b_2[b_3 - ln(T_0 / b_4)]^2
     * @f]
     *
     * @f[
     *  \Delta\lambda_0 = exp[a_1 + a_2/T_0] (exp[(a_3 + a_4/T_0^{3/4})\rho_0^{0.1}
     *                    + (\rho_0 / \rho_0,c - 1)\rho_0^{0.5}
     *                      (a_5 + a_6/T_0 + a_7/T_0^2) ] - 1)
     * @f]
     *
     * The constants a, b, and C are not related to any ones used in the previous
     * equations, their values are defined in Table I of @cite ely-hanley1983.
     *
     * @note The correlations for the reference fluid viscosity return a value of
     *       micro-grams/cm/s and the parts of the thermal conductivity that utilize
     *       the correlation fit parameters give values in units of mW/m/K. Appropriate
     *       conversions were made in the relevant functions.
     *       @see elyHanleyDiluteViscosity()
     *       @see elyHanleyReferenceThermalConductivity()
     */
    double thermalConductivity() override;

     /**
     * Returns the mixture high-pressure viscosity in Pa*s using the Lucas method.
     *
     * This uses the approach described in chapter 9-7. In this method, the mixture
     * viscosity at high pressure is computed using the pure-fluid high pressure
     * relation of Lucas with the difference being that mixture values of the
     * model parameters are used. These mixture values are computed using the mixing
     * rules described in @see computeMixtureParameters() .
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

    /**
     * Computes and stores the estimate of the critical properties for each species.
     *
     * This method sets the species composition vector to unity for species i and zero
     * for all other species, and then queries the thermo object for the critical
     * properties and stores them. It then resets the composition vector to the original
     * state. This method only needs to be called once, as the critical properties for
     * the pure species do not change.
     *
     */
    void initializeCriticalProperties();

    /**
     * Returns the stored value of the critical temperature.
     *
     * @param i  Species index
     */
    double Tcrit_i(size_t i);

    /**
     * Returns the stored value of the critical pressure.
     *
     * @param i  Species index
     */
    double Pcrit_i(size_t i);

    /**
     * Returns the stored value of the critical volume.
     *
     * @param i  Species index
     */
    double Vcrit_i(size_t i);

    /**
     * Returns the stored value of the critical compressibility.
     *
     * @param i  Species index
     */
    double Zcrit_i(size_t i);

    /**
     * Returns the viscosity of for the reference fluid of methane for low pressures.
     * Units are Pa*s
     *
     * Computes the temperature dependent (referred to as the dilute viscosity in the
     * reference) component only (eta_0) from the expression in Table III in
     * @cite ely-hanley1981
     *
     * @param T0
     */
    double elyHanleyDiluteViscosity(double T0);

    /**
     * Returns the thermal conductivity of for the reference fluid of methane
     * in units of W/m/K.
     *
     * Computes the temperature and density dependent reference fluid thermal conductivity
     * from the expression in Table I in @cite ely-hanley1983 .
     *
     * @param rho0
     * @param T0
     */
    double elyHanleyReferenceThermalConductivity(double rho0, double T0);

    /**
     * Computes the composition-dependent values of parameters that are needed for the
     * Lucas viscosity model.
     *
     * The equations for the mixing rules defined on page 9.23 of @cite poling2001 for
     * the Lucas method's composition dependent parameters. The primary mixing rules
     * are defined below, and the reduced properties are just the properties divided
     * by the pseudo-critical mixture properties defined below.
     *
     * @note Equation numbers are from @cite poling2001
     *
     * @f[
     *  T_{\text{c,m}} = \sum_i X_i T_{\text{c,i}}
     *
     *  \quad \text{( Equation 9-5.18)}
     * @f]
     *
     * Where @f$ T_{\text{c,i}} @f$ is the critical temperature of species i,
     * and @f$ X_i @f$ is the mole fraction of species i.
     *
     * @f[
     *  P_{\text{c,m}} = R T_{\text{c,m}} \frac{\sum_i X_i Z_{\text{c,i}}}{\sum_i X_i V_{\text{c,i}}}
     *
     *  \quad \text{(Equation 9-5.19)}
     * @f]
     *
     * Where @f$ Z_{\text{c,i}} @f$ is the critical compressibility of species i, and
     * @f$ V_{\text{c,i}} @f$ is the critical volume of species i.
     *
     * @f[
     *  M_m = \sum X_i M_i
     *
     *  \quad \text{(Equation 9-5.20)}
     * @f]
     *
     * Where @f$ M_i @f$ is the molecular weight of species i.
     *
     * @f[
     *  F_{P,m}^{\text{o}} = \sum X_i F_{P,i}^{\text{o}}
     *
     *  \quad \text{(Equation 9-5.21)}
     * @f]
     *
     * Where @f$ F_{P,i}^{\text{o}} @f$ is the low-pressure polarity correction
     * factor of species i from equation 9-4.18.
     *
     * @f[
     *  F_{Q,m}^{\text{o}} = \left ( \sum X_i F_{Q,i}^{\text{o}} \right ) A
     *
     *  \quad \text{(Equation 9-5.22)}
     * @f]
     *
     * Where @f$ F_{Q,i}^{\text{o}} @f$ is the low-pressure quantum correction factor
     * of species i from equation 9-4.19, and A is defined below.
     *
     * @f[
     *   A = 1 - 0.01 \left ( \frac{M_H}{M_L} \right )^{0.87}
     *
     *   \quad \text{(Equation 9-5.23)}
     * @f]
     *
     * For @f$ \frac{M_H}{M_L} > 9 @f$ and @f$ 0.05 < X_H < 0.7 @f$, otherwise A = 1.
     * In the above equation, $M_H$ and $M_L$ are the molecular weights of the
     * heaviest and lightest components in the mixture, and @f$ X_H @f$ is the mole
     * fraction of the heaviest component.
     *
     * While it isn't returned as a parameter, the species-specific reduced dipole
     * moment (@f$ \mu_r @f$) is used to compute the mixture polarity correction factor. It is defined
     * as:
     *
     * @f[
     *   \mu_r = 52.46 \frac{\mu^2 P_{\text{c,i}}}{T_{\text{c,i}}}
     *
     *   \quad \text{(Equation 9-4.17)}
     * @f]
     *
     */
    void computeMixtureParameters();

    /**
     * Returns the non-dimensional low-pressure mixture viscosity in using the Lucas
     * method.
     *
     * @f[
     *   \eta \xi = F_P^{\text{o}} F_Q^{\text{o}} [0.807 T_r^{0.618}
     *                                             - 0.357 e^{-0.449 T_r}
     *                                             + 0.340e^{-4.058 T_r} + 0.018]
     *
     *   \quad \text{(Equation 9-4.16)}
     * @f]
     *
     * This function is structured such that it can be used for pure species or
     * mixtures, with the only difference being the values that are passed to the
     * function (pure values versus mixture values).
     *
     * For the definition of the mixture rules, see @see computeMixtureParameters() .
     *
     * @param Tr Reduced temperature [unitless]
     * @param FP Polarity correction factor [unitless]
     * @param FQ Quantum correction factor [unitless]
     */
    double lowPressureNondimensionalViscosity(double Tr, double FP, double FQ);

    /**
     * Returns the non-dimensional high-pressure mixture viscosity in using the Lucas
     * method.
     *
     * @f[
     *   \eta \xi = Z_2 F_P F_Q
     *
     *   \quad \text{(Equation 9-6.12)}
     * @f]
     *
     * This returns the value of η*ξ (by multiplying both sides of 9-6.12 by ξ and
     * returning the right-side of the resulting equation).
     *
     * This function is structured such that it can be used for pure species or
     * mixtures, with the only difference being the values that are passed to the
     * function (pure values versus mixture values).
     *
     * For the definition of the mixture rules, see @see computeMixtureParameters() .
     *
     * @param Tr Reduced temperature [unitless]
     * @param Pr Reduced pressure [unitless]
     * @param FP_low Low-pressure polarity correction factor [unitless]
     * @param FQ_low Low-pressure quantum correction factor [unitless]
     * @param P_vap Vapor pressure [Pa]
     * @param P_crit Critical pressure [Pa]
     */
    double highPressureNondimensionalViscosity(double Tr, double Pr, double FP_low,
                                              double FQ_low, double P_vap,
                                              double P_crit);

    /**
     * Calculates quantum correction term of the Lucas method for a species based
     * on the reduced temperature(Tr) and molecular weight(MW), used in viscosity
     * calculation.
     *
     * @f[
     *  F_{Q}^{\text{o}} = 1.22 Q^{0.15} {1 + 0.00385[ (T_r - 12)^2 ]^{\frac{1}{MW}}
     *                     \text{sign} (T_r - 12 )}
     *
     *  \quad \text{(Equation 9-4.19)}
     * @f]
     *
     * @param Q  Species-specific constant
     * @param Tr  Reduced temperature [unitless]
     * @param MW  Molecular weight [kg/kmol]
     */
    double quantumCorrectionFactor(double Q, double Tr, double MW);

    /**
     * Returns the polarity correction term for a species based on reduced temperature,
     * reduced dipole moment, and critical compressibility. Used in the calculation of
     * viscosity.
     *
     * Calculates polarity correction term of the Lucas method for a species based
     * on the reduced temperature(Tr) and molecular weight(MW). Equation 9.4.18.
     *
     * @f[
     *  \begin{equation}
     *   F_P^0 =
     *   \begin{cases}
     *      1 & 0 \leq \mu_r < 0.022 \\
     *      1 + 30.55(0.292 - Z_c)^{1.72} & 0.022 \leq \mu_r < 0.075 \\
     *      1 + 30.55(0.292 - Z_c)^{1.72} \times 0.96 + 0.1(T_r - 0.7) & 0.075 \leq \mu_r
     *   \end{cases}
     *  \end{equation}
     * @f]
     *
     * @note The original description in Poling(2001) neglects to mention what happens
     * when the quantity raised to the 1.72 power goes negative. That is an undefined
     * operation that generates real and imaginary numbers. For now, only positive
     * values are allowed.
     *
     * @param mu_r  Species Reduced dipole moment
     * @param Tr  Reduced temperature
     * @param Z_c  Species Critical compressibility
     */
    double polarityCorrectionFactor(double mu_r, double Tr, double Z_c);


    private:
    //! Critical properties of each species.
    std::vector<double> m_Tcrit;
    std::vector<double> m_Pcrit;
    std::vector<double> m_Vcrit;
    std::vector<double> m_Zcrit;

    // Reference fluid properties, for methane (used by the thermalConductivity()
    // method)
    const double ref_MW = 16.04; // kg/kmol
    const double ref_Tc = 190.4; // K
    const double ref_Vc = 0.0986; // m^3/kmol
    const double ref_Zc = 0.288; // unitless
    const double ref_rhoc = 0.1628; // g/cm^3
    const double ref_acentric_factor = 0.011; // unitless

    // Parameters that are needed to calculate the viscosity using the Lucas method.
    double m_FQ_mix_o; // Quantum correction factor
    double m_FP_mix_o; // Polarity correction factor
    double m_Tr_mix; // Reduced temperature
    double m_Pr_mix; // Reduced pressure
    double m_Pc_mix; // Critical pressure
    double m_Tc_mix; // Critical temperature
    double m_MW_mix; // Molecular weight
    double m_P_vap_mix; // Vapor pressure
};

/**
 * Transport properties for high pressure gas mixtures using the Chung method for
 * viscosity and thermal conductivity.
 *
 * The implementation employs a method of corresponding states, using the Takahashi
 * @cite takahashi1975 approach for binary diffusion coefficients (using mixture
 * averaging rules for the mixture properties), and the Chung method for the viscosity
 * and thermal conductivity of a high-pressure gas mixture. All methods are described
 * in Poling et al. @cite poling2001 (viscosity in Ch. 9, thermal conductivity in
 * Ch. 10, and diffusion coefficients in Ch. 11).
 *
 * @note All equations that are cited in this implementation are from the 5th edition
 * of the book "The Properties of Gases and Liquids" by Poling, Prausnitz, and O'Connell.
 *
 * @ingroup tranprops
 */
class ChungHighPressureGasTransport : public MixTransport
{
protected:
    //! default constructor
    ChungHighPressureGasTransport()=default;

public:
    string transportModel() const override {
        return "high-pressure-chung";
    }

    void init(ThermoPhase* thermo, int mode=0, int log_level=0) override;

    /**
     * Returns the high-pressure mixture viscosity in Pa*s using the Chung method.
     *
     * Based on the high-pressure gas mixture viscosity model of Chung described in
     * chapter 9-7 of Poling. This method uses the pure species high-pressure viscosity
     * relation of Chung with the difference being that mixture values of the model
     * are computed using a set of mixing rules given by Chung
     * @see computeMixtureParameters() . The mixing rules are defined in section
     * 9-5 of @cite poling2001.
     *
     * Because this method is using the high-pressure viscosity model with mixture
     * parameters, @see highPressureViscosity() for details on the model.
     *
     */
    double viscosity() override;

    /**
    * Calculates the high-pressure mixture thermal conductivity using the Chung method.
    *
    * This method obtains appropriate mixture values of the parameters needed for the
    * Chung model and then calls the highPressureThermalConductivity() method to
    * obtain the mixture thermal conductivity.
    *
    * The mixture values of the pseudo-critical temperature and other model parameters
    * are calculated using the Chung mixing rules defined on page 9.25.
    * @see computeMixtureParameters() .
    *
    * The mixture value of the specific heat is computed using equation 10-6.6, which
    * is the mole fraction weighted sum of the pure species specific heats. This value
    * is not directly computed by the computeMixtureParameters() method.
    *
    * @f[
    *   C_{v,m} = \sum_i X_i C_{v,i}
    *
    *   \quad \text{(Equation 10-6.6)}
    * @f]
    *
    * Where @f$ C_{v,i} @f$ is the specific heat of species i, and @f$ X_i @f$ is the
    * mole fraction of species i.
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
     * Computes and stores the estimate of the critical properties for each species.
     *
     * This method sets the species composition vector to unity for species i and zero
     * for all other species, and then queries the thermo object for the critical
     * properties and stores them. It then resets the composition vector to the original
     * state. This method only needs to be called once, as the critical properties for
     * the pure species do not change.
     *
     */
    void initializeCriticalProperties();

     /**
     * Computes and stores pure-fluid specific properties each species.
     *
     *
     */
    void initializePureFluidProperties();

    /**
     * Returns the stored value of the critical temperature.
     *
     * @param i  Species index
     */
    double Tcrit_i(size_t i);

    /**
     * Returns the stored value of the critical pressure.
     *
     * @param i  Species index
     */
    double Pcrit_i(size_t i);

    /**
     * Returns the stored value of the critical volume.
     *
     * @param i  Species index
     */
    double Vcrit_i(size_t i);

    /**
     * Returns the stored value of the critical compressibility.
     *
     * @param i  Species index
     */
    double Zcrit_i(size_t i);

    /**
     * Computes the composition-dependent values of the parameters needed for
     * the Chung viscosity model.
     *
     * The equations for the mixing rules defined on page 9.25 for the Chung method's
     * composition dependent parameters. The primary mixing rules are defined below.
     *
     * @f[
     *  \sigma_m^3 = \sum_{i} \sum_{j} X_i X_j \sigma_{ij}^3
     *
     *  \quad \text{(Equation 9-5.25)}
     * @f]
     *
     * Where @f$ \sigma_{ij} @f$ is the molecular diameter
     *
     * @f[
     *  T_m^* = \frac{T}{\left( \frac{\epsilon}{k} \right )_m}
     *
     *  \quad \text{(Equation 9-5.26)}
     * @f]
     *
     * Where @f$ k @f$ is the Boltzmann constant and @f$ \epsilon @f$ is the minimum
     * of the pair-potential energy. In these equations, we do not need to worry about
     * what the values of @f$ \epsilon @f$ and @f$ k @f$ are.
     *
     * @f[
     *  \left( \frac{\epsilon}{k} \right)_m = \frac{\sum_{i} \sum_{j} X_i X_j
     *                                        \left( \frac{\epsilon_{ij}}{k} \right)
     *                                        \sigma_{ij}^3}{\sigma_m^3}
     *
     *  \quad \text{(Equation 9-5.27)}
     * @f]
     *
     * @f[
     *  MW_m = \left[ \frac{\sum_{i} \sum_{j} X_i X_j \left( \frac{\epsilon_{ij}}{k} \right)
     *         \sigma_{ij}^2 MW_{ij}^{\frac{1}{2}}}{\left( \frac{\epsilon}{k} \right)_m
     *         \sigma_m^2} \right]^2
     *
     *  \quad \text{(Equation 9-5.28)}
     * @f]
     *
     * Where MW is the molecular weight.
     *
     * @f[
     *   \omega_m = \frac{\sum_{i} \sum_{j} X_i X_j \omega_{ij} \sigma{ij}^3}{\sigma_m^3}
     *
     *   \quad \text{(Equation 9-5.29)}
     * @f]
     *
     * Where @f$ \omega @f$ is the acentric factor.
     *
     * @f[
     *   \mu_m^4 = \sigma_m^3 \sum_{i} \sum_{j} \left( \frac{X_i X_j \mu_i^2 \mu_j^2}
     *             {\sigma_{ij}^3} \right)
     *
     *   \quad \text{(Equation 9-5.30)}
     * @f]
     *
     * Where @f$ \mu @f$ is the dipole moment.
     *
     * @f[
     *  \kappa_m = \sum_{i} \sum_{j} X_i X_j \kappa_{ij}
     *
     *  \quad \text{(Equation 9-5.31)}
     * @f]
     *
     * Where @f$ \kappa @f$ is the association factor, which is used for highly polar
     * molecules. In this work, the value is assumed to be zero for all species.
     *
     * The combining rules for species-species values (subscripted with ij in the above
     * equations) are defined below.
     *
     * @f[
     *   \sigma_{i} = 0.809 V_{c,i}^{1/3}
     *
     *   \quad \text{(Equation 9-5.32)}
     * @f]
     *
     * Where @f$ V_{c,i} @f$ is the critical volume of species i.
     *
     * @f[
     *  \sigma_{ij} =  \xi_{ij} \left( \sigma_{i} \sigma_{j} \right)^{1/2}
     *
     *  \quad \text{(Equation 9-5.33)}
     * @f]
     *
     * Where @f$ \xi @f$ is a binary interaction parameter.
     *
     * @f[
     *  \left( \frac{\epsilon_i}{k} \right) = \frac{T_{c,i}}{1.2593}
     *
     *  \quad \text{(Equation 9-5.34)}
     * @f]
     *
     * @f[
     *  \frac{\epsilon_{ij}}{k} = \zeta_{ij} \left( \right)^{\frac{1}{2}}
     *
     *  \quad \text{(Equation 9-5.35)}
     * @f]
     *
     * Where @f$ \zeta @f$ is a binary interaction parameter.
     *
     * @f[
     *  \omega_{ij} = \frac{\omega_i + \omega_j}{2}
     *
     *  \quad \text{(Equation 9-5.37)}
     * @f]
     *
     * Where @f$ \omega @f$ is the acentric factor.
     *
     * @f[
     *  \kappa_{ij} = \left( \kappa_i \kappa_j \right)^{1/2}
     *
     *  \quad \text{(Equation 9-5.39)}
     * @f]
     *
     * Where @f$ \kappa @f$ is the association factor.
     *
     * @f[
     *  MW_{ij} = \frac{2 MW_i MW_j}{MW_i + MW_j}
     *
     *  \quad \text{(Equation 9-5.40)}
     * @f]
     *
     * @f$ \xi @f$ and @f$ \zeta @f$ are the binary interaction parameters, and are
     * assumed to be unity in this implementation, in keeping with the Chung method.
     *
     * The Chung viscosity correction factor is defined as:
     *
     * @f[
     *  F_{c,m} = 1 - 0.275 \omega_m + 0.059035 \mu_{r,m}^4 + \kappa_m
     *
     *  \quad \text{(Equation 9-5.41)}
     * @f]
     *
     * Where @f$ \omega_m @f$ is the mixture acentric factor, @f$ \mu_{r,m} @f$ is the
     * mixture reduced dipole moment, and @f$ \kappa_m @f$ is the mixture association
     * factor.
     *
     * The mixture reduced dipole moment is computed using:
     *
     * @f[
     *   \mu_{r,m} = \frac{131.3 \mu_m}{( V_{c,m} T_{c,m})^{\frac{1}{2}}}
     *
     *   \quad \text{(Equation 9-5.42)}
     * @f]
     *
     * Where @f$ V_{c,m} @f$  and @f$ T_{c,m} @f$ are computed using the following
     * equations.
     *
     * @f[
     *   V_{c,m} = \left( \frac{\sigma_m}{0.809} \right)^3
     *
     *   \quad \text{(Equation 9-5.43)}
     * @f]
     *
     * @f[
     *   T_{c,m} = 1.2593 \left( \frac{\epsilon}{k} \right)_m
     *
     *   \quad \text{(Equation 9-5.44)}
     * @f]
     *
     * In the equations, @f$ T_c @f$ must be in units of K, @f$ V_c @f$ must be in
     * units of cm^3/mol, and @f$ \mu @f$ must be in units of Debye.
     */
    void computeMixtureParameters();

    /**
     * Returns the low-pressure mixture viscosity in Pa*s using the Chung method.
     *
     * @f[
     *   \eta = 26.69 F_c \frac{(MW*T)^(\frac{1}{2})}{\sigma^2 \Omega}
     *
     *   \quad \text{(Equation 9-4.10)}
     * @f]
     *
     * T must be in units of K, MW must be units of kg/kmol, and @f$ \sigma @f$ must
     * be in units of Angstroms. The viscosity is computed in micropoise, but the
     * return value is in standard SI units (Pa*s).
     *
     * This function is structured such that it can be used for pure species or
     * mixtures, with the only difference being the values that are passed to the
     * function (pure values versus mixture values).
     *
     * @param T  Temperature [K]
     * @param T_star  Reduced temperature [unitless]
     * @param MW  Molecular weight [kg/kmol]
     * @param acentric_factor  Acentric factor [unitless]
     * @param mu_r  Dipole moment [Debye]
     * @param sigma  Lennard-Jones collision diameter [Angstroms]
     * @param kappa  Polar correction factor [unitless]
     */
    double lowPressureViscosity(double T, double T_star, double MW, double acentric_factor,
                                double mu_r, double sigma, double kappa);

    /**
     * Returns the high-pressure mixture viscosity in micropoise using the Chung
     * method.
     *
     * @f[
     *   \eta = \eta^* \frac{36.344 (M*T_c)^(\frac{1}{2})}{V^{\frac{2}{3}}}
     *
     *   \quad \text{(Equation 9-6.18)}
     * @f]
     *
     * where:
     *
     * @f[
     *   \begin{align*}
     *       \eta &= \text{viscosity, \mu P)} \\
     *       M &= \text{molecular wight, kg/kmol} \\
     *       T_c &= \text{critical temperature, K} \\
     *       V_c &= \text{critical molar volume, cm}^3 / \text{mol} \\
     *   \end{align*}
     * @f]
     *
     * and,
     *
     * @f[
     *  \eta^* = \frac{(T^*)^{\frac{1}{2}}}{\Omega_v} {F_c[(G_2)^{-1} + E_6 y]} + \eta^{**}
     *
     *  \quad \text{(Equation 9-6.19)}
     * @f]
     *
     * The values of @f$ T^* @f$ and @f$ F_c @f$ are defined as follows.
     *
     * @f[
     *   T^* = 1.2593 T_r
     *
     *  \quad \text{(Equation 9-4.9)}
     * @f]
     *
     * @f[
     *   F_c = 1 - 0.275 \omega + 0.059035 \mu_r^4 + \kappa
     *
     *   \quad \text{(Equation 9-4.11)}
     * @f]
     *
     * The value of @f$ \Omega_v @f$ is the viscosity collision integral evaluated at
     * the non-dimensional reduced temperature @f$ T^* @f$. @see neufeldCollisionIntegral() .
     *
     * This function is structured such that it can be used for pure species or
     * mixtures, with the only difference being the values that are passed to the
     * function (pure values versus mixture values).
     *
     * @param T_star  Reduced temperature [unitless]
     * @param MW  Molecular weight [kg/kmol]
     * @param rho  Density [mol/cm^3]
     * @param Vc  Critical volume [cm^3/mol]
     * @param Tc  Critical temperature [K]
     * @param acentric_factor  Acentric factor [unitless]
     * @param mu_r  Dipole moment [Debye]
     * @param kappa  Polar correction factor [unitless]
     */
    double highPressureViscosity(double T_star, double MW, double rho, double Vc,
        double Tc, double acentric_factor, double mu_r, double kappa);

    /**
     * Computes the high-pressure thermal conductivity using the Chung method in units
     * of W/m/K
     *
     * The Chung method for computing high-pressure thermal conductivity is described
     * on page 10.23 of @cite poling2001 .
     *
     *  @f[
     *    \lambda = \frac{31.2 \eta^0 \Psi}{M'} \left( G_2^{-1} + B_6 y \right)
     *              + q B_7 y^2 T_r^{1/2} G_2
     *
     *    \quad \text{(Equation 10-5.5)}
     * @f]
     *
     * where:
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
     *       V_c &= \text{critical molar volume, cm}^3/\text{mol} \\
     *   \end{align*}
     * @f]
     *
     * where,
     *
     * @f[
     *  y = \frac{V_c}{6V}
     *
     *  \quad \text{(Equation 10-5.6)}
     * @f]
     *
     * V is the molar volume of the fluid in cm^3/mol.
     *
     * @f[
     *  G_1 = \frac{1 - 0.5y}{(1-y)^3}
     *
     *  \quad \text{(Equation 10-5.7)}
     * @f]
     *
     * @f[
     *  G_2 = \frac{(B_1 / y)[1 - \text{exp}(-B_4 y)] + B_2 G_1 \text{exp}(B_5 y) + B_3 G_1}{B_1 B_4 + B_2 + B_3}
     *
     *  \quad \text{(Equation 10-5.8)}
     * @f]
     *
     * The coefficients @f$ B_1 @f$, through @f$ B_7 @f$  are functions of the
     * acentric factor, reduced dipole moment and the association factor.
     *
     * @f[
     *  B_i = a_i + b_i \omega + c_i \mu_r^4 + d_i \kappa
     *
     *  \quad \text{(Equation 10-5.9)}
     * @f]
     *
     * The constants in the above equation are from table 10-3 on page 10.23 of
     * @cite poling2001.
     *
     * The definition of the @f$ \Psi @f$ function is given by:
     *
     * @f[
     *  \Psi = 1 + \alpha {\frac{0.215 + 0.28288\alpha - 1.061\beta + 0.26665Z}{0.6366 + \beta Z + 1.061 \alpha \beta}}
     * @f]
     *
     * with,
     *
     * @f[
     *  \alpha = \frac{C_v}{R} - \frac{3}{2}
     * @f]
     *
     * @f[
     *  \beta = 0.7862 - 0.7109 \omega + 1.3168 \omega^2
     * @f]
     *
     * @f[
     *  Z = 2.0 + 10.5 T_r^2
     * @f]
     *
     * These functions are from page 10.12 of @cite poling2001 .
     *
     *
     * This method utilizes the low-pressure Chung viscosity as that is a required
     * parameter in the model, and thus calls the low pressure viscosity
     * implementation. This is why it requires parameters typically associated with
     * the viscosity calculation.
     *
     * This function is structured such that it can be used for pure species or
     * mixtures, with the only difference being the values that are passed to the
     * function (pure values versus mixture values).
     *
     * For mixtures, the mixture values of the input variables are computed using the
     * mixing rules of Chung, @see computeMixtureParameters() .
     *
     * @param T  Temperature [K]
     * @param T_star  Reduced temperature [unitless]
     * @param MW  Molecular weight [kg/kmol]
     * @param rho  Density [mol/cm^3]
     * @param Cv  Specific heat [J/kg/K]
     * @param Vc  Critical volume [cm^3/mol]
     * @param Tc  Critical temperature [K]
     * @param sigma  Lennard-Jones collision diameter [Angstroms]
     * @param acentric_factor  Acentric factor [unitless]
     * @param mu_r  Dipole moment [Debye]
     * @param kappa  Polar correction factor [unitless]
     * @return double
     */
    double highPressureThermalConductivity(double T, double T_star, double MW,
        double rho, double Cv, double Vc, double Tc, double sigma,
        double acentric_factor, double mu_r, double kappa);

    private:
    //! Critical properties of each species.
    std::vector<double> m_Tcrit;
    std::vector<double> m_Pcrit;
    std::vector<double> m_Vcrit;
    std::vector<double> m_Zcrit;

    //! Pure fluid properties of each species.
    std::vector<double> sigma_i;
    std::vector<double> epsilon_over_k_i;
    std::vector<double> MW_i;
    std::vector<double> acentric_factor_i;
    std::vector<double> kappa_i;

    // Parameters that are needed to calculate the viscosity using the Chung method.

    // Mixture critical properties used by the Chung viscosity model.
    double m_Vc_mix = 0;
    double m_Tc_mix = 0;

    // Values associated with the calculation of sigma and the molecular weight used in the Chung viscosity model.
    double m_sigma_mix = 0;
    double m_epsilon_over_k_mix = 0;
    double m_MW_mix = 0;

    // Values associated with the calculation of the Fc factor in the Chung viscosity model.
    double m_mu_mix = 0;
    double m_mu_r_mix = 0;
    double m_acentric_factor_mix = 0;
    double m_kappa_mix = 0;

};

}
#endif
