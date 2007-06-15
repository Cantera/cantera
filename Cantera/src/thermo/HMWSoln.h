/**
 *  @file HMWSoln.h
 *    Header file for Pitzer activity coefficient implementation 
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Id$
 */

#ifndef CT_HMWSOLN_H
#define CT_HMWSOLN_H

#include "MolalityVPSSTP.h"
#include "electrolytes.h"

namespace Cantera {


  /**
   * Major Parameters:
   *   The form of the Pitzer expression refers to the 
   *   form of the Gibbs free energy expression. The temperature
   *   dependence of the Pitzer coefficients are handled by
   *   another parameter.
   *
   * m_formPitzer = Form of the Pitzer expression
   *
   *  PITZERFORM_BASE = 0
   *
   *   Only one form is supported atm. This parameter is included for
   *   future expansion.
   *  
   */
#define PITZERFORM_BASE 0


  /*!
   * @name Temperature Dependence of the Pitzer Coefficients
   *
   *  Note, the temperature dependence of the
   * Gibbs free energy also depends on the temperature dependence
   * of the standard state and the temperature dependence of the
   * Debye-Huckel constant, which includes the dielectric constant
   * and the density. Therefore, this expression defines only part
   * of the temperature dependence for the mixture thermodynamic
   * functions.
   *
   *  PITZER_TEMP_CONSTANT
   *     All coefficients are considered constant wrt temperature
   *  PITZER_TEMP_LINEAR
   *     All coefficients are assumed to have a linear dependence
   *     wrt to temperature.
   *  PITZER_TEMP_COMPLEX1
   *     All coefficnets are assumed to have a complex functional
   *     based dependence wrt temperature;  See:
   *    (Silvester, Pitzer, J. Phys. Chem. 81, 19 1822 (1977)).
   *
   *       beta0 = q0 + q3(1/T - 1/Tr) + q4(ln(T/Tr)) +
   *               q1(T - Tr) + q2(T**2 - Tr**2)
   * 
   */
  //@{
#define PITZER_TEMP_CONSTANT   0
#define PITZER_TEMP_LINEAR     1
#define PITZER_TEMP_COMPLEX1   2
  //@}

  /*
   *  @name ways to calculate the value of A_Debye
   *
   *  These defines determine the way A_Debye is calculated
   */
 //@{
#define    A_DEBYE_CONST  0
#define    A_DEBYE_WATER  1
  //@}

  class WaterProps;
  class WaterPDSS;

  /**
   * Class %HMWSoln represents a dilute or concentrated liquid electrolyte
   * phase which obeys the Pitzer formulation for nonideality.
   *
   * As a prerequisite to the specification of thermodynamic quantities,
   * The concentrations of the ionic species are assumed to obey the 
   * electroneutrality condition. 
   *
   * <HR>
   * <H2> Specification of Species Standard %State Properties </H2>
   * <HR>
   *
   * The solvent is assumed to be liquid water. A real model for liquid
   * water (IAPWS 1995 formulation) is used as its standard state.
   * All standard state properties for the solvent are based on
   * this real model for water, and involve function calls
   * to the object that handles the real water model, #WaterPropsIAPWS.
   *
   * The standard states for solutes are on the unit molality basis. 
   * Therefore, in the documentation below, the normal \f$ o \f$ 
   * superscript is replaced with
   * the \f$ \triangle \f$ symbol. The reference state symbol is now
   *  \f$ \triangle, ref \f$.
   * 
   *  
   *  It is assumed that the reference state thermodynamics may be
   *  obtained by a pointer to a populated species thermodynamic property
   *  manager class (see ThermoPhase::m_spthermo). How to relate pressure
   *  changes to the reference state thermodynamics is resolved at this level.
   *
   *  For solutes that rely on ThermoPhase::m_spthermo, are assumed to
   *  have an incompressible standard state mechanical property.
   *  In other words, the molar volumes are independent of temperature
   *  and pressure.
   *
   *  For these incompressible,
   *  standard states, the molar internal energy is
   *  independent of pressure. Since the thermodynamic properties
   *  are specified by giving the standard-state enthalpy, the
   *  term \f$ P_0 \hat v\f$ is subtracted from the specified molar
   *  enthalpy to compute the molar internal energy. The entropy is
   *  assumed to be independent of the pressure.
   *
   * The enthalpy function is given by the following relation.
   *
   *       \f[
   *   \raggedright  h^\triangle_k(T,P) = h^{\triangle,ref}_k(T) 
   *         + \tilde{v}_k \left( P - P_{ref} \right) 
   *       \f]
   *
   * For an incompressible,
   * stoichiometric substance, the molar internal energy is
   * independent of pressure. Since the thermodynamic properties
   * are specified by giving the standard-state enthalpy, the
   * term \f$ P_{ref} \tilde v\f$ is subtracted from the specified reference molar
   * enthalpy to compute the molar internal energy.
   *
   *       \f[
   *            u^\triangle_k(T,P) = h^{\triangle,ref}_k(T) - P_{ref} \tilde{v}_k
   *       \f]
   *
   *
   * The solute standard state heat capacity and entropy are independent
   * of pressure. The solute standard state gibbs free energy is obtained
   * from the enthalpy and entropy functions.
   *
   * The vector Constituents::m_speciesSize[] is used to hold the
   * base values of species sizes. These are defined as the 
   * molar volumes of species at infinite dilution at 300 K and 1 atm
   * of water. m_speciesSize are calculated during the initialization of the
   * %HMWSoln object and are then not touched.
   *
   * The current model assumes that an incompressible molar volume for
   * all solutes. The molar volume for the water solvent, however,
   * is obtained from a pure water equation of state, waterSS.
   * Therefore, the water standard state varies with both T and P.
   * It is an error to request standard state water properties  at a T and P
   * where the water phase is not a stable phase, i.e., beyond its
   * spinodal curve.
   *   
   * <HR>
   * <H2> Specification of Solution Thermodynamic Properties </H2>
   * <HR>
   *
   * Chemical potentials
   * of the solutes,  \f$ \mu_k \f$, and the solvent, \f$ \mu_o \f$, which are based 
   * on the molality form, have the following general format:
   *
   * \f[
   *    \mu_k = \mu^{\triangle}_k(T,P) + R T ln(\gamma_k^{\triangle} \frac{m_k}{m^\triangle}) 
   * \f]
   * \f[
   *    \mu_o = \mu^o_o(T,P) + RT ln(a_o) 
   * \f]
   *
   * where \f$ \gamma_k^{\triangle} \f$ is the molality based activity coefficient for species
   * \f$k\f$.
   * 
   * Individual activity coefficients of ions can not be independently measured. Instead,
   * only binary pairs forming electroneutral solutions can be measured. This problem
   * leads to a redundancy in the evaluation of species standard state properties.
   * The redundancy issue is resolved by setting the standard state chemical potential
   * enthalpy, entropy, and volume for the hydrogen ion, H+, to zero, for every temperature
   * and pressure. After this convention is applied, all other standard state 
   * properties of ionic species contain meaningfull information.
   *
   *
   *  <H3> Ionic Strength </H3>
   *
   *  Most of the parameterizations within the model use the ionic strength 
   *  as a key variable. The ionic strength, \f$ I\f$ is defined as follows
   *
   *  \f[
   *    I = \frac{1}{2} \sum_k{m_k  z_k^2}
   *  \f]
   *
   *
   * \f$ m_k \f$ is the molality of the kth species. \f$ z_k \f$ is the charge
   * of the kth species. Note, the ionic strength is a defined units quantity.
   * The molality has defined units of gmol kg-1, and therefore the ionic
   * strength has units of sqrt( gmol kg-1).
   *
   * In some instances, from some authors, a different 
   * formulation is used for the ionic strength in the equations below. The different
   * formulation is due to the possibility of the existence of weak acids and how
   * association wrt to the weak acid equilibrium relation affects the calculation 
   * of the activity coefficients via the assumed value of the ionic strength.
   *
   * If we are to assume that the association reaction doesn't have an effect
   * on the ionic strength, then we will want to consider the associated weak
   * acid as in effect being fully dissociated, when we calculate an effective
   * value for the ionic strength. We will call this calculated value, the
   * stoichiometric ionic strength, \f$ I_s \f$, putting a subscript s to denote
   * it from the more straightforward calculation of \f$ I \f$.
   *
   *  \f[
   *    I_s = \frac{1}{2} \sum_k{m_k^s  z_k^2}
   *  \f]
   *
   *  Here, \f$ m_k^s \f$ is the value of the molalities calculated assuming that
   *  all weak acid-base pairs are in their fully dissociated states. This calculation may
   * be simplified by considering that the weakly associated acid may be made up of two
   * charged species, k1 and k2, each with their own charges, obeying the following relationship:
   *
   *   \f[
   *      z_k = z_{k1} +  z_{k2}
   *   \f]
   *  Then, we may only need to specify one charge value, say, \f$  z_{k1}\f$, 
   *  the cation charge number,
   *  in order to get both numbers, since we have already specified \f$ z_k \f$ 
   *  in the definition of original species.
   *  Then, the stoichiometric ionic strength may be calculated via the following formula.
   *
   *  \f[
   *    I_s = \frac{1}{2} \left(\sum_{k,ions}{m_k  z_k^2}+ 
   *               \sum_{k,weak_assoc}(m_k  z_{k1}^2 + m_k  z_{k2}^2) \right)
   *  \f]
   *
   *  The specification of which species are weakly associated acids is made in the input 
   *  file via the
   *  <TT> stoichIsMods </TT> XML block, where the charge for k1 is also specified. 
   *  An example is given below:
   * 
   * @code
   *          <stoichIsMods>
   *                NaCl(aq):-1.0
   *          </stoichIsMods>
   * @endcode
   *
   *
   *  Because we need the concept of a weakly associated acid in order to calculated 
   *  \f$ I_s \f$ we need to 
   *  catalog all species in the phase. This is done using the following categories:
   *
   *  -  <B>cEST_solvent</B>    :           Solvent species (neutral)
   *  -  <B>cEST_chargedSpecies</B>         Charged species (charged)
   *  -  <B>cEST_weakAcidAssociated</B>     Species which can break apart into charged species.
   *                                        It may or may not be charged.  These may or 
   *                                        may not be be included in the
   *                                        species solution vector.
   *  -  <B>cEST_strongAcidAssociated</B>   Species which always breaksapart into charged species.
   *                                        It may or may not be charged. Normally, these
   *                                        aren't included in the speciation vector.
   *  -  <B>cEST_polarNeutral </B>          Polar neutral species
   *  -  <B>cEST_nonpolarNeutral</B>        Non poloar neutral species
   *
   *  Polar and non-polar neutral species are differentiated, because some additions 
   *  to the activity 
   *  coefficient expressions distinguish between these two types of solutes. This is the so-called
   *  salt-out effect.
   *
   * The type of species is specified in the <TT>electrolyteSpeciesType</TT> XML block.
   * Note, this is not
   * considered a part of the specification of the standard state for the species, 
   * at this time. Therefore,
   * this information is put under the <TT>activityCoefficient</TT> XML block. An example 
   * is given below
   *
   * @code
   *         <electrolyteSpeciesType>
   *                H2L(L):solvent
   *                H+:chargedSpecies
   *                NaOH(aq):weakAcidAssociated
   *                NaCl(aq):strongAcidAssociated
   *                NH3(aq):polarNeutral
   *                O2(aq):nonpolarNeutral
   *         </electrolyteSpeciesType>
   * @endcode
   *
   *
   *  Much of the species electrolyte type information is infered from other information in the
   *  input file. For example, as species which is charged is given the "chargedSpecies" default
   *  category. A neutral solute species is put into the "nonpolarNeutral" category by default.
   *
   *
   *  <H3> Multicomponent Activity Coefficients for Solutes </H3>
   *
   *    In the formulas below the following conventions are used. The subscript <I>M</I> refers
   *    to a particular cation. The subscript X refers to a particular anion, whose
   *    activity is being currently evaluated. the subscript <I>a</I> refers to a summation
   *    over all anions in the solution, while the subscript <I>c</I> refers to a summation
   *    over all cations in the solutions.
   *
   *     The activity coefficient for a particular cation <I>M</I> is given by
   *
   *   \f[
   *      \ln(\gamma_M^\triangle) = -z_M^2(F) + \sum_a m_a \left( 2 B_{Ma} + Z C_{Ma} \right)
   *      + z_M   \left( \sum_a  \sum_c m_a m_c C_{ca} \right)
   *             + \sum_c m_c \left[ 2 \Phi_{Mc} + \sum_a m_a \psi_{Mca} \right]
   *             + \sum_{a < a'} \sum m_a m_{a'} \psi_{Ma{a'}}
   *             +  2 \sum_n m_n \lambda_{nM}
   *   \f]
   *
   *     The activity coefficient for a particular anion <I>X</I> is given by
   *
   *   \f[
   *      \ln(\gamma_X^\triangle) = -z_X^2(F) + \sum_a m_c \left( 2 B_{cX} + Z C_{cX} \right)
   *      + \left|z_X \right|  \left( \sum_a  \sum_c m_a m_c C_{ca} \right)
   *             + \sum_a m_a \left[ 2 \Phi_{Xa} + \sum_c m_c \psi_{cXa} \right]
   *             + \sum_{c < c'} \sum m_c m_{c'} \psi_{c{c'}X}
   *             +  2 \sum_n m_n \lambda_{nM}
   *   \f]
   *              where the function \f$ F \f$ is given by
   *
   *
   *   \f[
   *       F = - A_{\phi} \left[ \frac{\sqrt{I}}{1 + b \sqrt{I}} 
   *                 + \frac{2}{b} \ln{\left(1 + b\sqrt{I}\right)} \right]
   *                 + \sum_a \sum_c m_a m_c B'_{ca}
   *                 + \sum_{c < c'} \sum m_c m_{c'} \Phi'_{c{c'}}
   *                 + \sum_{a < a'} \sum m_a m_{a'} \Phi'_{a{a'}}
   *   \f]
   *
   *   where \f$ I\f$ is the ionic strength
   *
   *   \f[
   *       I = \frac{1}{2} \sum_k{m_k  z_k^2}
   *   \f]
   *
   *   and the function \f$ Z \f$ is given by
   *
   *   \f[
   *       Z = \sum_i m_i \left| z_i \right|
   *   \f]
   *
   *   In the above formulas, \f$ \Phi'_{c{c'}} \f$  and \f$  \Phi'_{a{a'}} \f$ are the
   *   ionic strength derivatives of \f$ \Phi_{c{c'}} \f$  and \f$  \Phi_{a{a'}} \f$,
   *   respectively.
   *
   *   The function \f$ B'_{MX} \f$ is defined as:
   *
   *   \f[
   *       B'_{MX} = \left( \frac{\beta^1_{MX} h(\alpha \sqrt{I})}{I}  \right) 
   *                 \left( \frac{\beta^2_{MX} h(\alpha \sqrt{I})}{I}  \right)
   *   \f]
   *
   *  where \f$ h(x) \f$ is defined as
   *
   *   \f[
   *       h(x) = g'(x) \frac{x}{2} = 
   *        \frac{2\left(1 - \left(1 + x + \frac{x^2}{2} \right)\exp(-x) \right)}{x^2}
   *   \f]
   *
   *   The activity coefficient for neutral species <I>N</I> is given by 
   *
   *   \f[
   *       \ln(\gamma_N^\triangle) = 2 \left( \sum_i m_i \lambda_{iN}\right)
   *   \f]
   *
   *
   *  <H3> Activity of the Water Solvent </H3>
   *
   *  The activity for the solvent water,\f$ a_o \f$, is not independent and must be 
   *  determined either from the Gibbs-Duhem relation or from taking the appropriate derivative
   *  of the same excess Gibbs free energy function as was used to formulate
   *  the solvent activity coefficients. Pitzer's description follows the later approach to
   *  derive a formula for the osmotic coefficient, \f$ \phi \f$.
   *
   *  \f[
   *       \phi - 1 = - \left( \frac{d\left(\frac{G^{ex}}{RT} \right)}{d(\tilde{M}_o n_o)}  \right)
   *                \frac{1}{\sum_{i \ne 0} m_i}
   *  \f]
   *
   *  The result is the following
   *
   *  \f[
   *     \phi - 1 =  
   *          \frac{2}{\sum_{i \ne 0} m_i}
   *           \left[
   *           \begin{array}{c}
   *        -  A_{\phi} \frac{I^{3/2}}{1 + b \sqrt{I}}  
   *        +   \sum_c  \sum_a m_c m_a \left( B^{\phi}_{ca} + Z C_{ca}\right) 
   *          \\
   *        +   \sum_{c < c'} \sum m_c m_{c'} \left[ \Phi^{\phi}_{c{c'}} + \sum_a m_a \Psi_{c{c'}a} \right]
   *        +   \sum_{a < a'} \sum m_a m_{a'} \left[ \Phi^{\phi}_{a{a'}} + \sum_c m_c \Psi_{a{a'}c} \right]
   *          \\
   *        + \sum_n \sum_c m_n m_c \lambda_{nc} +  \sum_n \sum_a m_n m_a \lambda_{na} 
   *        + \sum_{n < n'} \sum m_n m_{n'} \lambda_{n{n'}}
   *        + \frac{1}{2} \left( \sum_n m^2_n \lambda_{nn}\right)
   *          \end{array}
   *          \right]
   *  \f]
   *     
   *
   *
   * An example is given below.
   *
   * An example <TT> activityCoefficients </TT> XML block for this formulation is supplied below
   *
   *   * @code
   *  <activityCoefficients model="Beta_ij">
   *         <!-- A_Debye units = sqrt(kg/gmol) -->
   *         <A_Debye> 1.172576 </A_Debye>
   *         <!-- B_Debye units = sqrt(kg/gmol)/m   -->
   *         <B_Debye> 3.28640E9 </B_Debye>
   *         <ionicRadius default="3.042843"  units="Angstroms">
   *         </ionicRadius>
   *         <DHBetaMatrix>
   *               H+:Cl-:0.27
   *               Na+:Cl-:0.15
   *               Na+:OH-:0.06
   *         </DHBetaMatrix>
   *         <stoichIsMods>
   *                NaCl(aq):-1.0
   *         </stoichIsMods>
   *         <electrolyteSpeciesType>
   *                H+:chargedSpecies
   *                NaCl(aq):weakAcidAssociated
   *         </electrolyteSpeciesType>
   *  </activityCoefficients>
   * @endcode
   *
   *
   * <H3> Specification of the Debye Huckel Constants </H3>
   *
   *  In the equations above, the formulas for  \f$  A_{Debye} \f$ and \f$  B_{Debye} \f$ 
   *  are needed. The %DebyeHuckel object uses two methods for specifying these quantities.
   *  The default method is to assume that \f$  A_{Debye} \f$  is a constant, given
   *  in the initialization process, and storred in the
   *  member double, m_A_Debye. Optionally, a full water treatment may be employed that makes
   *  \f$ A_{Debye} \f$ a full function of <I>T</I> and <I>P</I>.
   *
   *   \f[
   *      A_{Debye} = \frac{F e B_{Debye}}{8 \pi \epsilon R T} {\left( C_o \tilde{M}_o \right)}^{1/2}
   *   \f]
   * where
   * 
   *  \f[
   *         B_{Debye} = \frac{F} {{(\frac{\epsilon R T}{2})}^{1/2}} 
   *  \f]
   *  Therefore:
   * \f[
   *   A_{Debye} = \frac{1}{8 \pi} 
   *                 {\left(\frac{2 N_a \rho_o}{1000}\right)}^{1/2}
   *                 {\left(\frac{N_a e^2}{\epsilon R T }\right)}^{3/2}
   * \f]
   *
   *            Units = sqrt(kg/gmol)
   *
   *
   *     where
   *      - \f$ N_a \f$ is Avrogadro's number
   *      - \f$ \rho_w \f$ is the density of water
   *      - \f$ e \f$ is the electronic charge
   *      - \f$ \epsilon = K \epsilon_o \f$ is the permitivity of water
   *           where \f$ K \f$ is the dielectric condstant of water,
   *           and  \f$ \epsilon_o \f$ is the permitivity of free space.
   *      - \f$ \rho_o \f$ is the density of the solvent in its standard state.
   *
   *            Nominal value at 298 K and 1 atm = 1.172576 (kg/gmol)<SUP>1/2</SUP>
   *                  based on:
   *                 -   \f$ \epsilon / \epsilon_0 \f$ = 78.54
   *                           (water at 25C)
   *                 -   \f$ \epsilon_0 \f$= 8.854187817E-12 C<SUP>2</SUP> N<SUP>-1</SUP> m<SUP>-2</SUP>
   *                 -   e = 1.60217653E-19 C
   *                 -   F = 9.6485309E7 C kmol<SUP>-1</SUP>
   *                 -   R = 8.314472E3 kg m<SUP>2</SUP> s<SUP>-2</SUP> kmol<SUP>-1</SUP> K<SUP>-1</SUP>
   *                 -   T = 298.15 K
   *                 -   B_Debye = 3.28640E9 (kg/gmol)<SUP>1/2</SUP> m<SUP>-1</SUP>
   *                 -   \f$N_a\f$ = 6.0221415E26 kmol<SUP>-1</SUP>
   *
   * An example of a fixed value implementation is given below.
   * @code
   *   <activityCoefficients model="Beta_ij">
   *         <!-- A_Debye units = sqrt(kg/gmol)  -->
   *         <A_Debye> 1.172576 </A_Debye>
   *         <!-- B_Debye units = sqrt(kg/gmol)/m  -->
   *         <B_Debye> 3.28640E9 </B_Debye>
   *   </activityCoefficients>
   * @endcode
   *
   * An example of a variable value implementation is given below.
   *
   * @code
   *   <activityCoefficients model="Beta_ij">
   *         <A_Debye model="water" /> 
   *         <!-- B_Debye units = sqrt(kg/gmol)/m  -->
   *         <B_Debye> 3.28640E9 </B_Debye>
   *   </activityCoefficients>
   * @endcode
   *
   * An example of a variable value implementation is given below.
   *
   * @code
   *   <activityCoefficients model="Beta_ij">
   *         <A_Debye model="water" /> 
   *         <!-- B_Debye units = sqrt(kg/gmol)/m  -->
   *         <B_Debye> 3.28640E9 </B_Debye>
   *   </activityCoefficients>
   * @endcode
   *
   *  Currently, \f$  B_{Debye} \f$ is a constant in the model, specified either by a default
   *  water value, or through the input file. This may have to be looked at, in the future.
   *
   * <HR>
   * <H2> %Application within %Kinetics Managers </H2>
   * <HR>
   *
   * For the time being, we have set the standard concentration for all species in
   * this phase equal to the default concentration of the solvent at 298 K and 1 atm. 
   * This means that the
   * kinetics operator essentially works on an activities basis, with units specified
   * as if it were on a concentration basis.
   *
   * For example, a bulk-phase binary reaction between liquid species j and k, producing
   * a new liquid species l would have the
   * following equation for its rate of progress variable, \f$ R^1 \f$, which has
   * units of kmol m-3 s-1.
   *
   *   \f[
   *    R^1 = k^1 C_j^a C_k^a =  k^1 (C_o a_j) (C_o a_k) 
   *   \f]
   * where
   *   \f[
   *      C_j^a = C_o a_j \quad and \quad C_k^a = C_o a_k
   *   \f]
   *   
   *  \f$ C_j^a \f$ is the activity concentration of species j, and 
   *  \f$ C_k^a \f$ is the activity concentration of species k. \f$ C_o \f$
   *  is the concentration of water at 298 K and 1 atm. \f$ a_j \f$ is
   *  the activity of species j at the current temperature and pressure
   *  and concentration of the liquid phase. \f$k^1 \f$ has units of m3 kmol-1 s-1.
   *
   *
   *
  *  The reverse rate constant can then be obtained from the law of microscopic reversibility
   * and the equilibrium expression for the system.
   *
   *   \f[
   *         \frac{a_j a_k}{ a_l} = K^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
   *   \f]
   *
   *  \f$  K^{o,1} \f$ is the dimensionless form of the equilibrium constant.
   *  
   *   \f[
   *    R^{-1} = k^{-1} C_l^a =  k^{-1} (C_o a_l)
   *   \f]
   *
   *  where
   *
   *    \f[
   *       k^{-1} =  k^1 K^{o,1} C_o
   *   \f]
   *
   *  \f$k^{-1} \f$ has units of s-1.
   * 
   *  Note, this treatment may be modified in the future, as events dictate.
   *
   * <HR>
   * <H2> Instantiation of the Class </H2>
   * <HR>
   *
   *   * The constructor for this phase is NOT located in the default ThermoFactory
   * for %Cantera. However, a new %DebyeHuckel object may be created by 
   * the following code snippets:
   *
   * @code
   *      DebyeHuckel *DH = new DebyeHuckel("DH_NaCl.xml", "NaCl_electrolyte");
   * @endcode
   *
   * or
   *
   * @code
   *    char iFile[80], file_ID[80];
   *    strcpy(iFile, "DH_NaCl.xml");
   *    sprintf(file_ID,"%s#NaCl_electrolyte", iFile);
   *    XML_Node *xm = get_XML_NameID("phase", file_ID, 0);
   *    DebyeHuckel *dh = new DebyeHuckel(*xm);
   * @endcode
   *
   * or by the following call to importPhase():
   *
   * @code
   *    char iFile[80], file_ID[80];
   *    strcpy(iFile, "DH_NaCl.xml");
   *    sprintf(file_ID,"%s#NaCl_electrolyte", iFile);
   *    XML_Node *xm = get_XML_NameID("phase", file_ID, 0);
   *    DebyeHuckel dhphase;
   *    importPhase(*xm, &dhphase);
   * @endcode
   *
   * <HR>
   * <H2> XML Example </H2>
   * <HR>
   *
   * The phase model name for this is called StoichSubstance. It must be supplied
   * as the model attribute of the thermo XML element entry.
   * Within the phase XML block,
   * the density of the phase must be specified. An example of an XML file
   * this phase is given below. 
   * 
   * @verbatim
   <phase id="NaCl_electrolyte" dim="3">
    <speciesArray datasrc="#species_waterSolution">
               H2O(L) Na+ Cl- H+ OH- NaCl(aq) NaOH(aq)
    </speciesArray>
    <state>
      <temperature units="K"> 300  </temperature>
      <pressure units="Pa">101325.0</pressure>
      <soluteMolalities>
             Na+:3.0
             Cl-:3.0
             H+:1.0499E-8
             OH-:1.3765E-6
             NaCl(aq):0.98492
             NaOH(aq):3.8836E-6
      </soluteMolalities>
    </state>
    <!-- thermo model identifies the inherited class
         from ThermoPhase that will handle the thermodynamics.
      -->
    <thermo model="DebyeHuckel">
       <standardConc model="solvent_volume" />
       <activityCoefficients model="Beta_ij">
                <!-- A_Debye units = sqrt(kg/gmol)  -->
                <A_Debye> 1.172576 </A_Debye>
                <!-- B_Debye units = sqrt(kg/gmol)/m   -->
                <B_Debye> 3.28640E9 </B_Debye>
                <ionicRadius default="3.042843"  units="Angstroms">
                </ionicRadius>
                <DHBetaMatrix>
                  H+:Cl-:0.27
                  Na+:Cl-:0.15
                  Na+:OH-:0.06
                </DHBetaMatrix>
                <stoichIsMods>
                   NaCl(aq):-1.0
                </stoichIsMods>
                <electrolyteSpeciesType>
                   H+:chargedSpecies
                   NaCl(aq):weakAcidAssociated
                </electrolyteSpeciesType>
       </activityCoefficients>
       <solvent> H2O(L) </solvent>
    </thermo>
    <elementArray datasrc="elements.xml"> O H Na Cl </elementArray>
  </phase> 
  @endverbatim
   *
   *
   *
   * @ingroup thermoprops
   *
   */
  class HMWSoln : public MolalityVPSSTP {

  public:
        
    //! Default Constructor 
    HMWSoln();

    //! Full constructor for setting up the entire ThermoPhase Object
    /*!
     * Working constructors
     *
     *  The two constructors below are the normal way
     *  the phase initializes itself. They are shells that call
     *  the routine initThermo(), with a reference to the
     *  XML database to get the info for the phase.
     *
     * @param inputFile Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     */
    HMWSoln(std::string inputFile, std::string id = "");

    //! Full constructor for creating the phase.
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase. 
     *                (default is the empty string)
     */
    HMWSoln(XML_Node& phaseRef, std::string id = "");


    //! Copy Constructor
    /*!
     * Copy constructor for the object. Constructed
     * object will be a clone of this object, but will
     * also own all of its data.
     * This is a wrapper around the assignment operator
     * 
     * @param right Object to be copied.
     */
    HMWSoln(const HMWSoln &right);

    //! Asignment operator
    /*!
     * Assignment operator for the object. Constructed
     * object will be a clone of this object, but will
     * also own all of its data.
     * 
     * @param right Object to be copied.
     */
    HMWSoln& operator=(const HMWSoln& right);

    
    //!  This is a special constructor, used to replicate test problems
    //!  during the initial verification of the object
    /*!  
     *
     *
     *  test problems:
     *  1 = NaCl problem - 5 species -
     *   the thermo is read in from an XML file
     *
     * speci   molality                        charge
     *  Cl-     6.0954          6.0997E+00      -1
     *  H+      1.0000E-08      2.1628E-09      1
     *  Na+     6.0954E+00      6.0997E+00      1
     *  OH-     7.5982E-07      1.3977E-06     -1
     *  HMW_params____beta0MX__beta1MX__beta2MX__CphiMX_____alphaMX__thetaij
     * 10
     * 1  2          0.1775  0.2945   0.0      0.00080    2.0      0.0
     * 1  3          0.0765  0.2664   0.0      0.00127    2.0      0.0
     * 1  4          0.0     0.0      0.0      0.0        0.0     -0.050
     * 2  3          0.0     0.0      0.0      0.0        0.0      0.036
     * 2  4          0.0     0.0      0.0      0.0        0.0      0.0
     * 3  4          0.0864  0.253    0.0      0.0044     2.0      0.0
     * Triplet_interaction_parameters_psiaa'_or_psicc'
     * 2
     * 1  2  3   -0.004
     * 1  3  4   -0.006
     *
     * @param testProb Hard -coded test problem to instantiate.
     *                 Current valid values are 1.
     */
    HMWSoln(int testProb);

    //! Destructor. 
    virtual ~HMWSoln();

    //! Duplicator from the ThermoPhase parent class
    /*!
     * Given a pointer to a ThermoPhase object, this function will
     * duplicate the ThermoPhase object and all underlying structures.
     * This is basically a wrapper around the copy constructor.
     *
     * @return returns a pointer to a ThermoPhase
     */
    ThermoPhase *duplMyselfAsThermoPhase() const;

    /**
     *   
     * @name  Utilities  
     * @{
     */

    /** 
     * Equation of state type flag. The base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Constants defined for this purpose are
     * listed in mix_defs.h.
     */
    virtual int eosType() const;

    /**
     * @} 
     * @name  Molar Thermodynamic Properties of the Solution --------------
     * @{
     */

    /// Molar enthalpy. Units: J/kmol. 
    /**
     * Molar enthalpy of the solution. Units: J/kmol.
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal enthalpy_mole() const;

    /**
     * Excess molar enthalpy of the solution from 
     * the mixing process. Units: J/ kmol.
     *
     * Note this is kmol of the total solution.
     */
    virtual doublereal relative_enthalpy() const;

    /**
     * Excess molar enthalpy of the solution from 
     * the mixing process on a molality basis.
     *  Units: J/ (kmol add salt).
     *
     * Note this is kmol of the guessed at salt composition
     */
    virtual doublereal relative_molal_enthalpy() const;


    /// Molar internal energy. Units: J/kmol. 
    /**
     * Molar internal energy of the solution. Units: J/kmol.
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal intEnergy_mole() const;

    /// Molar entropy. Units: J/kmol/K. 
    /**
     * Molar entropy of the solution. Units: J/kmol/K.
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     * \f[
     * \hat s(T, P, X_k) = \sum_k X_k \hat s^0_k(T)  
     *      - \hat R  \sum_k X_k log(X_k)
     * \f]
     * The reference-state pure-species entropies 
     * \f$ \hat s^0_k(T,p_{ref}) \f$ are computed by the
     *  species thermodynamic 
     * property manager. The pure species entropies are independent of 
     * temperature since the volume expansivities are equal to zero.
     * @see SpeciesThermo
     *
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal entropy_mole() const;

    /// Molar Gibbs function. Units: J/kmol. 
    /*
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal gibbs_mole() const;

    /// Molar heat capacity at constant pressure. Units: J/kmol/K. 
    /*
     *     
     */
    virtual doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K.
    /*
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal cv_mole() const;

    //@}
    /// @name Mechanical Equation of State Properties ---------------------
    //@{
    /**
     *   In this equation of state implementation, the density is a 
     *   function only of the mole fractions. Therefore, it can't be 
     *   an independent variable. Instead, the pressure is used as the
     *   independent variable. Functions which try to set the thermodynamic
     *   state by calling setDensity() may cause an exception to be
     *   thrown.  
     */

    /**
     * Pressure. Units: Pa.
     * For this incompressible system, we return the internally storred
     * independent value of the pressure.
     */ 
    virtual doublereal pressure() const;

    //! Set the internally storred pressure (Pa) at constant
    //! temperature and composition
    /*!
     *  This method sets the pressure within the object.
     *  The water model is a completely compressible model.
     *  Also, the dielectric constant is pressure dependent.
     *
     *  @param p input Pressure (Pa)
     *
     * @todo Implement a variable pressure capability
     */
    virtual void setPressure(doublereal p);

    /**
     * Calculate the density of the mixture using the partial 
     * molar volumes and mole fractions as input
     *
     * The formula for this is
     *
     * \f[ 
     * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}} 
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are
     * the molecular weights, and \f$V_k\f$ are the pure species
     * molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal
     * solution the partial molar volumes are equal to the pure
     * species molar volumes. We have additionally specified
     * in this class that the pure species molar volumes are
     * independent of temperature and pressure.
     *
     * NOTE: This is a non-virtual function, which is not a 
     *       member of the ThermoPhase base class. 
     */
    void calcDensity();

    //! Set the internally storred density (gm/m^3) of the phase.
    /*!
     * Overwritten setDensity() function is necessary because of
     * the underlying water model.
     *
     * @todo Now have a compressible ss equation for liquid water.
     *       Therefore, this phase is compressible. May still
     *       want to change the independent variable however. 
     *
     *  NOTE: This is an overwritten function from the State.h
     *        class
     * 
     * @param rho Input density (kg/m^3).
     */
    void setDensity(doublereal rho);


    //! Set the internally storred molar density (kmol/m^3) for the phase.
    /**
     * Overwritten setMolarDensity() function is necessary because of the
     * underlying water model.
     *
     * This function will now throw an error condition if the input
     * isn't exactly equal to the current molar density.
     *
     *  NOTE: This is a virtual function overwritten from the State.h
     *        class
     *
     * @param conc   Input molar density (kmol/m^3).
     */
    void setMolarDensity(doublereal conc);

    //! Set the temperature (K)
    /*!
     * Overwritten setTemperature(double) from State.h. This
     * function sets the temperature, and makes sure that
     * the value propagates to underlying objects, such as
     * the water standard state model.
     *
     * @todo Make State::setTemperature a virtual function
     *
     * @param temp Temperature in kelvin
     */
    virtual void setTemperature(doublereal temp);

    /**
     * The isothermal compressibility. Units: 1/Pa.
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const;

    /**
     * The thermal expansion coefficient. Units: 1/K.
     * The thermal expansion coefficient is defined as
     *
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal thermalExpansionCoeff() const;

    /**
     * @} 
     * @name Potential Energy
     * 
     * Species may have an additional potential energy due to the
     * presence of external gravitation or electric fields. These
     * methods allow specifying a potential energy for individual
     * species.
     * @{
     */



    /**
     * @}
     * @name Activities, Standard States,  and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
     * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T,P)\f$ is
     * the chemical potential at unit activity, which depends only
     * on temperature and the pressure.
     * Activity is assumed to be molality-based here.
     * @{
     */

    /**
     * This method returns an array of generalized concentrations
     * \f$ C_k\f$ that are defined such that 
     * \f$ a_k = C_k / C^0_k, \f$ where \f$ C^0_k \f$ 
     * is a standard concentration
     * defined below.  These generalized concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions. 
     *
     * @param c Array of generalized concentrations. The 
     *          units depend upon the implementation of the
     *          reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the activity (i.e., generalized) concentration for use
     *
     * For the time being, we will use the concentration of pure
     * solvent for the the standard concentration of all species.
     * This has the effect of making mass-action reaction rates
     * based on the molality of species proportional to the
     * molality of the species.
     *
     * @param k Optional parameter indicating the species. The default
     *         is to assume this refers to species 0.
     * @return 
     *   Returns the standard Concentration in units of 
     *   m<SUP>3</SUP> kmol<SUP>-1</SUP>.
     *
     * @param k Species index
     */
    virtual doublereal standardConcentration(int k=0) const;

    
    //! Returns the natural logarithm of the standard 
    //! concentration of the kth species
    /*!
     * @param k Species index
     */
    virtual doublereal logStandardConc(int k=0) const;

   
    //! Returns the units of the standard and generalized concentrations.
    /*!
     * Note they have the same units, as their
     * ratio is defined to be equal to the activity of the kth
     * species in the solution, which is unitless.
     *
     * This routine is used in print out applications where the
     * units are needed. Usually, MKS units are assumed throughout
     * the program and in the XML input files.
     *
     * The base %ThermoPhase class assigns the default quantities
     * of (kmol/m3) for all species.
     * Inherited classes are responsible for overriding the default 
     * values if necessary.
     *
     * @param uA Output vector containing the units
     *  uA[0] = kmol units - default  = 1
     *  uA[1] = m    units - default  = -nDim(), the number of spatial
     *                                dimensions in the Phase class.
     *  uA[2] = kg   units - default  = 0;
     *  uA[3] = Pa(pressure) units - default = 0;
     *  uA[4] = Temperature units - default = 0;
     *  uA[5] = time units - default = 0
     * @param k species index. Defaults to 0.
     * @param sizeUA output int containing the size of the vector.
     *        Currently, this is equal to 6.
     */
    virtual void getUnitsStandardConc(double *uA, int k = 0,
				      int sizeUA = 6) const;

    //! Get the array of non-dimensional activities at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     *
     * We resolve this function at this level by calling
     * on the activityConcentration function. However, 
     * derived classes may want to override this default
     * implementation.
     *
     * (note solvent is on molar scale).
     *
     * @param ac  Output vector of activities. Length: m_kk.
     */
    virtual void getActivities(doublereal* ac) const;

   
    //! Get the array of non-dimensional molality-based 
    //! activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     *  note solvent is on molar scale. The solvent molar
     *  based activity coefficient is returned.
     *
     * @param acMolality Vector of Molality-based activity coefficients
     *                   Length: m_kk
     */
    virtual void 
    getMolalityActivityCoefficients(doublereal* acMolality) const;

    //@}
    /// @name  Partial Molar Properties of the Solution -----------------
    //@{

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     *
     * This function returns a vector of chemical potentials of the 
     * species in solution.
     *
     * \f[
     *    \mu_k = \mu^{\triangle}_k(T,P) + R T ln(\gamma_k^{\triangle} m_k)
     * \f] 
     *
     * @param mu  Output vector of species chemical 
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const;
   
    //! Returns an array of partial molar enthalpies for the species
    //! in the mixture. Units (J/kmol)
    /*!
     * For this phase, the partial molar enthalpies are equal to the
     * standard state enthalpies modified by the derivative of the
     * molality-based activity coefficent wrt temperature
     *
     *  \f[
     * \bar h_k(T,P) = h^{\triangle}_k(T,P) - R T^2 \frac{d \ln(\gamma_k^\triangle)}{dT}
     * \f]
     * The solvent partial molar enthalpy is equal to 
     *  \f[
     * \bar h_o(T,P) = h^{o}_o(T,P) - R T^2 \frac{d \ln(a_o}{dT}
     * \f]
     *
     *
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

  
    //! Returns an array of partial molar entropies of the species in the
    //! solution. Units: J/kmol/K.
    /**
     * Maxwell's equations provide an insight in how to calculate this
     *   (p.215 Smith and Van Ness)
     *
     *      d(chemPot_i)/dT = -sbar_i
     *      
     * For this phase, the partial molar entropies are equal to the
     * SS species entropies plus the ideal solution contribution
     * plus complicated functions of the
     * temperature derivative of the activity coefficents.
     *
     *  \f[
     *     \bar s_k(T,P) =  \hat s^0_k(T) - R log(M0 * molality[k])
     * \f]
     * \f[
     *      \bar s_solvent(T,P) =  \hat s^0_solvent(T) 
     *                  - R ((xmolSolvent - 1.0) / xmolSolvent)
     * \f]
     *
     *  @param sbar    Output vector of species partial molar entropies.
     *                 Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;
     
    //! Return an array of partial molar volumes for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  For this solution, the partial molar volumes are functions
     *  of the pressure derivatives of the activity coefficients.
     *
     *  @param vbar   Output vector of speciar partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //! Return an array of partial molar heat capacities for the
    //! species in the mixture.  Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat 
     *                capacities at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;


    //@}

    /// @name  Properties of the Standard State of the Species
    //          in the Solution --
    //@{

     
    //! Get the array of chemical potentials at unit activity for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current
     * temperature and pressure of the solution
     *
     *  Get the standard state chemical potentials of the species.
     *  This is the array of chemical potentials at unit activity 
     *  \f$ \mu^0_k(T,P) \f$.
     *  Activity is molality based in this object.
     *  We define these here as the chemical potentials of the pure
     *  species at the temperature and pressure of the solution.
     *  This function is used in the evaluation of the 
     *  equilibrium constant Kc. Therefore, Kc will also depend
     *  on T and P. This is the norm for liquid and solid systems.
     *
     *  units = J / kmol
     *
     * @param mu      Output vector of chemical potentials. 
     *                Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu) const;

   
    //! Get the nondimensional Gibbs functions for the species
    //! in their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     *  The standard states of the solutes are on the unit molality basis.
     *  \f[
     *  \mu^{\triangle}_k(T,P) = \mu^{\triangle,ref}_k(T) + (P - P_{ref}) * V_k
     * \f]
     *
     *  where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
     * \f$ \mu^{\triangle,ref}_k(T)\f$ is the chemical potential of pure
     * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
     *
     * A real water model is used. Therefore, \f$ \mu^{o}_0(T,P) \f$ is a 
     * complicated function of temperature and pressure.
     *
     * @param grt  Output vector of nondimensional standard state gibbs free energies
     *             Length: m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    //! Get the Gibbs functions for the standard
    //! state of the species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     *  The standard states are on the unit molality basis.
     * Units are Joules/kmol
     * @param gpure  Output vector of  standard state gibbs free energies
     *               Length: m_kk.
     */
    virtual void getPureGibbs(doublereal* gpure) const;

 
    //! Get the nondimensional Enthalpy functions for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     *  The standard states are on the unit molality basis.
     * We assume an incompressible constant partial molar
     * volume for the solutes.
     *
     *  \f[
     *      h^{\triangle}_k(T,P) = h^{\triangle,ref}_k(T) + (P - P_{ref}) * V_k
     *  \f]
     *
     * where \f$V_k\f$ is the molar volume of SS species <I>k</I>.
     * \f$ h^{ref}_k(T)\f$ is the enthalpy of the SS
     * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
     *
     * The solvent water enthalpy is obtained from a pure water
     * equation of state model.
     *
     * @param hrt      Output vector of nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    //! Get the array of nondimensional Entropy functions for the
    //! standard state species at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     *
     *  The standard states are on the unit molality basis.
     *
     *  \f[
     *      s^{\triangle}_k(T,P) = s^{\triangle,ref}_k(T) 
     *  \f]
     *
     * Note, this is equal to the reference state entropies
     * due to the zero volume expansivity:
     * i.e., (dS/dp)_T = (dV/dT)_P = 0.0
     *
     * The solvent water entropy is obtained from a pure water
     * equation of state model.
     *
     * @param sr   Output vector of  nondimensional standard state entropies.
     *             Length: m_kk. The solvent water is species 0, always.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    //! Get the nondimensional Heat Capacities at constant
    //! pressure for the species standard states
    //! at the current <I>T</I> and <I>P</I> of the solution
    /*!
     *  The standard states are on the unit molality basis.
     * For the solutes:
     * \f[
     *  Cp^\triangle_k(T,P) = Cp^{\triangle,ref}_k(T)
     * \f]
     *
     * \f$ Cp^{ref}_k(T)\f$ is the constant pressure heat capacity
     * of species <I>k</I> at the reference pressure, \f$p_{ref}\f$.
     *
     * The solute heat capacity is obtained from a pure water
     * equation of state model, so it depends on T and P.
     *
     * @param cpr Vector of length m_kk, which on return cpr[k]
     *           will contain the nondimensional 
     *           constant pressure heat capacity for species k.
     */
    virtual void getCp_R(doublereal* cpr) const;

 
    //!  Get the molar volumes of the species standard states at the current
    //!  <I>T</I> and <I>P</I> of the solution.
    /*!
     * The current model assumes that an incompressible molar volume for
     * all solutes. The molar volume for the water solvent, however,
     * is obtained from a pure water equation of state, waterSS.
     * Therefore, the water standard state varies with both T and P.
     * It is an error to request the water molar volume at a T and P
     * where the water phase is not stable phase.
     *
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk. The solvent water is species 0, always.
     */
    virtual void getStandardVolumes(doublereal *vol) const;

    //!  Returns the vector of nondimensional
    //!  Gibbs Free Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param grt     Output vector containing the nondimensional reference state
     *                Gibbs Free energies.  Length: m_kk.
     */
    virtual void getGibbs_RT_ref(doublereal *grt) const;

    //!  Returns the vector of nondimensional
    //!  enthalpies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param hrt     Output vector containing the nondimensional reference state enthalpies
     *                Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal *hrt) const;

    /*!
     *  Returns the vector of nondimensional
     *  entropies of the reference state at the current temperature
     *  of the solution and the reference pressure for each species.
     *
     * @param er      Output vector containing the nondimensional reference state
     *                entropies.  Length: m_kk.
       */
    virtual void getEntropy_R_ref(doublereal *er) const;

    //!  Returns the vector of nondimensional
    //!  constant pressure heat capacities of the reference state
    //!  at the current temperature of the solution
    //!  and reference pressure for each species.
    /*!
     *
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void getCp_R_ref(doublereal *cprt) const;

    //!  Get the molar volumes of the species reference states at the current
    //!  <I>T</I> and <I>P_ref</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal *vol) const;

  protected:

    //! Updates the standard state thermodynamic functions at the current T and P of the solution.
    /*!
     * @internal
     *
     * This function gets called for every call to a public function in this
     * class. It checks to see whether the temperature or pressure has changed and
     * thus whether the ss thermodynamics functions must be recalculated.
     *
     * @param pres  Pressure at which to evaluate the standard states.
     *              The default, indicated by a -1.0, is to use the current pressure
     */                      
    virtual void _updateStandardStateThermo(doublereal pres = -1.0) const;

    //@}
    /// @name Thermodynamic Values for the Species Reference States ---
    //@{


    ///////////////////////////////////////////////////////
    //
    //  The methods below are not virtual, and should not
    //  be overloaded.
    //
    //////////////////////////////////////////////////////

    /**
     * @name Specific Properties
     * @{
     */


    /**
     * @name Setting the State
     *
     * These methods set all or part of the thermodynamic
     * state.
     * @{
     */

    //@}

    /**
     * @name Chemical Equilibrium
     * Chemical equilibrium.
     * @{
     */
  public:

    //!This method is used by the ChemEquil equilibrium solver.
    /*!
     * It sets the state such that the chemical potentials satisfy
     * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) \f] where 
     * \f$ \lambda_m \f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     *
     * @param lambda_RT Input vector of dimensionless element potentials
     *                  The length is equal to nElements().
     */ 
    virtual void setToEquilState(const doublereal* lambda_RT) {
      err("setToEquilState");
    }

    //@}



    //! Set the equation of state parameters
    /*!
     * @internal
     *  The number and meaning of these depends on the subclass. 
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     */
    virtual void setParameters(int n, doublereal* c);

    //! Get the equation of state parameters in a vector
    /*!
     * @internal
     * The number and meaning of these depends on the subclass. 
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     */
    virtual void getParameters(int &n, doublereal * const c) const;

    /**
     * Set equation of state parameter values from XML
     * entries. This method is called by function importPhase in
     * file importCTML.cpp when processing a phase definition in
     * an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. 
     *   
     * @param eosdata An XML_Node object corresponding to
     * the "thermo" entry for this phase in the input file.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);
 
    //---------------------------------------------------------
    /// @name Critical state properties.
    /// These methods are only implemented by some subclasses.
        
    //@{
        
    /// Critical temperature (K).
    virtual doublereal critTemperature() const {
      err("critTemperature"); return -1.0;
    }
        
    /// Critical pressure (Pa).
    virtual doublereal critPressure() const {
      err("critPressure"); return -1.0;
    }
        
    /// Critical density (kg/m3).
    virtual doublereal critDensity() const {
      err("critDensity"); return -1.0;
    }                
        
    //@}
        
    /// @name Saturation properties.
    /// These methods are only implemented by subclasses that 
    /// implement full liquid-vapor equations of state.
    ///
    virtual doublereal satTemperature(doublereal p) const {
      err("satTemperature"); return -1.0;
    }
          
    //! Get the saturation pressure for a given temperature. 
    /*!
     * Note the limitations of this function. Stability considerations
     * concernting multiphase equilibrium are ignored in this 
     * calculation. Therefore, the call is made directly to the SS of 
     * water underneath. The object is put back into its original
     * state at the end of the call.
     *
     * @todo This is probably not implemented correctly. The stability
     *       of the salt should be added into this calculation. The
     *       underlying water model may be called to get the stability
     *       of the pure water solution, if needed.
     *
     * @param T  Temperature (kelvin)
     */
    virtual doublereal satPressure(doublereal T) const;
        
    virtual doublereal vaporFraction() const {
      err("vaprFraction"); return -1.0;
    }
        
    virtual void setState_Tsat(doublereal t, doublereal x) {
      err("setState_sat"); 
    }

    virtual void setState_Psat(doublereal p, doublereal x) {
      err("setState_sat"); 
    }

    //@}


    /*
     *  -------------- Utilities -------------------------------
     */

    /**
     * Return a reference to the species thermodynamic property
     * manager. 
     *
     * @todo This method will fail if no species thermo
     * manager has been installed.
     */
    SpeciesThermo& speciesThermo() { return *m_spthermo; }

    //! Initialization of a HMWSoln phase using an xml file
    /*!
     * This routine is a precursor to initThermo(XML_Node*)
     * routine, which does most of the work.
     *
     * @param inputFile XML file containing the description of the
     *        phase
     *
     * @param id  Optional parameter identifying the name of the
     *            phase. If none is given, the first XML
     *            phase element will be used.
     */
    void constructPhaseFile(std::string inputFile, std::string id);

    //!   Import and initialize a HMWSoln phase 
    //!   specification in an XML tree into the current object.
    /*!
     *   Here we read an XML description of the phase.
     *   We import descriptions of the elements that make up the
     *   species in a phase.
     *   We import information about the species, including their
     *   reference state thermodynamic polynomials. We then freeze
     *   the state of the species.
     *
     *   Then, we read the species molar volumes from the xml 
     *   tree to finish the initialization.
     *
     * @param phaseNode This object must be the phase node of a
     *             complete XML tree
     *             description of the phase, including all of the
     *             species data. In other words while "phase" must
     *             point to an XML phase object, it must have
     *             sibling nodes "speciesData" that describe
     *             the species in the phase.
     *
     * @param id   ID of the phase. If nonnull, a check is done
     *             to see if phaseNode is pointing to the phase
     *             with the correct id. 
     */
    void constructPhaseXML(XML_Node& phaseNode, std::string id);

    /**
     * @internal Initialize. This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase.
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();

    /*
     * initThermoXML()                 (virtual from ThermoPhase)
     *
     *  This gets called from importPhase(). It processes the XML file
     *  after the species are set up. This is the main routine for
     *  reading in activity coefficient parameters.
     *
     * @param phaseNode This object must be the phase node of a
     *             complete XML tree
     *             description of the phase, including all of the
     *             species data. In other words while "phase" must
     *             point to an XML phase object, it must have
     *             sibling nodes "speciesData" that describe
     *             the species in the phase.
     * @param id   ID of the phase. If nonnull, a check is done
     *             to see if phaseNode is pointing to the phase
     *             with the correct id.
     */
    virtual void initThermoXML(XML_Node& phaseNode, std::string id);

    /**
     * Report the molar volume of species k
     *
     * units - \f$ m^3 kmol^-1 \f$
     *
     * @param k species index
     *
     * @deprecated
     *   The getPartialMolarVolumes() expression is more precise.
     */
    double speciesMolarVolume(int k) const;

    /**
     * Value of the Debye Huckel constant as a function of temperature
     * and pressure.
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     *
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calcualtion
     *                    or -1 to indicate the current pressure
     */
    virtual double A_Debye_TP(double temperature = -1.0, 
			      double pressure = -1.0) const;

    /**
     * Value of the derivative of the Debye Huckel constant with 
     * respect to temperature as a function of temperature
     * and pressure.
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     *
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calcualtion
     *                    or -1 to indicate the current pressure
     */
    virtual double dA_DebyedT_TP(double temperature = -1.0, 
				 double pressure = -1.0) const;

    /**
     * Value of the derivative of the Debye Huckel constant with 
     * respect to pressure, as a function of temperature
     * and pressure.
     *
     *      A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *  Units = sqrt(kg/gmol)
     *
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calcualtion
     *                    or -1 to indicate the current pressure
     */
    virtual double dA_DebyedP_TP(double temperature = -1.0, 
				 double pressure = -1.0) const;

    /**
     * Return Pitzer's definition of A_L. This is basically the
     * derivative of the A_phi multiplied by 4 R T**2
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *            dA_phidT = d(A_Debye)/dT / 3.0
     *            A_L = dA_phidT * (4 * R * T * T)
     *
     *            Units = sqrt(kg/gmol) (RT)
     * 
     *
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calcualtion
     *                    or -1 to indicate the current pressure
     */
    double ADebye_L(double temperature = -1.0,
		    double pressure = -1.0) const;


    /**
     * Return Pitzer's definition of A_J. This is basically the
     * temperature derivative of A_L, and the second derivative
     * of A_phi
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *            dA_phidT = d(A_Debye)/dT / 3.0
     *            A_J = 2 A_L/T + 4 * R * T * T * d2(A_phi)/dT2
     *
     *            Units = sqrt(kg/gmol) (R)
     *
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calcualtion
     *                    or -1 to indicate the current pressure
     */
    double ADebye_J(double temperature = -1.0,
		    double pressure = -1.0) const;
    /**
     * Return Pitzer's definition of A_V. This is the
     * derivative wrt pressure of A_phi multiplied by - 4 R T
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *            dA_phidT = d(A_Debye)/dP / 3.0
     *            A_V = - dA_phidP * (4 * R * T)
     *
     *            Units = sqrt(kg/gmol) (RT) / Pascal
     *
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calcualtion
     *                    or -1 to indicate the current pressure
     * 
     */
    double ADebye_V(double temperature = -1.0,
		    double pressure = -1.0) const;
    
    //! Value of the 2nd derivative of the Debye Huckel constant with 
    //! respect to temperature as a function of temperature
    //! and pressure.
    /*!
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     *
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calcualtion
     *                    or -1 to indicate the current pressure
     */
    virtual double d2A_DebyedT2_TP(double temperature = -1.0, 
				   double pressure = -1.0) const;

    //! Reports the ionic radius of the kth species
    /*!
     *  @param k Species index 
     */
    double AionicRadius(int k = 0) const;

    /**
     *
     * formPitzer():
     *
     *  Returns the form of the Pitzer parameterization used
     */
    int formPitzer() const { return m_formPitzer; }

    /**
     * Print out all of the input coefficients.
     */
    void printCoeffs () const;


    //@}
         
  protected:

    /**
     * This is the form of the Pitzer parameterization
     * used in this model.
     * The options are described at the top of this document,
     * and in the general documentation.
     * The list is repeated here:
     *
     * PITZERFORM_BASE  = 0    (only one supported atm)
     *
     */
    int m_formPitzer;

    /**
     * This is the form of the temperature dependence of Pitzer
     * parameterization used in the model.
     *
     *       PITZER_TEMP_CONSTANT   0
     *       PITZER_TEMP_LINEAR     1
     *       PITZER_TEMP_COMPLEX1   2
     */
    int m_formPitzerTemp;

    /**
     * Format for the generalized concentration:
     *
     *  0 = unity
     *  1 = molar_volume
     *  2 = solvent_volume    (default)
     *
     * The generalized concentrations can have three different forms
     * depending on the value of the member attribute m_formGC, which
     * is supplied in the constructor.
     *                          <TABLE>
     *  <TR><TD> m_formGC </TD><TD> GeneralizedConc </TD><TD> StandardConc </TD></TR>
     *  <TR><TD> 0        </TD><TD> X_k             </TD><TD> 1.0          </TD></TR>
     *  <TR><TD> 1        </TD><TD> X_k / V_k       </TD><TD> 1.0 / V_k    </TD></TR>
     *  <TR><TD> 2        </TD><TD> X_k / V_N       </TD><TD> 1.0 / V_N    </TD></TR>
     *                         </TABLE>
     *
     * The value and form of the generalized concentration will affect
     * reaction rate constants involving species in this phase.
     *
     * (HKM Note: Using option #1 may lead to spurious results and
     *  has been included only with warnings. The reason is that it
     *  molar volumes of electrolytes may often be negative. The
     *  molar volume of H+ is defined to be zero too. Either options
     *  0 or 2 are the appropriate choice. Option 0 leads to 
     *  bulk reaction rate constants which have units of s-1.
     *  Option 2 leads to bulk reaction rate constants for 
     *  bimolecular rxns which have units of m-3 kmol-1 s-1.)
     */
    int m_formGC;

    //! Vector containing the electrolyte species type
    /*!
     * The possible types are:
     *  - solvent
     *  - Charged Species
     *  - weakAcidAssociated
     *  - strongAcidAssociated
     *  - polarNeutral
     *  - nonpolarNeutral  .
     */
    vector_int  m_electrolyteSpeciesType;

    /**
     * Species molar volumes \f$ m^3 kmol^-1 \f$
     *  -> m_speciesSize in Constituents.h
     */
    //array_fp m_speciesMolarVolume;

    /**
     *  a_k = Size of the ionic species in the DH formulation
     *        units = meters
     */
    array_fp m_Aionic;

    /**
     * Current value of the ionic strength on the molality scale
     * Associated Salts, if present in the mechanism,
     * don't contribute to the value of the ionic strength
     * in this version of the Ionic strength.
     */
    mutable double m_IionicMolality;

    /**
     * Maximum value of the ionic strength allowed in the
     * calculation of the activity coefficients.
     */
    double m_maxIionicStrength;

    /**
     * Reference Temperature for the Pitzer formulations.
     */
    double m_TempPitzerRef;

  protected:
    /**
     * Stoichiometric ionic strength on the molality scale.
     * This differs from m_IionicMolality in the sense that
     * associated salts are treated as unassociated salts,
     * when calculating the Ionic strength by this method.
     */
    mutable double m_IionicMolalityStoich;

  public:
    /**
     * Form of the constant outside the Debye-Huckel term
     * called A. It's normally a function of temperature 
     * and pressure. However, it can be set from the
     * input file in order to aid in numerical comparisons.
     * Acceptable forms:
     *
     *       A_DEBYE_CONST  0
     *       A_DEBYE_WATER  1
     *
     * The A_DEBYE_WATER form may be used for water solvents
     * with needs to cover varying temperatures and pressures.
     * Note, the dielectric constant of water is a relatively
     * strong function of T, and its variability must be
     * accounted for,
     */
    int m_form_A_Debye;

  protected:
    /**
     * A_Debye -> this expression appears on the top of the
     *            ln actCoeff term in the general Debye-Huckel
     *            expression
     *            It depends on temperature. And, therefore,
     *            most be recalculated whenever T or P changes.
     *            This variable is a local copy of the calculation.
     *            
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *                 where B_Debye = F / sqrt(epsilon R T/2) 
     *                                 (dw/1000)^(1/2)
     *
     *            A_Debye = (1/ (8 Pi)) (2 Pi * Na * dw/1000)^(1/2)
     *                       (e * e / (epsilon * kb * T))^(3/2) 
     *
     *            Units = sqrt(kg/gmol)
     *
     *            Nominal value = 1.172576 sqrt(kg/gmol)
     *                  based on:
     *                    epsilon/epsilon_0 = 78.54
     *                           (water at 25C)
     *                    epsilon_0 = 8.854187817E-12 C2 N-1 m-2
     *                    e = 1.60217653 E-19 C
     *                    F = 9.6485309E7 C kmol-1
     *                    R = 8.314472E3 kg m2 s-2 kmol-1 K-1
     *                    T = 298.15 K
     *                    B_Debye = 3.28640E9 sqrt(kg/gmol)/m
     *                    dw = C_0 * M_0 (density of water) (kg/m3)
     *                       = 1.0E3 at 25C
     */
    mutable double m_A_Debye;

    
    //!  Water standard state calculator
    /*!
     *  derived from the equation of state for water. 
     */
    WaterPDSS *m_waterSS;

    //! density of standard-state water
    /*!
     * internal temporary variable
     */
    double m_densWaterSS;

    /**
     *  Pointer to the water property calculator
     */
    WaterProps *m_waterProps;

    /**
     * Vector containing the species reference exp(-G/RT) functions
     * at T = m_tlast
     */
    mutable vector_fp      m_expg0_RT;

    /**
     * Vector of potential energies for the species.
     */
    mutable vector_fp      m_pe;

    /**
     * Temporary array used in equilibrium calculations
     */
    mutable vector_fp      m_pp;

    /**
     * vector of size m_kk, used as a temporary holding area.
     */
    mutable vector_fp      m_tmpV;

    /**
     * Stoichiometric species charge -> This is for calculations
     * of the ionic strength which ignore ion-ion pairing into
     * neutral molecules. The Stoichiometric species charge is the
     * charge of one of the ion that would occur if the species broke
     * into two charged ion pairs.
     *  NaCl ->   m_speciesCharge_Stoich = -1;
     *  HSO4- -> H+ + SO42-              = -2
     *      -> The other charge is calculated.
     * For species that aren't ion pairs, its equal to the
     * m_speciesCharge[] value.
     */
    vector_fp  m_speciesCharge_Stoich;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Beta0_ij[i][j] is the value of the Beta0 coefficient
     *  for the ij salt. It will be nonzero iff i and j are
     *  both charged and have opposite sign. The array is also
     *  symmetric.
     *     counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     */
    mutable vector_fp  m_Beta0MX_ij;

    //! Derivative of Beta0_ij[i][j] wrt T
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp  m_Beta0MX_ij_L;

    //! Derivative of Beta0_ij[i][j] wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp  m_Beta0MX_ij_LL;

    //! Derivative of Beta0_ij[i][j] wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp  m_Beta0MX_ij_P;

    //! Array of coefficients for Beta0, a variable in Pitzer's papers
    /*!
     *  column index is counterIJ
     *  m_Beta0MX_ij_coeff.ptrColumn(counterIJ) is a double* containing
     *  the vector of coefficients for the counterIJ interaction.
     */
    mutable Array2D    m_Beta0MX_ij_coeff;

    /*!
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Beta1_ij[i][j] is the value of the Beta1 coefficient
     *  for the ij salt. It will be nonzero iff i and j are
     *  both charged and have opposite sign. The array is also
     *  symmetric. 
     *     counterIJ where counterIJ = m_counterIJ[i][j]
     *     is used to access this array.
     */
    mutable vector_fp m_Beta1MX_ij;

    //! Derivative of Beta1_ij[i][j] wrt T
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Beta1MX_ij_L;

    //! Derivative of Beta1_ij[i][j] wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Beta1MX_ij_LL;

    //! Derivative of Beta1_ij[i][j] wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Beta1MX_ij_P;

    //! Array of coefficients for Beta1, a variable in Pitzer's papers
    /*!
     * column index is counterIJ
     *  m_Beta1MX_ij_coeff.ptrColumn(counterIJ) is a double* containing
     *  the vector of coefficients for the counterIJ interaction.
     */
    mutable Array2D   m_Beta1MX_ij_coeff;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Beta2_ij[i][j] is the value of the Beta2 coefficient
     *  for the ij salt. It will be nonzero iff i and j are
     *  both charged and have opposite sign, and i and j
     *  both have charges of 2 or more. The array is also
     *  symmetric.
     *  counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     */
    vector_fp m_Beta2MX_ij;

    //! Derivative of Beta2_ij[i][j] wrt T
    /*!
     *  vector index is counterIJ
     */
    vector_fp m_Beta2MX_ij_L;

    //! Derivative of Beta2_ij[i][j] wrt TT
    /*!
     *  vector index is counterIJ
     */
    vector_fp m_Beta2MX_ij_LL;

    //! Derivative of Beta2_ij[i][j] wrt P
    /*!
     *  vector index is counterIJ
     */
    vector_fp m_Beta2MX_ij_P;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Alpha1MX_ij[i][j] is the value of the alpha1 coefficient
     *  for the ij interaction. It will be nonzero iff i and j are
     *  both charged and have opposite sign, and i and j
     *  both have charges of 2 or more. The array is also
     *  symmetric.
     *  counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     */
    vector_fp m_Alpha1MX_ij;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  CphiMX_ij[i][j] is the value of the Cphi coefficient
     *  for the ij interaction. It will be nonzero iff i and j are
     *  both charged and have opposite sign, and i and j
     *  both have charges of 2 or more. The array is also
     *  symmetric.
     *  counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     */
    mutable vector_fp m_CphiMX_ij;

    //! Derivative of Cphi_ij[i][j] wrt T
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_CphiMX_ij_L;

    //! Derivative of Cphi_ij[i][j] wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_CphiMX_ij_LL;

    //! Derivative of Cphi_ij[i][j] wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_CphiMX_ij_P;

    //! Array of coefficients for Beta1, a variable in Pitzer's papers
    /*!
     *  column index is counterIJ
     *  m_CphiMX_ij_coeff.ptrColumn(counterIJ) is a double* containing
     *  the vector of coefficients for the counterIJ interaction.
     */
    mutable Array2D   m_CphiMX_ij_coeff;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Theta_ij[i][j] is the value of the theta coefficient
     *  for the ij interaction. It will be nonzero for charged
     *  ions with the same sign. It is symmetric.
     *  counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     *
     *  HKM Recent Pitzer papers have used a functional form
     *      for Theta_ij, which depends on the ionic strength.
     */
    vector_fp m_Theta_ij;

    //! Derivative of Theta_ij[i][j] wrt T
    /*!
     *  vector index is counterIJ
     */
    vector_fp m_Theta_ij_L;

    //! Derivative of Theta_ij[i][j] wrt TT
    /*!
     *  vector index is counterIJ
     */
    vector_fp m_Theta_ij_LL;

    //! Derivative of Theta_ij[i][j] wrt P
    /*!
     *  vector index is counterIJ
     */
    vector_fp m_Theta_ij_P;

    /**
     * Array of 3D data sed in the Pitzer/HMW formulation.
     * Psi_ijk[n] is the value of the psi coefficient for the
     * ijk interaction where
     *
     *   n = k + j * m_kk + i * m_kk * m_kk;
     *
     * It is potentially nonzero everywhere. 
     * The first two coordinates are symmetric wrt cations,
     * and the last two coordinates are symmetric wrt anions.
     */
    vector_fp m_Psi_ijk;

    //! Derivitive of Psi_ijk[n] wrt T
    /*! 
     *  see m_Psi_ijk for reference on the indexing into this variable. 
     */
    vector_fp m_Psi_ijk_L;

    //! Derivitive of Psi_ijk[n] wrt TT
    /*! 
     *  see m_Psi_ijk for reference on the indexing into this variable. 
     */
    vector_fp m_Psi_ijk_LL;

    //! Derivitive of Psi_ijk[n] wrt P
    /*! 
     *  see m_Psi_ijk for reference on the indexing into this variable. 
     */
    vector_fp m_Psi_ijk_P;

    //! Lambda coefficient for the ij interaction
    /*!
     * Array of 2D data used in the Pitzer/HMW formulation.
     * Lambda_ij[i][j] represents the lambda coefficient for the
     * ij interaction. This is a general interaction representing
     * neutral species. The neutral species occupy the first
     * index, i.e., i. The charged species occupy the j coordinate.
     * neutral, neutral interactions are also included here.
     */
    Array2D   m_Lambda_ij;

    //! Derivative of  Lambda_ij[i][j] wrt T. see m_Lambda_ij
    Array2D   m_Lambda_ij_L;

    //! Derivative of  Lambda_ij[i][j] wrt TT
    Array2D   m_Lambda_ij_LL;

    //! Derivative of  Lambda_ij[i][j] wrt P
    Array2D   m_Lambda_ij_P;

    
    //!  Logarithm of the activity coefficients on the molality
    //!  scale.
    /*!
     *       mutable because we change this if the composition
     *       or temperature or pressure changes.
     *
     *  index is the species index
     */
    mutable vector_fp m_lnActCoeffMolal;

    //!  Derivative of the Logarithm of the activity coefficients on the molality
    //!  scale wrt T
    /*!
     *  index is the species index
     */
    mutable vector_fp m_dlnActCoeffMolaldT;

    //!  Derivative of the Logarithm of the activity coefficients on the molality
    //!  scale wrt TT
    /*!
     *  index is the species index
     */
    mutable vector_fp m_d2lnActCoeffMolaldT2;

    //!  Derivative of the Logarithm of the activity coefficients on the molality
    //!  scale wrt P
    /*!
     *  index is the species index
     */
    mutable vector_fp m_dlnActCoeffMolaldP;

    /*
     * -------- Temporary Variables Used in the Activity Coeff Calc
     */

    
    //! a counter variable for keeping track of symmetric binary
    //! interactions amongst the solute species.
    /*!
     * n = m_kk*i + j
     * m_CounterIJ[n] = counterIJ
     */
    mutable array_int m_CounterIJ;

    /**
     *  This is elambda, MEC
     */
    mutable double elambda[17];

    /**
     *  This is elambda1, MEC
     */
    mutable double elambda1[17];

    /**
     *  Various temporary arrays used in the calculation of
     *  the Pitzer activity coefficents.
     *  The subscript, L, denotes the same quantity's derivative
     *  wrt temperature
     */

    //! This is the value of g(x) in Pitzer's papers
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_gfunc_IJ;

    //! hfunc, was called gprime in Pitzer's paper. However,
    //! it's not the derivative of gfunc(x), so I renamed it.
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_hfunc_IJ;

    //! Intermediate variable called BMX in Pitzer's paper
    //! This is the basic cation - anion interaction
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BMX_IJ;

    //! Derivative of BMX_IJ wrt T
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BMX_IJ_L;

    //! Derivative of BMX_IJ wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BMX_IJ_LL;

    //! Derivative of BMX_IJ wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BMX_IJ_P;

    //! Intermediate variable called BprimeMX in Pitzer's paper
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BprimeMX_IJ;

    //! Derivative of BprimeMX wrt T 
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BprimeMX_IJ_L;

    //! Derivative of BprimeMX wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BprimeMX_IJ_LL;

    //! Derivative of BprimeMX wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BprimeMX_IJ_P;

    //! Intermediate variable called BphiMX in Pitzer's paper
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BphiMX_IJ;

    //! Derivative of BphiMX_IJ wrt T
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BphiMX_IJ_L;

    //! Derivative of BphiMX_IJ wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BphiMX_IJ_LL;

    //! Derivative of BphiMX_IJ wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_BphiMX_IJ_P;

    //! Intermediate variable called Phi in Pitzer's paper
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Phi_IJ;

    //! Derivative of m_Phi_IJ wrt T
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Phi_IJ_L;

    //! Derivative of m_Phi_IJ wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Phi_IJ_LL;

    //! Derivative of m_Phi_IJ wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Phi_IJ_P;

    //! Intermediate variable called Phiprime in Pitzer's paper
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Phiprime_IJ;

    //! Intermediate variable called PhiPhi in Pitzer's paper
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_PhiPhi_IJ;

    //! Derivative of m_PhiPhi_IJ wrt T
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_PhiPhi_IJ_L;

    //! Derivative of m_PhiPhi_IJ wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_PhiPhi_IJ_LL;

    //! Derivative of m_PhiPhi_IJ wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_PhiPhi_IJ_P;

    //! Intermediate variable called CMX in Pitzer's paper
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_CMX_IJ;

    //! Derivative of m_CMX_IJ wrt T
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_CMX_IJ_L;

    //! Derivative of m_CMX_IJ wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_CMX_IJ_LL;

    //! Derivative of m_CMX_IJ wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_CMX_IJ_P;

    //! Intermediate storage of the activity coefficient itself
    /*!
     *  vector index is the species index
     */
    mutable vector_fp m_gamma;

  private:

    //! Local error routine
    /*!
     * @param msg print out a message and error exit
     */
    doublereal err(std::string msg) const;

    //!  Initialize all of the species - dependent lengths in the object
    void initLengths();

    /*
     * This function will be called to update the internally storred
     * natural logarithm of the molality activity coefficients 
     */
    void s_update_lnMolalityActCoeff() const;

  public:

    //!  Calculates the temperature derivative of the
    //!  natural logarithm of the molality activity coefficients.
    /*!
     * Public function makes sure that all dependent data is
     * up to date, before calling a private function
     */
    void s_Pitzer_dlnMolalityActCoeff_dT() const;

    //!  Calculates the Pressure derivative of the
    //!  natural logarithm of the molality activity coefficients.
    /*!
     * Public function makes sure that all dependent data is
     * up to date, before calling a private function
     */
    void s_Pitzer_dlnMolalityActCoeff_dP() const;

  private:
   
     //! This function calculates the temperature derivative of the
     //! natural logarithm of the molality activity coefficients.
     /*!
      * Private function does the work
      */
    void s_update_dlnMolalityActCoeff_dT() const;

    /**
     * This function calcultes the temperature second derivative
     * of the natural logarithm of the molality activity 
     * coefficients.
     */
    void s_update_d2lnMolalityActCoeff_dT2() const;
    /**
     * This function calculates the pressure derivative of the
     * natural logarithm of the molality activity coefficients.
     */
    void s_update_dlnMolalityActCoeff_dP() const;

    //! Calculates the Pitzer coefficients' dependence on the temperature.
    /*!
     *  It will also calculate the temperature
     *  derivatives of the coefficients, as they are important
     *  in the calculation of the latent heats and the
     *  heat capacities of the mixtures.
     *
     *   @param doDerivs If >= 1, then the routine will calculate
     *                  the first derivative. If >= 2, the 
     *                  routine will calculate the first and second
     *                  temperature derivative.
     *                  default = 2
     */
    void s_updatePitzerCoeffWRTemp(int doDerivs = 2) const;

    /**
     * This function does the main pitzer coefficient 
     * calculation
     */
    void s_updatePitzerSublnMolalityActCoeff() const;

    
    //! Calculate the lambda interactions.
    /*!
     *
     * Calculate E-lambda terms for charge combinations of like sign,
     *   using method of Pitzer (1975).
     *
     * @param is Ionic strength
     */
    void calc_lambdas(double is) const;

    /**
     *  Calculate etheta and etheta_prime
     *
     *  This interaction will be nonzero for species with the 
     *  same charge. this routine is not to be called for 
     *  neutral species; it core dumps or error exits.
     *
     * MEC implementation routine.
     *
     *  @param z1 charge of the first molecule
     *  @param z2 charge of the second molecule
     *  @param etheta return pointer containing etheta
     *  @param etheta_prime Return pointer containing etheta_prime.
     *
     *  This routine uses the internal variables, 
     *   elambda[] and elambda1[].
     *
     *  There is no prohibition against calling 
     *
     */
    void calc_thetas(int z1, int z2,
		     double *etheta, double *etheta_prime) const;

    //! Set up a counter variable for keeping track of symmetric binary
    //! interactactions amongst the solute species.
    /*!
     * The purpose of this is to squeeze the ij parameters into a
     * compressed single counter.
     *
     * n = m_kk*i + j 
     * m_Counter[n] = counter
     */
    void counterIJ_setup() const;

    //! Process an XML node called "binarySaltParameters"
    /*!
     * This node contains all of the parameters necessary to describe
     * the Pitzer model for that particular binary salt.
     * This function reads the XML file and writes the coefficients
     * it finds to an internal data structures.
     *
     * @param BinSalt  reference to the XML_Node named binarySaltParameters
     *                 containing the
     *                 anion - cation interaction
     */
    void readXMLBinarySalt(XML_Node &BinSalt);

    //! Process an XML node called "thetaAnion"
    /*!
     * This node contains all of the parameters necessary to describe
     * the binary interactions between two anions.
     *
     * @param BinSalt  reference to the XML_Node named thetaAnion
     *                 containing the
     *                 anion - anion interaction
     */
    void readXMLThetaAnion(XML_Node &BinSalt);

    //! Process an XML node called "thetaCation"
    /*!
     * This node contains all of the parameters necessary to describe
     * the binary interactions between two cations.
     *
     * @param BinSalt  reference to the XML_Node named thetaCation
     *                 containing the
     *                 cation - cation interaction
     */
    void readXMLThetaCation(XML_Node &BinSalt);

    //! Process an XML node called "psiCommonAnion"
    /*!
     * This node contains all of the parameters necessary to describe
     * the ternary interactions between one anion and two cations.
     *
     * @param BinSalt  reference to the XML_Node named psiCommonAnion
     *                 containing the
     *                 anion - cation1 - cation2 interaction
     */
    void readXMLPsiCommonAnion(XML_Node &BinSalt);

    //! Process an XML node called "psiCommonCation"
    /*!
     * This node contains all of the parameters necessary to describe
     * the ternary interactions between one cation and two anions.
     *
     * @param BinSalt  reference to the XML_Node named psiCommonCation
     *                 containing the
     *                 cation - anion1 - anion2 interaction
     */
    void readXMLPsiCommonCation(XML_Node &BinSalt);

    //! Process an XML node called "lambdaNeutral"
    /*!
     * This node contains all of the parameters necessary to describe
     * the binary interactions between one neutral species and
     * any other species (neutral or otherwise) in the mechanism.
     *
     * @param BinSalt  reference to the XML_Node named lambdaNeutral
     *                 containing multiple
     *                 Neutral - species interactions
     */
    void readXMLLambdaNeutral(XML_Node &BinSalt);

    //! utility function to assign an integer value from a string
    //! for the ElectrolyteSpeciesType field.
    /*!
     *  @param estString string name of the electrolyte species type
     */
    static int interp_est(std::string estString);

  public:
    /*!
     * Turn on copious debug printing when this
     * is true and DEBUG_MODE is turned on.
     */
    mutable int m_debugCalc;

    //! Return int specifying the amount of debug printing
    /*!
     *  This will return 0 if DEBUG_MODE is not turned on
     */
    int debugPrinting();
  };

}
        
#endif

