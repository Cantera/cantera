/**
 *  @file HMWSoln.h
 *    Headers for the %HMWSoln ThermoPhase object, which models concentrated
 *    electrolyte solutions
 *    (see \ref thermoprops and \link Cantera::HMWSoln HMWSoln \endlink) .
 *
 * Class %HMWSoln represents a concentrated liquid electrolyte phase which
 * obeys the Pitzer formulation for nonideality using molality-based
 * standard states.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_HMWSOLN_H
#define CT_HMWSOLN_H

#include "MolalityVPSSTP.h"
#include "electrolytes.h"

namespace Cantera
{

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
 *     All coefficients are assumed to have a complex functional
 *     based dependence wrt temperature;  See:
 *    (Silvester, Pitzer, J. Phys. Chem. 81, 19 1822 (1977)).
 *
 *       beta0 = q0 + q3(1/T - 1/Tr) + q4(ln(T/Tr)) +
 *               q1(T - Tr) + q2(T**2 - Tr**2)
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
class PDSS_Water;

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
 * to the object that handles the real water model, #Cantera::WaterPropsIAPWS.
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
 *          h^\triangle_k(T,P) = h^{\triangle,ref}_k(T)
 *              + \tilde{v}_k \left( P - P_{ref} \right)
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
 * The vector Phase::m_speciesSize[] is used to hold the
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
 * properties of ionic species contain meaningful information.
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
 * strength has units of sqrt( gmol kg<SUP>-1</SUP>).
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
 *  -  <B>cEST_strongAcidAssociated</B>   Species which always breaks apart into charged species.
 *                                        It may or may not be charged. Normally, these
 *                                        aren't included in the speciation vector.
 *  -  <B>cEST_polarNeutral </B>          Polar neutral species
 *  -  <B>cEST_nonpolarNeutral</B>        Non polar neutral species
 *
 *  Polar and non-polar neutral species are differentiated, because some additions
 *  to the activity
 *  coefficient expressions distinguish between these two types of solutes.
 *  This is the so-called salt-out effect.
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
 *  Much of the species electrolyte type information is inferred from other information in the
 *  input file. For example, as species which is charged is given the "chargedSpecies" default
 *  category. A neutral solute species is put into the "nonpolarNeutral" category by default.
 *
 *
 *  <H3> Specification of the Excess Gibbs Free Energy </H3>
 *
 *  Pitzer's formulation may best be represented as a specification of the excess gibbs
 *  free energy, \f$ G^{ex} \f$, defined as the deviation of the total gibbs free energy from
 *  that of an ideal molal solution.
 *     \f[
 *         G = G^{id} + G^{ex}
 *     \f]
 *
 *  The ideal molal solution contribution, not equal to an ideal solution contribution
 *  and in fact containing a singularity at the zero solvent mole fraction limit, is
 *  given below.
 *     \f[
 *         G^{id} = n_o \mu^o_o + \sum_{k\ne o} n_k \mu_k^{\triangle}
 *               + \tilde{M}_o n_o ( RT (\sum{m_i(\ln(m_i)-1)}))
 *     \f]
 *
 *  From the excess Gibbs free energy formulation, the activity coefficient expression
 *  and the osmotic coefficient expression for the solvent may be defined, by
 *  taking the appropriate derivatives. Using this approach guarantees that the
 *  entire system will obey the Gibbs-Duhem relations.
 *
 *  Pitzer employs the following general expression for the excess Gibbs free energy
 *
 *  \f[
 *    \begin{array}{cclc}
 *     \frac{G^{ex}}{\tilde{M}_o  n_o RT} &= &
 *          \left( \frac{4A_{Debye}I}{3b} \right) \ln(1 + b \sqrt{I})
 *        +   2 \sum_c \sum_a m_c m_a B_{ca}
 *        +     \sum_c \sum_a m_c m_a Z C_{ca}
 *          \\&&
 *        +   \sum_{c < c'} \sum m_c m_{c'} \left[ 2 \Phi_{c{c'}} + \sum_a m_a \Psi_{c{c'}a} \right]
 *        +   \sum_{a < a'} \sum m_a m_{a'} \left[ 2 \Phi_{a{a'}} + \sum_c m_c \Psi_{a{a'}c} \right]
 *          \\&&
 *        + 2 \sum_n \sum_c m_n m_c \lambda_{nc} + 2 \sum_n \sum_a m_n m_a \lambda_{na}
 *        + 2 \sum_{n < n'} \sum m_n m_{n'} \lambda_{n{n'}}
 *        +  \sum_n m^2_n \lambda_{nn}
 *    \end{array}
 *  \f]
 *
 *  <I>a</I> is a subscript over all anions, <I>c</I> is a subscript extending over all
 *  cations, and  <I>i</I> is a subscript that extends over all anions and cations.
 *  <I>n</I> is a subscript that extends only over neutral solute molecules.
 *  The second line contains cross terms where cations affect
 *  cations and/or cation/anion pairs,
 *  and anions affect anions or cation/anion pairs. Note part of the coefficients,
 *  \f$ \Phi_{c{c'}} \f$ and  \f$ \Phi_{a{a'}} \f$  stem from the theory
 *  of unsymmetrical mixing of electrolytes with different charges. This
 *  theory depends on the total ionic strength of the solution, and therefore,
 *  \f$ \Phi_{c{c'}} \f$ and  \f$ \Phi_{a{a'}} \f$  will depend on <I>I</I>, the
 *  ionic strength.  \f$ B_{ca}\f$ is a strong function of the
 *  total ionic strength,  <I>I</I>,
 *  of the electrolyte. The rest of the coefficients are assumed to be independent of the
 *  molalities or ionic strengths. However, all coefficients are potentially functions
 *  of the temperature and pressure of the solution.
 *
 *  <I>A</I> is the Debye-Huckel constant. Its specification is described in its own
 *  section below.
 *
 *  \f$ I\f$ is the ionic strength of the solution, and is given by:
 *
 *  \f[
 *      I = \frac{1}{2} \sum_k{m_k  z_k^2}
 *  \f]
 *
 *  In contrast to several other Debye-Huckel implementations (see \ref DebyeHuckel), the
 *  parameter \f$ b\f$ in the above equation is a constant that
 *  does not vary with respect to ion identity. This is an important simplification
 *  as it avoids troubles with satisfaction of the Gibbs-Duhem analysis.
 *
 *  The function \f$ Z \f$ is given by
 *
 *  \f[
 *      Z = \sum_i m_i \left| z_i \right|
 *  \f]
 *
 *  The value of \f$ B_{ca}\f$ is given by the following function
 *
 *  \f[
 *      B_{ca} = \beta^{(0)}_{ca} + \beta^{(1)}_{ca} g(\alpha^{(1)}_{ca} \sqrt{I})
 *             + \beta^{(2)}_{ca} g(\alpha^{(2)}_{ca} \sqrt{I})
 *  \f]
 *
 *  where
 *
 *  \f[
 *      g(x) = 2 \frac{(1 - (1 + x)\exp[-x])}{x^2}
 *  \f]
 *
 *  The formulation for \f$ B_{ca}\f$ combined with the formulation of the
 *  Debye-Huckel term in the eqn. for the excess Gibbs free energy stems
 *  essentially from an empirical fit to the ionic strength dependent data
 *  based over a wide sampling of binary electrolyte systems. \f$ C_{ca} \f$,
 *  \f$ \lambda_{nc} \f$, \f$ \lambda_{na} \f$, \f$ \lambda_{nn} \f$,
 *  \f$ \Psi_{c{c'}a} \f$, \f$ \Psi_{a{a'}c} \f$ are experimentally derived
 *  coefficients that may have pressure and/or temperature dependencies.
 *
 *  The \f$ \Phi_{c{c'}} \f$ and \f$ \Phi_{a{a'}} \f$ formulations are
 *  slightly more complicated. \f$ b \f$ is a universal
 *  constant defined to be equal to \f$ 1.2\ kg^{1/2}\ gmol^{-1/2} \f$. The exponential
 *  coefficient \f$ \alpha^{(1)}_{ca} \f$ is usually
 *  fixed at \f$ \alpha^{(1)}_{ca} = 2.0\ kg^{1/2} gmol^{-1/2}\f$
 *  except for 2-2 electrolytes, while other parameters were fit to experimental
 *  data. For 2-2 electrolytes, \f$ \alpha^{(1)}_{ca} = 1.4\ kg^{1/2}\ gmol^{-1/2}\f$
 *  is used in combination with either \f$ \alpha^{(2)}_{ca} = 12\ kg^{1/2}\ gmol^{-1/2}\f$
 *  or \f$ \alpha^{(2)}_{ca} = k A_\psi \f$, where <I>k</I> is a constant. For electrolytes other
 *  than 2-2 electrolytes the \f$ \beta^{(2)}_{ca} g(\alpha^{(2)}_{ca} \sqrt{I}) \f$  term
 *  is not used in the fitting procedure; it is only used for divalent metal
 *  solfates and other high-valence electrolytes which exhibit significant
 *  association at low ionic strengths.
 *
 *  The \f$ \beta^{(0)}_{ca} \f$,  \f$ \beta^{(1)}_{ca}\f$,  \f$ \beta^{(2)}_{ca} \f$,
 *  and \f$ C_{ca} \f$ binary coefficients are referred to as ion-interaction or
 *  Pitzer parameters. These Pitzer parameters may vary with temperature and pressure
 *  but they do not depend on the ionic strength. Their values and temperature
 *  derivatives of their values have been tabulated for a range of electrolytes
 *
 *  The \f$ \Phi_{c{c'}} \f$ and \f$ \Phi_{a{a'}} \f$ contributions, which
 *  capture cation-cation and anion-anion interactions, also have an
 *  ionic strength dependence.
 *
 *  Ternary contributions \f$ \Psi_{c{c'}a} \f$ and \f$ \Psi_{a{a'}c} \f$
 *  have been measured also for some systems. The success of the Pitzer
 *  method lies in its ability to model nonlinear activity coefficients
 *  of complex multicomponent systems with just binary and minor
 *  ternary contributions, which can be independently measured in
 *  binary or ternary subsystems.
 *
 *
 *  <H3> Multicomponent Activity Coefficients for Solutes </H3>
 *
 *  The formulas for activity coefficients of solutes may be obtained by taking the
 *  following derivative of the excess Gibbs Free Energy formulation described above:
 *
 *  \f[
 *     \ln(\gamma_k^\triangle) = \frac{d\left( \frac{G^{ex}}{M_o n_o RT} \right)}{d(m_k)}\Bigg|_{n_i}
 *  \f]
 *
 *  In the formulas below the following conventions are used. The subscript <I>M</I> refers
 *  to a particular cation. The subscript X refers to a particular anion, whose
 *  activity is being currently evaluated. the subscript <I>a</I> refers to a summation
 *  over all anions in the solution, while the subscript <I>c</I> refers to a summation
 *  over all cations in the solutions.
 *
 *  The activity coefficient for a particular cation <I>M</I> is given by
 *
 *  \f[
 *      \ln(\gamma_M^\triangle) = -z_M^2(F) + \sum_a m_a \left( 2 B_{Ma} + Z C_{Ma} \right)
 *      + z_M   \left( \sum_a  \sum_c m_a m_c C_{ca} \right)
 *             + \sum_c m_c \left[ 2 \Phi_{Mc} + \sum_a m_a \Psi_{Mca} \right]
 *             + \sum_{a < a'} \sum m_a m_{a'} \Psi_{Ma{a'}}
 *             +  2 \sum_n m_n \lambda_{nM}
 *  \f]
 *
 *  The activity coefficient for a particular anion <I>X</I> is given by
 *
 *  \f[
 *      \ln(\gamma_X^\triangle) = -z_X^2(F) + \sum_a m_c \left( 2 B_{cX} + Z C_{cX} \right)
 *      + \left|z_X \right|  \left( \sum_a  \sum_c m_a m_c C_{ca} \right)
 *             + \sum_a m_a \left[ 2 \Phi_{Xa} + \sum_c m_c \Psi_{cXa} \right]
 *             + \sum_{c < c'} \sum m_c m_{c'} \Psi_{c{c'}X}
 *             +  2 \sum_n m_n \lambda_{nM}
 *  \f]
 *  where the function \f$ F \f$ is given by
 *
 *  \f[
 *       F = - A_{\phi} \left[ \frac{\sqrt{I}}{1 + b \sqrt{I}}
 *                 + \frac{2}{b} \ln{\left(1 + b\sqrt{I}\right)} \right]
 *                 + \sum_a \sum_c m_a m_c B'_{ca}
 *                 + \sum_{c < c'} \sum m_c m_{c'} \Phi'_{c{c'}}
 *                 + \sum_{a < a'} \sum m_a m_{a'} \Phi'_{a{a'}}
 *  \f]
 *
 *  We have employed the definition of \f$ A_{\phi} \f$, also used by Pitzer
 *  which is equal to
 *
 *   \f[
 *     A_{\phi} = \frac{A_{Debye}}{3}
 *   \f]
 *
 *  In the above formulas, \f$ \Phi'_{c{c'}} \f$  and \f$ \Phi'_{a{a'}} \f$ are the
 *  ionic strength derivatives of \f$ \Phi_{c{c'}} \f$  and \f$  \Phi_{a{a'}} \f$,
 *  respectively.
 *
 *  The function \f$ B'_{MX} \f$ is defined as:
 *
 *  \f[
 *       B'_{MX} = \left( \frac{\beta^{(1)}_{MX} h(\alpha^{(1)}_{MX} \sqrt{I})}{I}  \right)
 *                 \left( \frac{\beta^{(2)}_{MX} h(\alpha^{(2)}_{MX} \sqrt{I})}{I}  \right)
 *  \f]
 *
 *  where \f$ h(x) \f$ is defined as
 *
 *   \f[
 *       h(x) = g'(x) \frac{x}{2} =
 *        \frac{2\left(1 - \left(1 + x + \frac{x^2}{2} \right)\exp(-x) \right)}{x^2}
 *   \f]
 *
 *  The activity coefficient for neutral species <I>N</I> is given by
 *
 *  \f[
 *      \ln(\gamma_N^\triangle) = 2 \left( \sum_i m_i \lambda_{iN}\right)
 *  \f]
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
 *  The osmotic coefficient may be related to the water activity by the following relation:
 *
 *  \f[
 *       \phi = - \frac{1}{\tilde{M}_o \sum_{i \neq o} m_i} \ln(a_o)
 *           = - \frac{n_o}{\sum_{i \neq o}n_i} \ln(a_o)
 *  \f]
 *
 *
 *  The result is the following
 *
 *  \f[
 *    \begin{array}{ccclc}
 *      \phi - 1 &= &
 *          \frac{2}{\sum_{i \ne 0} m_i}
 *           \bigg[ &
 *        -  A_{\phi} \frac{I^{3/2}}{1 + b \sqrt{I}}
 *        +   \sum_c  \sum_a m_c m_a \left( B^{\phi}_{ca} + Z C_{ca}\right)
 *          \\&&&
 *        +   \sum_{c < c'} \sum m_c m_{c'} \left[ \Phi^{\phi}_{c{c'}} + \sum_a m_a \Psi_{c{c'}a} \right]
 *        +   \sum_{a < a'} \sum m_a m_{a'} \left[ \Phi^{\phi}_{a{a'}} + \sum_c m_c \Psi_{a{a'}c} \right]
 *          \\&&&
 *        + \sum_n \sum_c m_n m_c \lambda_{nc} +  \sum_n \sum_a m_n m_a \lambda_{na}
 *        + \sum_{n < n'} \sum m_n m_{n'} \lambda_{n{n'}}
 *        + \frac{1}{2} \left( \sum_n m^2_n \lambda_{nn}\right)
 *          \bigg]
 *    \end{array}
 *  \f]
 *
 *  It can be shown that the expression
 *
 *  \f[
 *     B^{\phi}_{ca} = \beta^{(0)}_{ca} + \beta^{(1)}_{ca} \exp{(- \alpha^{(1)}_{ca} \sqrt{I})}
 *             + \beta^{(2)}_{ca} \exp{(- \alpha^{(2)}_{ca} \sqrt{I} )}
 *  \f]
 *
 *  is consistent with the expression \f$ B_{ca} \f$ in the \f$ G^{ex} \f$ expression
 *  after carrying out the derivative wrt \f$ m_M \f$.
 *
 *  Also taking into account that  \f$ {\Phi}_{c{c'}} \f$ and
 *  \f$ {\Phi}_{a{a'}} \f$ has an ionic strength dependence.
 *
 *  \f[
 *    \Phi^{\phi}_{c{c'}} = {\Phi}_{c{c'}} + I \frac{d{\Phi}_{c{c'}}}{dI}
 *  \f]
 *
 *  \f[
 *    \Phi^{\phi}_{a{a'}} = \Phi_{a{a'}} + I \frac{d\Phi_{a{a'}}}{dI}
 *  \f]
 *
 *
 *  <H3> Temperature and Pressure Dependence of the Pitzer Parameters </H3>
 *
 *  In general most of the coefficients introduced in the previous section may
 *  have a temperature and pressure dependence. The temperature and pressure
 *  dependence of these coefficients strongly influence the value of the
 *  excess Enthalpy and excess Volumes of Pitzer solutions. Therefore, these
 *  are readily measurable quantities.
 *  HMWSoln provides several
 *  different methods for putting these dependencies into the coefficients.
 *  HMWSoln has an implementation described by Silverter and Pitzer (1977),
 *  which was used to fit experimental data for NaCl over an extensive range,
 *  below the critical temperature of water.
 *  They found a temperature functional form for fitting the 3 following
 *  coefficients that describe the Pitzer parameterization for a single salt
 *  to be adequate to describe how the excess gibbs free energy values for
 *  the binary salt changes with respect to temperature.
 *  The following functional form
 *  was used to fit the temperature dependence of the Pitzer Coefficients
 *  for each cation - anion pair, M X.
 *
 *  \f[
 *      \beta^{(0)}_{MX} = q^{b0}_0
 *                       + q^{b0}_1 \left( T - T_r \right)
 *                       + q^{b0}_2 \left( T^2 - T_r^2 \right)
 *                       + q^{b0}_3 \left( \frac{1}{T} - \frac{1}{T_r}\right)
 *                       + q^{b0}_4 \ln \left( \frac{T}{T_r} \right)
 *  \f]
 *  \f[
 *      \beta^{(1)}_{MX} = q^{b1}_0  + q^{b1}_1 \left( T - T_r \right)
 *                       + q^{b1}_{2} \left( T^2 - T_r^2 \right)
 *  \f]
 *  \f[
 *     C^{\phi}_{MX} = q^{Cphi}_0
 *                   + q^{Cphi}_1 \left( T - T_r \right)
 *                   + q^{Cphi}_2 \left( T^2 - T_r^2 \right)
 *                   + q^{Cphi}_3 \left( \frac{1}{T} - \frac{1}{T_r}\right)
 *                   + q^{Cphi}_4 \ln \left( \frac{T}{T_r} \right)
 *  \f]
 *
 *  where
 *
 *  \f[
 *      C^{\phi}_{MX} =  2 {\left| z_M z_X \right|}^{1/2} C_{MX}
 *  \f]
 *
 *  In later papers, Pitzer has added additional temperature dependencies
 *  to all of the other remaining second and third order virial coefficients.
 *  Some of these dependencies are justified and motivated by theory.
 *  Therefore,
 *  a formalism wherein all of the coefficients in the base theory have
 *  temperature dependencies associated with them has been implemented
 *  within the  %HMWSoln object. Much of the formalism, however,
 *  has been unexercised.
 *
 *  In the %HMWSoln object, the temperature dependence of the Pitzer
 *  parameters are specified in the following way.
 *
 *    - PIZTER_TEMP_CONSTANT      - string name "CONSTANT"
 *      - Assumes that all coefficients are independent of temperature
 *        and pressure
 *    - PIZTER_TEMP_COMPLEX1     - string name "COMPLEX" or "COMPLEX1"
 *      - Uses the full temperature dependence for the
 *        \f$\beta^{(0)}_{MX} \f$ (5 coeffs),
 *        the  \f$\beta^{(1)}_{MX} \f$ (3 coeffs),
 *        and \f$ C^{\phi}_{MX} \f$ (5 coeffs) parameters described above.
 *    - PITZER_TEMP_LINEAR        - string name "LINEAR"
 *      - Uses just the temperature dependence for the
 *        \f$\beta^{(0)}_{MX} \f$, the \f$\beta^{(1)}_{MX} \f$,
 *        and \f$ C^{\phi}_{MX} \f$ coefficients described above.
 *        There are 2 coefficients for each term.
 *
 *  The temperature dependence is specified in an attributes field
 *  in the <TT>  activityCoefficients </TT>  XML block,
 *  called <TT> TempModel </TT>. Permissible values for that
 *  attribute are <TT> CONSTANT, COMPLEX1</TT>, and <TT> LINEAR.</TT>
 *
 *  The specification of the binary interaction between a cation and
 *  an anion is given by the coefficients, \f$ B_{MX}\f$ and
 *  \f$ C_{MX}\f$
 *  The specification of \f$ B_{MX}\f$ is a function of
 *  \f$\beta^{(0)}_{MX} \f$, \f$\beta^{(1)}_{MX} \f$,
 *  \f$\beta^{(2)}_{MX} \f$,  \f$\alpha^{(1)}_{MX} \f$, and
 *  \f$\alpha^{(2)}_{MX} \f$.
 *  \f$ C_{MX}\f$ is calculated from \f$C^{\phi}_{MX} \f$
 *  from the formula above.
 *  All of the underlying coefficients are specified in the
 *  XML element block <TT> binarySaltParameters </TT>, which
 *  has the attribute <TT> cation </TT> and <TT> anion </TT>
 *  to identify the interaction. XML elements named
 *  <TT> beta0, beta1, beta2, Cphi, Alpha1, Alpha2 </TT>
 *  within each <TT> binarySaltParameters </TT> block
 *  specify the parameters. Within each of these blocks
 *  multiple parameters describing temperature or pressure
 *  dependence are serially listed in the order that they
 *  appear in the equation in this document. An example of
 *  the <TT> beta0 </TT> block that fits the <TT> COMPLEX1 </TT> temperature
 *  dependence given above is
 *
 * @code
    <binarySaltParameters cation="Na+" anion="OH-">
      <beta0> q0, q1, q2, q3, q4  </beta0>
    </binarySaltParameters>
   @endcode
 *
 *  The parameters for \f$ \beta^{(0)}\f$ fit the following equation:
 *
 *  \f[
 *        \beta^{(0)} = q_0^{{\beta}0} + q_1^{{\beta}0} \left( T - T_r \right)
 *               + q_2^{{\beta}0} \left( T^2 - T_r^2 \right)
 *               + q_3^{{\beta}0} \left( \frac{1}{T} - \frac{1}{T_r} \right)
 *               + q_4^{{\beta}0} \ln \left( \frac{T}{T_r} \right)
 *  \f]
 *
 *  This same <TT> COMPLEX1 </TT> temperature
 *  dependence given above is used for the following parameters:
 *  \f$ \beta^{(0)}_{MX} \f$, \f$ \beta^{(1)}_{MX} \f$,
 *  \f$ \beta^{(2)}_{MX} \f$, \f$ \Theta_{cc'} \f$,  \f$\Theta_{aa'} \f$,
 *  \f$ \Psi_{c{c'}a} \f$ and \f$ \Psi_{ca{a'}} \f$.
 *
 *
 *  <H3> Like-Charged Binary Ion Parameters and the Mixing Parameters </H3>
 *
 *  The previous section contained the functions, \f$ \Phi_{c{c'}} \f$,
 *  \f$ \Phi_{a{a'}} \f$ and their derivatives wrt the
 *  ionic strength, \f$ \Phi'_{c{c'}} \f$ and \f$ \Phi'_{a{a'}} \f$.
 *  Part of these terms come from theory.
 *
 *  Since like charged ions repel each other and are generally
 *  not near each other, the virial coefficients for same-charged ions
 *  are small. However, Pitzer doesn't ignore these in his
 *  formulation. Relatively larger and longer range terms between
 *  like-charged ions exist however, which appear only for
 *  unsymmetrical mixing of same-sign charged ions with different
 *  charges. \f$ \Phi_{ij} \f$, where \f$ ij \f$ is either \f$ a{a'} \f$
 *  or \f$ c{c'} \f$ is given by
 *
 *  \f[
 *      {\Phi}_{ij} = \Theta_{ij} + \,^E \Theta_{ij}(I)
 *  \f]
 *
 *  \f$ \Theta_{ij} \f$ is the small virial coefficient expansion term.
 *  Dependent in general on temperature and pressure, its ionic
 *  strength dependence is ignored in Pitzer's approach.
 *  \f$ \,^E\Theta_{ij}(I) \f$ accounts for the electrostatic
 *  unsymmetrical mixing effects and is dependent only on the
 *  charges of the ions i, j, the total ionic strength and on
 *  the dielectric constant and density of the solvent.
 *  This seems to be a relatively well-documented part of the theory.
 *  They theory below comes from Pitzer summation (Pitzer) in the
 *  appendix. It's also mentioned in Bethke's book (Bethke), and
 *  the equations are summarized in Harvie & Weare (1980).
 *  Within the code, \f$ \,^E\Theta_{ij}(I) \f$ is evaluated according
 *  to the algorithm described in Appendix B [Pitzer] as
 *
 *  \f[
 *     \,^E\Theta_{ij}(I) = \left( \frac{z_i z_j}{4I} \right)
 *        \left( J(x_{ij})  - \frac{1}{2} J(x_{ii})
 *                          - \frac{1}{2} J(x_{jj})  \right)
 *  \f]
 *
 *  where \f$ x_{ij} = 6 z_i z_j A_{\phi} \sqrt{I} \f$ and
 *
 *   \f[
 *      J(x) = \frac{1}{x} \int_0^{\infty}{\left( 1 + q +
 *             \frac{1}{2} q^2 - e^q \right) y^2 dy}
 *   \f]
 *
 *  and \f$ q = - (\frac{x}{y}) e^{-y} \f$. \f$ J(x) \f$ is evaluated by
 *  numerical integration.
 *
 *  The \f$  \Theta_{ij} \f$ term is a constant that is specified
 *  by the XML element <TT> thetaCation </TT> and
 *  <TT> thetaAnion </TT>,  which
 *  has the attribute <TT> cation1 </TT>, <TT> cation2 </TT> and
 *  <TT> anion1 </TT>, <TT> anion2 </TT> respectively
 *  to identify the interaction. No temperature or
 *  pressure dependence of this parameter is currently allowed.
 *  An example of the block is presented below.
 *
 * @code
     <thetaCation cation1="Na+" cation2="H+">
                <Theta> 0.036 </Theta>
     </thetaCation>
   @endcode
 *
 *
 *  <H3> Ternary Pitzer Parameters </H3>
 *
 *  The \f$  \Psi_{c{c'}a} \f$ and \f$  \Psi_{ca{a'}} \f$ terms
 *  represent ternary interactions between two cations and
 *  an anion and two anions and a cation, respectively.
 *  In Pitzer's implementation these terms are usually small
 *  in absolute size. Currently these parameters do not have
 *  any dependence on temperature, pressure, or ionic strength.
 *
 *  Their values are input using the XML element
 *  <TT> psiCommonCation </TT> and  <TT> psiCommonAnion </TT>.
 *  The species id's are specified in attribute fields in
 *  the XML element. The fields  <TT>cation</TT>,
 *  <TT> anion1</TT>, and <TT> anion2</TT>
 *  are used for  <TT>psiCommonCation</TT>. The fields <TT> anion</TT>,
 *  <TT>cation1</TT> and  <TT>cation2</TT> are used for
 *  <TT> psiCommonAnion</TT>. An example block is given below.
 *  The <TT> Theta </TT> field below is a duplicate of the
 *  <TT> thetaAnion </TT> field mentioned above. The two fields
 *  are input into the same block for convenience, and because
 *  their data are highly correlated, in practice.
 *  It is an error for the
 *  two blocks to specify different information about
 *  thetaAnion (or thetaCation) in different blocks. It's
 *  ok to specify duplicate but consistent information
 *  in multiple blocks.
 *
 * @code
     <psiCommonCation cation="Na+" anion1="Cl-" anion2="OH-">
           <Theta> -0.05 </Theta>
           <Psi> -0.006 </Psi>
     </psiCommonCation>
    @endcode
 *
 *  <H3> Treatment of Neutral Species </H3>
 *
 *  Binary virial-coefficient-like interactions between two neutral
 *  species may be specified in the \f$ \lambda_{mn} \f$ terms
 *  that appear in the formulas above.
 *  Currently these interactions are independent of temperature,
 *  pressure, and ionic strength. Also, currently, the neutrality
 *  of the species are not checked. Therefore, this interaction
 *  may involve charged species in the solution as well.
 *  The identity of the species is specified by the
 *  <TT>species1</TT> and <TT>species2</TT> attributes to the XML
 *  <TT>lambdaNeutral</TT> node. These terms are symmetrical;
 *  <TT>species1</TT> and <TT>species2</TT> may be reversed and
 *  the term will be the same. An example is given below.
 *
 * @code
     <lambdaNeutral species1="CO2" species2="CH4">
           <lambda> 0.05 </lambda>
     </lambdaNeutral>
   @endcode
 *
 *  <H3> Example of the Specification of Parameters for the Activity
 *       Coefficients </H3>
 *
 * An example is given below.
 *
 * An example <TT> activityCoefficients </TT> XML block for this
 * formulation is supplied below
 *
 * @verbatim
  <activityCoefficients model="Pitzer" TempModel="complex1">
              <!-- Pitzer Coefficients
                   These coefficients are from Pitzer's main
                   paper, in his book.
                -->
              <A_Debye model="water" />
              <ionicRadius default="3.042843"  units="Angstroms">
              </ionicRadius>
              <binarySaltParameters cation="Na+" anion="Cl-">
                <beta0> 0.0765, 0.008946, -3.3158E-6,
                        -777.03, -4.4706
                </beta0>
                <beta1> 0.2664, 6.1608E-5, 1.0715E-6, 0.0, 0.0 </beta1>
                <beta2> 0.0,   0.0, 0.0, 0.0, 0.0  </beta2>
                <Cphi> 0.00127, -4.655E-5, 0.0,
                       33.317, 0.09421
                </Cphi>
                <Alpha1> 2.0 </Alpha1>
              </binarySaltParameters>

              <binarySaltParameters cation="H+" anion="Cl-">
                <beta0> 0.1775, 0.0, 0.0, 0.0, 0.0 </beta0>
                <beta1> 0.2945, 0.0, 0.0, 0.0, 0.0 </beta1>
                <beta2> 0.0,    0.0, 0.0, 0.0, 0.0 </beta2>
                <Cphi> 0.0008, 0.0, 0.0, 0.0, 0.0 </Cphi>
                <Alpha1> 2.0 </Alpha1>
              </binarySaltParameters>

              <binarySaltParameters cation="Na+" anion="OH-">
                <beta0> 0.0864, 0.0, 0.0, 0.0, 0.0 </beta0>
                <beta1> 0.253,  0.0, 0.0  0.0, 0.0 </beta1>
                <beta2> 0.0     0.0, 0.0, 0.0, 0.0 </beta2>
                <Cphi> 0.0044,  0.0, 0.0, 0.0, 0.0 </Cphi>
                <Alpha1> 2.0 </Alpha1>
              </binarySaltParameters>

              <thetaAnion anion1="Cl-" anion2="OH-">
                <Theta> -0.05,  0.0, 0.0, 0.0, 0.0 </Theta>
              </thetaAnion>

              <psiCommonCation cation="Na+" anion1="Cl-" anion2="OH-">
                <Theta> -0.05,  0.0, 0.0, 0.0, 0.0 </Theta>
                <Psi> -0.006 </Psi>
              </psiCommonCation>

              <thetaCation cation1="Na+" cation2="H+">
                <Theta> 0.036,  0.0, 0.0, 0.0, 0.0 </Theta>
              </thetaCation>

              <psiCommonAnion anion="Cl-" cation1="Na+" cation2="H+">
                <Theta> 0.036,  0.0, 0.0, 0.0, 0.0 </Theta>
                <Psi> -0.004 </Psi>
              </psiCommonAnion>

     </activityCoefficients>
  @endverbatim
 *
 *
 * <H3> Specification of the Debye-Huckel Constant </H3>
 *
 *  In the equations above, the formula for  \f$  A_{Debye} \f$
 *  is needed. The %HMWSoln object uses two methods for specifying these quantities.
 *  The default method is to assume that \f$  A_{Debye} \f$  is a constant, given
 *  in the initialization process, and stored in the
 *  member double, m_A_Debye. Optionally, a full water treatment may be employed that makes
 *  \f$ A_{Debye} \f$ a full function of <I>T</I> and <I>P</I> and creates nontrivial entries for
 *  the excess heat capacity, enthalpy, and excess volumes of solution.
 *
 *  \f[
 *      A_{Debye} = \frac{F e B_{Debye}}{8 \pi \epsilon R T} {\left( C_o \tilde{M}_o \right)}^{1/2}
 *  \f]
 *  where
 *
 *  \f[
 *      B_{Debye} = \frac{F} {{(\frac{\epsilon R T}{2})}^{1/2}}
 *  \f]
 *  Therefore:
 *  \f[
 *      A_{Debye} = \frac{1}{8 \pi}
 *                  {\left(\frac{2 N_a \rho_o}{1000}\right)}^{1/2}
 *                  {\left(\frac{N_a e^2}{\epsilon R T }\right)}^{3/2}
 *  \f]
 *
 *  Units = sqrt(kg/gmol)
 *
 *  where
 *  - \f$ N_a \f$ is Avogadro's number
 *  - \f$ \rho_w \f$ is the density of water
 *  - \f$ e \f$ is the electronic charge
 *  - \f$ \epsilon = K \epsilon_o \f$ is the permittivity of water
 *  - \f$ K \f$ is the dielectric constant of water,
 *  - \f$ \epsilon_o \f$ is the permittivity of free space.
 *  - \f$ \rho_o \f$ is the density of the solvent in its standard state.
 *
 *  Nominal value at 298 K and 1 atm = 1.172576 (kg/gmol)<SUP>1/2</SUP>
 *  based on:
 *  - \f$ \epsilon / \epsilon_0 \f$ = 78.54 (water at 25C)
 *  - T = 298.15 K
 *  - B_Debye = 3.28640E9 (kg/gmol)<SUP>1/2</SUP> m<SUP>-1</SUP>
 *
 *  An example of a fixed value implementation is given below.
 * @code
 *   <activityCoefficients model="Pitzer">
 *         <!-- A_Debye units = sqrt(kg/gmol)  -->
 *         <A_Debye> 1.172576 </A_Debye>
 *         <!-- object description continues -->
 *   </activityCoefficients>
 * @endcode
 *
 * An example of a variable value implementation within the %HMWSoln object is given below.
 * The model attribute, "water", triggers the full implementation.
 *
 * @code
 *   <activityCoefficients model="Pitzer">
 *         <!-- A_Debye units = sqrt(kg/gmol)  -->
 *         <A_Debye model="water" />
 *         <!-- object description continues -->
 *   </activityCoefficients>
 * @endcode
 *
 *
 * <H3>  Temperature and Pressure Dependence of the Activity Coefficients </H3>
 *
 *  Temperature dependence of the activity coefficients leads to nonzero terms
 *  for the excess enthalpy and entropy of solution. This means that the
 *  partial molar enthalpies, entropies, and heat capacities are all
 *  non-trivial to compute. The following formulas are used.
 *
 *  The partial molar enthalpy, \f$ \bar s_k(T,P) \f$:
 *
 *  \f[
 * \bar h_k(T,P) = h^{\triangle}_k(T,P)
 *               - R T^2 \frac{d \ln(\gamma_k^\triangle)}{dT}
 * \f]
 * The solvent partial molar enthalpy is equal to
 *  \f[
 * \bar h_o(T,P) = h^{o}_o(T,P) - R T^2 \frac{d \ln(a_o)}{dT}
 *       = h^{o}_o(T,P)
 *       + R T^2 (\sum_{k \neq o} m_k)  \tilde{M_o} (\frac{d \phi}{dT})
 * \f]
 *
 *  The partial molar entropy, \f$ \bar s_k(T,P) \f$:
 *
 *  \f[
 *     \bar s_k(T,P) =  s^{\triangle}_k(T,P)
 *             - R \ln( \gamma^{\triangle}_k \frac{m_k}{m^{\triangle}}))
 *                    - R T \frac{d \ln(\gamma^{\triangle}_k) }{dT}
 * \f]
 * \f[
 *      \bar s_o(T,P) = s^o_o(T,P) - R \ln(a_o)
 *                    - R T \frac{d \ln(a_o)}{dT}
 * \f]
 *
 * The partial molar heat capacity, \f$ C_{p,k}(T,P)\f$:
 *
 *  \f[
 *     \bar C_{p,k}(T,P) =  C^{\triangle}_{p,k}(T,P)
 *             - 2 R T \frac{d \ln( \gamma^{\triangle}_k)}{dT}
 *                    - R T^2 \frac{d^2 \ln(\gamma^{\triangle}_k) }{{dT}^2}
 * \f]
 * \f[
 *      \bar C_{p,o}(T,P) = C^o_{p,o}(T,P)
 *                   - 2 R T \frac{d \ln(a_o)}{dT}
 *                    - R T^2 \frac{d^2 \ln(a_o)}{{dT}^2}
 * \f]
 *
 *  The pressure dependence of the activity coefficients leads to non-zero terms
 *  for the excess Volume of the solution.
 *  Therefore, the partial molar volumes are functions
 *  of the pressure derivatives of the activity coefficients.
 * \f[
 *     \bar V_k(T,P) =  V^{\triangle}_k(T,P)
 *                    + R T \frac{d \ln(\gamma^{\triangle}_k) }{dP}
 * \f]
 * \f[
 *      \bar V_o(T,P) = V^o_o(T,P)
 *                    + R T \frac{d \ln(a_o)}{dP}
 * \f]
 *
 *  The majority of work for these functions take place in the internal
 *  routines that calculate the first and second derivatives of the log
 *  of the activity coefficients wrt temperature,
 *  s_update_dlnMolalityActCoeff_dT(), s_update_d2lnMolalityActCoeff_dT2(),
 *  and the first
 *  derivative of the log activity coefficients wrt pressure,
 *  s_update_dlnMolalityActCoeff_dP().
 *
 * <HR>
 * <H2> %Application within %Kinetics Managers </H2>
 * <HR>
 *
 * For the time being, we have set the standard concentration for all solute
 * species in
 * this phase equal to the default concentration of the solvent at the system temperature
 * and pressure multiplied by Mnaught (kg solvent / gmol solvent). The solvent
 * standard concentration is just equal to its standard state concentration.
 *
 *
 * This means that the
 * kinetics operator essentially works on an generalized concentration basis (kmol / m3),
 * with units for the kinetic rate constant specified
 * as if all reactants (solvent or solute) are on a concentration basis (kmol /m3).
 * The concentration will be modified by the activity coefficients.
 *
 * For example, a bulk-phase binary reaction between liquid solute species
 * <I>j</I> and <I>k</I>, producing
 * a new liquid solute species <I>l</I> would have the
 * following equation for its rate of progress variable, \f$ R^1 \f$, which has
 * units of kmol m-3 s-1.
 *
 *   \f[
 *    R^1 = k^1 C_j^a C_k^a =  k^1 (C^o_o \tilde{M}_o a_j) (C^o_o \tilde{M}_o a_k)
 *   \f]
 *
 * where
 *
 *   \f[
 *      C_j^a = C^o_o \tilde{M}_o a_j \quad and \quad C_k^a = C^o_o \tilde{M}_o a_k
 *   \f]
 *
 *  \f$ C_j^a \f$ is the activity concentration of species <I>j</I>, and
 *  \f$ C_k^a \f$ is the activity concentration of species <I>k</I>. \f$ C^o_o \f$
 *  is the concentration of water at 298 K and 1 atm. \f$ \tilde{M}_o \f$
 *  has units of kg solvent per gmol solvent and is equal to
 *
 * \f[
 *     \tilde{M}_o = \frac{M_o}{1000}
 * \f]
 *
 * \f$ a_j \f$ is
 *  the activity of species <I>j</I> at the current temperature and pressure
 *  and concentration of the liquid phase is given by the molality based
 *  activity coefficient multiplied by the molality of the jth species.
 *
 * \f[
 *      a_j  =  \gamma_j^\triangle m_j = \gamma_j^\triangle \frac{n_j}{\tilde{M}_o n_o}
 * \f]
 *
 * \f$k^1 \f$ has units of m<SUP>3</SUP>  kmol<SUP>-1</SUP> s<SUP>-1</SUP>.
 *
 *  Therefore the generalized activity concentration of a solute species has the following form
 *
 *  \f[
 *      C_j^a = C^o_o \frac{\gamma_j^\triangle n_j}{n_o}
 *  \f]
 *
 *  The generalized activity concentration of the solvent has the same units, but it's a simpler form
 *
 *  \f[
 *      C_o^a = C^o_o a_o
 *  \f]
 *
 *  The reverse rate constant can then be obtained from the law of microscopic reversibility
 *  and the equilibrium expression for the system.
 *
 *   \f[
 *      \frac{a_j a_k}{ a_l} = K^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
 *   \f]
 *
 *  \f$ K^{o,1} \f$ is the dimensionless form of the equilibrium constant.
 *
 *   \f[
 *       R^{-1} = k^{-1} C_l^a =  k^{-1} (C_o  \tilde{M}_o a_l)
 *   \f]
 *
 *  where
 *
 *   \f[
 *       k^{-1} =  k^1 K^{o,1} C_o \tilde{M}_o
 *   \f]
 *
 *  \f$ k^{-1} \f$ has units of s<SUP>-1</SUP>.
 *
 *  Note, this treatment may be modified in the future, as events dictate.
 *
 * <HR>
 * <H2> Instantiation of the Class </H2>
 * <HR>
 *
 * The constructor for this phase is now located in the default ThermoFactory
 * for %Cantera. The following code snippet may be used to initialize the phase
 * using the default construction technique within %Cantera.
 *
 * @code
 *      ThermoPhase *HMW = newPhase("HMW_NaCl.xml", "NaCl_electrolyte");
 * @endcode
 *
 *
 * A new %HMWSoln object may be created by  the following code snippets:
 *
 * @code
 *      HMWSoln *HMW = new HMWSoln("HMW_NaCl.xml", "NaCl_electrolyte");
 * @endcode
 *
 * or
 *
 * @code
 *    XML_Node *xm = get_XML_NameID("phase", "HMW_NaCl.xml#NaCl_electrolyte", 0);
 *    HMWSoln *dh = new HMWSoln(*xm);
 * @endcode
 *
 * or by the following call to importPhase():
 *
 * @code
 *    XML_Node *xm = get_XML_NameID("phase", "HMW_NaCl.xml#NaCl_electrolyte", 0);
 *    HMWSoln dhphase;
 *    importPhase(*xm, &dhphase);
 * @endcode
 *
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
             H2O(L) Na+ Cl- H+ OH-
  </speciesArray>
  <state>
    <temperature units="K"> 300  </temperature>
    <pressure units="Pa">101325.0</pressure>
    <soluteMolalities>
           Na+:3.0
           Cl-:3.0
           H+:1.0499E-8
           OH-:1.3765E-6
    </soluteMolalities>
  </state>
  <!-- thermo model identifies the inherited class
       from ThermoPhase that will handle the thermodynamics.
    -->
  <thermo model="HMW">
     <standardConc model="solvent_volume" />
   <activityCoefficients model="Pitzer" TempModel="complex1">
              <!-- Pitzer Coefficients
                   These coefficients are from Pitzer's main
                   paper, in his book.
                -->
              <A_Debye model="water" />
              <ionicRadius default="3.042843"  units="Angstroms">
              </ionicRadius>
              <binarySaltParameters cation="Na+" anion="Cl-">
                <beta0> 0.0765, 0.008946, -3.3158E-6,
                        -777.03, -4.4706
                </beta0>
                <beta1> 0.2664, 6.1608E-5, 1.0715E-6 </beta1>
                <beta2> 0.0    </beta2>
                <Cphi> 0.00127, -4.655E-5, 0.0,
                       33.317, 0.09421
                </Cphi>
                <Alpha1> 2.0 </Alpha1>
              </binarySaltParameters>

              <binarySaltParameters cation="H+" anion="Cl-">
                <beta0> 0.1775, 0.0, 0.0, 0.0, 0.0</beta0>
                <beta1> 0.2945, 0.0, 0.0 </beta1>
                <beta2> 0.0    </beta2>
                <Cphi> 0.0008, 0.0, 0.0, 0.0, 0.0 </Cphi>
                <Alpha1> 2.0 </Alpha1>
              </binarySaltParameters>

              <binarySaltParameters cation="Na+" anion="OH-">
                <beta0> 0.0864, 0.0, 0.0, 0.0, 0.0 </beta0>
                <beta1> 0.253, 0.0, 0.0 </beta1>
                <beta2> 0.0    </beta2>
                <Cphi> 0.0044, 0.0, 0.0, 0.0, 0.0 </Cphi>
                <Alpha1> 2.0 </Alpha1>
              </binarySaltParameters>

              <thetaAnion anion1="Cl-" anion2="OH-">
                <Theta> -0.05 </Theta>
              </thetaAnion>

              <psiCommonCation cation="Na+" anion1="Cl-" anion2="OH-">
                <Theta> -0.05 </Theta>
                <Psi> -0.006 </Psi>
              </psiCommonCation>

              <thetaCation cation1="Na+" cation2="H+">
                <Theta> 0.036 </Theta>
              </thetaCation>

              <psiCommonAnion anion="Cl-" cation1="Na+" cation2="H+">
                <Theta> 0.036 </Theta>
                <Psi> -0.004 </Psi>
              </psiCommonAnion>

     </activityCoefficients>

     <solvent> H2O(L) </solvent>
  </thermo>
  <elementArray datasrc="elements.xml"> O H Na Cl </elementArray>
  <kinetics model="none" >
  </kinetics>
</phase>
@endverbatim
 *
 *
 *
 * @ingroup thermoprops
 *
 */
class HMWSoln : public MolalityVPSSTP
{

public:

    //! Default Constructor
    HMWSoln();

    //! Construct and initialize an HMWSoln ThermoPhase object
    //! directly from an ASCII input file
    /*!
     *  This constructor is a shell that calls the routine initThermo(), with
     *  a reference to the XML database to get the info for the phase.
     *
     * @param inputFile Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     */
    HMWSoln(const std::string& inputFile, const std::string& id = "");

    //! Construct and initialize an HMWSoln ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     */
    HMWSoln(XML_Node& phaseRef, const std::string& id = "");

    //! Copy Constructor
    /*!
     * Copy constructor for the object. Constructed
     * object will be a clone of this object, but will
     * also own all of its data.
     * This is a wrapper around the assignment operator
     *
     * @param right Object to be copied.
     */
    HMWSoln(const HMWSoln& right);

    //! Assignment operator
    /*!
     * Assignment operator for the object. Constructed
     * object will be a clone of this object, but will
     * also own all of its data.
     *
     * @param right Object to be copied.
     */
    HMWSoln& operator=(const HMWSoln& right);

    //! This is a special constructor, used to replicate test problems
    //! during the initial verification of the object
    /*!
     *  test problems:
     *  1 = NaCl problem - 5 species - the thermo is read in from an XML file
     *
     *      speci   molality                        charge
     *      Cl-     6.0954          6.0997E+00      -1
     *      H+      1.0000E-08      2.1628E-09      1
     *      Na+     6.0954E+00      6.0997E+00      1
     *      OH-     7.5982E-07      1.3977E-06     -1
     *
     *      HMW_params____beta0MX__beta1MX__beta2MX__CphiMX_____alphaMX__thetaij
     *      10
     *      1  2          0.1775  0.2945   0.0      0.00080    2.0      0.0
     *      1  3          0.0765  0.2664   0.0      0.00127    2.0      0.0
     *      1  4          0.0     0.0      0.0      0.0        0.0     -0.050
     *      2  3          0.0     0.0      0.0      0.0        0.0      0.036
     *      2  4          0.0     0.0      0.0      0.0        0.0      0.0
     *      3  4          0.0864  0.253    0.0      0.0044     2.0      0.0
     *
     *      Triplet_interaction_parameters_psiaa'_or_psicc'
     *      2
     *      1  2  3   -0.004
     *      1  3  4   -0.006
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
    ThermoPhase* duplMyselfAsThermoPhase() const;

    //! Import, construct, and initialize a HMWSoln phase
    //! specification from an XML tree into the current object.
    /*
     * This routine is a precursor to constructPhaseXML(XML_Node*)
     * routine, which does most of the work.
     *
     * @param inputfile XML file containing the description of the phase
     *
     * @param id  Optional parameter identifying the name of the
     *            phase. If none is given, the first XML
     *            phase element will be used.
     */
    void constructPhaseFile(std::string inputFile, std::string id);

    //! Import and initialize a HMWSoln phase specification in an XML tree
    //! into the current object.
    /*!
     * Here we read an XML description of the phase. We import descriptions of
     * the elements that make up the species in a phase. We import information
     * about the species, including their reference state thermodynamic
     * polynomials. We then freeze the state of the species.
     *
     * Then, we read the species molar volumes from the xml tree to finish the
     * initialization.
     *
     * @param phaseNode This object must be the phase node of a complete XML tree
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

    //! @name  Utilities
    //! @{

    virtual int eosType() const;

    //! @}
    //! @name  Molar Thermodynamic Properties of the Solution
    //! @{

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
    virtual doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K.
    /*
     *      (HKM -> Bump up to Parent object)
     */
    virtual doublereal cv_mole() const;

    //!@}
    //! @name Mechanical Equation of State Properties
    /*!
     * In this equation of state implementation, the density is a function
     * only of the mole fractions. Therefore, it can't be an independent
     * variable. Instead, the pressure is used as the independent variable.
     * Functions which try to set the thermodynamic state by calling
     * setDensity() may cause an exception to be thrown.
     */
    //!@{

    /**
     * Pressure. Units: Pa.
     * For this incompressible system, we return the internally stored
     * independent value of the pressure.
     */
    virtual doublereal pressure() const;

    //! Set the internally stored pressure (Pa) at constant
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

protected:
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

public:
    //! Returns the current value of the density
    /*!
     *  @return value of the density. Units: kg/m^3
     */
    virtual doublereal density() const;

    //! Set the internally stored density (kg/m^3) of the phase.
    /*!
     * Overwritten setDensity() function is necessary because the
     * density is not an independent variable.
     *
     * This function will now throw an error condition.
     *
     * Note, in general, setting the phase density is now a nonlinear
     * calculation. P and T are the fundamental variables. This
     * routine should be revamped to do the nonlinear problem.
     *
     * @todo May have to adjust the strategy here to make
     * the eos for these materials slightly compressible, in order
     * to create a condition where the density is a function of
     * the pressure.
     * @todo Now have a compressible ss equation for liquid water.
     *       Therefore, this phase is compressible. May still
     *       want to change the independent variable however.
     *
     * @param rho Input density (kg/m^3).
     */
    void setDensity(const doublereal rho);

    //! Set the internally stored molar density (kmol/m^3) for the phase.
    /**
     * Overwritten setMolarDensity() function is necessary because of the
     * underlying water model.
     *
     * This function will now throw an error condition if the input
     * isn't exactly equal to the current molar density.
     *
     * @param conc   Input molar density (kmol/m^3).
     */
    void setMolarDensity(const doublereal conc);

    //! Set the temperature (K)
    /*!
     * This function sets the temperature, and makes sure that
     * the value propagates to underlying objects, such as
     * the water standard state model.
     *
     * @todo Make Phase::setTemperature a virtual function
     *
     * @param temp Temperature in kelvin
     */
    virtual void setTemperature(const doublereal temp);

    //! Set the temperature (K) and pressure (Pa)
    /*!
     *  Set the temperature and pressure.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     */
    virtual void setState_TP(doublereal t, doublereal p);

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

    //! This method returns an array of generalized activity concentrations
    /*!
     * The generalized activity concentrations, \f$ C_k^a\f$,  are defined such that
     * \f$ a_k = C^a_k / C^0_k, \f$ where \f$ C^0_k \f$
     * is a standard concentration
     * defined below.  These generalized concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions.
     *
     *  The generalized activity concentration of a solute species has the following form
     *
     *  \f[
     *      C_j^a = C^o_o \frac{\gamma_j^\triangle n_j}{n_o}
     *  \f]
     *
     *  The generalized activity concentration of the solvent has the same units, but it's a simpler form
     *
     *  \f[
     *      C_o^a = C^o_o a_o
     *  \f]
     *
     * @param c Array of generalized concentrations. The
     *          units are kmol m-3 for both the solvent and the solute species
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the activity (i.e., generalized) concentration for use
     *
     * We have set the standard concentration for all solute  species in
     * this phase equal to the default concentration of the solvent at the system temperature
     * and pressure multiplied by Mnaught (kg solvent / gmol solvent). The solvent
     * standard concentration is just equal to its standard state concentration.
     *
     *   \f[
     *      C_j^0 = C^o_o \tilde{M}_o \quad and  C_o^0 = C^o_o
     *   \f]
     *
     * The consequence of this is that the standard concentrations have unequal units
     * between the solvent and the solute. However, both the solvent and the solute
     * activity concentrations will have the same units of kmol kg<SUP>-3</SUP>.
     *
     * This means that the
     * kinetics operator essentially works on an generalized concentration basis (kmol / m3),
     * with units for the kinetic rate constant specified
     * as if all reactants (solvent or solute) are on a concentration basis (kmol /m3).
     * The concentration will be modified by the activity coefficients.
     *
     * For example, a bulk-phase binary reaction between liquid solute species
     * <I>j</I> and <I>k</I>, producing
     * a new liquid solute species <I>l</I> would have the
     * following equation for its rate of progress variable, \f$ R^1 \f$, which has
     * units of kmol m-3 s-1.
     *
     *   \f[
     *    R^1 = k^1 C_j^a C_k^a =  k^1 (C^o_o \tilde{M}_o a_j) (C^o_o \tilde{M}_o a_k)
     *   \f]
     *
     * where
     *
     *   \f[
     *      C_j^a = C^o_o \tilde{M}_o a_j \quad and \quad C_k^a = C^o_o \tilde{M}_o a_k
     *   \f]
     *
     * \f$ C_j^a \f$ is the activity concentration of species <I>j</I>, and
     * \f$ C_k^a \f$ is the activity concentration of species <I>k</I>. \f$ C^o_o \f$
     * is the concentration of water at 298 K and 1 atm. \f$ \tilde{M}_o \f$
     * has units of kg solvent per gmol solvent and is equal to
     *
     * \f[
     *     \tilde{M}_o = \frac{M_o}{1000}
     * \f]
     *
     * \f$ a_j \f$ is
     *  the activity of species <I>j</I> at the current temperature and pressure
     *  and concentration of the liquid phase is given by the molality based
     *  activity coefficient multiplied by the molality of the jth species.
     *
     * \f[
     *      a_j  =  \gamma_j^\triangle m_j = \gamma_j^\triangle \frac{n_j}{\tilde{M}_o n_o}
     * \f]
     *
     * \f$k^1 \f$ has units of m<SUP>3</SUP>  kmol<SUP>-1</SUP> s<SUP>-1</SUP>.
     *
     *  Therefore the generalized activity concentration of a solute species has the following form
     *
     *  \f[
     *      C_j^a = C^o_o \frac{\gamma_j^\triangle n_j}{n_o}
     *  \f]
     *
     *  The generalized activity concentration of the solvent has the same units, but it's a simpler form
     *
     *  \f[
     *      C_o^a = C^o_o a_o
     *  \f]
     *
     * @param k Optional parameter indicating the species. The default
     *         is to assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of
     *   m<SUP>3</SUP> kmol<SUP>-1</SUP>.
     *
     * @param k Species index
     */
    virtual doublereal standardConcentration(size_t k=0) const;

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
    virtual void getUnitsStandardConc(double* uA, int k = 0,
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

    //! @}
    //! @name  Partial Molar Properties of the Solution
    //! @{

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
     * molality-based activity coefficient wrt temperature
     *
     *  \f[
     * \bar h_k(T,P) = h^{\triangle}_k(T,P)
     *               - R T^2 \frac{d \ln(\gamma_k^\triangle)}{dT}
     * \f]
     * The solvent partial molar enthalpy is equal to
     *  \f[
     * \bar h_o(T,P) = h^{o}_o(T,P) - R T^2 \frac{d \ln(a_o)}{dT}
     *       = h^{o}_o(T,P)
     *       + R T^2 (\sum_{k \neq o} m_k)  \tilde{M_o} (\frac{d \phi}{dT})
     * \f]
     *
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;


    //! Returns an array of partial molar entropies of the species in the
    //! solution. Units: J/kmol/K.
    /*!
     * Maxwell's equations provide an answer for how calculate this
     *   (p.215 Smith and Van Ness)
     *
     *      d(chemPot_i)/dT = -sbar_i
     *
     * For this phase, the partial molar entropies are equal to the
     * SS species entropies plus the ideal solution contribution
     * plus complicated functions of the
     * temperature derivative of the activity coefficients.
     *
     *  \f[
     *     \bar s_k(T,P) =  s^{\triangle}_k(T,P)
     *             - R \ln( \gamma^{\triangle}_k \frac{m_k}{m^{\triangle}}))
     *                    - R T \frac{d \ln(\gamma^{\triangle}_k) }{dT}
     * \f]
     * \f[
     *      \bar s_o(T,P) = s^o_o(T,P) - R \ln(a_o)
     *                    - R T \frac{d \ln(a_o)}{dT}
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
     *  \f[
     *     \bar V_k(T,P)  = V^{\triangle}_k(T,P)
     *                    + R T \frac{d \ln(\gamma^{\triangle}_k) }{dP}
     * \f]
     * \f[
     *      \bar V_o(T,P) = V^o_o(T,P)
     *                    + R T \frac{d \ln(a_o)}{dP}
     * \f]
     *
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //! Return an array of partial molar heat capacities for the
    //! species in the mixture.  Units: J/kmol/K
    /*!
     *  The following formulas are implemented within the code.
     *
     *  \f[
     *     \bar C_{p,k}(T,P) =  C^{\triangle}_{p,k}(T,P)
     *             - 2 R T \frac{d \ln( \gamma^{\triangle}_k)}{dT}
     *                    - R T^2 \frac{d^2 \ln(\gamma^{\triangle}_k) }{{dT}^2}
     * \f]
     * \f[
     *      \bar C_{p,o}(T,P) = C^o_{p,o}(T,P)
     *                   - 2 R T \frac{d \ln(a_o)}{dT}
     *                    - R T^2 \frac{d^2 \ln(a_o)}{{dT}^2}
     * \f]
     *
     * @param cpbar   Output vector of species partial molar heat
     *                capacities at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

public:
    //! @}
    //! @name Chemical Equilibrium
    //! @{

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
        updateStandardStateThermo();
        throw NotImplementedError("HMWSoln::setToEquilState");
    }

    //@}

    //! Get the saturation pressure for a given temperature.
    /*!
     * Note the limitations of this function. Stability considerations
     * concerning multiphase equilibrium are ignored in this
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
    virtual doublereal satPressure(doublereal T);

    /*
     *  -------------- Utilities -------------------------------
     */

    //! Internal initialization required after all species have
    //! been added
    /*!
     * @internal Initialize. This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase().
     */
    virtual void initThermo();

    //! Initialize the phase parameters from an XML file.
    /*!
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
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Value of the Debye Huckel constant as a function of temperature
    //! and pressure.
    /*!
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     *
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calculation
     *                    or -1 to indicate the current pressure
     */
    virtual double A_Debye_TP(double temperature = -1.0,
                              double pressure = -1.0) const;

    //! Value of the derivative of the Debye Huckel constant with
    //! respect to temperature as a function of temperature
    //! and pressure.
    /*!
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     *
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calculation
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
     * @param pressure    Pressure of the derivative calculation
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
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calculation
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
     * @param pressure    Pressure of the derivative calculation
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
     * @param pressure    Pressure of the derivative calculation
     *                    or -1 to indicate the current pressure
     */
    double ADebye_V(double temperature = -1.0,
                    double pressure = -1.0) const;

    //! Value of the 2nd derivative of the Debye Huckel constant with
    //! respect to temperature as a function of temperature
    //! and pressure.
    /*!
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     *
     * @param temperature  Temperature of the derivative calculation
     *                     or -1 to indicate the current temperature
     *
     * @param pressure    Pressure of the derivative calculation
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
     * formPitzer():
     *
     *  Returns the form of the Pitzer parameterization used
     */
    int formPitzer() const {
        return m_formPitzer;
    }

    //! Print out all of the input Pitzer coefficients.
    void printCoeffs() const;

    //!  Get the array of unscaled non-dimensional molality based
    //!  activity coefficients at the current solution temperature,
    //!  pressure, and solution concentration.
    /*!
     *  See Denbigh p. 278 for a thorough discussion. This class must be overwritten in
     *  classes which derive from %MolalityVPSSTP. This function takes over from the
     *  molar-based activity coefficient calculation, getActivityCoefficients(), in
     *  derived classes.
     *
     * @param acMolality Output vector containing the molality based activity coefficients.
     *                   length: m_kk.
     */
    void getUnscaledMolalityActivityCoefficients(doublereal* acMolality) const;

private:

    //! Apply the current phScale to a set of activity Coefficients
    /*!
     *  See the Eq3/6 Manual for a thorough discussion.
     */
    void s_updateScaling_pHScaling() const;

    //!  Apply the current phScale to a set of derivatives of the activity Coefficients
    //!  wrt temperature
    /*!
     *  See the Eq3/6 Manual for a thorough discussion of the need
     */
    void s_updateScaling_pHScaling_dT() const;

    //!  Apply the current phScale to a set of 2nd derivatives of the activity Coefficients
    //!  wrt temperature
    /*!
     *  See the Eq3/6 Manual for a thorough discussion of the need
     */
    void s_updateScaling_pHScaling_dT2() const;

    //!  Apply the current phScale to a set of derivatives of the activity Coefficients
    //!  wrt pressure
    /*!
     *  See the Eq3/6 Manual for a thorough discussion of the need
     */
    void s_updateScaling_pHScaling_dP() const;

    //!  Calculate the Chlorine activity coefficient on the NBS scale
    /*!
     *  We assume here that the m_IionicMolality variable is up to date.
     */
    doublereal s_NBS_CLM_lnMolalityActCoeff() const;

    //!  Calculate the temperature derivative of the Chlorine activity coefficient
    //!  on the NBS scale
    /*!
     *  We assume here that the m_IionicMolality variable is up to date.
     */
    doublereal s_NBS_CLM_dlnMolalityActCoeff_dT() const;

    //!  Calculate the second temperature derivative of the Chlorine activity coefficient
    //!  on the NBS scale
    /*!
     *  We assume here that the m_IionicMolality variable is up to date.
     */
    doublereal s_NBS_CLM_d2lnMolalityActCoeff_dT2() const;

    //!  Calculate the pressure derivative of the Chlorine activity coefficient
    /*!
     *  We assume here that the m_IionicMolality variable is up to date.
     */
    doublereal s_NBS_CLM_dlnMolalityActCoeff_dP() const;

    //@}

private:

    /**
     * This is the form of the Pitzer parameterization
     * used in this model.
     * The options are described at the top of this document,
     * and in the general documentation.
     * The list is repeated here:
     *
     * PITZERFORM_BASE  = 0    (only one supported atm)
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
     *  a_k = Size of the ionic species in the DH formulation
     *        units = meters
     */
    vector_fp m_Aionic;

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

    //! Reference Temperature for the Pitzer formulations.
    double m_TempPitzerRef;

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

private:
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
     *            A_Debye = (1/ (8 Pi)) (2 Na * dw/1000)^(1/2)
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
    PDSS* m_waterSS;

    //! density of standard-state water
    /*!
     * internal temporary variable
     */
    double m_densWaterSS;

    //! Pointer to the water property calculator
    WaterProps* m_waterProps;

    //! Temporary array used in equilibrium calculations
    mutable vector_fp      m_pp;

    //! vector of size m_kk, used as a temporary holding area.
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
    mutable vector_fp m_Beta2MX_ij;

    //! Derivative of Beta2_ij[i][j] wrt T
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Beta2MX_ij_L;

    //! Derivative of Beta2_ij[i][j] wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Beta2MX_ij_LL;

    //! Derivative of Beta2_ij[i][j] wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Beta2MX_ij_P;

    //! Array of coefficients for Beta2, a variable in Pitzer's papers
    /*!
     * column index is counterIJ
     *  m_Beta2MX_ij_coeff.ptrColumn(counterIJ) is a double* containing
     *  the vector of coefficients for the counterIJ interaction.
     *  This was added for the YMP database version of the code since it
     *  contains temperature-dependent parameters for some 2-2 electrolytes.
     */
    mutable Array2D   m_Beta2MX_ij_coeff;

    /**
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Alpha1MX_ij[i][j] is the value of the alpha1 coefficient
     *  for the ij interaction. It will be nonzero iff i and j are
     *  both charged and have opposite sign.
     *  It is symmetric wrt i, j.
     *  counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     */
    vector_fp m_Alpha1MX_ij;

    /**
      *  Array of 2D data used in the Pitzer/HMW formulation.
      *  Alpha2MX_ij[i][j] is the value of the alpha2 coefficient
      *  for the ij interaction. It will be nonzero iff i and j are
      *  both charged and have opposite sign, and i and j
      *  both have charges of 2 or more, usually.
      *  It is symmetric wrt i, j.
      *  counterIJ, where counterIJ = m_counterIJ[i][j],
      *  is used to access this array.
      */
    vector_fp m_Alpha2MX_ij;

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

    //! Array of coefficients for CphiMX, a parameter in the activity
    //! coefficient formulation
    /*!
     *  Column index is counterIJ
     *  m_CphiMX_ij_coeff.ptrColumn(counterIJ) is a double* containing
     *  the vector of coefficients for the counterIJ interaction.
     */
    mutable Array2D m_CphiMX_ij_coeff;

    //! Array of 2D data for Theta_ij[i][j] in the Pitzer/HMW formulation.
    /*!
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
    mutable vector_fp m_Theta_ij;

    //! Derivative of Theta_ij[i][j] wrt T
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Theta_ij_L;

    //! Derivative of Theta_ij[i][j] wrt TT
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Theta_ij_LL;

    //! Derivative of Theta_ij[i][j] wrt P
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_Theta_ij_P;

    //! Array of coefficients for Theta_ij[i][j] in the Pitzer/HMW formulation.
    /*!
     *  Theta_ij[i][j] is the value of the theta coefficient
     *  for the ij interaction. It will be nonzero for charged
     *  ions with the same sign. It is symmetric.
     *  Column index is counterIJ.
     *  counterIJ where counterIJ = m_counterIJ[i][j]
     *  is used to access this array.
     *
     *  m_Theta_ij_coeff.ptrColumn(counterIJ) is a double* containing
     *  the vector of coefficients for the counterIJ interaction.
     */
    Array2D m_Theta_ij_coeff;

    //! Array of 3D data used in the Pitzer/HMW formulation.
    /*!
     * Psi_ijk[n] is the value of the psi coefficient for the
     * ijk interaction where
     *
     *   n = k + j * m_kk + i * m_kk * m_kk;
     *
     * It is potentially nonzero everywhere.
     * The first two coordinates are symmetric wrt cations,
     * and the last two coordinates are symmetric wrt anions.
     */
    mutable vector_fp m_Psi_ijk;

    //! Derivative of Psi_ijk[n] wrt T
    /*!
     *  see m_Psi_ijk for reference on the indexing into this variable.
     */
    mutable vector_fp m_Psi_ijk_L;

    //! Derivative of Psi_ijk[n] wrt TT
    /*!
     *  see m_Psi_ijk for reference on the indexing into this variable.
     */
    mutable vector_fp m_Psi_ijk_LL;

    //! Derivative of Psi_ijk[n] wrt P
    /*!
     *  see m_Psi_ijk for reference on the indexing into this variable.
     */
    mutable vector_fp m_Psi_ijk_P;

    //! Array of coefficients for Psi_ijk[n] in the Pitzer/HMW formulation.
    /*!
     * Psi_ijk[n] is the value of the psi coefficient for the
     * ijk interaction where
     *
     *   n = k + j * m_kk + i * m_kk * m_kk;
     *
     * It is potentially nonzero everywhere.
     * The first two coordinates are symmetric wrt cations,
     * and the last two coordinates are symmetric wrt anions.
     *
     *
     *  m_Psi_ijk_coeff.ptrColumn(n) is a double* containing
     *  the vector of coefficients for the n interaction.
     */
    Array2D m_Psi_ijk_coeff;

    //! Lambda coefficient for the ij interaction
    /*!
     * Array of 2D data used in the Pitzer/HMW formulation.
     * Lambda_nj[n][j] represents the lambda coefficient for the
     * ij interaction. This is a general interaction representing
     * neutral species. The neutral species occupy the first
     * index, i.e., n. The charged species occupy the j coordinate.
     * neutral, neutral interactions are also included here.
     */
    mutable Array2D   m_Lambda_nj;

    //! Derivative of  Lambda_nj[i][j] wrt T. see m_Lambda_ij
    mutable Array2D   m_Lambda_nj_L;

    //! Derivative of  Lambda_nj[i][j] wrt TT
    mutable Array2D   m_Lambda_nj_LL;

    //! Derivative of  Lambda_nj[i][j] wrt P
    mutable Array2D   m_Lambda_nj_P;

    //! Array of coefficients for Lambda_nj[i][j] in the Pitzer/HMW formulation.
    /*!
     *  Lambda_ij[i][j] is the value of the theta coefficient
     *  for the ij interaction.
     *  Array of 2D data used in the Pitzer/HMW formulation.
     *  Lambda_ij[i][j] represents the lambda coefficient for the
     *  ij interaction. This is a general interaction representing
     *  neutral species. The neutral species occupy the first
     *  index, i.e., i. The charged species occupy the j coordinate.
     *  Neutral, neutral interactions are also included here.
     *
     *      n = j + m_kk * i
     *
     *  m_Lambda_ij_coeff.ptrColumn(n) is a double* containing
     *  the vector of coefficients for the (i,j) interaction.
     */
    Array2D  m_Lambda_nj_coeff;

    //! Mu coefficient for the self-ternary neutral coefficient
    /*!
     * Array of 2D data used in the Pitzer/HMW formulation.
     * Mu_nnn[i] represents the Mu coefficient for the
     * nnn interaction. This is a general interaction representing
     * neutral species interacting with itself.
     */
    mutable vector_fp  m_Mu_nnn;

    //! Mu coefficient temperature derivative for the self-ternary neutral coefficient
    /*!
     * Array of 2D data used in the Pitzer/HMW formulation.
     * Mu_nnn_L[i] represents the Mu coefficient temperature derivative for the
     * nnn interaction. This is a general interaction representing
     * neutral species interacting with itself.
     */
    mutable vector_fp  m_Mu_nnn_L;

    //! Mu coefficient 2nd temperature derivative for the self-ternary neutral coefficient
    /*!
     * Array of 2D data used in the Pitzer/HMW formulation.
     * Mu_nnn_L[i] represents the Mu coefficient 2nd temperature derivative for the
     * nnn interaction. This is a general interaction representing
     * neutral species interacting with itself.
     */
    mutable vector_fp  m_Mu_nnn_LL;

    //! Mu coefficient pressure derivative for the self-ternary neutral coefficient
    /*!
     * Array of 2D data used in the Pitzer/HMW formulation.
     * Mu_nnn_L[i] represents the Mu coefficient pressure derivative for the
     * nnn interaction. This is a general interaction representing
     * neutral species interacting with itself.
     */
    mutable vector_fp  m_Mu_nnn_P;

    //! Array of coefficients form_Mu_nnn term
    Array2D   m_Mu_nnn_coeff;

    //!  Logarithm of the activity coefficients on the molality
    //!  scale.
    /*!
     *   mutable because we change this if the composition
     *   or temperature or pressure changes.
     *
     *  index is the species index
     */
    mutable vector_fp m_lnActCoeffMolal_Scaled;

    //!  Logarithm of the activity coefficients on the molality
    //!  scale.
    /*!
     *   mutable because we change this if the composition
     *   or temperature or pressure changes.
     *
     *  index is the species index
     */
    mutable vector_fp m_lnActCoeffMolal_Unscaled;

    //!  Derivative of the Logarithm of the activity coefficients on the molality
    //!  scale wrt T
    /*!
     *  index is the species index
     */
    mutable vector_fp m_dlnActCoeffMolaldT_Scaled;

    //!  Derivative of the Logarithm of the activity coefficients on the molality
    //!  scale wrt T
    /*!
     *  index is the species index
     */
    mutable vector_fp m_dlnActCoeffMolaldT_Unscaled;

    //!  Derivative of the Logarithm of the activity coefficients on the molality
    //!  scale wrt TT
    /*!
     *  index is the species index
     */
    mutable vector_fp m_d2lnActCoeffMolaldT2_Scaled;

    //!  Derivative of the Logarithm of the activity coefficients on the molality
    //!  scale wrt TT
    /*!
     *  index is the species index
     */
    mutable vector_fp m_d2lnActCoeffMolaldT2_Unscaled;

    //!  Derivative of the Logarithm of the activity coefficients on the
    //!  molality scale wrt P
    /*!
     *  index is the species index
     */
    mutable vector_fp m_dlnActCoeffMolaldP_Scaled;

    //!  Derivative of the Logarithm of the activity coefficients on the
    //!  molality scale wrt P
    /*!
     *  index is the species index
     */
    mutable vector_fp m_dlnActCoeffMolaldP_Unscaled;

    /*
     * -------- Temporary Variables Used in the Activity Coeff Calc
     */

    //! Cropped and modified values of the molalities used in activity
    //! coefficient calculations
    mutable vector_fp m_molalitiesCropped;

    //! Boolean indicating whether the molalities are cropped or are modified
    mutable bool m_molalitiesAreCropped;

    //! a counter variable for keeping track of symmetric binary
    //! interactions amongst the solute species.
    /*!
     * n = m_kk*i + j
     * m_CounterIJ[n] = counterIJ
     */
    mutable vector_int m_CounterIJ;

    //! This is elambda, MEC
    mutable double elambda[17];

    //! This is elambda1, MEC
    mutable double elambda1[17];

    /**
     *  Various temporary arrays used in the calculation of
     *  the Pitzer activity coefficients.
     *  The subscript, L, denotes the same quantity's derivative
     *  wrt temperature
     */

    //! This is the value of g(x) in Pitzer's papers
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_gfunc_IJ;

    //! This is the value of g2(x2) in Pitzer's papers
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_g2func_IJ;

    //! hfunc, was called gprime in Pitzer's paper. However,
    //! it's not the derivative of gfunc(x), so I renamed it.
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_hfunc_IJ;

    //! hfunc2, was called gprime in Pitzer's paper. However,
    //! it's not the derivative of gfunc(x), so I renamed it.
    /*!
     *  vector index is counterIJ
     */
    mutable vector_fp m_h2func_IJ;

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
    mutable vector_fp m_gamma_tmp;

    //! Logarithm of the molal activity coefficients
    /*!
     *   Normally these are all one. However, stability schemes will change that
     */
    mutable vector_fp      IMS_lnActCoeffMolal_;

    //! IMS Cutoff type
    int IMS_typeCutoff_;

    //! value of the solute mole fraction that centers the cutoff polynomials
    //! for the cutoff =1 process;
    doublereal IMS_X_o_cutoff_;

    //! gamma_o value for the cutoff process at the zero solvent point
    doublereal IMS_gamma_o_min_;

    //! gamma_k minimum for the cutoff process at the zero solvent point
    doublereal IMS_gamma_k_min_;

    //! Parameter in the polyExp cutoff treatment having to do with rate of exp decay
    doublereal IMS_cCut_;

    //! Parameter in the polyExp cutoff treatment
    /*!
     *  This is the slope of the f function at the zero solvent point
     *  Default value is 0.6
     */
    doublereal IMS_slopefCut_;

    //! Parameter in the polyExp cutoff treatment
    /*!
     *  This is the slope of the g function at the zero solvent point
     *  Default value is 0.0
     */
    doublereal IMS_slopegCut_;

    //! @name Parameters in the polyExp cutoff treatment having to do with rate of exp decay
    //! @{
    doublereal IMS_dfCut_;
    doublereal IMS_efCut_;
    doublereal IMS_afCut_;
    doublereal IMS_bfCut_;
    doublereal IMS_dgCut_;
    doublereal IMS_egCut_;
    doublereal IMS_agCut_;
    doublereal IMS_bgCut_;
    //! @}

    //! value of the solvent mole fraction that centers the cutoff polynomials
    //! for the cutoff =1 process;
    doublereal MC_X_o_cutoff_;

    //! gamma_o value for the cutoff process at the zero solvent point
    doublereal MC_X_o_min_;

    //! Parameter in the Molality Exp cutoff treatment
    /*!
     *  This is the slope of the p function at the zero solvent point
     *  Default value is 0.0
     */
    doublereal MC_slopepCut_;

    //! @name Parameters in the Molality Exp cutoff treatment
    //! @{
    doublereal MC_dpCut_;
    doublereal MC_epCut_;
    doublereal MC_apCut_;
    doublereal MC_bpCut_;
    doublereal MC_cpCut_;
    doublereal CROP_ln_gamma_o_min;
    doublereal CROP_ln_gamma_o_max;
    doublereal CROP_ln_gamma_k_min;
    doublereal CROP_ln_gamma_k_max;

    //! This is a boolean-type vector indicating whether
    //! a species's activity coefficient is in the cropped regime
    /*!
     *  * 0 = Not in cropped regime
     *  * 1 = In a transition regime where it is altered but there
     *        still may be a temperature or pressure dependence
     *  * 2 = In a cropped regime where there is no temperature
     *        or pressure dependence
     */
    mutable std::vector<int> CROP_speciesCropped_;
    //! @}

    //!  Initialize all of the species-dependent lengths in the object
    void initLengths();

    //! Apply the current phScale to a set of activity Coefficients or
    //! activities
    /*!
     *  See the Eq3/6 Manual for a thorough discussion.
     *
     * @param acMolality input/Output vector containing the molality based
     *                   activity coefficients. length: m_kk.
     */
    virtual void applyphScale(doublereal* acMolality) const;

private:
    /*
     * This function will be called to update the internally stored
     * natural logarithm of the molality activity coefficients
     */
    void s_update_lnMolalityActCoeff() const;

    //! This function calculates the temperature derivative of the
    //! natural logarithm of the molality activity coefficients.
    /*!
     * This function does all of the direct work. The solvent activity
     * coefficient is on the molality scale. It's derivative is too.
     */
    void s_update_dlnMolalityActCoeff_dT() const;

    /**
     * This function calculates the temperature second derivative
     * of the natural logarithm of the molality activity
     * coefficients.
     */
    void s_update_d2lnMolalityActCoeff_dT2() const;

    /**
     * This function calculates the pressure derivative of the
     * natural logarithm of the molality activity coefficients.
     *
     * Assumes that the activity coefficients are current.
     */
    void s_update_dlnMolalityActCoeff_dP() const;

    //! This function will be called to update the internally stored
    //! natural logarithm of the molality activity coefficients
    /*
     * Normally they are all one. However, sometimes they are not,
     * due to stability schemes
     *
     *    gamma_k_molar =  gamma_k_molal / Xmol_solvent
     *
     *    gamma_o_molar = gamma_o_molal
     */
    void s_updateIMS_lnMolalityActCoeff() const;

private:
    //! Calculate the Pitzer portion of the activity coefficients.
    /**
     *  This is the main routine in the whole module. It calculates the
     *  molality based activity coefficients for the solutes, and
     *  the activity of water.
     */
    void s_updatePitzer_lnMolalityActCoeff() const;

    //!  Calculates the temperature derivative of the
    //!  natural logarithm of the molality activity coefficients.
    /*!
     * Public function makes sure that all dependent data is
     * up to date, before calling a private function
     */
    void s_updatePitzer_dlnMolalityActCoeff_dT() const;

    /**
     * This function calculates the temperature second derivative of the
     * natural logarithm of the molality activity coefficients.
     *
     * It is assumed that the Pitzer activity coefficient and first derivative
     * routine are called immediately preceding the call to this routine.
     */
    void s_updatePitzer_d2lnMolalityActCoeff_dT2() const;

    //!  Calculates the Pressure derivative of the
    //!  natural logarithm of the molality activity coefficients.
    /*!
     * It is assumed that the Pitzer activity coefficient and first derivative
     * routine are called immediately preceding the calling of this routine.
     */
    void s_updatePitzer_dlnMolalityActCoeff_dP() const;

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
    void s_updatePitzer_CoeffWRTemp(int doDerivs = 2) const;

    //! Calculate the lambda interactions.
    /*!
     *
     * Calculate E-lambda terms for charge combinations of like sign, using
     * method of Pitzer (1975). This implementation is based on Bethke,
     * Appendix 2.
     *
     * @param is Ionic strength
     */
    void calc_lambdas(double is) const;
    mutable doublereal m_last_is;

    /**
     *  Calculate etheta and etheta_prime
     *
     * This interaction accounts for the mixing effects of like-signed ions
     * with different charges. This interaction will be nonzero for species
     * with the same charge. this routine is not to be called for neutral
     * species; it core dumps or error exits.
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
                     double* etheta, double* etheta_prime) const;

    //! Set up a counter variable for keeping track of symmetric binary
    //! interactions amongst the solute species.
    /*!
     * The purpose of this is to squeeze the ij parameters into a
     * compressed single counter.
     *
     * n = m_kk*i + j
     * m_Counter[n] = counter
     */
    void counterIJ_setup() const;

    //! Calculate the cropped molalities
    /*!
     * This is an internal routine that calculates values
     * of m_molalitiesCropped from m_molalities
     */
    void calcMolalitiesCropped() const;

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
    void readXMLBinarySalt(XML_Node& BinSalt);

    //! Process an XML node called "thetaAnion"
    /*!
     * This node contains all of the parameters necessary to describe
     * the binary interactions between two anions.
     *
     * @param BinSalt  reference to the XML_Node named thetaAnion
     *                 containing the
     *                 anion - anion interaction
     */
    void readXMLThetaAnion(XML_Node& BinSalt);

    //! Process an XML node called "thetaCation"
    /*!
     * This node contains all of the parameters necessary to describe
     * the binary interactions between two cations.
     *
     * @param BinSalt  reference to the XML_Node named thetaCation
     *                 containing the
     *                 cation - cation interaction
     */
    void readXMLThetaCation(XML_Node& BinSalt);

    //! Process an XML node called "psiCommonAnion"
    /*!
     * This node contains all of the parameters necessary to describe
     * the ternary interactions between one anion and two cations.
     *
     * @param BinSalt  reference to the XML_Node named psiCommonAnion
     *                 containing the
     *                 anion - cation1 - cation2 interaction
     */
    void readXMLPsiCommonAnion(XML_Node& BinSalt);

    //! Process an XML node called "psiCommonCation"
    /*!
     * This node contains all of the parameters necessary to describe
     * the ternary interactions between one cation and two anions.
     *
     * @param BinSalt  reference to the XML_Node named psiCommonCation
     *                 containing the
     *                 cation - anion1 - anion2 interaction
     */
    void readXMLPsiCommonCation(XML_Node& BinSalt);

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
    void readXMLLambdaNeutral(XML_Node& BinSalt);

    //! Process an XML node called "MunnnNeutral"
    /*!
     * This node contains all of the parameters necessary to describe
     * the self-ternary interactions for one neutral species.
     *
     * @param BinSalt  reference to the XML_Node named Munnn
     *                 containing the self-ternary interaction
     */
    void readXMLMunnnNeutral(XML_Node& BinSalt);

    //! Process an XML node called "zetaCation"
    /*!
     * This node contains all of the parameters necessary to describe
     * the ternary interactions between one neutral, one cation, and one anion.
     *
     * @param BinSalt  reference to the XML_Node named psiCommonCation
     *                 containing the
     *                 neutral - cation - anion interaction
     */
    void readXMLZetaCation(const XML_Node& BinSalt);

    //! Process an XML node called "croppingCoefficients"
    //! for the cropping coefficients values
    /*!
     * @param acNode Activity Coefficient XML Node
     */
    void readXMLCroppingCoefficients(const XML_Node& acNode);

    //! Precalculate the IMS Cutoff parameters for typeCutoff = 2
    void  calcIMSCutoffParams_();

    //! Calculate molality cut-off parameters
    void  calcMCCutoffParams_();

    //! Utility function to assign an integer value from a string
    //! for the ElectrolyteSpeciesType field.
    /*!
     *  @param estString string name of the electrolyte species type
     */
    static int interp_est(const std::string& estString);

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
