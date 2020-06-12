/**
 *  @file DebyeHuckel.h
 *    Headers for the DebyeHuckel ThermoPhase object, which models dilute
 *    electrolyte solutions
 *    (see \ref thermoprops and \link Cantera::DebyeHuckel DebyeHuckel \endlink) .
 *
 * Class DebyeHuckel represents a dilute liquid electrolyte phase which
 * obeys the Debye Huckel formulation for nonideality.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DEBYEHUCKEL_H
#define CT_DEBYEHUCKEL_H

#include "MolalityVPSSTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{

/*!
 * @name Formats for the Activity Coefficients
 *
 *   These are possible formats for the molality-based activity coefficients.
 */
//@{
#define DHFORM_DILUTE_LIMIT 0
#define DHFORM_BDOT_AK 1
#define DHFORM_BDOT_ACOMMON 2
#define DHFORM_BETAIJ 3
#define DHFORM_PITZER_BETAIJ 4
//@}
/*!
 *  @name  Acceptable ways to calculate the value of A_Debye
 */
//@{
#define A_DEBYE_CONST 0
#define A_DEBYE_WATER 1
//@}

class WaterProps;
class PDSS_Water;

/**
 * @ingroup thermoprops
 *
 * Class DebyeHuckel represents a dilute liquid electrolyte phase which obeys
 * the Debye Huckel formulation for nonideality.
 *
 * The concentrations of the ionic species are assumed to obey the
 * electroneutrality condition.
 *
 * ## Specification of Species Standard State Properties
 *
 * The standard states are on the unit molality basis. Therefore, in the
 * documentation below, the normal \f$ o \f$ superscript is replaced with the
 * \f$ \triangle \f$ symbol. The reference state symbol is now
 *  \f$ \triangle, ref \f$.
 *
 * It is assumed that the reference state thermodynamics may be obtained by a
 * pointer to a populated species thermodynamic property manager class (see
 * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
 * state thermodynamics is resolved at this level.
 *
 * For an incompressible, stoichiometric substance, the molar internal energy is
 * independent of pressure. Since the thermodynamic properties are specified by
 * giving the standard-state enthalpy, the term \f$ P_0 \hat v\f$ is subtracted
 * from the specified molar enthalpy to compute the molar internal energy. The
 * entropy is assumed to be independent of the pressure.
 *
 * The enthalpy function is given by the following relation.
 *
 *       \f[
 *         h^\triangle_k(T,P) = h^{\triangle,ref}_k(T)
 *         + \tilde v \left( P - P_{ref} \right)
 *       \f]
 *
 * For an incompressible, stoichiometric substance, the molar internal energy is
 * independent of pressure. Since the thermodynamic properties are specified by
 * giving the standard-state enthalpy, the term \f$ P_{ref} \tilde v\f$ is
 * subtracted from the specified reference molar enthalpy to compute the molar
 * internal energy.
 *
 *       \f[
 *            u^\triangle_k(T,P) = h^{\triangle,ref}_k(T) - P_{ref} \tilde v
 *       \f]
 *
 * The standard state heat capacity and entropy are independent of pressure. The
 * standard state Gibbs free energy is obtained from the enthalpy and entropy
 * functions.
 *
 * The current model assumes that an incompressible molar volume for all
 * solutes. The molar volume for the water solvent, however, is obtained from a
 * pure water equation of state, waterSS. Therefore, the water standard state
 * varies with both T and P. It is an error to request standard state water
 * properties at a T and P where the water phase is not a stable phase, i.e.,
 * beyond its spinodal curve.
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * Chemical potentials of the solutes, \f$ \mu_k \f$, and the solvent, \f$ \mu_o
 * \f$, which are based on the molality form, have the following general format:
 *
 * \f[
 *    \mu_k = \mu^{\triangle}_k(T,P) + R T ln(\gamma_k^{\triangle} \frac{m_k}{m^\triangle})
 * \f]
 * \f[
 *    \mu_o = \mu^o_o(T,P) + RT ln(a_o)
 * \f]
 *
 * where \f$ \gamma_k^{\triangle} \f$ is the molality based activity coefficient
 * for species \f$k\f$.
 *
 * Individual activity coefficients of ions can not be independently measured.
 * Instead, only binary pairs forming electroneutral solutions can be measured.
 *
 * ### Ionic Strength
 *
 * Most of the parameterizations within the model use the ionic strength as a
 * key variable. The ionic strength, \f$ I\f$ is defined as follows
 *
 *  \f[
 *    I = \frac{1}{2} \sum_k{m_k  z_k^2}
 *  \f]
 *
 * \f$ m_k \f$ is the molality of the kth species. \f$ z_k \f$ is the charge of
 * the kth species. Note, the ionic strength is a defined units quantity. The
 * molality has defined units of gmol kg-1, and therefore the ionic strength has
 * units of sqrt( gmol kg-1).
 *
 * In some instances, from some authors, a different formulation is used for the
 * ionic strength in the equations below. The different formulation is due to
 * the possibility of the existence of weak acids and how association wrt to the
 * weak acid equilibrium relation affects the calculation of the activity
 * coefficients via the assumed value of the ionic strength.
 *
 * If we are to assume that the association reaction doesn't have an effect on
 * the ionic strength, then we will want to consider the associated weak acid as
 * in effect being fully dissociated, when we calculate an effective value for
 * the ionic strength. We will call this calculated value, the stoichiometric
 * ionic strength, \f$ I_s \f$, putting a subscript s to denote it from the more
 * straightforward calculation of \f$ I \f$.
 *
 *  \f[
 *    I_s = \frac{1}{2} \sum_k{m_k^s  z_k^2}
 *  \f]
 *
 * Here, \f$ m_k^s \f$ is the value of the molalities calculated assuming that
 * all weak acid-base pairs are in their fully dissociated states. This
 * calculation may be simplified by considering that the weakly associated acid
 * may be made up of two charged species, k1 and k2, each with their own
 * charges, obeying the following relationship:
 *
 *   \f[
 *      z_k = z_{k1} +  z_{k2}
 *   \f]
 * Then, we may only need to specify one charge value, say, \f$  z_{k1}\f$, the
 * cation charge number, in order to get both numbers, since we have already
 * specified \f$ z_k \f$ in the definition of original species. Then, the
 * stoichiometric ionic strength may be calculated via the following formula.
 *
 *  \f[
 *    I_s = \frac{1}{2} \left(\sum_{k,ions}{m_k  z_k^2}+
 *               \sum_{k,weak_assoc}(m_k  z_{k1}^2 + m_k  z_{k2}^2) \right)
 *  \f]
 *
 * The specification of which species are weakly associated acids is made in the
 * input file via the `stoichIsMods` XML block, where the charge for k1 is also
 * specified. An example is given below:
 *
 * @code
 *          <stoichIsMods>
 *                NaCl(aq):-1.0
 *          </stoichIsMods>
 * @endcode
 *
 * Because we need the concept of a weakly associated acid in order to calculate
 * \f$ I_s \f$ we need to catalog all species in the phase. This is done using
 * the following categories:
 *
 *  -  `cEST_solvent`                Solvent species (neutral)
 *  -  `cEST_chargedSpecies`         Charged species (charged)
 *  -  `cEST_weakAcidAssociated`     Species which can break apart into charged species.
 *                                   It may or may not be charged.  These may or
 *                                   may not be be included in the
 *                                   species solution vector.
 *  -  `cEST_strongAcidAssociated`   Species which always breaks apart into charged species.
 *                                   It may or may not be charged. Normally, these aren't included
 *                                   in the speciation vector.
 *  -  `cEST_polarNeutral`           Polar neutral species
 *  -  `cEST_nonpolarNeutral`        Non polar neutral species
 *
 * Polar and non-polar neutral species are differentiated, because some
 * additions to the activity coefficient expressions distinguish between these
 * two types of solutes. This is the so-called salt-out effect.
 *
 * The type of species is specified in the `electrolyteSpeciesType` XML block.
 * Note, this is not considered a part of the specification of the standard
 * state for the species, at this time. Therefore, this information is put under
 * the `activityCoefficient` XML block. An example is given below
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
 * Much of the species electrolyte type information is inferred from other
 * information in the input file. For example, as species which is charged is
 * given the "chargedSpecies" default category. A neutral solute species is put
 * into the "nonpolarNeutral" category by default.
 *
 * The specification of solute activity coefficients depends on the model
 * assumed for the Debye-Huckel term. The model is set by the internal parameter
 * #m_formDH. We will now describe each category in its own section.
 *
 * ### Debye-Huckel Dilute Limit
 *
 * DHFORM_DILUTE_LIMIT = 0
 *
 * This form assumes a dilute limit to DH, and is mainly for informational purposes:
 * \f[
 *     \ln(\gamma_k^\triangle) = - z_k^2 A_{Debye} \sqrt{I}
 * \f]
 * where \f$ I\f$ is the ionic strength
 * \f[
 *   I = \frac{1}{2} \sum_k{m_k  z_k^2}
 * \f]
 *
 * The activity for the solvent water,\f$ a_o \f$, is not independent and must
 * be determined from the Gibbs-Duhem relation.
 *
 * \f[
 *      \ln(a_o) = \frac{X_o - 1.0}{X_o} + \frac{ 2 A_{Debye} \tilde{M}_o}{3} (I)^{3/2}
 * \f]
 *
 * ### Bdot Formulation
 *
 *    DHFORM_BDOT_AK       = 1
 *
 * This form assumes Bethke's format for the Debye Huckel activity coefficient:
 *
 * \f[
 *    \ln(\gamma_k^\triangle) = -z_k^2 \frac{A_{Debye} \sqrt{I}}{ 1 + B_{Debye}  a_k \sqrt{I}}
 *                        + \log(10) B^{dot}_k  I
 * \f]
 *
 * Note, this particular form where \f$ a_k \f$ can differ in multielectrolyte
 * solutions has problems with respect to a Gibbs-Duhem analysis. However, we
 * include it here because there is a lot of data fit to it.
 *
 * The activity for the solvent water,\f$ a_o \f$, is not independent and must
 * be determined from the Gibbs-Duhem relation. Here, we use:
 *
 *  \f[
 *       \ln(a_o) = \frac{X_o - 1.0}{X_o}
 *        + \frac{ 2 A_{Debye} \tilde{M}_o}{3} (I)^{1/2}
 *                        \left[ \sum_k{\frac{1}{2} m_k z_k^2 \sigma( B_{Debye} a_k \sqrt{I} ) } \right]
 *                        - \frac{\log(10)}{2} \tilde{M}_o I \sum_k{ B^{dot}_k m_k}
 *  \f]
 *    where
 *  \f[
 *     \sigma (y) = \frac{3}{y^3} \left[ (1+y) - 2 \ln(1 + y) - \frac{1}{1+y} \right]
 *  \f]
 *
 * Additionally, Helgeson's formulation for the water activity is offered as an
 * alternative.
 *
 * ### Bdot Formulation with Uniform Size Parameter in the Denominator
 *
 * DHFORM_BDOT_AUNIFORM = 2
 *
 * This form assumes Bethke's format for the Debye-Huckel activity coefficient
 *
 * \f[
 *    \ln(\gamma_k^\triangle) = -z_k^2 \frac{A_{Debye} \sqrt{I}}{ 1 + B_{Debye}  a \sqrt{I}}
 *                        + \log(10) B^{dot}_k  I
 * \f]
 *
 * The value of a is determined at the beginning of the calculation, and not changed.
 *
 *  \f[
 *       \ln(a_o) = \frac{X_o - 1.0}{X_o}
 *        + \frac{ 2 A_{Debye} \tilde{M}_o}{3} (I)^{3/2} \sigma( B_{Debye} a \sqrt{I} )
 *                        - \frac{\log(10)}{2} \tilde{M}_o I \sum_k{ B^{dot}_k m_k}
 *  \f]
 *
 * ### Beta_IJ formulation
 *
 *     DHFORM_BETAIJ        = 3
 *
 * This form assumes a linear expansion in a virial coefficient form. It is used
 * extensively in the book by Newmann, "Electrochemistry Systems", and is the
 * beginning of more complex treatments for stronger electrolytes, fom Pitzer
 * and from Harvey, Moller, and Weire.
 *
 * \f[
 *    \ln(\gamma_k^\triangle) = -z_k^2 \frac{A_{Debye} \sqrt{I}}{ 1 + B_{Debye}  a \sqrt{I}}
 *                         + 2 \sum_j \beta_{j,k} m_j
 * \f]
 *
 * In the current treatment the binary interaction coefficients, \f$
 * \beta_{j,k}\f$, are independent of temperature and pressure.
 *
 * \f[
 *    \ln(a_o) = \frac{X_o - 1.0}{X_o}
 *     + \frac{ 2 A_{Debye} \tilde{M}_o}{3} (I)^{3/2} \sigma( B_{Debye} a \sqrt{I} )
 *     -  \tilde{M}_o  \sum_j \sum_k \beta_{j,k} m_j m_k
 * \f]
 *
 * In this formulation the ionic radius, \f$ a \f$, is a constant. This must be
 * supplied to the model, in an <DFN> ionicRadius </DFN> XML block.
 *
 * The \f$ \beta_{j,k} \f$ parameters are binary interaction parameters. They
 * are supplied to the object in an `DHBetaMatrix` XML block. There are in
 * principle \f$ N (N-1) /2 \f$ different, symmetric interaction parameters,
 * where \f$ N \f$ are the number of solute species in the mechanism. An example
 * is given below.
 *
 * An example `activityCoefficients` XML block for this formulation is supplied
 * below
 *
 * @code
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
 * ### Pitzer Beta_IJ formulation
 *
 *     DHFORM_PITZER_BETAIJ  = 4
 *
 * This form assumes an activity coefficient formulation consistent with a
 * truncated form of Pitzer's formulation. Pitzer's formulation is equivalent to
 * the formulations above in the dilute limit, where rigorous theory may be
 * applied.
 *
 * \f[
 *    \ln(\gamma_k^\triangle) = -z_k^2 \frac{A_{Debye}}{3} \frac{\sqrt{I}}{ 1 + B_{Debye}  a \sqrt{I}}
 *       -2 z_k^2 \frac{A_{Debye}}{3}  \frac{\ln(1 + B_{Debye}  a  \sqrt{I})}{ B_{Debye}  a}
 *                         + 2 \sum_j \beta_{j,k} m_j
 * \f]
 * \f[
 *       \ln(a_o) = \frac{X_o - 1.0}{X_o}
 *        + \frac{ 2 A_{Debye} \tilde{M}_o}{3} \frac{(I)^{3/2} }{1 +  B_{Debye}  a \sqrt{I} }
 *        -  \tilde{M}_o  \sum_j \sum_k \beta_{j,k} m_j m_k
 * \f]
 *
 * ### Specification of the Debye Huckel Constants
 *
 * In the equations above, the formulas for  \f$  A_{Debye} \f$ and \f$
 * B_{Debye} \f$ are needed. The DebyeHuckel object uses two methods for
 * specifying these quantities. The default method is to assume that \f$
 * A_{Debye} \f$  is a constant, given in the initialization process, and stored
 * in the member double, m_A_Debye. Optionally, a full water treatment may be
 * employed that makes \f$ A_{Debye} \f$ a full function of *T* and *P*.
 *
 * \f[
 *      A_{Debye} = \frac{F e B_{Debye}}{8 \pi \epsilon R T} {\left( C_o \tilde{M}_o \right)}^{1/2}
 * \f]
 * where
 *  \f[
 *         B_{Debye} = \frac{F} {{(\frac{\epsilon R T}{2})}^{1/2}}
 *  \f]
 * Therefore:
 * \f[
 *   A_{Debye} = \frac{1}{8 \pi}
 *                 {\left(\frac{2 N_a \rho_o}{1000}\right)}^{1/2}
 *                 {\left(\frac{N_a e^2}{\epsilon R T }\right)}^{3/2}
 * \f]
 * where
 *   - \f$ N_a \f$ is Avogadro's number
 *   - \f$ \rho_w \f$ is the density of water
 *   - \f$ e \f$ is the electronic charge
 *   - \f$ \epsilon = K \epsilon_o \f$ is the permittivity of water
 *   - \f$ K \f$ is the dielectric constant of water
 *   - \f$ \epsilon_o \f$ is the permittivity of free space
 *   - \f$ \rho_o \f$ is the density of the solvent in its standard state.
 *
 * Nominal value at 298 K and 1 atm = 1.172576 (kg/gmol)^(1/2) based on:
 *   - \f$ \epsilon / \epsilon_0 \f$ = 78.54 (water at 25C)
 *   - T = 298.15 K
 *   - B_Debye = 3.28640E9 (kg/gmol)^(1/2) / m
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
 * @code
 *   <activityCoefficients model="Beta_ij">
 *         <A_Debye model="water" />
 *         <!-- B_Debye units = sqrt(kg/gmol)/m  -->
 *         <B_Debye> 3.28640E9 </B_Debye>
 *   </activityCoefficients>
 * @endcode
 *
 * Currently, \f$  B_{Debye} \f$ is a constant in the model, specified either by
 * a default water value, or through the input file. This may have to be looked
 * at, in the future.
 *
 * ## %Application within Kinetics Managers
 *
 * For the time being, we have set the standard concentration for all species in
 * this phase equal to the default concentration of the solvent at 298 K and 1
 * atm. This means that the kinetics operator essentially works on an activities
 * basis, with units specified as if it were on a concentration basis.
 *
 * For example, a bulk-phase binary reaction between liquid species j and k,
 * producing a new liquid species l would have the following equation for its
 * rate of progress variable, \f$ R^1 \f$, which has units of kmol m-3 s-1.
 *
 * \f[
 *    R^1 = k^1 C_j^a C_k^a =  k^1 (C_o a_j) (C_o a_k)
 * \f]
 * where
 * \f[
 *      C_j^a = C_o a_j \quad and \quad C_k^a = C_o a_k
 * \f]
 *
 * \f$ C_j^a \f$ is the activity concentration of species j, and
 * \f$ C_k^a \f$ is the activity concentration of species k. \f$ C_o \f$
 * is the concentration of water at 298 K and 1 atm. \f$ a_j \f$ is the activity
 * of species j at the current temperature and pressure and concentration of the
 * liquid phase. \f$k^1 \f$ has units of m3 kmol-1 s-1.
 *
 * The reverse rate constant can then be obtained from the law of microscopic
 * reversibility and the equilibrium expression for the system.
 *
 * \f[
 *     \frac{a_j a_k}{ a_l} = K^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
 * \f]
 *
 * \f$  K^{o,1} \f$ is the dimensionless form of the equilibrium constant.
 *
 * \f[
 *    R^{-1} = k^{-1} C_l^a =  k^{-1} (C_o a_l)
 * \f]
 * where
 * \f[
 *     k^{-1} =  k^1 K^{o,1} C_o
 * \f]
 *
 * \f$k^{-1} \f$ has units of s-1.
 *
 * Note, this treatment may be modified in the future, as events dictate.
 *
 * ## Instantiation of the Class
 *
 * The constructor for this phase is NOT located in the default ThermoFactory
 * for %Cantera. However, a new DebyeHuckel object may be created by
 * the following code snippets:
 *
 * @code
 *      DebyeHuckel *DH = new DebyeHuckel("DH_NaCl.xml", "NaCl_electrolyte");
 * @endcode
 *
 * or
 *
 * @code
 *    XML_Node *xm = get_XML_NameID("phase", "DH_NaCl.xml#NaCl_electrolyte", 0);
 *    DebyeHuckel *dh = new DebyeHuckel(*xm);
 * @endcode
 *
 * or by the following call to importPhase():
 *
 * @code
 *    XML_Node *xm = get_XML_NameID("phase", "DH_NaCl.xml#NaCl_electrolyte", 0);
 *    DebyeHuckel dhphase;
 *    importPhase(*xm, &dhphase);
 * @endcode
 *
 * ## XML Example
 *
 * The phase model name for this is called StoichSubstance. It must be supplied
 * as the model attribute of the thermo XML element entry. Within the phase XML
 * block, the density of the phase must be specified. An example of an XML file
 * this phase is given below.
 *
 * @code
 * <phase id="NaCl_electrolyte" dim="3">
 *  <speciesArray datasrc="#species_waterSolution">
 *             H2O(L) Na+ Cl- H+ OH- NaCl(aq) NaOH(aq)
 *  </speciesArray>
 *  <state>
 *    <temperature units="K"> 300  </temperature>
 *    <pressure units="Pa">101325.0</pressure>
 *    <soluteMolalities>
 *           Na+:3.0
 *           Cl-:3.0
 *           H+:1.0499E-8
 *           OH-:1.3765E-6
 *           NaCl(aq):0.98492
 *           NaOH(aq):3.8836E-6
 *    </soluteMolalities>
 *  </state>
 *  <!-- thermo model identifies the inherited class
 *       from ThermoPhase that will handle the thermodynamics.
 *    -->
 *  <thermo model="DebyeHuckel">
 *     <standardConc model="solvent_volume" />
 *     <activityCoefficients model="Beta_ij">
 *              <!-- A_Debye units = sqrt(kg/gmol)  -->
 *              <A_Debye> 1.172576 </A_Debye>
 *              <!-- B_Debye units = sqrt(kg/gmol)/m   -->
 *              <B_Debye> 3.28640E9 </B_Debye>
 *              <ionicRadius default="3.042843"  units="Angstroms">
 *              </ionicRadius>
 *              <DHBetaMatrix>
 *                H+:Cl-:0.27
 *                Na+:Cl-:0.15
 *                Na+:OH-:0.06
 *              </DHBetaMatrix>
 *              <stoichIsMods>
 *                 NaCl(aq):-1.0
 *              </stoichIsMods>
 *              <electrolyteSpeciesType>
 *                 H+:chargedSpecies
 *                 NaCl(aq):weakAcidAssociated
 *              </electrolyteSpeciesType>
 *     </activityCoefficients>
 *     <solvent> H2O(L) </solvent>
 *  </thermo>
 *  <elementArray datasrc="elements.xml"> O H Na Cl </elementArray>
 * </phase>
 * @endcode
 */
class DebyeHuckel : public MolalityVPSSTP
{
public:
    //! Default Constructor
    DebyeHuckel();

    virtual ~DebyeHuckel();

    //! Full constructor for creating the phase.
    /*!
     *  @param inputFile  File name containing the definition of the phase
     *  @param id         id attribute containing the name of the phase.
     */
    DebyeHuckel(const std::string& inputFile, const std::string& id = "");

    //! Full constructor for creating the phase.
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id       id attribute containing the name of the phase.
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    DebyeHuckel(XML_Node& phaseRef, const std::string& id = "");

    //! @name  Utilities
    //! @{

    virtual std::string type() const {
        return "DebyeHuckel";
    }

    //! @}
    //! @name  Molar Thermodynamic Properties of the Solution
    //! @{

    virtual doublereal enthalpy_mole() const;

    /// Molar entropy. Units: J/kmol/K.
    /**
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
     * @see MultiSpeciesThermo
     */
    virtual doublereal entropy_mole() const;

    virtual doublereal gibbs_mole() const;
    virtual doublereal cp_mole() const;

    //@}
    /** @name Mechanical Equation of State Properties
     //@{
     * In this equation of state implementation, the density is a function only
     * of the mole fractions. Therefore, it can't be an independent variable.
     * Instead, the pressure is used as the independent variable. Functions
     * which try to set the thermodynamic state by calling setDensity() will
     * cause an exception to be thrown.
     */

protected:
    virtual void calcDensity();

public:
    //! Set the internally stored density (gm/m^3) of the phase.
    /*!
     * Overridden setDensity() function is necessary because the density is not
     * an independent variable.
     *
     * This function will now throw an error condition
     *
     * @internal May have to adjust the strategy here to make the eos for these
     *     materials slightly compressible, in order to create a condition where
     *     the density is a function of the pressure.
     *
     * This function will now throw an error condition if the input isn't
     * exactly equal to the current density.
     *
     * @todo Now have a compressible ss equation for liquid water. Therefore,
     *       this phase is compressible. May still want to change the
     *       independent variable however.
     *
     * @param rho Input density (kg/m^3).
     * @deprecated Functionality merged with base function after Cantera 2.5.
     *             (superseded by isCompressible check in Phase::setDensity)
     */
    virtual void setDensity(const doublereal rho);

    //! Set the internally stored molar density (kmol/m^3) of the phase.
    /**
     * Overridden setMolarDensity() function is necessary because the density
     * is not an independent variable.
     *
     * This function will now throw an error condition if the input isn't
     * exactly equal to the current molar density.
     *
     * @param conc   Input molar density (kmol/m^3).
     * @deprecated Functionality merged with base function after Cantera 2.5.
     *             (superseded by isCompressible check in Phase::setDensity)
     */
    virtual void setMolarDensity(const doublereal conc);

    /**
     * @}
     * @name Activities, Standard States, and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is related to the
     * chemical potential by \f[ \mu_k = \mu_k^0(T) + \hat R T \log a_k. \f] The
     * quantity \f$\mu_k^0(T,P)\f$ is the chemical potential at unit activity,
     * which depends only on temperature and the pressure. Activity is assumed
     * to be molality-based here.
     * @{
     */

    virtual void getActivityConcentrations(doublereal* c) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize the activity
     * (i.e., generalized) concentration in kinetics calculations.
     *
     * For the time being, we will use the concentration of pure solvent for the
     * the standard concentration of all species. This has the effect of making
     * reaction rates based on the molality of species proportional to the
     * molality of the species.
     *
     * @param k Optional parameter indicating the species. The default is to
     *         assume this refers to species 0.
     * @return the standard Concentration in units of m^3/kmol
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    //! Get the array of non-dimensional activities at the current solution
    //! temperature, pressure, and solution concentration.
    /*!
     * (note solvent activity coefficient is on molar scale).
     *
     * @param ac  Output vector of activities. Length: m_kk.
     */
    virtual void getActivities(doublereal* ac) const;

    //! Get the array of non-dimensional molality-based activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     *  note solvent is on molar scale. The solvent molar based activity
     *  coefficient is returned.
     *
     *  Note, most of the work is done in an internal private routine
     *
     * @param acMolality Vector of Molality-based activity coefficients
     *                   Length: m_kk
     */
    virtual void getMolalityActivityCoefficients(doublereal* acMolality) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     *
     * This function returns a vector of chemical potentials of the species in
     * solution.
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
     * \bar h_k(T,P) = h^{\triangle}_k(T,P) - R T^2 \frac{d \ln(\gamma_k^\triangle)}{dT}
     * \f]
     * The solvent partial molar enthalpy is equal to
     *  \f[
     * \bar h_o(T,P) = h^{o}_o(T,P) - R T^2 \frac{d \ln(a_o}{dT}
     * \f]
     *
     * The temperature dependence of the activity coefficients currently
     * only occurs through the temperature dependence of the Debye constant.
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
     * \f[
     *    \frac{d\mu_i}{dT} = -\bar{s}_i
     * \f]
     *
     * For this phase, the partial molar entropies are equal to the SS species
     * entropies plus the ideal solution contribution:
     * \f[
     *     \bar s_k(T,P) =  \hat s^0_k(T) - R log(M0 * molality[k])
     * \f]
     * \f[
     *     \bar s_{solvent}(T,P) =  \hat s^0_{solvent}(T)
     *                 - R ((xmolSolvent - 1.0) / xmolSolvent)
     * \f]
     *
     * The reference-state pure-species entropies,\f$ \hat s^0_k(T) \f$, at the
     * reference pressure, \f$ P_{ref} \f$, are computed by the species
     * thermodynamic property manager. They are polynomial functions of
     * temperature.
     * @see MultiSpeciesThermo
     *
     *  @param sbar    Output vector of species partial molar entropies.
     *                 Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    virtual void getPartialMolarCp(doublereal* cpbar) const;

    //! Return an array of partial molar volumes for the species in the mixture.
    //! Units: m^3/kmol.
    /*!
     * For this solution, the partial molar volumes are normally equal to the
     * constant species molar volumes, except when the activity coefficients
     * depend on pressure.
     *
     * The general relation is
     *
     *       vbar_i = d(chemPot_i)/dP at const T, n
     *              = V0_i + d(Gex)/dP)_T,M
     *              = V0_i + RT d(lnActCoeffi)dP _T,M
     *
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //@}

    /*
     *  -------------- Utilities -------------------------------
     */

    virtual bool addSpecies(shared_ptr<Species> spec);
    virtual void initThermo();
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Return the Debye Huckel constant as a function of temperature
    //! and pressure (Units = sqrt(kg/gmol))
    /*!
     *  The default is to assume that it is constant, given in the
     *  initialization process, and stored in the member double, m_A_Debye.
     *  Optionally, a full water treatment may be employed that makes
     *  \f$ A_{Debye} \f$ a full function of T and P.
     *
     *   \f[
     *      A_{Debye} = \frac{F e B_{Debye}}{8 \pi \epsilon R T} {\left( C_o \tilde{M}_o \right)}^{1/2}
     *   \f]
     *  where
     *  \f[
     *         B_{Debye} = \frac{F} {{(\frac{\epsilon R T}{2})}^{1/2}}
     *  \f]
     *  Therefore:
     *  \f[
     *   A_{Debye} = \frac{1}{8 \pi}
     *                 {\left(\frac{2 N_a \rho_o}{1000}\right)}^{1/2}
     *                 {\left(\frac{N_a e^2}{\epsilon R T }\right)}^{3/2}
     *  \f]
     *
     *  where
     *  - Units = sqrt(kg/gmol)
     *  - \f$ N_a \f$ is Avogadro's number
     *  - \f$ \rho_w \f$ is the density of water
     *  - \f$ e \f$ is the electronic charge
     *  - \f$ \epsilon = K \epsilon_o \f$ is the permittivity of water
     *  - \f$ K \f$ is the dielectric constant of water,
     *  - \f$ \epsilon_o \f$ is the permittivity of free space.
     *  - \f$ \rho_o \f$ is the density of the solvent in its standard state.
     *
     *  Nominal value at 298 K and 1 atm = 1.172576 (kg/gmol)^(1/2)
     *  based on:
     *    - \f$ \epsilon / \epsilon_0 \f$ = 78.54 (water at 25C)
     *    - T = 298.15 K
     *    - B_Debye = 3.28640E9 (kg/gmol)^(1/2)/m
     *
     * @param temperature Temperature in kelvin. Defaults to -1, in which
     *                    case the   temperature of the phase is assumed.
     * @param pressure Pressure (Pa). Defaults to -1, in which
     *                    case the pressure of the phase is assumed.
     */
    virtual double A_Debye_TP(double temperature = -1.0,
                              double pressure = -1.0) const;

    //! Value of the derivative of the Debye Huckel constant with
    //! respect to temperature.
    /*!
     * This is a function of temperature and pressure. See A_Debye_TP() for
     * a definition of \f$ A_{Debye} \f$.
     *
     * Units = sqrt(kg/gmol) K-1
     *
     * @param temperature Temperature in kelvin. Defaults to -1, in which
     *                    case the   temperature of the phase is assumed.
     * @param pressure Pressure (Pa). Defaults to -1, in which
     *                    case the pressure of the phase is assumed.
     */
    virtual double dA_DebyedT_TP(double temperature = -1.0,
                                 double pressure = -1.0) const;

    //! Value of the 2nd derivative of the Debye Huckel constant with
    //! respect to temperature as a function of temperature and pressure.
    /*!
     * This is a function of temperature and pressure. See A_Debye_TP() for
     * a definition of \f$ A_{Debye} \f$.
     *
     * Units = sqrt(kg/gmol) K-2
     *
     * @param temperature Temperature in kelvin. Defaults to -1, in which
     *                    case the   temperature of the phase is assumed.
     * @param pressure Pressure (Pa). Defaults to -1, in which
     *                    case the pressure of the phase is assumed.
     */
    virtual double d2A_DebyedT2_TP(double temperature = -1.0,
                                   double pressure = -1.0) const;

    //! Value of the derivative of the Debye Huckel constant with
    //! respect to pressure, as a function of temperature and pressure.
    /*!
     * This is a function of temperature and pressure. See A_Debye_TP() for
     * a definition of \f$ A_{Debye} \f$.
     *
     * Units = sqrt(kg/gmol) Pa-1
     *
     * @param temperature Temperature in kelvin. Defaults to -1, in which
     *                    case the   temperature of the phase is assumed.
     * @param pressure Pressure (Pa). Defaults to -1, in which
     *                    case the pressure of the phase is assumed.
     */
    virtual double dA_DebyedP_TP(double temperature = -1.0,
                                 double pressure = -1.0) const;

    //! Reports the ionic radius of the kth species
    /*!
     * @param k  species index.
     */
    double AionicRadius(int k = 0) const;

    //! Set the DebyeHuckel parameterization form. Must be one of
    //! 'dilute-limit', 'B-dot-with-variable-a', 'B-dot-with-common-a',
    //! 'beta_ij', or 'Pitzer-with-beta_ij'.
    void setDebyeHuckelModel(const std::string& form);

    //! Returns the form of the Debye-Huckel parameterization used
    int formDH() const {
        return m_formDH;
    }

    //! Set the A_Debye parameter. If a negative value is provided, enables
    //! calculation of A_Debye using the detailed water equation of state.
    void setA_Debye(double A);

    void setB_Debye(double B) { m_B_Debye = B; }
    void setB_dot(double bdot);
    void setMaxIonicStrength(double Imax) { m_maxIionicStrength = Imax; }
    void useHelgesonFixedForm(bool mode=true) { m_useHelgesonFixedForm = mode; }

    //! Set the default ionic radius [m] for each species
    void setDefaultIonicRadius(double value);

    //! Set the value for the beta interaction between species sp1 and sp2.
    void setBeta(const std::string& sp1, const std::string& sp2, double value);

    //! Returns a reference to M_Beta_ij
    Array2D& get_Beta_ij() {
        return m_Beta_ij;
    }

private:
    //! Static function that implements the non-polar species salt-out
    //! modifications.
    /*!
     *  Returns the calculated activity coefficients.
     *
     * @param IionicMolality Value of the ionic molality (sqrt(gmol/kg))
     */
    static double _nonpolarActCoeff(double IionicMolality);

    //! Formula for the osmotic coefficient that occurs in the GWB.
    /*!
     *  It is originally from Helgeson for a variable NaCl brine. It's to be
     *  used with extreme caution.
     */
    double _osmoticCoeffHelgesonFixedForm() const;

    //! Formula for the log of the water activity that occurs in the GWB.
    /*!
     *  It is originally from Helgeson for a variable NaCl brine. It's to be
     *  used with extreme caution.
     */
    double _lnactivityWaterHelgesonFixedForm() const;
    //@}

protected:
    //! form of the Debye-Huckel parameterization used in the model.
    /*!
     * The options are described at the top of this document,
     * and in the general documentation.
     * The list is repeated here:
     *
     * DHFORM_DILUTE_LIMIT  = 0       (default)
     * DHFORM_BDOT_AK       = 1
     * DHFORM_BDOT_AUNIFORM = 2
     * DHFORM_BETAIJ        = 3
     * DHFORM_PITZER_BETAIJ = 4
     */
    int m_formDH;

    //! Vector containing the electrolyte species type
    /*!
     * The possible types are:
     *  - solvent
     *  - Charged Species
     *  - weakAcidAssociated
     *  - strongAcidAssociated
     *  - polarNeutral
     *  - nonpolarNeutral
     *  .
     */
    vector_int m_electrolyteSpeciesType;

    //! a_k = Size of the ionic species in the DH formulation. units = meters
    vector_fp m_Aionic;

    //! Current value of the ionic strength on the molality scale
    mutable double m_IionicMolality;

    //! Maximum value of the ionic strength allowed in the calculation of the
    //! activity coefficients.
    double m_maxIionicStrength;

public:
    //! If true, then the fixed for of Helgeson's activity for water is used
    //! instead of the rigorous form obtained from Gibbs-Duhem relation. This
    //! should be used with caution, and is really only included as a validation
    //! exercise.
    bool m_useHelgesonFixedForm;
protected:
    //! Stoichiometric ionic strength on the molality scale
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
    //! Current value of the Debye Constant, A_Debye
    /**
     * A_Debye -> this expression appears on the top of the ln actCoeff term in
     *            the general Debye-Huckel expression It depends on temperature
     *            and pressure.
     *
     *            A_Debye = (F e B_Debye) / (8 Pi epsilon R T)
     *
     *            Units = sqrt(kg/gmol)
     *
     *            Nominal value(298K, atm) = 1.172576 sqrt(kg/gmol)
     *                  based on:
     *                    epsilon/epsilon_0 = 78.54
     *                           (water at 25C)
     *                    T = 298.15 K
     *                    B_Debye = 3.28640E9 sqrt(kg/gmol)/m
     *
     *            note in Pitzer's nomenclature, A_phi = A_Debye/3.0
     */
    mutable double m_A_Debye;

    //! Current value of the constant that appears in the denominator
    /**
     * B_Debye -> this expression appears on the bottom of the ln actCoeff term
     *            in the general Debye-Huckel expression It depends on
     *            temperature
     *
     *            B_Bebye = F / sqrt( epsilon R T / 2 )
     *
     *            Units = sqrt(kg/gmol) / m
     *
     *            Nominal value = 3.28640E9 sqrt(kg/gmol) / m
     *                  based on:
     *                    epsilon/epsilon_0 = 78.54
     *                           (water at 25C)
     *                    T = 298.15 K
     */
    double m_B_Debye;

    //! Array of B_Dot values
    /**
     *  This expression is an extension of the Debye-Huckel expression used
     *  in some formulations to extend DH to higher molalities. B_dot is
     *  specific to the major ionic pair.
     */
    vector_fp m_B_Dot;

    //! Pointer to the Water standard state object
    /*!
     *  derived from the equation of state for water.
     */
    PDSS_Water* m_waterSS;

    //! Storage for the density of water's standard state
    /*!
     * Density depends on temperature and pressure.
     */
    double m_densWaterSS;

    //! Pointer to the water property calculator
    std::unique_ptr<WaterProps> m_waterProps;

    //! vector of size m_kk, used as a temporary holding area.
    mutable vector_fp m_tmpV;

    /**
     * Stoichiometric species charge -> This is for calculations
     * of the ionic strength which ignore ion-ion pairing into
     * neutral molecules. The Stoichiometric species charge is the
     * charge of one of the ion that would occur if the species broke
     * into two charged ion pairs.
     *  NaCl ->   m_speciesCharge_Stoich = -1;
     *  HSO4- -> H+ + SO42-              = -2
     *      -> The other charge is calculated.
     * For species that aren't ion pairs, it's equal to the
     * m_speciesCharge[] value.
     */
    vector_fp m_speciesCharge_Stoich;

    /**
     *  Array of 2D data used in the DHFORM_BETAIJ formulation
     *  Beta_ij.value(i,j) is the coefficient of the jth species
     *  for the specification of the chemical potential of the ith
     *  species.
     */
    Array2D m_Beta_ij;

    //! Logarithm of the activity coefficients on the molality scale.
    /*!
     *  mutable because we change this if the composition or temperature or
     *  pressure changes.
     */
    mutable vector_fp m_lnActCoeffMolal;

    //! Derivative of log act coeff wrt T
    mutable vector_fp m_dlnActCoeffMolaldT;

    //! 2nd Derivative of log act coeff wrt T
    mutable vector_fp m_d2lnActCoeffMolaldT2;

    //! Derivative of log act coeff wrt P
    mutable vector_fp m_dlnActCoeffMolaldP;

private:
    //! Calculate the log activity coefficients
    /*!
     * This function updates the internally stored natural logarithm of the
     * molality activity coefficients. This is the main routine for
     * implementing the activity coefficient formulation.
     */
    void s_update_lnMolalityActCoeff() const;

    //! Calculation of temperature derivative of activity coefficient
    /*!
     * Using internally stored values, this function calculates the temperature
     * derivative of the logarithm of the activity coefficient for all species
     * in the mechanism.
     *
     * We assume that the activity coefficients are current in this routine. The
     * solvent activity coefficient is on the molality scale. Its derivative is
     * too.
     */
    void s_update_dlnMolalityActCoeff_dT() const;

    //! Calculate the temperature 2nd derivative of the activity coefficient
    /*!
     * Using internally stored values, this function calculates the temperature
     * 2nd derivative of the logarithm of the activity coefficient for all
     * species in the mechanism.
     *
     * We assume that the activity coefficients are current in this routine.
     * Solvent activity coefficient is on the molality scale. Its derivatives
     * are too.
     */
    void s_update_d2lnMolalityActCoeff_dT2() const;

    //! Calculate the pressure derivative of the activity coefficient
    /*!
     * Using internally stored values, this function calculates the pressure
     * derivative of the logarithm of the activity coefficient for all species
     * in the mechanism.
     *
     * We assume that the activity coefficients, molalities, and A_Debye are
     * current. Solvent activity coefficient is on the molality scale. Its
     * derivatives are too.
     */
    void s_update_dlnMolalityActCoeff_dP() const;
};

}

#endif
