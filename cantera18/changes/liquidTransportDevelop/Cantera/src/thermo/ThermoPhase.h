/**
 *  @file ThermoPhase.h
 * Header file for class ThermoPhase, the base class for phases with
 * thermodynamic properties, and the text for the Module thermoprops
 * (see \ref thermoprops and class \link Cantera::ThermoPhase
 * ThermoPhase\endlink).
 */

/*
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#ifndef CT_THERMOPHASE_H
#define CT_THERMOPHASE_H

#include "Phase.h"


namespace Cantera {

  /*!
   * @name CONSTANTS - Specification of the Molality conventention
   */
  //@{
  //! Standard state uses the molar convention
  const int    cAC_CONVENTION_MOLAR    = 0;
  //! Standard state uses the molality convention
  const int    cAC_CONVENTION_MOLALITY = 1;
  //@}

  /*!
   * @name CONSTANTS - Specification of the SS conventention
   */
  //@{
  //! Standard state uses the molar convention
  const int cSS_CONVENTION_TEMPERATURE = 0;
  //! Standard state uses the molality convention
  const int cSS_CONVENTION_VPSS = 1;
  //! Standard state thermodynamics is obtained from slave %ThermoPhase objects
  const int cSS_CONVENTION_SLAVE = 2;
  //@}


  class XML_Node;

  /**
   * @defgroup thermoprops Thermodynamic Properties
   *
   *
   * These classes are used to compute the thermodynamic properties of
   * phases of matter. The main base class for describing thermodynamic 
   * properties of phases within %Cantera is called ThermoPhase. %ThermoPhase
   * is a large class that describes the interface within %Cantera to Thermodynamic
   * functions for a phase. 
   *
   * The calculation of thermodynamic functions within %ThermoPhase is 
   * broken down roughly into two or more steps. First, the standard state 
   * properties
   * of all of the species are calculated at the current temperature and at
   * either
   * the current pressure or at a reference pressure. If the calculation is
   * carried out at a reference pressure instead of at the current pressure
   * the calculation is called a "reference state properties" calculation,
   * just to make the distinction (even though it may be considered to be
   * a fixed-pressure standard-state calculation). The next step is to 
   * adjust the reference state calculation to the current pressure. The 
   * thermodynamic
   * functions then are considered to be at the standard state of each species.
   * Lastly the mixing contributions are added to arrive at the thermodynamic
   * functions for the solution.
   *
   * The %ThermoPhase class provides interfaces to thermodynamic properties 
   * calculated for 
   * the reference state of each species, the standard state values for 
   * each species, the thermodynamic functions for solution values, both
   * on a per mole of solution basis (i.e., enthalpy_mole()), on a per kg of
   * solution basis,  and on a
   * partial molar basis for each species (i.e.,  
   * getPartialMolarEnthalpies(double *hbar)).
   * At each level, functions for the enthalpy, entropy, Gibbs free energy,
   * internal energy, and volume are provided. So, 5 levels (reference state,
   * standard state, partial molar, per mole of solution, and per mass of 
   * solution) 
   * and 5 functions multiplied together makes 25 possible functions. That's 
   * why %ThermoPhase is such a large class.
   *
   *   <H3>
   *      Categorizing the Different %ThermoPhase Objects
   *   </H3>
   *
   *   ThermoPhase objects may be catelogged into four general bins.
   *
   *   The first type are those whose underlying species have a reference state associated
   *   with them. The reference state describes the thermodynamic functions for a
   *   species at a single reference pressure, \f$p_0\f$. The thermodynamic functions
   *   are specified via derived objects of the SpeciesThermoInterpType object class, and usually
   *   consist of polynomials in temperature such as the NASA polynomial or the SHOMATE
   *   polynomial.  Calculators for these
   *   reference states, which manage the calculation for all of the species
   *   in a phase, are all derived from the virtual base class SimpleThermo. Calculators
   *   are needed because the actual calculation of the reference state thermodynamics
   *   has been shown to be relatively expensive. A great deal of work has gone
   *   into devising efficient schemes for calculating the thermodynamic polynomials
   *   of a set of species in a phase, in particular gas species in ideal gas phases
   *   whose reference state thermodynamics is specified by NASA polynomials.
   * 
   *   The reference state thermodynamics combined with the mixing rules and
   *   an assumption about the pressure dependence yields the thermodynamic functions for
   *   the phase. 
   *   Expressions involving the specification of the fugacities of species would fall into 
   *   this category of %ThermoPhase objects. Note, however, that at this time, we do not
   *   have any nontrivial examples of these types of phases.
   *   In general, the independent variables that completely describe the state of the
   *   system  for this class are temperature, the
   *   phase density, and \f$ N - 1 \f$ species mole or mass fractions. 
   *   Additionally, if the
   *   phase involves charged species, the phase electric potential is an added independent variable.
   *   Examples of the first class of %ThermoPhase functions, which includes the
   *   IdealGasPhase object, the most commonly used object with %Cantera, are given below.
   *
   *    - IdealGasPhase       in IdealGasPhase.h
   *    - StoichSubstance     in StoichSubstance.h
   *    - SurfPhase           in SurfPhase.h
   *    - EdgePhase           in EdgePhase.h
   *    - LatticePhase        in LatticePhase.h
   *    - LatticeSolidPhase   in LatticeSolidPhase.h
   *    - ConstDensityThermo  in ConstDensityThermo.h
   *    - PureFluidPhase      in PureFluidPhase.h
   *    - IdealSolidSolnPhase in IdealSolidSolnPhase.h
   *    - VPStandardStateTP   in VPStandardStateTP.h
   *
   *   The second class of objects are actually all derivatives of the VPStandardState
   *   class listed above. These classes assume that there exists a standard state
   *   for each species in the phase, where the Thermodynamic functions are specified
   *   as a function of temperature and pressure.  Standard state objects for each
   *   species are all derived from the PDSS virtual base class. Calculators for these
   *   standard state, which coordinate the calculation for all of the species
   *   in a phase, are all derived from the virtual base class VPSSMgr.
   *   In turn, these standard states may employ reference state calculation to
   *   aid in their calculations. And the VPSSMgr calculators may also employ
   *   SimpleThermo calculators to help in calculating the properties for all of the
   *   species in a phase. However, there are some PDSS objects which do not employ
   *   reference state calculations. An example of this is real equation of state for
   *   liquid water used within the calculation of brine thermodynamcis. 
   *   In general, the independent variables that completely describe the state of the
   *   system  for this class are temperature, the
   *   phase pressure, and N - 1 species mole or mass fractions or molalities. 
   *    The standard state thermodynamics combined with the mixing rules yields
   *   the thermodynamic functions for the phase. Mixing rules are given in terms
   *   of specifying the molar-base activity coefficients or activities.
   *   Lists of phases which belong to this group are given below
   *
   *    - IdealSolnGasVPSS  in IdealSolnGasVPSS.h 
   *    - MolalityVPSSTP    in MolalityVPSSTP.h 
   *
   *   Note, the ideal gas and ideal solution approximations are lumped together
   *   in the class IdealSolnGasVPSS, because at this level they look alike having
   *   the same mixing rules with respect to the specification of the excess 
   *   thermodynamic properties.
   *
   *   The third class of objects are actually all derivatives of the MolalityVPSSTP
   *   object. They assume that the standard states are temperature and 
   *   pressure dependent. But, they also assume that the standard states are 
   *   molality-based. In other words they assume that the standard state of the solute
   *   species are in a pseudo state of 1 molality but at infinite dilution.
   *   A solvent must be specified in these calculations. The solvent is assumed
   *   to be species zero, and its standard state is the pure solvent state.
   *   Lists of phases which belong to this group are:
   *
   *   - DebyeHuckel     in DebyeHuckel.h
   *   - IdealMolalSoln  in IdealMolalSoln.h  
   *   - HMWSoln         in HMWSoln.h
   *  
   *   The fourth class of %ThermoPhase objects are stoichiometric phases.
   *   Stoichiometric phases are phases which consist of one and only one
   *   species. The class  SingleSpeciesTP is the base class for these
   *   substances. Within the class, the general %ThermoPhase interface is
   *   dumbed down so that phases consisting of one species may be
   *   succinctly described.
   *   These phases may have PDSS classes or SimpleThermo calculators associated
   *   with them.
   *   In general, the independent variables that completely describe the state of the
   *   system  for this class are temperature and either the
   *   phase density or the phase pressure.
   *   Lists of classes in this group are given below.
   *
   *   - StoichSubstanceSSTP  in StoichSubstanceSSTP.h
   *   - WaterSSTP            in WaterSSTP.h
   *
   *   The reader may note that there are duplications in functionality in the
   *   above lists. This is true. And, it's used for the internal verification of
   *   capabilities within %Cantera's unit tests.
   *
   *
   * <H3>
   * Setting the %State of the phase
   * </H3>
   *
   *   Typically, the way the ThermoPhase object works is that there are a set 
   *   of functions that set the state of the phase via setting the internal 
   *   independent variables. Then, there are another set of functions that
   *   query the thermodynamic functions evalulated at the current %State of the
   *   phase. Internally, most of the intermediate work generally occurs at the
   *   point where the internal state of the system is set and not at the time
   *   when individual thermodynamic functions are queried (though the actual
   *   breakdown in work is dependent on the individual derived ThermoPhase object).
   *   Therefore, for efficiency, the user should lump together queries of thermodynamic functions
   *   after setting the state. Moreover, in setting the state, if the
   *   density is the independent variable, the following order should be
   *   used:
   *
   *      - Set the temperature
   *      - Set the mole or mass fractions or set the molalities
   *      - set the pressure. 
   *
   *   For classes which inherit from VPStandardStateTP, the above order may
   *   be used, or the following order may be used. It's not important.
   *
   *      - Set the temperature
   *      - Set the pressure
   *      - Set the mole or mass fractions or set the molalities
   *
   *   The following functions are used to set the state:
   *
   *     <TABLE>
   *      <TR>
   *        <TD> \link ThermoPhase::setState_TPX() setState_TPX()\endlink </TD>
   *        <TD>    Sets the temperature, mole fractions and then the pressure
   *                of the phase. </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link ThermoPhase::setState_TPY() setState_TPY()\endlink    </TD>
   *        <TD>    Set the temperature, mass fractions and then the pressure
   *                of the phase. </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link MolalityVPSSTP::setState_TPM() setState_TPM()\endlink     </TD>
   *        <TD>    Set the temperature, solute molalities, and then the
   *                pressure of the phase. Only from %ThermoPhase objects which
   *                inherit from MolalityVPSSTP
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link ThermoPhase::setState_TP() setState_TP()\endlink     </TD>
   *        <TD>    Set the temperature, and then the pressure
   *                of the phase. The mole fractions are assumed fixed. 
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link ThermoPhase::setState_PX() setState_PX()\endlink     </TD>
   *        <TD>    Set the mole fractions and then the pressure
   *                         of the phase. The temperature is assumed fixed.
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link ThermoPhase::setState_PY() setState_PY()\endlink     </TD>
   *        <TD>  Set the mass fractions and then the pressure
   *                         of the phase. The temperature is assumed fixed.
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link ThermoPhase::setState_HP() setState_HP()\endlink     </TD>
   *        <TD> Set the total specific enthalpy and the pressure
   *                         of the phase using an iterative process.
   *                         The mole fractions are assumed fixed
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link ThermoPhase::setState_UV() setState_UV()\endlink     </TD>
   *        <TD>  Set the total specific internal energy and the pressure
   *                         of the phase using an iterative process.
   *                         The mole fractions are assumed fixed.
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link ThermoPhase::setState_SP() setState_SP()\endlink     </TD>
   *        <TD>  Set the total specific internal energy and the pressure
   *                         of the phase using an iterative process.
   *                         The mole fractions are assumed fixed.
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link ThermoPhase::setState_HPX() setState_HP()\endlink     </TD>
   *        <TD> Set the total specific enthalpy and the pressure
   *                         of the phase using an iterative process.
   *                         The mole fractions are assumed fixed
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link ThermoPhase::setState_UVX() setState_UV()\endlink     </TD>
   *        <TD>  Set the total specific internal energy and the pressure
   *                         of the phase using an iterative process.
   *                         The mole fractions are assumed fixed.
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link State::setConcentrations() setConcentrations()\endlink     </TD>
   *        <TD> Set the concentrations of all the species in the
   *                         phase. Note this implicitly specifies the pressure and
   *                         density of the phase. The temperature is assumed fixed.
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link State::setDensity() setDensity()\endlink     </TD>
   *        <TD>  Set the total density of the phase. The temperature and
   *                         mole fractions are assumed fixed. Note this implicity
   *                         sets the pressure of the phase.
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link State::setTemperature() setTemperature()\endlink     </TD>
   *        <TD> Set the temperature of the phase. The density and
   *                         the mole fractions of the phase are fixed.
   *        </TD>
   *      </TR>
   *      <TR>
   *        <TD> \link ThermoPhase::setToEquilState() setToEquilState()\endlink     </TD>
   *        <TD>  Sets the mole fractions of the phase to their 
   *                         equilibrium values assuming fixed temperature and
   *                         total density.
   *        </TD>
   *      </TR>
   *    </TABLE>
   *
   *
   * 
   *  Some of the functions, like setState_TPX() have multiple forms depending upon
   *  the format for how the species compositions are set. 
   * 
   *
   *  Molar Basis vs. Molality Basis
   *
   * <H3>
   * Mechanical properties
   * </H3>
   *
   *  The %ThermoPhase object specifies the mechanical equation of state of the
   *  phase. Functions which are defined at the %ThermoPhase level to give the
   *  user more information about the mechanical properties are:
   *
   *       - ThermoPhase::pressure()
   *       - ThermoPhase::isothermalCompressibility()
   *       - ThermoPhase::thermalExpansionCoeff()
   *       .
   *
   * <H3>
   * Treatment of the %Phase Potential and the electrochemical potential of a species
   * </H3>
   *
   *  The electrochemical potential of species k in a phase p, \f$ \zeta_k \f$,
   *  is related to the chemical potential via
   *  the following equation, 
   * 
   *       \f[
   *            \zeta_{k}(T,P) = \mu_{k}(T,P) + z_k \phi_p 
   *       \f]
   *
   *   where  \f$ \nu_k \f$ is the charge of species k, and \f$ \phi_p \f$ is
   *   the electric potential of phase p.
   *      
   *  The potential  \f$ \phi_p \f$ is tracked and internally storred within
   *  the base %ThermoPhase object. It constitutes a specification of the
   *  internal state of the phase; it's the third state variable, the first 
   *  two being temperature and density (or, pressure, for incompressible
   *  equations of state). It may be set with the function,
   *  ThermoPhase::setElectricPotential(),
   *  and may be queried with the function ThermoPhase::electricPotential().
   *  
   *  Note, the overall electrochemical potential of a phase may not be
   *  changed by the potential because many phases enforce charge
   *  neutrality:
   *
   *       \f[
   *            0 = \sum_k z_k X_k
   *       \f]
   *
   *  Whether charge neutrality is necessary for a phase is also specified
   *  within the ThermoPhase object, by the function call
   *  ThermoPhase::chargeNeutralityNecessary(). Note, that it is not
   *  necessary for the IdealGas phase, currently. However, it is
   *  necessary for liquid phases such as Cantera::DebyeHuckel and
   *  Cantera::HMWSoln for the proper specification of the chemical potentials.
   *
   *
   *  This equation, when applied to the \f$ \zeta_k \f$ equation described
   *  above, results in a zero net change in the effective Gibbs free
   *  energy of the phase. However, specific charged species in the phase
   *  may increase or decrease their electochemical potentials, which will
   *  have an effect on interfacial reactions involving charged species,
   *  when there is a potential drop between phases. This effect is used
   *  within the Cantera::InterfaceKinetics and Cantera::EdgeKinetics kinetics 
   *  objects classes.
   *  
   *
   *  Other internal state variables, that track the treatment of other
   *  potential energy contributions, by adding contributions to the
   *  chemical potential to create an effective chemical potential, 
   *  may be added at a later time.
   *
   *  <H3>
   *   Specification of Activities and Activity Conventions
   *  </H3>
   * 
   *
   * The activity \f$a_k\f$ and activity coefficient \f$ \gamma_k \f$ of a 
   * species in solution is related to the chemical potential by 
   *
   * \f[ 
   *    \mu_k = \mu_k^0(T,P) + \hat R T \log a_k.= \mu_k^0(T,P) + \hat R T \log x_k \gamma_k
   * \f] 
   *
   * The quantity \f$\mu_k^0(T,P)\f$ is
   * the standard chemical potential at unit activity, 
   * which depends on the temperature and pressure, 
   * but not on the composition. The
   * activity is dimensionless. Within liquid electrolytes its common to use a 
   * molality convention, where solute species employ the molality-based
   * activity coefficients:
   *
   * \f[
   *  \mu_k =  \mu_k^\triangle(T,P) + R T ln(a_k^{\triangle}) = 
   *            \mu_k^\triangle(T,P) + R T ln(\gamma_k^{\triangle} \frac{m_k}{m^\triangle}) 
   * \f]
   *
   * And, the solvent employs the following convention
   * \f[
   *    \mu_o = \mu^o_o(T,P) + RT ln(a_o) 
   * \f]
   *
   * where \f$ a_o \f$ is often redefined in terms of the osmotic coefficient \f$ \phi \f$.
   *
   *   \f[
   *       \phi = \frac{- ln(a_o)}{\tilde{M}_o \sum_{i \ne o} m_i}
   *   \f]
   *
   *  %ThermoPhase classes which employ the molality based convention are all derived
   *  from the MolalityVPSSTP class. See the class description for further information
   *  on its capabilities. 
   *
   *  The activity convention used by a %ThermoPhase object 
   *  may be queried via the ThermoPhase::activityConvention() function. A zero means molar based,
   *  while a one means molality based.
   *
   *  The function ThermoPhase::getActivities() returns a vector of activities. Whether these are
   *  molar-based or molality-based depends on the value of activityConvention().
   *
   *  The function getActivityCoefficients() always returns molar-based activity
   *  coefficients regardless of the activity convention used. The function
   *  MolalityVPSSTP::getMolalityActivityCoefficients() returns molality
   *  based activity coefficients for those ThermoPhase objects derived
   *  from the MolalityVPSSTP class. The function MolalityVPSSTP::osmoticCoefficient()
   *  returns the osmotic coefficient.
   *
   *  <H3>
   *   Activity Concentrations: Relationship of %ThermoPhase to %Kinetics Expressions
   *  </H3>
   *
   *   %Cantera can handle both thermodynamics and kinetics mechanisms. Reversible
   *   kinetics
   *   mechanisms within %Cantera must be compatible with thermodynamics in the
   *   sense that at equilibrium, or at infinite times, the concentrations
   *   of species must conform to thermodynamics. This means that for every
   *   valid reversible kinetics reaction in a mechanism, it must be reducible to 
   *   an expression involving the ratio of the product activity to
   *   the reactant activities being equal to the exponential of the
   *   dimensionless standard state gibbs free energies of reaction.
   *   Irreversible kinetics reactions do not have this requirement; however,
   *   their usage can yield unexpected and inconsistent results in many
   *   situations.
   *   The actual units used in a kinetics expression depend
   *   on the context or the relative field of study. For example, in 
   *   gas phase kinetics, species in kinetics expressions are expressed in
   *   terms of concentrations, i.e., gmol cm-3. In solid phase studies,
   *   however, kinetics is usually expressed in terms of unitless activities,
   *   which most often equate to solid phase mole fractions. In order to 
   *   accomodate variability here, %Cantera has come up with the idea
   *   of activity concentrations, \f$ C^a_k \f$. Activity concentrations are the expressions
   *   used directly in kinetics expressions.    
   *   These activity (or generalized) concentrations are used
   *   by kinetics manager classes to compute the forward and
   *   reverse rates of elementary reactions. Note that they may
   *   or may not have units of concentration --- they might be
   *   partial pressures, mole fractions, or surface coverages,
   *   The activity concentrations for species <I>k</I>, \f$ C^a_k \f$, are 
   *   related to the activity for species, k, \f$ a_k \f$,
   *   via the following expression:
   *
   *   \f[
   *       a_k = C^a_k / C^0_k
   *   \f]
   *
   *  \f$ C^0_k \f$ are called standard concentrations. They serve as multiplicative factors
   *  bewteen the activities and the generalized concentrations. Standard concentrations
   *  may be different for each species. They may depend on both the temperature
   *  and the pressure. However, they may not depend
   *  on the composition of the phase. For example, for the IdealGasPhase object
   *  the standard concentration is defined as
   *
   *  \f[
   *     C^0_k = P/ R T 
   *  \f]
   *  
   *  In many solid phase kinetics problems,
   *
   *   \f[
   *     C^0_k = 1.0 ,
   *  \f]
   *
   *  is employed making the units for activity concentrations in solids unitless.
   *  
   *  %ThermoPhase member functions dealing with this concept include 
   *  ThermoPhase::getActivityConcentrations() , which provides a vector of the current
   *  activity concentrations. The function, ThermoPhase::standardConcentration(int k=0) returns
   *  the standard concentration of the kth species. The function,
   *  ThermoPhase::logStandardConc(int k=0), returns the natural log of the kth standard
   *  concentration. The function  ThermoPhase::getUnitsStandardConc() returns a vector of 
   *  doubles, specifying the MKS units of the standard concentration of the
   *  kth species.
   *
   *
   *  <H3>
   *  Initialization of %ThermoPhase Objects within %Cantera
   *  </H3>
   *
   *  Instantiation of %ThermoPhase properties occurs by reading and
   *  processing the XML data contained within an ctxml data file.
   *  First a call to  newPhase(std::string file, std::string id) or
   *  newPhase(XML_Node &phase)
   *  is made. The arguments serve to specify the
   *  XML data structure containing the phase information.
   *
   *  Within newPhase() a determination of what type of %ThermoPhase object should be 
   *  used is made. This is done within the routine ThermoFactory::newThermoPhase(std::string model)
   *  or related routines.
   *  Once the correct %ThermoPhase derived object is selected and instantiated with a 
   *  bare constructor, the
   *  function Cantera::importPhase() is called with the  %ThermoPhase derived object as
   *  one of its arguments.
   *
   *  Within importPhase(), a decision is made as to what type of
   *  standard state, i.e.,
   *  either a reference state (just T dependent)  or a standard state
   *  (both P and T dependent), is to be used to calculate the
   *  standard state properties of the species within the phase.
   *  If only a reference state is needed 
   *  then a call to 
   *  \link #newSpeciesThermoMgr(std::vector<XML_Node*> spData_nodes, SpeciesThermoFactory* f=0, bool opt=false) newSpeciesThermoMgr()\endlink
   *  is made in order
   *  pick a manager, i.e., a derivative of the SpeciesThermo 
   *  object, to use.
   *
   *  If a temperature and pressure dependent standard state is needed 
   *  then a call to VPSSMgrFactory::newVPSSMgr()
   *  is made in order
   *  pick a manager, i.e., a derivative of the VPSSMgr 
   *  object, to use. Along with the VPSSMgr designation comes a 
   *  determination of whether there is an accompanying  SpeciesThermo 
   *  and what type of SpeciesThermo object to use in the
   *  VPSSMgr calculations.
   *  
   *  Once these determinations are made, the %ThermoPhase object is
   *  ready to start reading in the species information, which includes
   *  all of the available standard state information about the
   *  species. this is done within the routine installSpecies().
   *  
   *  Within installSpecies(), most of the common steps for adding a
   *  species are carried out. The element stoichiometry is read
   *  and elements are added as needed to the list of elements
   *  kept with the ThermoPhase object. The charge of the species
   *  is read in. The species is added into the list
   *  of species kept within the ThermoPhase object. Lastly, the
   *  standard state thermodynamics for the species is read in.
   *  For reference states, the routine, SpeciesThermoFactory::installThermoForSpecies(),
   *  is used to read in the data. Essentially, this routine is a
   *  factory routine for picking the correct subroutine to 
   *  call to read the XML data from the input file and install the
   *  correct  SpeciesThermoInterpType object into the SpeciesThermo object. 
   *
   *  Within installSpecies(), for standard states, the routine, 
   *  SpeciesThermoFactory::installVPThermoForSpecies() is
   *  called. However, this is just a shell routine for calling
   *  the VPSSMgr's derived VPSSMgr::createInstallPDSS() routine.
   *  Within the  VPSSMgr::createInstallPDSS() routine of the derived VPSSMgr's
   *  object, the XML data from the input file is read and the
   *  calculations for the species standard state is installed.
   *  Additionally, the  derived PDSS object is created and installed
   *  into the VPStandardStateTP list containing all of the PDSS objects
   *  for that phase.
   *
   *  Now that all of the species standard states are read in and
   *  installed into the ThermoPhase object, control once again
   *  is returned to the importPhase() function.  Two derived functions
   *  are then called. The first one, ThermoPhase::initThermo(), is called. In this
   *  routine, all internal arrays within the %ThermoPhase object are
   *  dimensioned according to the number of elements and species.
   *  Then, the function ThermoPhase::initThermoXML() is called.
   *  This function is tasked with reading in all of the thermodynamic
   *  function information specific to the calculation of the
   *  phase information. This includes all of the information about
   *  the activity coefficient calculation.
   *
   *  After the ThermoPhase::initThermoXML() is finished, the
   *  ThermoPhase routine is ready to receive requests for 
   *  thermodynamic property information.
   *  
   *
   *  There is an alternative way to instantiate %ThermoPhase objects that
   *  is applicable to a significant proportion of %ThermoPhase classes.
   *  The phase may be instantiated via a constructor that invokes the
   *  XML data structure wherein the phase information is to be read directly.
   *  In this case, the call to newPhase() and the call to 
   *  ThermoFactory::newThermoPhase(std::string model)
   *  is not made. However, soon after that, the call to importPhase() is
   *  made and thereafter instantiation follows the initialization course described
   *  previously in order to avoid as much duplicate code as possible.
   *  This alternative way to  instantiate %ThermoPhase objects has the 
   *  advantage of working well with hard-coded situations. And, it
   *  works well also with situations where new %ThermoPhase classes 
   *  are being developed and haven't yet made their way into the
   *  factory routines.
   *
   *  <H3>
   *  Adding Additional Thermodynamics Models
   *  </H3>
   *
   *  In general, factory routines throw specific errors when encountering
   *  unknown thermodynamics models in XML files. All of the error classes
   *  derive from the class, CanteraError.
   *  The newVPSSMgr() routines throws the UnknownVPSSMgr class error when
   *  they encounter an unknown string in the XML input file specifying the
   *  VPSSMgr class to use.
   *
   *  Many of the important member functions in factory routines are
   *  virtual classes. This means that a user may write their own
   *  factory classes which inherit from the base %Cantera factory classes
   *  to provide additional %ThermoPhase classes.
   * 
   *
   * @see newPhase(std::string file, std::string id) Description for how to
   *               read ThermoPhases from XML files.
   * @see newPhase(XML_Node &phase) How to call the Factory routine to create 
   *          and initialize %ThermoPhase objects.
   * @ingroup phases
   */

  
  //!   Base class for a phase with thermodynamic properties. 
  /*!
   * Class %ThermoPhase is the base class for the family of classes
   * that represent phases of matter of any type. It defines a
   * common public interface, and implements a few methods. Most of
   * the methods, however, are declared virtual and are meant to be
   * overloaded in derived classes.  The standard way used
   * throughout Cantera to compute properties of phases of matter is
   * through pointers of type ThermoPhase* that point to objects of
   * subclasses of ThermoPhase.
   * 
   * Class %ThermoPhase extends class Phase by adding methods to compute 
   * thermodynamic
   * properties in addition to the ones (temperature, density,
   * composition) that class Phase provides. The distinction is that
   * the methods declared in ThermoPhase require knowing the
   * particular equation of state of the phase of interest, while
   * those of class Phase do not, since they only involve data values
   * stored within the object.
   *
   * Instances of subclasses of %ThermoPhase should be created using
   * the factory class ThermoFactory, not by calling the constructor
   * directly. This allows new classes to be used with the various
   * Cantera language interfaces.
   * 
   * To implement a new equation of state, derive a class from
   * ThermoPhase and overload the virtual methods in
   * ThermoPhase. Methods that are not needed can be left
   * unimplimented, which will cause an exception to be thrown if it
   * is called.
   *
   * Relationship with the kinetics operator:
   *
   * Describe activity coefficients.
   *
   * Describe K_a, K_p, and K_c, These are three different equilibrium
   * constants. 
   *
   *   K_a is the calculation of the equilibrium constant from the
   *   standard state Gibbs free energy values. It is by definition
   *   dimensionless.
   *
   *   K_p is the calculation of the equilibrium constant from the
   *   reference state gibbs free energy values. It is by definition
   *   dimensionless. The pressure dependence is handled entirely
   *   on the rhs of the equilibrium expression.
   *
   *   K_c is the equilibrium constant calculated from the
   *   activity concentrations. The dimensions depend on the number
   *   of products and reactants.
   *    
   *
   * The kinetics manager requires the calculation of K_c for the
   * calculation of the reverse rate constant
   * 
   *
   * @ingroup thermoprops
   * @ingroup phases
   */
  class ThermoPhase : public Phase {

    public:

    //! Constructor. Note that ThermoPhase is meant to be used as
    //! a base class, so this constructor should not be called
    //! explicitly.
    ThermoPhase();

    //! Destructor. Deletes the species thermo manager.
    virtual ~ThermoPhase();
   
    //!Copy Constructor for the %ThermoPhase object. 
    /*!
     * @param right  ThermoPhase to be copied
     */
    ThermoPhase(const ThermoPhase &right);
	
    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %ThermoPhase object to be copied into the
     *                 current one. 
     */
    ThermoPhase& operator=(const ThermoPhase &right);

     //! Duplication routine for objects which inherit from 
     //!  ThermoPhase.
     /*!
     *  This virtual routine can be used to duplicate %ThermoPhase objects
     *  inherited from %ThermoPhase even if the application only has
     *  a pointer to %ThermoPhase to work with.
     * 
     *  These routines are basically wrappers around the derived copy
     *  constructor.
     */
    virtual ThermoPhase *duplMyselfAsThermoPhase() const;
    
    /**
     *   
     * @name  Information Methods  
     * @{
     */
        
    //! Equation of state type flag.
    /*!
     *  The base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Constants defined for this purpose are
     * listed in mix_defs.h.
     */
    virtual int eosType() const { return 0; }
    
    /**
     * Returns the reference pressure in Pa. This function is a wrapper
     * that calls the species thermo refPressure function.
     */
    virtual doublereal refPressure() const {
      return m_spthermo->refPressure();
    }

        
    //! Minimum temperature for which the thermodynamic data for the species 
    //! or phase are valid.
    /*!
     * If no argument is supplied, the
     * value returned will be the lowest temperature at which the
     * data for \e all species are valid. Otherwise, the value
     * will be only for species \a k. This function is a wrapper
     * that calls the species thermo minTemp function.
     *
     * @param k index of the species. Default is -1, which will return the max of the min value
     *          over all species.
     */
    virtual doublereal minTemp(int k = -1) const {
      return m_spthermo->minTemp(k);
    }
        
#ifdef H298MODIFY_CAPABILITY

    //! Report the 298 K Heat of Formation of the standard state of one species (J kmol-1)
    /*!
     *   The 298K Heat of Formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param k    species index
     *   @return     Returns the current value of the Heat of Formation at 298K and 1 bar
     */
    doublereal Hf298SS(const int k) const {
      return (m_spthermo->reportOneHf298(k));
    }

    //! Modify the value of the 298 K Heat of Formation of one species in the phase (J kmol-1)
    /*!
     *   The 298K heat of formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param  k           Species k
     *   @param  Hf298New    Specify the new value of the Heat of Formation at 298K and 1 bar                      
     */
    virtual void modifyOneHf298SS(const int k, const doublereal Hf298New) {
      m_spthermo->modifyOneHf298(k, Hf298New);
    }

#else

    //! Report the 298 K Heat of Formation of the standard state of one species (J kmol-1)
    /*!
     *   The 298K Heat of Formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param k    species index
     *   @return     Returns the current value of the Heat of Formation at 298K and 1 bar
     */
    doublereal Hf298SS(const int k) const {
       return err("Hf298SS - H298MODIFY_CAPABILITY not compiled in"); 
    }

    //! Modify the value of the 298 K Heat of Formation of one species in the phase (J kmol-1)
    /*!
     *   The 298K heat of formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param  k           Species k
     *   @param  Hf298New    Specify the new value of the Heat of Formation at 298K and 1 bar                      
     */
    virtual void modifyOneHf298SS(const int k, const doublereal Hf298New) {
      (void) err("Hf298SS - H298MODIFY_CAPABILITY not compiled in"); 
    }
#endif

    //! Maximum temperature for which the thermodynamic data for the species 
    //! are valid. 
    /*!
     * If no argument is supplied, the
     * value returned will be the highest temperature at which the
     * data for \e all species are valid. Otherwise, the value
     * will be only for species \a k. This function is a wrapper
     * that calls the species thermo maxTemp function.
     *
     * @param k index of the species. Default is -1, which will return the min of the max value
     *          over all species.
     */
    virtual doublereal maxTemp(int k = -1) const {
      return m_spthermo->maxTemp(k);
    }

    //! Returns the chargeNeutralityNecessity boolean
    /*!
     * Some phases must have zero net charge in order for their thermodynamics functions to be valid.
     * If this is so, then the value returned from this function is true.
     * If this is not the case, then this is false. Now, ideal gases have this parameter set to false,
     * while solution with  molality-based activity coefficients have this parameter set to true. 
     */
    bool chargeNeutralityNecessary() const {
      return m_chargeNeutralityNecessary;
    }
      
    /**
     * @} 
     * @name  Molar Thermodynamic Properties of the Solution
     * @{
     */

    /// Molar enthalpy. Units: J/kmol. 
    virtual doublereal enthalpy_mole() const {
      return err("enthalpy_mole");
    }

    /// Molar internal energy. Units: J/kmol. 
    virtual doublereal intEnergy_mole() const {
      return enthalpy_mole() - pressure()* molarVolume();
    }

    /// Molar entropy. Units: J/kmol/K. 
    virtual doublereal entropy_mole() const {
      return err("entropy_mole");
    }

    /// Molar Gibbs function. Units: J/kmol. 
    virtual doublereal gibbs_mole() const {
      return enthalpy_mole() - temperature()*entropy_mole();
    }

    /// Molar heat capacity at constant pressure. Units: J/kmol/K. 
    virtual doublereal cp_mole() const {
      return err("cp_mole");
    }

    /// Molar heat capacity at constant volume. Units: J/kmol/K. 
    virtual doublereal cv_mole() const {
      return err("cv_mole");
    }


    /**
     * @}
     * @name Mechanical Properties
     * @{
     */
       
    //! Return the thermodynamic pressure (Pa).
    /*!
     *  This method must be overloaded in derived classes. Since the
     *  mass density, temperature, and mass fractions are stored,
     *  this method should use these values to implement the
     *  mechanical equation of state \f$ P(T, \rho, Y_1, \dots,
     *  Y_K) \f$.
     */
    virtual doublereal pressure() const {
      return err("pressure");
    }

    //! Set the internally storred pressure (Pa) at constant
    //! temperature and composition
    /*!
     *   This method must be reimplemented in derived classes, where it
     *   may involve the solution of a nonlinear equation. Within %Cantera,
     *   the independent variable is the density. Therefore, this function
     *   solves for the density that will yield the desired input pressure.
     *   The temperature and composition iare held constant during this process.
     *
     *  This base class function will print an error, if not overwritten.
     *
     *  @param p input Pressure (Pa)
     */
    virtual void setPressure(doublereal p) {
      err("setPressure");
    }
      
    //! Returns  the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     *  or
     * \f[
     * \kappa_T = \frac{1}{\rho}\left(\frac{\partial \rho}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const {
      err("isothermalCompressibility"); return -1.0;
    }

    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    /*!      
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal thermalExpansionCoeff() const {
      err("thermalExpansionCoeff()"); return -1.0;
    }

    /// @deprecated
    virtual void updateDensity() {
      deprecatedMethod("ThermoPhase","updateDensity","");
    }

    /**
     * @} 
     * @name Electric Potential
     * 
     * The phase may be at some non-zero electrical
     * potential. These methods set or get the value of the
     * electric potential.
     */
    //@{
      
    //! Set the electric potential of this phase (V).
    /*!
     * This is used by classes InterfaceKinetics and EdgeKinetics to
     * compute the rates of charge-transfer reactions, and in computing
     * the electrochemical potentials of the species.
     * 
     * Each phase may have its own electric potential.
     *
     * @param v Input value of the electric potential in Volts
     */
    void setElectricPotential(doublereal v) {
      m_phi = v;
    }

    //! Returns the electric potential of this phase (V).
    /*!
     *  Units are Volts (which are Joules/coulomb)
     */ 
    doublereal electricPotential() const { return m_phi; }

    /**
     * @}
     * @name Activities, Standard States, and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is related
     * to the chemical potential by \f[ \mu_k = \mu_k^0(T,P) +
     * \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T,P)\f$ is
     * the standard chemical potential at unit activity, 
     * which depends on  temperature and pressure, 
     * but not on composition. The
     * activity is dimensionless.
     * @{
     */

	
    //! This method returns the convention used in specification
    //! of the activities, of which there are currently two, molar-
    //! and molality-based conventions.
    /*!
     * Currently, there are two activity conventions:
     *  - Molar-based activities
     *       %Unit activity of species at either a hypothetical pure
     *       solution of the species or at a hypothetical
     *       pure ideal solution at infinite dilution
     *   cAC_CONVENTION_MOLAR 0
     *      - default
     *  
     *  - Molality-based acvtivities
     *       (unit activity of solutes at a hypothetical 1 molal
     *        solution referenced to infinite dilution at all
     *        pressures and temperatures).
     *   cAC_CONVENTION_MOLALITY 1
     */
    virtual int activityConvention() const;
    
    //! This method returns the convention used in specification
    //! of the standard state, of which there are currently two,
    //! temperature based, and variable pressure based.
    /*!
     * Currently, there are two standard state conventions:
     *  - Temperature-based activities
     *   cSS_CONVENTION_TEMPERATURE 0
     *      - default
     *
     *  -  Variable Pressure and Temperature -based activities
     *   cSS_CONVENTION_VPSS 1
     *
     *  -  Thermodynamics is set via slave ThermoPhase objects with
     *     nothing being carried out at this %ThermoPhase object level
     *   cSS_CONVENTION_SLAVE 2
     */
    virtual int standardStateConvention() const;
        
    //! This method returns an array of generalized concentrations
    /*!
     * \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k /
     * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration
     * defined below and \f$ a_k \f$ are activities used in the
     * thermodynamic functions.  These activity (or generalized)
     * concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions. Note that they may
     * or may not have units of concentration --- they might be
     * partial pressures, mole fractions, or surface coverages,
     * for example.
     *
     * @param c Output array of generalized concentrations. The 
     *           units depend upon the implementation of the
     *           reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const {
      err("getActivityConcentrations");
    }

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the activity (i.e., generalized) concentration. In many cases, this quantity
     * will be the same for all species in a phase - for example,
     * for an ideal gas \f$ C^0_k = P/\hat R T \f$. For this
     * reason, this method returns a single value, instead of an
     * array.  However, for phases in which the standard
     * concentration is species-specific (e.g. surface species of
     * different sizes), this method may be called with an
     * optional parameter indicating the species.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return 
     *   Returns the standard concentration. The units are by definition
     *   dependent on the ThermoPhase and kinetics manager representation.
     */
    virtual doublereal standardConcentration(int k=0) const {
      err("standardConcentration");
      return -1.0;
    }  
        
    //! Natural logarithm of the standard concentration of the kth species.
    /*!
     * @param k    index of the species (defaults to zero)
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
     * Note, for molality based formulations, this returns the 
     * molality based activities.
     *
     * We resolve this function at this level by calling
     * on the activityConcentration function. However, 
     * derived classes may want to override this default
     * implementation.
     *
     * @param a   Output vector of activities. Length: m_kk.
     */
    virtual void getActivities(doublereal* a) const;
    
    //! Get the array of non-dimensional molar-based activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const {
      if (m_kk == 1) {
	ac[0] = 1.0;
	} else {
	  err("getActivityCoefficients");
	}
    }

    //! Get the array of non-dimensional molar-based ln activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param lnac Output vector of ln activity coefficients. Length: m_kk.
     */
    virtual void getLnActivityCoefficients(doublereal * const lnac) const;
      
    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{
      
    /**
     * Get the array of non-dimensional species chemical potentials
     * These are partial molar Gibbs free energies.
     * \f$ \mu_k / \hat R T \f$.
     * Units: unitless
     *
     * @param mu  Output vector of dimensionless chemical potentials.
     *            Length: m_kk.
     */
    virtual void getChemPotentials_RT(doublereal* mu) const {
      err("getChemPotentials_RT");
    }
      
     
    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the 
     * species in solution at the current temperature, pressure
     * and mole fraction of the solution.
     *
     * @param mu  Output vector of species chemical 
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const {
      err("getChemPotentials");
    }
    
    //!  Get the species electrochemical potentials. 
    /*!
     *  These are partial molar quantities.  This method adds a term \f$ F z_k
     *  \phi_p \f$ to each chemical potential.  
     *  The electrochemical potential of species k in a phase p, \f$ \zeta_k \f$,
     *  is related to the chemical potential via
     *  the following equation, 
     *
     *       \f[
     *            \zeta_{k}(T,P) = \mu_{k}(T,P) + F z_k \phi_p 
     *       \f]
     *
     * @param mu  Output vector of species electrochemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    void getElectrochemPotentials(doublereal* mu) const {
      getChemPotentials(mu);
      double ve = Faraday * electricPotential();
      for (int k = 0; k < m_kk; k++) {
	mu[k] += ve*charge(k);
      }
    }
    
    //! Returns an array of partial molar enthalpies for the species
    //! in the mixture. Units (J/kmol)
    /*!
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const {
      err("getPartialMolarEnthalpies");
    }
      
    //! Returns an array of partial molar entropies of the species in the
    //! solution. Units: J/kmol/K.
    /*!
     * @param sbar    Output vector of species partial molar entropies.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const {
      err("getPartialMolarEntropies");
    }

    //! Return an array of partial molar internal energies for the 
    //! species in the mixture.  Units: J/kmol.
    /*!
     * @param ubar    Output vector of speciar partial molar internal energies.
     *                Length = m_kk. units are J/kmol.
     */
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const {
      err("getPartialMolarIntEnergies");
    }
    
    //! Return an array of partial molar heat capacities for the
    //! species in the mixture.  Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat 
     *                capacities at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const {
      err("getPartialMolarCp");
    }
        
    //! Return an array of partial molar volumes for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  @param vbar   Output vector of speciar partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const {
      err("getPartialMolarVolumes");
    }

    //! Return an array of derivatives of partial molar volumes wrt temperature for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  The derivative is at constant pressure
     *
     *  @param d_vbar_dT   Output vector of derivatives of species partial molar volumes wrt T.
     *                     Length = m_kk. units are m^3/kmol/K.
     */
    virtual void getdPartialMolarVolumes_dT(doublereal* d_vbar_dT) const {
      err("getdPartialMolarVolumes_dT");
    }

    //! Return an array of derivatives of partial molar volumes wrt pressure  for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  The derivative is at constant temperature
     *
     *  @param d_vbar_dP   Output vector of derivatives of species partial molar volumes wrt P.
     *                     Length = m_kk. units are m^3/kmol/Pa.
     */
    virtual void getdPartialMolarVolumes_dP(doublereal* d_vbar_dP) const {
      err("getdPartialMolarVolumes_dP");
    }

    //@}
    /// @name Properties of the Standard State of the Species in the Solution 
    //@{
    
    //! Get the array of chemical potentials at unit activity for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current
     * temperature and pressure of the solution
     *
     * @param mu      Output vector of chemical potentials. 
     *                Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu) const {
      err("getStandardChemPotentials");
      }
       
    //! Get the nondimensional Enthalpy functions for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param hrt      Output vector of  nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const {
      err("getEnthalpy_RT");
    }

    //! Get the array of nondimensional Entropy functions for the
    //! standard state species at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param sr   Output vector of  nondimensional standard state entropies.
     *             Length: m_kk.
     */
    virtual void getEntropy_R(doublereal* sr) const {
      err("getEntropy_R");
    }

    //! Get the nondimensional Gibbs functions for the species
    //! in their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param grt  Output vector of nondimensional standard state gibbs free energies
     *             Length: m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const {
      err("getGibbs_RT");
    }
      
    //! Get the Gibbs functions for the standard
    //! state of the species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * Units are Joules/kmol
     * @param gpure  Output vector of  standard state gibbs free energies
     *               Length: m_kk.
     */
    virtual void getPureGibbs(doublereal* gpure) const {
      err("getPureGibbs");
    }
      
    //!  Returns the vector of nondimensional Internal Energies  of the standard
    //!  state species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * @param urt  output vector of nondimensional standard state internal energies
     *             of the species. Length: m_kk. 
     */
    virtual void getIntEnergy_RT(doublereal *urt) const {
      err("getIntEnergy_RT");
    }
      
    //! Get the nondimensional Heat Capacities at constant
    //! pressure for the species standard states
    //! at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * @param cpr   Output vector of nondimensional standard state heat capacities
     *              Length: m_kk.
     */
    virtual void getCp_R(doublereal* cpr) const {
      err("getCp_R");
    }
      
    //!  Get the molar volumes of the species standard states at the current
    //!  <I>T</I> and <I>P</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes(doublereal *vol) const {
      err("getStandardVolumes");
    }

    //!  Get the derivative of the molar volumes of the species standard states wrt temperature at the current
    //!  <I>T</I> and <I>P</I> of the solution.
    /*!  
     *  The derivative is at constant pressure
     *   units = m^3 / kmol / K
     *
     * @param d_vol_dT Output vector containing derivatives of standard state volumes wrt T
     *                 Length: m_kk.
     */
    virtual void getdStandardVolumes_dT(doublereal *d_vol_dT) const {
      err("getdStandardVolumes_dT");
    }

    //!  Get the derivative molar volumes of the species standard states wrt pressure at the current
    //!  <I>T</I> and <I>P</I> of the solution.
    /*!
     *  The derivative is at constant temperature.
     * units = m^3 / kmol / Pa
     *
     * @param d_vol_dP  Output vector containing the derivative of standard state volumes wrt P.
     *                  Length: m_kk.
     */
    virtual void getdStandardVolumes_dP(doublereal *d_vol_dP) const {
      err("getdStandardVolumes_dP");
    }

    //@}
    /// @name Thermodynamic Values for the Species Reference States 
    //@{

    //!  Returns the vector of nondimensional
    //!  enthalpies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     *  This base function will throw a CanteraException unless
     *  it is overwritten in a derived class.
     *
     * @param hrt     Output vector containing the nondimensional reference state 
     *                enthalpies
     *                Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal *hrt) const {
      err("getEnthalpy_RT_ref");
    }
     
    //!  Returns the vector of nondimensional
    //!  Gibbs Free Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param grt     Output vector containing the nondimensional reference state 
     *                Gibbs Free energies.  Length: m_kk.
     */
    virtual void getGibbs_RT_ref(doublereal *grt) const {
      err("getGibbs_RT_ref");
    }
                   
    //!  Returns the vector of the
    //!  gibbs function of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     *  units = J/kmol
     *
     * @param g       Output vector containing the  reference state 
     *                Gibbs Free energies.  Length: m_kk. Units: J/kmol.
     */
    virtual void getGibbs_ref(doublereal *g) const {
      err("getGibbs_ref");
    }
      
    //!  Returns the vector of nondimensional
    //!  entropies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    /*!
     * @param er      Output vector containing the nondimensional reference state 
     *                entropies.  Length: m_kk.
     */
    virtual void getEntropy_R_ref(doublereal *er) const {
      err("getEntropy_R_ref");
    }
    
    //! Returns the vector of nondimensional
    //!  internal Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    /*!
     * @param urt    Output vector of nondimensional reference state
     *               internal energies of the species.
     *               Length: m_kk
     */
    virtual void getIntEnergy_RT_ref(doublereal *urt) const {
      err("getIntEnergy_RT_ref");
    }
    
    //!  Returns the vector of nondimensional
    //!  constant pressure heat capacities of the reference state
    //!  at the current temperature of the solution
    //!  and reference pressure for each species.
    /*!
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void getCp_R_ref(doublereal *cprt) const {
      err("getCp_R_ref()");
    }
     
    //!  Get the molar volumes of the species reference states at the current
    //!  <I>T</I> and <I>P_ref</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal *vol) const {
      err("getStandardVolumes_ref");
    }

    //! Sets the reference composition
    /*!
     *  @param x   Mole fraction vector to set the reference composition to.
     *             If this is zero, then the reference mole fraction
     *             is set to the current mole fraction vector.
     */
    virtual void setReferenceComposition(const doublereal * const x);

    //! Gets the reference composition
    /*!
     *  The reference mole fraction is a safe mole fraction.
     *  @param x   Mole fraction vector containing the reference composition.
     */
    virtual void getReferenceComposition(doublereal * const x) const;
   
    //
    //  The methods below are not virtual, and should not
    //  be overloaded.
    //
   
    //@}
    /**
     * @name Specific Properties
     * @{
     */

    /**
     * Specific enthalpy. Units: J/kg. 
     */
    doublereal enthalpy_mass() const {
      return enthalpy_mole()/meanMolecularWeight();
    }

    /**
     * Specific internal energy. Units: J/kg. 
     */
    doublereal intEnergy_mass() const {
      return intEnergy_mole()/meanMolecularWeight();
    }

    /**
     * Specific entropy. Units: J/kg/K. 
     */
    doublereal entropy_mass() const {
      return entropy_mole()/meanMolecularWeight();
    }

    /**
     * Specific Gibbs function. Units: J/kg. 
     */
    doublereal gibbs_mass() const {
      return gibbs_mole()/meanMolecularWeight();
    }

    /**
     * Specific heat at constant pressure. Units: J/kg/K. 
     */
    doublereal cp_mass() const {
      return cp_mole()/meanMolecularWeight();
    }

    /**
     * Specific heat at constant volume. Units: J/kg/K. 
     */
    doublereal cv_mass() const {
      return cv_mole()/meanMolecularWeight();
    }
    //@}

    //! Return the Gas Constant multiplied by the current temperature
    /*!
     *  The units are Joules kmol-1
     */
    doublereal _RT() const {
      return temperature() * GasConstant;
    }

    /**
     * @name Setting the State
     *
     * These methods set all or part of the thermodynamic
     * state.
     * @{
     */

    //! Set the temperature (K), pressure (Pa), and mole fractions.
    /*!
     * Note, the mole fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param x    Vector of mole fractions.
     *             Length is equal to m_kk.
     */
    virtual void setState_TPX(doublereal t, doublereal p, const doublereal* x);

    //! Set the temperature (K), pressure (Pa), and mole fractions.  
    /*!
     * Note, the mole fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param x    Composition map of mole fractions. Species not in
     *             the composition map are assumed to have zero mole fraction
     */
    void setState_TPX(doublereal t, doublereal p, compositionMap& x);

    //! Set the temperature (K), pressure (Pa), and mole fractions.  
    /*!
     * Note, the mole fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param x    String containing a composition map of the mole fractions. Species not in
     *             the composition map are assumed to have zero mole fraction
     */
    void setState_TPX(doublereal t, doublereal p, const std::string& x);

    //! Set the internally storred temperature (K), pressure (Pa), and mass fractions of the phase.
    /*!
     * Note, the mass fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param y    Vector of mass fractions.
     *             Length is equal to m_kk.
     */
    void setState_TPY(doublereal t, doublereal p, const doublereal* y);

    //! Set the internally storred temperature (K), pressure (Pa), and mass fractions of the phase
    /*!
     * Note, the mass fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param y    Composition map of mass fractions. Species not in
     *             the composition map are assumed to have zero mass fraction
     */
    void setState_TPY(doublereal t, doublereal p, compositionMap& y);
        
    //! Set the internally storred temperature (K), pressure (Pa), and mass fractions of the phase
    /*!
     * Note, the mass fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param y    String containing a composition map of the mass fractions. Species not in
     *             the composition map are assumed to have zero mass fraction
     */
    void setState_TPY(doublereal t, doublereal p, const std::string& y);

    //! Set the temperature (K) and pressure (Pa)
    /*!
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     */
    void setState_TP(doublereal t, doublereal p);

    //! Set the pressure (Pa) and mole fractions. 
    /*!
     * Note, the mole fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param p    Pressure (Pa)
     * @param x    Vector of mole fractions.
     *             Length is equal to m_kk.
     */
    void setState_PX(doublereal p, doublereal* x);

    //! Set the internally storred pressure (Pa) and mass fractions. 
    /*!
     * Note, the temperature is held constant during this operation.
     * Note, the mass fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param p    Pressure (Pa)
     * @param y    Vector of mass fractions.
     *             Length is equal to m_kk.
     */
    void setState_PY(doublereal p, doublereal* y);

    //! Set the internally storred specific enthalpy (J/kg) and pressure (Pa) of the phase.
    /*!
     * @param h     Specific enthalpy (J/kg)
     * @param p     Pressure (Pa)
     * @param tol   Optional parameter setting the tolerance of the
     *              calculation. Defaults to 1.0E-4
     */
    virtual void setState_HP(doublereal h, doublereal p, doublereal tol = 1.e-4);

    //! Set the specific internal energy (J/kg) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific internal energy and specific volume have the value of the input parameters.
     *
     * @param u    specific internal energy (J/kg)
     * @param v    specific volume (m^3/kg).
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation. Defaults to 1.0E-4
     */
    virtual void setState_UV(doublereal u, doublereal v, doublereal tol = 1.e-4);

    //! Set the internally stored mole fraction, specific enthalpy (J/kg) and pressure (Pa) of the phase.
    /*!
     * @param h     Specific enthalpy (J/kg)
     * @param p     Pressure (Pa)
     * @param x    Vector of mole fractions.
     *             Length is equal to m_kk.
     * @param tol   Optional parameter setting the tolerance of the
     *              calculation. Defaults to 1.0E-4
     */
    virtual void setState_HPX(doublereal h, doublereal p, doublereal* x, doublereal tol = 1.e-4);

    //! Set the specific mole fraction, internal energy (J/kg) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific internal energy and specific volume have the value of the input parameters.
     *
     * @param u    specific internal energy (J/kg)
     * @param v    specific volume (m^3/kg).
     * @param x    Vector of mole fractions.
     *             Length is equal to m_kk.
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation. Defaults to 1.0E-4
     */
    virtual void setState_UVX(doublereal u, doublereal v, doublereal* x, doublereal tol = 1.e-4);

    //! Set the internally stored mass fraction, specific enthalpy (J/kg) and pressure (Pa) of the phase.
    /*!
     * @param h     Specific enthalpy (J/kg)
     * @param p     Pressure (Pa)
     * @param y    Vector of mass fractions.
     *             Length is equal to m_kk.
     * @param tol   Optional parameter setting the tolerance of the
     *              calculation. Defaults to 1.0E-4
     */
    virtual void setState_HPY(doublereal h, doublereal p, doublereal* y, doublereal tol = 1.e-4);

    //! Set the specific mass fraction, internal energy (J/kg) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific internal energy and specific volume have the value of the input parameters.
     *
     * @param u    specific internal energy (J/kg)
     * @param v    specific volume (m^3/kg).
     * @param y    Vector of mass fractions.
     *             Length is equal to m_kk.
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation. Defaults to 1.0E-4
     */
    virtual void setState_UVY(doublereal u, doublereal v, doublereal* y, doublereal tol = 1.e-4);

  private:

    //! Carry out work in HP and UV calculations.
    /*!
     * @param h     Specific enthalpy or internal energy (J/kg)
     * @param p     Pressure (Pa) or specific volume (m^3/kg)
     * @param tol   Optional parameter setting the tolerance of the
     *              calculation. Defaults to 1.0E-4
     * @param doUV  True if solving for UV, false for HP.
     */
    void setState_HPorUV(doublereal h, doublereal p, 
			 doublereal tol = 1.e-4, bool doUV = false);

  public:

    //! Set the specific entropy (J/kg/K) and pressure (Pa).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific entropy and the pressure have the value of the input parameters.
     *
     * @param s    specific entropy (J/kg/K)
     * @param p    specific pressure (Pa).
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation. Defaults to 1.0E-4
     */
    virtual void setState_SP(doublereal s, doublereal p, doublereal tol = 1.e-4);

    //! Set the specific entropy (J/kg/K) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific entropy and specific volume have the value of the input parameters.
     *
     * @param s    specific entropy (J/kg/K)
     * @param v    specific volume (m^3/kg).
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation. Defaults to 1.0E-4
     */
    virtual void setState_SV(doublereal s, doublereal v, doublereal tol = 1.e-4);

  private:

    //! Carry out work in SP and SV calculations.
    /*!
     * @param s     Specific entropy (J/kg)
     * @param p     Pressure (Pa) or specific volume (m^3/kg)
     * @param tol   Optional parameter setting the tolerance of the
     *              calculation. Defaults to 1.0E-4
     * @param doSV  True if solving for SV, false for SP.
     */
    void setState_SPorSV(doublereal s, doublereal p, 
			 doublereal tol = 1.e-4, bool doSV = false);

  public:

    //@}
      
    /**
     * @name Chemical Equilibrium
     * Chemical equilibrium.
     * @{
     */
      
      
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

    //! Stores the element potentials in the ThermoPhase object
    /*!
     * Called by function 'equilibrate' in ChemEquil.h to transfer
     * the element potentials to this object after every successful
     *  equilibration routine.
     * The element potentials are storred in their dimensionless
     * forms, calculated by dividing by RT.
     *
     *    @param lambda Input vector containing the element potentials.
     *           Length = nElements. Units are Joules/kmol.
     */
    void setElementPotentials(const vector_fp& lambda);
  

    //!  Returns the element potentials storred in the ThermoPhase object
    /*!
     * Returns the storred element potentials.
     * The element potentials are retrieved from their storred
     * dimensionless forms by multiplying by RT.
     * @param lambda Output vector containing the element potentials.
     *        Length = nElements. Units are Joules/kmol.
     * @return bool indicating whether thare are any valid storred element
     *         potentials. The calling routine should check this 
     *         bool. In the case that there aren't any, lambda is not
     *         touched.
     */
    bool getElementPotentials(doublereal* lambda) const;

    //@}

        
    //---------------------------------------------------------
    /// @name Critical State Properties.
    /// These methods are only implemented by some subclasses, and may 
    /// be moved out of ThermoPhase at a later date.
        
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

    /** @name Saturation Properties.
     *
     * These methods are only implemented by subclasses that 
     * implement full liquid-vapor equations of state. They may be
     * moved out of %ThermoPhase at a later date.
     */
    //@{

    //! Return the saturation temperature given the pressure
    /*!
     * @param p Pressure (Pa)
     */
    virtual doublereal satTemperature(doublereal p) const {
      err("satTemperature"); return -1.0;
    }

    //! Return the saturation pressure given the temperature
    /*!
     * @param t Temperature (Kelvin)
     */
    virtual doublereal satPressure(doublereal t) const {
      err("satPressure"); return -1.0;
    }

    //! Return the fraction of vapor at the current conditions
    virtual doublereal vaporFraction() const {
      err("vaprFraction"); return -1.0;
    }
        
    //! Set the state to a saturated system at a particular temperature
    /*!
     * @param t  Temperature (kelvin)
     * @param x  Fraction of vapor
     */
    virtual void setState_Tsat(doublereal t, doublereal x) {
      err("setState_sat"); 
    }

    //! Set the state to a saturated system at a particular pressure
    /*!
     * @param p  Pressure (Pa)
     * @param x  Fraction of vapor
     */
    virtual void setState_Psat(doublereal p, doublereal x) {
      err("setState_sat"); 
    }
       
    //@}


    //! @name Initialization Methods - For Internal Use (%ThermoPhase)
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an 
     * input file. They are not normally used in application programs.
     * To see how they are used, 
     * see files importCTML.cpp and  ThermoFactory.cpp.
     */
    //@{

    //! Store a reference pointer to the XML tree containing the species data
    //! for this phase. 
    /*!
     *   This is used to access data needed to construct transport manager later.
     *   @internal
     *
     * @param k      Species index
     * @param data   Pointer to the XML_Node data containing
     *               information about the species in the phase.
     */
    void saveSpeciesData(const int k, const XML_Node* const data);
      
    //!  Return a pointer to the vector of XML nodes containing the species
    //!  data for this phase.
    const std::vector<const XML_Node *> & speciesData() const;
      
    //!  Install a species thermodynamic property manager. 
    /*!
     * The species thermodynamic property manager
     * computes properties of the pure species for use in
     * constructing solution properties. It is meant for internal
     * use, and some classes derived from ThermoPhase may not use
     * any species thermodynamic property manager. This method is
     * called by function importPhase() in importCTML.cpp.
     *
     * @param spthermo input pointer to the species thermodynamic property
     *                 manager.
     *
     *  @internal
     */
    void setSpeciesThermo(SpeciesThermo* spthermo);
        
    //! Return a changeable reference to the calculation manager
    //! for species reference-state thermodynamic properties
    /*!
     *
     * @param k   Speices id. The default is -1, meaning return the default
     *
     * @internal
     */
    virtual SpeciesThermo& speciesThermo(int k = -1);

    /**
     * @internal
     * Initialization of a ThermoPhase object using an
     * ctml file.
     *
     *   This routine is a precursor to initThermoXML(XML_Node*)
     *   routine, which does most of the work.
     *   Here we read extra information about the XML description
     *   of a phase. Regular information about elements and species
     *   and their reference state thermodynamic information
     *   have already been read at this point.
     *   For example, we do not need to call this function for
     *   ideal gas equations of state.
     *
     * @param inputFile XML file containing the description of the
     *        phase
     *
     * @param id  Optional parameter identifying the name of the
     *            phase. If none is given, the first XML
     *            phase element encountered will be used.
     */
    virtual void initThermoFile(std::string inputFile, std::string id);


    //!Import and initialize a ThermoPhase object  using an XML tree.
    /*!
     * @internal
     *
     *   Here we read extra information about the XML description
     *   of a phase. Regular information about elements and species
     *   and their reference state thermodynamic information
     *   have already been read at this point.
     *   For example, we do not need to call this function for
     *   ideal gas equations of state. This function is called from importPhase() 
     *   after the elements and the species are initialized with 
     *   default ideal solution level data.
     *
     *   The default implementation in ThermoPhase calls the
     *   virtual function initThermo() and then sets the "state" of the
     *   phase by looking for an XML element named "state", and then
     *   interpreting its contents by calling the virtual function
     *   setStateFromXML().
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
    
    //! Initialize the ThermoPhase object after all species have been set up
    /*!
     * @internal Initialize.
     *
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called from ThermoPhase::initThermoXML(),
     * which is called from importPhase(),
     * just prior to returning from function importPhase().
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();
    
    //! Add in species from Slave phases
    /*!
     *  This hook is used for  cSS_CONVENTION_SLAVE phases
     *
     *  @param phaseNode   XML Element for the phase
     */
    virtual void installSlavePhases(Cantera::XML_Node* phaseNode);

    // The following methods are used by the clib interface
    // library, and should not be used by application programs.

    /*!
     * @internal 
     * Index number.  This method can be used to identify the
     * location of a phase object in a list, and is used by the
     * interface library (clib) routines for this purpose.
     */
    int index() const { return m_index; }


    /**
     * @internal Set the index number. The Cantera interface
     * library uses this method to set the index number to the
     * location of the pointer to this object in the pointer array
     * it maintains. Using this method for any other purpose will
     * lead to unpredictable results if used in conjunction with
     * the interface library.
     *
     * @param m  Input the index number.
     */ 
    void setIndex(int m) { m_index = m; }
    

    //! Set the equation of state parameters
    /*!
     * @internal
     *  The number and meaning of these depends on the subclass. 
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     */
    virtual void setParameters(int n, doublereal* const c) {}


    //! Get the equation of state parameters in a vector
    /*!
     * @internal
     * The number and meaning of these depends on the subclass. 
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     */
    virtual void getParameters(int &n, doublereal * const c) const {}

      
    //! Set equation of state parameter values from XML entries.
    /*!
     *
     * This method is called by function importPhase() in
     * file importCTML.cpp when processing a phase definition in
     * an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. Note, this method is called before the phase is
     * initialzed with elements and/or species.
     *   
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata) {}
      
    
    //! Set the initial state of the phase to the conditions 
    //! specified in the state XML element.
    /*!
     *
     * This method sets the temperature, pressure, and mole 
     * fraction vector to a set default value.
     *
     * @param state AN XML_Node object corresponding to
     *              the "state" entry for this phase in the
     *              input file.
     */
    virtual void setStateFromXML(const XML_Node& state);

    /**
     * @} 
     * @name  Derivatives of Thermodynamic Variables needed for Applications
     * @{
     */

    //! Get the change in activity coefficients wrt changes in state (temp, mole fraction, etc) along
    //! a line in parameter space or along a line in physical space
    /*!
     *
     * @param dTds           Input of temperature change along the path
     * @param dXds           Input vector of changes in mole fraction along the path. length = m_kk
     *                       Along the path length it must be the case that the mole fractions sum to one.
     * @param dlnActCoeffds  Output vector of the directional derivatives of the 
     *                       log Activity Coefficients along the path. length = m_kk
     *                       units are 1/units(s). if s is a physical coordinate then the units are 1/m.
     */
    virtual void getdlnActCoeffds(const doublereal dTds, const doublereal * const dXds,
				  doublereal *dlnActCoeffds) const {
      err("getdlnActCoeffds");
    }

    //! Get the array of ln mole fraction derivatives of the log activity coefficients - diagonal component only
    /*!
     * This function is a virtual method.  For ideal mixtures 
     * (unity activity coefficients), this can return zero.  
     * Implementations should take the derivative of the 
     * logarithm of the activity coefficient with respect to the 
     * logarithm of the mole fraction variable 
     * that represents the standard state.  
     * This quantity is to be used in conjunction with derivatives of 
     * that mole fraction variable when the derivative of the chemical 
     * potential is taken.  
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnX_diag    Output vector of derivatives of the 
     *                                log Activity Coefficients wrt the mole fractions. length = m_kk
     */
    virtual void getdlnActCoeffdlnX_diag(doublereal *dlnActCoeffdlnX_diag) const {
      err("getdlnActCoeffdlnX_diag");
    }

    //! Get the array of log species mole number derivatives of the log activity coefficients
    /*!
     *  This function is a virtual method. 
     *  For ideal mixtures  (unity activity coefficients), this can return zero.  
     *  Implementations should take the derivative of the 
     *  logarithm of the activity coefficient with respect to the 
     *  logarithm of the concentration-like variable (i.e. moles)
     *  that represents the standard state.  
     *  This quantity is to be used in conjunction with derivatives of 
     *  that species mole number variable when the derivative of the chemical 
     *  potential is taken.  
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnN_diag    Output vector of derivatives of the 
     *                                log Activity Coefficients. length = m_kk
     */
    virtual void getdlnActCoeffdlnN_diag(doublereal *dlnActCoeffdlnN_diag) const {
      err("getdlnActCoeffdlnN_diag");
    }

    //! Get the array of derivatives of the log activity coefficients with respect to the log of the species mole numbers
    /*!
     * Implementations should take the derivative of the logarithm of the activity coefficient with respect to a
     * species log mole number (with all other species mole numbers held constant). The default treatment in the
     * %ThermoPhase object is to set this vector to zero.
     * 
     *  units = 1 / kmol
     *
     *  dlnActCoeffdlnN[ ld * k  + m]  will contain the derivative of log act_coeff for the <I>m</I><SUP>th</SUP> 
     *                               species with respect to the number of moles of the <I>k</I><SUP>th</SUP> species.
     *
     * \f[
     *        \frac{d \ln(\gamma_m) }{d \ln( n_k ) }\Bigg|_{n_i}
     * \f]
     *
     * @param ld               Number of rows in the matrix
     * @param dlnActCoeffdlnN    Output vector of derivatives of the 
     *                           log Activity Coefficients. length = m_kk * m_kk        
     */
    virtual void getdlnActCoeffdlnN(const int ld, doublereal * const dlnActCoeffdlnN);

    virtual void getdlnActCoeffdlnN_numderiv(const int ld, doublereal * const dlnActCoeffdlnN);

    /**
     * @} 
     * @name Printing
     * @{
     */

    //! returns a summary of the state of the phase as a string
    /*!
     * @param show_thermo If true, extra information is printed out
     *                    about the thermodynamic state of the system.
     */
    virtual std::string report(bool show_thermo = true) const;

    //! returns a summary of the state of the phase to a comma separated file
    /*!
     * @param csvFile     ofstream file to print comma separated data for
     *                    the phase
     */
    virtual void reportCSV(std::ofstream& csvFile) const;

    //@}

  protected:

    //! Pointer to the calculation manager for species
    //! reference-state thermodynamic properties
    /*!
     *   This class is called when the reference-state thermodynamic properties
     *   of all the species in the phase needs to be evaluated.
     */
    SpeciesThermo* m_spthermo;

    //! Vector of pointers to the species databases.
    /*!
     * This is used to access data needed to
     * construct the transport manager and other properties
     * later in the initialization process.
     * We create a copy of the XML_Node data read in here. Therefore, we own this
     * data.
     */
    std::vector<const XML_Node *> m_speciesData;

    //! Index number of the phase
    /*!
     * The Cantera interface library uses this member to set the index number to the
     * location of the pointer to this object in the pointer array of ThermoPhase's
     * it maintains. Using this member for any other purpose will
     * lead to unpredictable results if used in conjunction with
     * the interface library.
     */
    int m_index;

    //! Storred value of the electric potential for this phase
    /*!
     * Units are Volts
     */
    doublereal m_phi;

    /// Vector of element potentials.
    ///    -> length equal to number of elements
    vector_fp m_lambdaRRT;

    //! Boolean indicating whether there is a valid set of saved element potentials 
    //! for this phase
    bool m_hasElementPotentials;

    //! Boolean indicating whether a charge neutrality condition is a necessity
    /*!
     * Note, the charge neutrality condition is not a necessity for ideal gas phases. There may
     * be a net charge in those phases, because the NASA polynomials for ionized species 
     * in Ideal gases take this condition into account. 
     * However, liquid phases usually require charge neutrality in order for their derived
     * thermodynamics to be valid.
     */ 
    bool m_chargeNeutralityNecessary;

    //! Contains the standard state convention
    int m_ssConvention;

    //! Reference Mole Fraction Composition
    /*!
     *  Occasionally, the need arises to find a safe mole fraction vector to initialize
     *  the object to. This contains such a vector.
     *  The algorithm will pick up the mole fraction vector that is applied from 
     *  the state xml file in the input file 
     */
    std::vector<doublereal> xMol_Ref;

  private:

    //! Error function that gets called for unhandled cases
    /*!
     * @param msg String containing the message.
     */
    doublereal err(std::string msg) const;

  };

  //! typedef for the ThermoPhase class
  typedef ThermoPhase thermophase_t;

  //! typedef for the ThermoPhase class
  typedef ThermoPhase thermo_t;

  //! Format a summary of the mixture state for output.
  /*!
   * @param th  ThermoPhase object to create a report about
   * @param show_thermo Boolean indicating whether the thermo functions
   *                    of the phase should be written out
   *
   * @return  Returns a string containing the report
   */

  std::string report(const ThermoPhase& th, const bool show_thermo = true);

 
}
        
#endif

