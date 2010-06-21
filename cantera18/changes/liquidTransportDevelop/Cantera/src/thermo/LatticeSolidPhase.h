/**
 *  @file LatticeSolidPhase.h
 *  Header for a simple thermodynamics model of a bulk solid phase
 *  derived from ThermoPhase,
 *  assuming an ideal solution model based on a lattice of solid atoms
 *  (see \ref thermoprops and class \link Cantera::LatticeSolidPhase LatticeSolidPhase\endlink).

 */

/*  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2005 California Institute of Technology
 *
 */

#ifndef CT_LATTICESOLID_H
#define CT_LATTICESOLID_H

#include "config.h"

#ifdef WITH_LATTICE_SOLID

#include "ct_defs.h"

#include "mix_defs.h"
#include "ThermoPhase.h"
#include "SpeciesThermo.h"
#include "LatticePhase.h"
#include "utilities.h"



namespace Cantera {

  //! A phase that is comprised of an additive combination of other lattice phases
  /*!
   *  This is the main way Cantera describes semiconductors and other solid phases.
   *  This Thermophase object calculates its properties as a sum over other LatticePhase objects. Each of the %LatticePhase
   *  objects is a ThermoPhase object by itself.
   *
   *  The sum over the LatticePhase objects is carried out by weighting each LatticePhase object
   *  value with the molarDensity of the LatticePhase. Then the resulting quantity is divided by
   *  the molar density of the total compound. The LatticeSolidPhase object therefore only contains a 
   *  listing of the number of Lattice Phases
   *  that comprises the solid and it contains a value for the molar density of the entire mixture.
   *
   *  Let's take FeS2 as an example, which may be thought of as a combination of two lattices: Fe and S lattice.
   *  The Fe sublattice has a molar density of 1 gmol cm-3. The S sublattice has a molar density of 2 gmol cm-3.
   *  We then define the LatticeSolidPhase object as having a nominal composition of FeS2, and having a
   *  molar density of 1 gmol cm-3.  All quantities pertaining to the FeS2 compound will be have weights
   *  associated with the sublattices. The Fe sublattice will have a weight of 1.0 associated with it. The
   *  S sublattice will have a weight of 2.0 associated with it. 
   *
   *  Currently, the molar density is set to a constant.
   *
   *  The results from this LatticeSolidPhase model reduces to the LatticePhase model when there is one
   *  lattice phase and the molar densities of the sublattice and the molar density within the LatticeSolidPhase
   *  have the same values.
   *
   *  The mole fraction vector has been redefined within the LatticeSolidPhase object. The mole fractions sum
   *  to one within each of the individual lattice phases. The routine getMoleFraction() and setMoleFraction()
   *  have been redefined to use this convention.
   *      
   */
  class LatticeSolidPhase : public ThermoPhase {

  public:

    //! Base empty constructor
    LatticeSolidPhase();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    LatticeSolidPhase(const LatticeSolidPhase &right);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    LatticeSolidPhase& operator=(const LatticeSolidPhase& right);

    //! Destructor
    virtual ~LatticeSolidPhase();

    //! Duplication function
    /*!
     * This virtual function is used to create a duplicate of the
     * current phase. It's used to duplicate the phase when given
     * a ThermoPhase pointer to the phase.
     *
     * @return It returns a ThermoPhase pointer.
     */
    ThermoPhase *duplMyselfAsThermoPhase() const;

    //! Equation of state type flag.
    /*!
     *  Redefine this to return cLatticeSolid, listed in mix_defs.h.
     */
    virtual int eosType() const { return cLatticeSolid; }

    //! This method returns the convention used in specification
    //! of the standard state, of which there are currently two,
    //! temperature based, and variable pressure based.
    /*!
     *  All of the thermo is determined by slave %ThermoPhase routines.
     */
    virtual int standardStateConvention() const {
      return cSS_CONVENTION_SLAVE; 
    }

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * The molar enthalpy is determined by the following formula, where \f$ C_n \f$ is the
     * lattice molar density of the nth lattice, and \f$ C_T \f$ is the molar density
     * of the solid compound.
     *
     * \f[
     *   \tilde h(T,P) = \frac{\sum_n C_n \tilde h_n(T,P) }{C_T},
     * \f]
     *
     * \f$ \tilde h_n(T,P) \f$ is the enthalpy of the n<SUP>th</SUP> lattice. 
     *
     *  units J/kmol
     */
    virtual doublereal enthalpy_mole() const;

  
    //! Return the Molar Internal Energy. Units: J/kmol.
    /*!
     * The molar internal energy is determined by the following formula, where \f$ C_n \f$ is the
     * lattice molar density of the nth lattice, and \f$ C_T \f$ is the molar density
     * of the solid compound.
     *
     * \f[
     *   \tilde u(T,P) = \frac{\sum_n C_n \tilde u_n(T,P) }{C_T},
     * \f]
     *
     * \f$ \tilde u_n(T,P) \f$ is the internal energy of the n<SUP>th</SUP> lattice. 
     *
     *  units J/kmol
     */
    virtual doublereal intEnergy_mole() const;

    //! Return the Molar Entropy. Units: J/kmol/K.
    /*!
     * The molar entropy is determined by the following formula, where \f$ C_n \f$ is the
     * lattice molar density of the nth lattice, and \f$ C_T \f$ is the molar density
     * of the solid compound.
     *
     * \f[
     *   \tilde s(T,P) = \frac{\sum_n C_n \tilde s_n(T,P) }{C_T},
     * \f]
     *
     * \f$ \tilde s_n(T,P) \f$ is the molar entropy of the n<SUP>th</SUP> lattice. 
     *
     *  units J/kmol/K
     */
    virtual doublereal entropy_mole() const;

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * The molar enthalpy is determined by the following formula, where \f$ C_n \f$ is the
     * lattice molar density of the nth lattice, and \f$ C_T \f$ is the molar density
     * of the solid compound.
     *
     * \f[
     *   \tilde h(T,P) = \frac{\sum_n C_n \tilde h_n(T,P) }{C_T},
     * \f]
     *
     * \f$ \tilde h_n(T,P) \f$ is the enthalpy of the n<SUP>th</SUP> lattice. 
     *
     *  units J/kmol
     */
    virtual doublereal gibbs_mole() const;

    //! Return the constant pressure heat capacity. Units: J/kmol/K
    /*!
     * The molar constant pressure heat capacity is determined by the following formula, where \f$ C_n \f$ is the
     * lattice molar density of the nth lattice, and \f$ C_T \f$ is the molar density
     * of the solid compound.
     *
     * \f[
     *   \tilde c_{p,n}(T,P) = \frac{\sum_n C_n \tilde c_{p,n}(T,P) }{C_T},
     * \f]
     *
     * \f$ \tilde c_{p,n}(T,P) \f$ is the heat capacity of the n<SUP>th</SUP> lattice. 
     *
     *  units J/kmol/K
     */
    virtual doublereal cp_mole() const;

    //! Return the constant volume heat capacity. Units: J/kmol/K
    /*!
     * The molar constant volume heat capacity is determined by the following formula, where \f$ C_n \f$ is the
     * lattice molar density of the nth lattice, and \f$ C_T \f$ is the molar density
     * of the solid compound.
     *
     * \f[
     *   \tilde c_{v,n}(T,P) = \frac{\sum_n C_n \tilde c_{v,n}(T,P) }{C_T},
     * \f]
     *
     * \f$ \tilde c_{v,n}(T,P) \f$ is the heat capacity of the n<SUP>th</SUP> lattice. 
     *
     *  units J/kmol/K
     */
    virtual doublereal cv_mole() const {
      return cp_mole();
    }

    //! Report the Pressure. Units: Pa.
    /*!
     *  This method simply returns the storred pressure value.
     */
    virtual doublereal pressure() const {
      return m_press;
    }

    //! Set the pressure at constant temperature. Units: Pa.
    /*!
     *
     * @param p Pressure (units - Pa)
     */
    virtual void setPressure(doublereal p) {
      m_press = p;
      setMolarDensity(m_molar_density);
    }
    virtual void setMoleFractions(const doublereal *x);
    virtual void getMoleFractions(doublereal *x) const;

    doublereal moleFraction(const int k) const {
      return err("not implemented");
    }
    

    void getMassFractions(doublereal* const y) const {
      err("not implemented");
    }
    doublereal massFraction(const int k) const {
      return err("not implemented");
    }

    virtual void setMassFractions(const doublereal *y) {
      err("not implemented");
    }
    virtual void setMassFractions_NoNorm(const doublereal* const y) {
      err("not implemented");
    }

    void getConcentrations(doublereal* const c) const {
      err("not implemented");
    }

    doublereal concentration(int k) const {
      return err("not implemented");
    }

    virtual void setConcentrations(const doublereal* const conc) {
      err("not implemented");
    }


    virtual void getActivityConcentrations(doublereal* c) const;

    virtual void getActivityCoefficients(doublereal* ac) const;

    virtual void getChemPotentials(doublereal* mu) const;
    virtual void getStandardChemPotentials(doublereal* mu0) const;
    virtual doublereal standardConcentration(int k=0) const;
    virtual doublereal logStandardConc(int k=0) const;

    virtual void initThermo();

    virtual void setParametersFromXML(const XML_Node& eosdata);

    void setLatticeMoleFractions(int n, std::string x);


#ifdef H298MODIFY_CAPABILITY
  
    //! Modify the value of the 298 K Heat of Formation of one species in the phase (J kmol-1)
    /*!
     *   The 298K heat of formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param  k           Species k
     *   @param  HF298New    Specify the new value of the Heat of Formation at 298K and 1 bar                      
     */
    virtual void modifyOneHf298SS(const int k, const doublereal Hf298New) {
       m_spthermo->modifyOneHf298(k, Hf298New);
       m_tlast += 0.0001234;
    }
#endif

  private:
    //! error routine
    /*!
     *  @param msg Message
     *
     *  @return nothing
     */
    doublereal err(std::string msg) const;

  protected:
        
    int m_mm;
    int m_kk;
    mutable doublereal        m_tlast;
    doublereal                m_press;
    doublereal                m_molar_density;


    int                       m_nlattice;
    std::vector<LatticePhase*>     m_lattice;
    mutable vector_fp                 m_x;

  private:

    void _updateThermo() const;
  };
}
        
#endif // #ifdef WITH_LATTICE_SOLID

#endif
