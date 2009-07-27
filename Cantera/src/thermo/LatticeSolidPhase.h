/**
 *  @file LatticeSolidPhase.h
 *  Header for a simple thermodynamics model of a bulk solid phase
 *  derived from ThermoPhase,
 *  assuming an ideal solution model based on a lattice of solid atoms
 *  (see \ref thermoprops and class \link Cantera::LatticeSolidPhase LatticeSolidPhase\endlink).

 */

/*  $Author: hkmoffa $
 *  $Date: 2008/12/13 01:59:49 $
 *  $Revision: 1.5 $
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
#include "utilities.h"



namespace Cantera {

  class LatticePhase;

  //! Additive combination of lattice phases
  /*!
   *
   */
  class LatticeSolidPhase : public ThermoPhase  {

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

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * For an ideal solution,
     * \f[
     * \hat h(T,P) = \sum_k X_k \hat h^0_k(T),
     * \f]
     * and is a function only of temperature.
     * The standard-state pure-species Enthalpies
     * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic
     * property manager.
     *
     * \see SpeciesThermo
     */
    virtual doublereal enthalpy_mole() const;

    virtual doublereal intEnergy_mole() const;

    virtual doublereal entropy_mole() const;

    virtual doublereal gibbs_mole() const;

    virtual doublereal cp_mole() const;

    virtual doublereal cv_mole() const {
      return cp_mole();
    }

    virtual doublereal pressure() const {
      return m_press;
    }

    virtual void setPressure(doublereal p) {
      m_press = p;
      setMolarDensity(m_molar_density);
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
  
    virtual void modifyOneHf298SS(const int k, const doublereal Hf298New) {
       m_spthermo->modifyOneHf298(k, Hf298New);
       m_tlast += 0.0001234;
    }
#endif

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
        
#endif
#endif
