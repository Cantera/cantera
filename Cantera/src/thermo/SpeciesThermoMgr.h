/**
 *  @file SpeciesThermoMgr.h
 *  This file contains descriptions of templated subclasses of 
 *  the virtual base class, SpeciesThermo, which 
 *  include  SpeciesThermoDuo and SpeciesThermo1
 *  (see \ref mgrsrefcalc and classes 
 *   \link Cantera::SpeciesThermoDuo SpeciesThermoDuo\endlink and
 *  \link Cantera::SpeciesThermo1 SpeciesThermo1\endlink)
 *
 * $Author: hkmoffa $
 * $Revision: 1.12 $
 * $Date: 2009/02/11 20:03:08 $
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_SPECIESTHERMO_MGR_H
#define CT_SPECIESTHERMO_MGR_H

#include "ct_defs.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "SpeciesThermo.h"
#include <map>

namespace Cantera {

  //! Invokes the 'updateProperties' method of all objects in the list.
  /*!
   *  This templated function has one template, InputIter. It should 
   *  point to a class such as one that inherits from the virtual
   *  base class, SpeciesThermoInterpType, which has 
   *  an updateProperties(T, Cp_R, h_RT, s)R) function
   *
   * @param begin Beginning iterator
   * @param end   end iterator
   * @param T       Temperature (Kelvin)
   * @param cp_R    Vector of Dimensionless heat capacities.
   *                (length m_kk).
   * @param h_RT    Vector of Dimensionless enthalpies.
   *                (length m_kk).
   * @param s_R     Vector of Dimensionless entropies.
   *                (length m_kk).
   *
   * @ingroup mgrsrefcalc
   */
  template<class InputIter>
  inline void _updateAll(InputIter begin, 
			 InputIter end,
			 doublereal T,  
			 vector_fp&  cp_R, 
			 vector_fp&  h_RT, 
			 vector_fp&  s_R)
  {
    for (; begin != end; ++begin)
      begin->updateProperties(T, cp_R, h_RT, s_R);
  }

  //! Iterates through a list of objects which implement a method
  //! 'minTemp()', and returns the largest 'minTemp' value.
  /*!
   *  This templated function has one template, InputIter. It should 
   *  point to a class such as one that inherits from either
   *  SpeciesThermoInterpType or SpeciesThermo, which have a minTemp() function
   *
   * @param begin Beginning iterator
   * @param end   end iterator
   *
   * @ingroup mgrsrefcalc
   */
  template<class InputIter>
  doublereal _minTemp(InputIter begin, InputIter end) {
    doublereal _minT = 0.0;
    for (; begin != end; ++begin)
      _minT = fmaxx(_minT, begin->minTemp());
    return _minT;
  }

  //! Iterates through a list of objects which implement a method
  //! 'maxTemp()', and returns the smallest 'maxTemp' value.
  /*!
   *  This templated function has one template, InputIter. It should 
   *  point to a class such as one that inherits from either
   *  SpeciesThermoInterpType or SpeciesThermo which have a minTemp() function
   *
   * @param begin Beginning iterator
   * @param end   end iterator
   *
   * @ingroup mgrsrefcalc
   */
  template<class _InputIter>
  doublereal _maxTemp(_InputIter begin, _InputIter end) {
    doublereal _maxT = 1.e10;
    for (; begin != end; ++begin)
      _maxT = fminn(_maxT, begin->maxTemp());
    return _maxT;
  }

  ///////////////////////  Exceptions //////////////////////////////

  //! Exception thrown if species reference pressures don't match.
  /*!
   * @ingroup mgrsrefcalc
   */
  class RefPressureMismatch : public CanteraError {
  public:
    //! constructor
    /*!
     * @param proc   name of the procecdure
     * @param prnew   reference pressure
     * @param prold   old reference pressure
     */
    RefPressureMismatch(std::string proc, doublereal prnew, 
			doublereal prold) : CanteraError(proc, 
							 "Species reference pressure ("
							 + fp2str(prnew) + ") does not match previously-defined "
							 + "reference pressure (" + fp2str(prold) + ")") {}
    //! destructor
    virtual ~RefPressureMismatch() throw() {}
  };

  //! Unknown species thermo manager string error
  /*!
   * @ingroup mgrsrefcalc
   */
  class UnknownSpeciesThermo : public CanteraError {
  public:

    //! constructor
    /*!
     * @param proc   name of the procecdure
     * @param type   unknown type
     */
    UnknownSpeciesThermo(std::string proc, int type) :
      CanteraError(proc, "Specified species parameterization type (" + int2str(type) 
		   + ") does not match any known type.") {}

    //! Alternate constructor
    /*!
     * @param proc   name of the procecdure
     * @param stype   String name for the unknown type
     */
    UnknownSpeciesThermo(std::string proc, std::string stype) :
      CanteraError(proc, "Specified species parameterization type (" + stype 
		   + ") does not match any known type.") {}
    //! destructor
    virtual ~UnknownSpeciesThermo() throw() {}
  };


  /**
   *  This species thermo manager requires that all species have one
   *  of two parameterizations.
   *
   * Note this seems to be a slow way to do things, and it may be on its way out.
   *
   * @ingroup mgrsrefcalc
   */
  template<class T1, class T2>
  class SpeciesThermoDuo : public SpeciesThermo {
	
  public:
    //! Constructor
    SpeciesThermoDuo();

    //! Destructor
    virtual ~SpeciesThermoDuo();

    //! copy constructor
    /*!
     * @param right Object to be copied
     */
    SpeciesThermoDuo(const SpeciesThermoDuo &right);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    SpeciesThermoDuo& operator=(const SpeciesThermoDuo &right);

    //! Duplication routine for objects which inherit from 
    //! %SpeciesThermo
    /*!
     *  This virtual routine can be used to duplicate %SpeciesThermo  objects
     *  inherited from %SpeciesThermo even if the application only has
     *  a pointer to %SpeciesThermo to work with.
     *  ->commented out because we first need to add copy constructors
     *   and assignment operators to all of the derived classes.
     */
    virtual SpeciesThermo *duplMyselfAsSpeciesThermo() const;

    /**
     * install a new species thermodynamic property
     * parameterization for one species.  
     *
     * @param name      Name of the species
     * @param sp        The 'update' method will update the property 
     *                  values for this species 
     *                  at position i index in the property arrays.  
     * @param type      int flag specifying the type of parameterization to be
     *                 installed. 
     * @param c        vector of coefficients for the parameterization. 
     *                 This vector is simply passed through to the
     *                 parameterization constructor.
     * @param minTemp  minimum temperature for which this parameterization
     *                 is valid.
     * @param maxTemp  maximum temperature for which this parameterization
     *                 is valid.
     * @param refPressure standard-state pressure for this 
     *                    parameterization. 
     * @see speciesThermoTypes.h 
     */
    virtual void install(std::string name, int sp, int type,
			 const doublereal* c,
			 doublereal minTemp, doublereal maxTemp, 
			 doublereal refPressure);

    //! Install a new species thermodynamic property
    //! parameterization for one species.
    /*!
     * @param stit_ptr Pointer to the SpeciesThermoInterpType object
     *          This will set up the thermo for one species
     */
    virtual void install_STIT(SpeciesThermoInterpType *stit_ptr);

    //! Compute the reference-state properties for all species.
    /*!
     * Given temperature T in K, this method updates the values of
     * the non-dimensional heat capacity at constant pressure,
     * enthalpy, and entropy, at the reference pressure, Pref
     * of each of the standard states.
     *
     * @param t       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    virtual void update(doublereal t, doublereal* cp_R,
			doublereal* h_RT, doublereal* s_R) const;

    //! Minimum temperature.
    /*!
     * If no argument is supplied, this
     * method returns the minimum temperature for which \e all
     * parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the minimum
     * temperature for species k in the phase.
     *
     * @param k    Species index
     */ 
    virtual doublereal minTemp(int k = -1) const;

    //! Maximum temperature.
    /*!
     * If no argument is supplied, this
     * method returns the maximum temperature for which \e all
     * parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the maximum
     * temperature for parameterization k.
     *
     * @param k index for parameterization k
     */
    virtual doublereal maxTemp(int k = -1) const;
         
    /**
     * The reference-state pressure for species k.
     *
     * returns the reference state pressure in Pascals for
     * species k. If k is left out of the argument list,
     * it returns the reference state pressure for the first
     * species.
     * Note that some SpeciesThermo implementations, such
     * as those for ideal gases, require that all species
     * in the same phase have the same reference state pressures.
     *
     * @param k index for parameterization k
     */
    virtual doublereal refPressure(int k = -1) const;

    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     *
     * @param k  Species index
     */
    virtual int reportType(int k) const;

    /*!
     * This utility function reports back the type of 
     * parameterization and all of the parameters for the 
     * species, index.
     *
     * @param index     Species index
     * @param type      Integer type of the standard type
     * @param c         Vector of coefficients used to set the
     *                  parameters for the standard state.
     * @param minTemp   output - Minimum temperature
     * @param maxTemp   output - Maximum temperature
     * @param refPressure output - reference pressure (Pa).
     *
     */
    virtual void reportParams(int index, int &type, 
			      doublereal * const c, 
			      doublereal &minTemp, 
			      doublereal &maxTemp, 
			      doublereal &refPressure) const;

    //! Modify parameters for the standard state
    /*!
     * @param index Species index
     * @param c     Vector of coefficients used to set the
     *              parameters for the standard state.
     */
    virtual void modifyParams(int index, doublereal *c);


#ifdef H298MODIFY_CAPABILITY
  
    virtual doublereal reportOneHf298(int k) const {
      throw CanteraError("reportHF298", "unimplemented");
    }

    virtual void modifyOneHf298(const int k, const doublereal Hf298New) {
      throw CanteraError("reportHF298", "unimplemented");
    }
  
#endif

  private:

    //! Thermo Type 1
    T1 m_thermo1;
    //! Thermo Type 2
    T2 m_thermo2;
    //! Reference pressure
    doublereal m_p0;
    //! map from species to type
    std::map<int, int> speciesToType;
  };
  
  //!  This species thermo manager requires that all species have the
  //!  same parameterization.
  /*!
   *
   *  This is a templated class. The first template is called SPM. SPM
   *  is an object that calculates the thermo for one species. This
   *  class contains a vector of SPM's, one for each
   *  species. Together, the vector of SPM's is itself a SpeciesThermo
   *  class.
   *
   * @todo The form of the template class, SPM, is basically
   *       unspecified. it needs to be nailed down to a specific
   *       form. One way to do this is with a virtual base class
   *       formulation. Note, that the specification could be that it
   *       inherits from the class SpeciesThermo, itself.
   *
   * @deprecated Note this is currently unused and it may be on its way out.
   *
   * @ingroup mgrsrefcalc
   */
  template<class SPM>
  class SpeciesThermo1 : public SpeciesThermo {
	
  public:
    //! base constructor
    SpeciesThermo1();
    //! destructor
    virtual ~SpeciesThermo1();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    SpeciesThermo1(const SpeciesThermo1 &right);

    //! Asignment Operator
    /*!
     * @param right Object to be copied
     */
    SpeciesThermo1 & operator=(const SpeciesThermo1 &right);

    //! Duplication routine for objects which inherit from 
    //! %SpeciesThermo
    /*!
     *  This virtual routine can be used to duplicate %SpeciesThermo  objects
     *  inherited from %SpeciesThermo even if the application only has
     *  a pointer to %SpeciesThermo to work with.
     *  ->commented out because we first need to add copy constructors
     *   and assignment operators to all of the derived classes.
     */
    virtual SpeciesThermo *duplMyselfAsSpeciesThermo() const;
    
    //! Install one species into this Species Thermo Manager
    /*!
     * @param name       Name of the species
     * @param sp         Species index
     * @param type       species type in terms of an int
     * @param c          Parameters for the species thermo
     */
    virtual void install(std::string name, int sp, int type,
			 const vector_fp& c);

    //! update the object, because the temperature changed
    /*!
     * @param t          temperature(Kelvin)
     * @param cp_R       vector of dimensionless heat capacity
     * @param h_RT       vector of dimensionless enthalpy
     * @param s_R        vector of dimensionless entropy
     */
    virtual void update(doublereal t, vector_fp& cp_R, 
			vector_fp& h_RT, vector_fp& s_R) const;

    //! update the object for one species, because the temperature changed
    /*!
     * @param k          species index
     * @param t          temperature(Kelvin)
     * @param cp_R       vector of dimensionless heat capacity
     * @param h_RT       vector of dimensionless enthalpy
     * @param s_R        vector of dimensionless entropy
     */
    virtual void update_one(int k, doublereal t, vector_fp& cp_R, 
			    vector_fp& h_RT, vector_fp& s_R) const;

    //! Returns the minimum temperature
    /*!
     * @param k species index. Defaults to -1.
     */
    virtual doublereal minTemp(int k = -1) const;

    //! Returns the maximum temperature
    /*!
     * @param k species index. Defaults to -1.
     */
    virtual doublereal maxTemp(int k = -1) const;
                    
    //! returns the reference pressure
    /*!
     * @param k species index. Defaults to -1.
     */
    virtual doublereal refPressure(int k = -1) const;

    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     *  Note, all parameterizations are the same, by definition, here
     *
     * @param k  Species index
     */
    virtual int reportType(int k) const;

    /*!
     * This utility function reports back the type of 
     * parameterization and all of the parameters for the 
     * species, index.
     *
     * @param index     Species index
     * @param type      Integer type of the standard type
     * @param c         Vector of coefficients used to set the
     *                  parameters for the standard state.
     * @param minTemp   output - Minimum temperature
     * @param maxTemp   output - Maximum temperature
     * @param refPressure output - reference pressure (Pa).
     */
    virtual void reportParams(int index, int &type, 
			      doublereal * const c, 
			      doublereal &minTemp, 
			      doublereal &maxTemp, 
			      doublereal &refPressure) const;

    //! Modify parameters for the standard state
    /*!
     * @param index Species index
     * @param c     Vector of coefficients used to set the
     *              parameters for the standard state.
     */
    virtual void modifyParams(int index, doublereal *c);

#ifdef H298MODIFY_CAPABILITY

    virtual doublereal reportOneHf298(int k) const {
      throw CanteraError("reportHF298", "unimplemented");
    }

    virtual void modifyOneHf298(const int k, const doublereal Hf298New) {
      throw CanteraError("reportHF298", "unimplemented");
    }

#endif
  private:
    //! Vector of SPM objects. There are m_kk of them
    std::vector<SPM> m_thermo;

    //! Reference pressure (Pa)
    doublereal m_pref;
  };
 
  // ------------------------- cpp part of file -------------------------------------

  // Definitions for the SpeciesThermoDuo<T1,T2> templated class

  template<class T1, class T2>
  SpeciesThermoDuo<T1, T2>::SpeciesThermoDuo() 
  {
  }
 
  template<class T1, class T2>
  SpeciesThermoDuo<T1, T2>::~SpeciesThermoDuo()
  {
  }

  template<class T1, class T2>
  SpeciesThermoDuo<T1, T2>::SpeciesThermoDuo(const SpeciesThermoDuo &right) {
    *this = operator=(right);
  }

  template<class T1, class T2>
  SpeciesThermoDuo<T1, T2> & 
  SpeciesThermoDuo<T1, T2>::operator=(const SpeciesThermoDuo &right) {
    if (&right == this) return *this;
    
    m_thermo1 = right.m_thermo1;
    m_thermo2 = right.m_thermo2;
    m_p0      = right.m_p0;
    speciesToType = right.speciesToType;
    
    return *this;
  }

  template<class T1, class T2>
  SpeciesThermo *
  SpeciesThermoDuo<T1, T2>::duplMyselfAsSpeciesThermo() const {
    SpeciesThermoDuo<T1,T2> *nt = new SpeciesThermoDuo<T1,T2>(*this);
    return (SpeciesThermo *) nt;
  }

  template<class T1, class T2>
  void 
  SpeciesThermoDuo<T1, T2>::install(std::string name, int sp, int type,
				    const doublereal* c,
				    doublereal minTemp,
				    doublereal maxTemp, 
				    doublereal refPressure) {
    m_p0 = refPressure;
    if (type == m_thermo1.ID) {
      m_thermo1.install(name, sp, 0, c, minTemp, maxTemp,
			refPressure);
      speciesToType[sp] = m_thermo1.ID;
    } else if (type == m_thermo2.ID) {
      m_thermo2.install(name, sp, 0, c, minTemp, maxTemp,
			refPressure);
      speciesToType[sp] = m_thermo2.ID;
    } else {
      throw UnknownSpeciesThermo("SpeciesThermoDuo:install",type);
    }
  }

  template<class T1, class T2>
  void 
  SpeciesThermoDuo<T1, T2>::install_STIT(SpeciesThermoInterpType *stit_ptr) {
    throw CanteraError("install_STIT", "not implemented");
  }

  template<class T1, class T2>
  void 
  SpeciesThermoDuo<T1, T2>::update(doublereal t, doublereal* cp_R,
				   doublereal* h_RT, doublereal* s_R) const {
    m_thermo1.update(t, cp_R, h_RT, s_R);
    m_thermo2.update(t, cp_R, h_RT, s_R);
  }

  template<class T1, class T2>
  doublereal 
  SpeciesThermoDuo<T1, T2>::minTemp(int k) const {
    doublereal tm1 = m_thermo1.minTemp();
    doublereal tm2 = m_thermo2.minTemp();
    return (tm1 < tm2 ? tm2 : tm1);
  }
  
  template<class T1, class T2>
  doublereal
  SpeciesThermoDuo<T1, T2>::maxTemp(int k) const {
    doublereal tm1 = m_thermo1.maxTemp();
    doublereal tm2 = m_thermo2.maxTemp();
    return (tm1 < tm2 ? tm1 : tm2);
  }
  
  template<class T1, class T2>
  doublereal
  SpeciesThermoDuo<T1, T2>::refPressure(int k) const {
    return m_p0;
  }
  
  template<class T1, class T2>
  int 
  SpeciesThermoDuo<T1, T2>::reportType(int k) const {
    std::map<int, int>::const_iterator p = speciesToType.find(k);
    if (p != speciesToType.end()) {
      const int type = p->second;
      return type;
    } 
    return -1;
  }
 
  template<class T1, class T2>
  void 
  SpeciesThermoDuo<T1, T2>::reportParams(int index, int &type, 
					 doublereal * const c, 
					 doublereal &minTemp, 
					 doublereal &maxTemp, 
					 doublereal &refPressure) const {
    int ctype = reportType(index);
    if (ctype == m_thermo1.ID) {
      m_thermo1.reportParams(index, type, c, minTemp, maxTemp, 
			     refPressure);
    } else if (ctype == m_thermo2.ID) {
      m_thermo2.reportParams(index, type, c, minTemp, maxTemp, 
			     refPressure);
    } else {
      throw CanteraError("  ", "confused");
    }
  }

  template<class T1, class T2>
  void 
  SpeciesThermoDuo<T1, T2>::modifyParams(int index, doublereal *c) {
    int ctype = reportType(index);
    if (ctype == m_thermo1.ID) {
      m_thermo1.modifyParams(index, c);
    } else if (ctype == m_thermo2.ID) {
      m_thermo2.modifyParams(index,  c);
    } else {
      throw CanteraError("modifyParams", "confused");
    }
  }

  // Definitions for the SpeciesTherm1 templated class

  template<class SPM>
  SpeciesThermo1<SPM>::SpeciesThermo1() : 
    m_pref(0.0) 
  {
  }

  template<class SPM>
  SpeciesThermo1<SPM>::~SpeciesThermo1()
  {
  }

  template<class SPM>
  SpeciesThermo1<SPM>::SpeciesThermo1(const SpeciesThermo1 &right) :
    m_pref(0.0)
  {
    *this = operator=(right);
  }

  template<class SPM>
  SpeciesThermo1<SPM> & 
  SpeciesThermo1<SPM>::operator=(const SpeciesThermo1 &right)
  {
    if (&right == this) return *this;  
    m_thermo = right.m_thermo;
    m_pref   = right.m_pref;
    return *this;
  }

  template<class SPM>
  SpeciesThermo *
  SpeciesThermo1<SPM>::duplMyselfAsSpeciesThermo() const {
    SpeciesThermo1<SPM> *nt = new SpeciesThermo1<SPM>(*this);
    return (SpeciesThermo *) nt;
  }
  
  template<class SPM>
  void 
  SpeciesThermo1<SPM>::install(std::string name, int sp, int type, const vector_fp& c)
  {
    m_thermo.push_back(SPM(sp, c));
    if (m_pref) {
      if (m_thermo.begin()->refPressure() != m_pref) {
	throw RefPressureMismatch("SpeciesThermo1:install",
				  refPressure(), m_pref);
      }
    }
    else  m_pref = m_thermo.begin()->refPressure();
  }

  template<class SPM>
  inline void
  SpeciesThermo1<SPM>::update(doublereal t, vector_fp& cp_R, 
			      vector_fp& h_RT, vector_fp& s_R) const {
    _updateAll(m_thermo.begin(), m_thermo.end(), t, cp_R, h_RT, s_R);
  }
  
  template<class SPM>
  void 
  SpeciesThermo1<SPM>::update_one(int k, doublereal t, vector_fp& cp_R, 
				  vector_fp& h_RT, vector_fp& s_R) const {
    m_thermo[k]->update(t, cp_R, h_RT, s_R);
  }
  
  template<class SPM>
  doublereal 
  SpeciesThermo1<SPM>::minTemp(int k) const {
    if (k < 0)
      return _minTemp(m_thermo.begin(), m_thermo.end());
    else
      return m_thermo[k].minTemp();
  }

  template<class SPM>
  doublereal 
  SpeciesThermo1<SPM>::maxTemp(int k) const {
    if (k < 0)
      return _maxTemp(m_thermo.begin(), m_thermo.end());
    else
      return m_thermo[k].maxTemp();
  }

  template<class SPM>
  doublereal 
  SpeciesThermo1<SPM>::refPressure(int k) const {
    return m_pref;
  }

  template<class SPM>
  int 
  SpeciesThermo1<SPM>::reportType(int k) const {
    return m_thermo[k]->reportType(-1);	  
  }

  template<class SPM>
  void 
  SpeciesThermo1<SPM>::reportParams(int index, int &type, 
				    doublereal * const c, 
				    doublereal &minTemp, 
				    doublereal &maxTemp, 
				    doublereal &refPressure) const {
    m_thermo[index]->reportParameters(index, type, c, minTemp, maxTemp, refPressure);
  }

  template<class SPM>
  void
  SpeciesThermo1<SPM>::modifyParams(int index, doublereal *c) {
    m_thermo[index]->modifyParameters(index, c);
  }
}

#endif
