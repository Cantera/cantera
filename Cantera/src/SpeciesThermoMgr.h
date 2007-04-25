/**
 *  @file SpeciesThermoMgr.h
 *  This file contains descriptions of templated subclasses of 
 *  the virtual base class, SpeciesThermo, which 
 *  include  SpeciesThermoDuo and SpeciesThermo1
 *  (see \ref spthermo and classes 
 *   \link Cantera::SpeciesThermoDuo SpeciesThermoDuo\endlink and
 *  \link Cantera::SpeciesThermo1 SpeciesThermo1\endlink)
 *
 * $Author$
 * $Revision$
 * $Date$
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
   * @ingroup spthermo
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
   * @ingroup spthermo
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
   * @ingroup spthermo
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
   * @ingroup spthermo
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
    virtual ~RefPressureMismatch() {}
  };

  //! Unknown species thermo manager string error
  /*!
   * @ingroup spthermo
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
    //! destructor
    virtual ~UnknownSpeciesThermo() {}
  };



  /**
   *  This species thermo manager requires that all species have one
   *  of two parameterizations.
   *
   * Note this seems to be a slow way to do things, and it may be on its way out.
   *
   * @ingroup spthermo
   */
  template<class T1, class T2>
  class SpeciesThermoDuo : public SpeciesThermo {
	
  public:
    //! Constructor
    SpeciesThermoDuo() {}
    //! Destructor
    virtual ~SpeciesThermoDuo(){}

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
			doublereal* h_RT, doublereal* s_R) const {
      m_thermo1.update(t, cp_R, h_RT, s_R);
      m_thermo2.update(t, cp_R, h_RT, s_R);
    }

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
    virtual doublereal minTemp(int k = -1) const {
      doublereal tm1 = m_thermo1.minTemp();
      doublereal tm2 = m_thermo2.minTemp();
      return (tm1 < tm2 ? tm2 : tm1);
    }

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
    virtual doublereal maxTemp(int k = -1) const {
      doublereal tm1 = m_thermo1.maxTemp();
      doublereal tm2 = m_thermo2.maxTemp();
      return (tm1 < tm2 ? tm1 : tm2);
    }               
         
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
    virtual doublereal refPressure(int k = -1) const {
      return m_p0;
    }

    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     *
     * @param k  Species index
     */
    virtual int reportType(int k) const {
      std::map<int, int>::const_iterator p = speciesToType.find(k);
      if (p != speciesToType.end()) {
	const int type = p->second;
	return type;
      } 
      return -1;
    }

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

    //! Modify parameters for the standard state
    /*!
     * @param index Species index
     * @param c     Vector of coefficients used to set the
     *              parameters for the standard state.
     */
    virtual void modifyParams(int index, doublereal *c) {
      int ctype = reportType(index);
      if (ctype == m_thermo1.ID) {
	m_thermo1.modifyParams(index, c);
      } else if (ctype == m_thermo2.ID) {
	m_thermo2.modifyParams(index,  c);
      } else {
	throw CanteraError("modifyParams", "confused");
      }
    }


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
   * @ingroup spthermo
   */
  template<class SPM>
  class SpeciesThermo1 : public SpeciesThermo {
	
  public:
    //! base constructor
    SpeciesThermo1() : m_pref(0.0) {}
    //! destructor
    virtual ~SpeciesThermo1(){}
    
    //! Install one species into this Species Thermo Manager
    /*!
     * @param name       Name of the species
     * @param sp         Species index
     * @param type       species type in terms of an int
     * @param c          Parameters for the species thermo
     */
    virtual void install(std::string name, int sp, int type, const vector_fp& c) {
      m_thermo.push_back(SPM(sp, c));
      if (m_pref) {
	if (m_thermo.begin()->refPressure() != m_pref) {
	  throw RefPressureMismatch("SpeciesThermo1:install",
				    refPressure(), m_pref);
	}
      }
      else  m_pref = m_thermo.begin()->refPressure();
    }

    //! update the object, because the temperature changed
    /*!
     * @param t          temperature(Kelvin)
     * @param cp_R       vector of dimensionless heat capacity
     * @param h_RT       vector of dimensionless enthalpy
     * @param s_R        vector of dimensionless entropy
     */
    virtual void update(doublereal t, vector_fp& cp_R, 
			vector_fp& h_RT, vector_fp& s_R) const {
      _updateAll(m_thermo.begin(),m_thermo.end(), 
		 t, cp_R, h_RT, s_R);
    }

    //! update the object for one species, because the temperature changed
    /*!
     * @param k          species index
     * @param t          temperature(Kelvin)
     * @param cp_R       vector of dimensionless heat capacity
     * @param h_RT       vector of dimensionless enthalpy
     * @param s_R        vector of dimensionless entropy
     */
    virtual void update_one(int k, doublereal t, vector_fp& cp_R, 
			    vector_fp& h_RT, vector_fp& s_R) const {
      m_thermo[k]->update(t, cp_R, h_RT, s_R);
    }

    //! returns the minimum temperature
    /*!
     * @param k species index. Defaults to -1.
     */
    virtual doublereal minTemp(int k = -1) const {
      if (k < 0)
	return _minTemp(m_thermo.begin(), m_thermo.end());
      else
	return m_thermo[k].minTemp();
    }

    //! returns the maximum temperature
    /*!
     * @param k species index. Defaults to -1.
     */
    virtual doublereal maxTemp(int k = -1) const {
      if (k < 0)
	return _maxTemp(m_thermo.begin(), m_thermo.end());
      else
	return m_thermo[k].maxTemp();
    }   
                    
    //! returns the reference pressure
    /*!
     * @param k species index. Defaults to -1.
     */
    virtual doublereal refPressure(int k = -1) const {
      return m_pref;
    }

    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     *  Note, all parameterizations are the same, by definition, here
     *
     * @param k  Species index
     */
    virtual int reportType(int k) const {
      return m_thermo[k]->reportType(-1);	  
    }

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
			      doublereal &refPressure) const {
      m_thermo[index]->reportParameters(index, type, c, minTemp, maxTemp, refPressure);
    }

    //! Modify parameters for the standard state
    /*!
     * @param index Species index
     * @param c     Vector of coefficients used to set the
     *              parameters for the standard state.
     */
    virtual void modifyParams(int index, doublereal *c) {
      m_thermo[index]->modifyParameters(index, c);
    }

  private:
    //! Vector of SPM objects. There are m_kk of them
    std::vector<SPM> m_thermo;
    //! Reference pressure (Pa)
    doublereal m_pref;
  };
  //#endif

}

#endif
