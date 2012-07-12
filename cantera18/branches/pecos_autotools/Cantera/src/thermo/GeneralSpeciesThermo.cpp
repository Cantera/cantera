/**
 *  @file GeneralSpeciesThermo.cpp
 *  Declarations for a completely general species thermodynamic property
 *  manager for a phase (see \ref spthermo and
 * \link Cantera::GeneralSpeciesThermo GeneralSpeciesThermo\endlink).
 */
/*
 * $Id$
 */
// Copyright 2001-2004  California Institute of Technology
                
#include "GeneralSpeciesThermo.h"
#include "NasaPoly1.h"
#include "NasaPoly2.h"
#include "ShomatePoly.h"
#include "ConstCpPoly.h"
#include "Mu0Poly.h"
#include "StatMech.h"
#ifdef WITH_ADSORBATE
#include "AdsorbateThermo.h"
#endif

#include "SpeciesThermoFactory.h"
#include <iostream>

using namespace std;
                                                           
namespace Cantera {


    /*
     * Constructors
     */
    GeneralSpeciesThermo::GeneralSpeciesThermo() :
	SpeciesThermo(),
	m_tlow_max(0.0), 
	m_thigh_min(1.0E30),
	m_p0(OneAtm),
        m_kk(0) 
    {
	m_tlow_max = 0.0;
	m_thigh_min = 1.0E30;
    }
 
    GeneralSpeciesThermo::
    GeneralSpeciesThermo(const GeneralSpeciesThermo &b) :
	m_tlow_max(b.m_tlow_max), 
	m_thigh_min(b.m_thigh_min),
	m_kk(b.m_kk) 
    {
      m_sp.resize(m_kk, 0);
      for (int k = 0; k < m_kk; k++) {
	SpeciesThermoInterpType *bk = b.m_sp[k];
	if (bk) {
	  m_sp[k] = bk->duplMyselfAsSpeciesThermoInterpType();
	}
      }
    }

    GeneralSpeciesThermo& 
    GeneralSpeciesThermo::operator=(const GeneralSpeciesThermo &b) {
	if (&b != this) {
	  m_tlow_max = b.m_tlow_max;
	  m_thigh_min = b.m_thigh_min;
	
	  for (int k = 0; k < m_kk; k++) {
	    SpeciesThermoInterpType *sp = m_sp[k];
	    if (sp) {
	      delete sp;
	      m_sp[k] = 0;
	    }
	  }
	  m_kk = b.m_kk;
	  m_sp.resize(m_kk, 0);
	  for (int k = 0; k < m_kk; k++) {
	    SpeciesThermoInterpType *bk = b.m_sp[k];
	    if (bk) {
	      m_sp[k] = bk->duplMyselfAsSpeciesThermoInterpType();
	    }
	  }
	}
	return *this;
    }

    GeneralSpeciesThermo::~GeneralSpeciesThermo() {
	for (int k = 0; k < m_kk; k++) {
	  SpeciesThermoInterpType *sp = m_sp[k];
	  if (sp) {
	    delete sp;
	    m_sp[k] = 0;
	  }
	}
    }


    SpeciesThermo *
    GeneralSpeciesThermo::duplMyselfAsSpeciesThermo() const {
	GeneralSpeciesThermo *gsth = new GeneralSpeciesThermo(*this);
	return (SpeciesThermo *) gsth;
    }

    
    /*
     * Install parameterization for a species.
     * @param index    Species index
     * @param type     parameterization type
     * @param c        coefficients. The meaning of these depends on 
     *                 the parameterization.
     */
    void GeneralSpeciesThermo::install(std::string name,
				       int index,
				       int type, 
				       const doublereal* c, 
				       doublereal minTemp,
				       doublereal maxTemp,
				       doublereal refPressure) {
	/*
	 * Resize the arrays if necessary, filling the empty
	 * slots with the zero pointer.
	 */
      
      if(minTemp <= 0.0)
	{
	  throw CanteraError("Error in GeneralSpeciesThermo.cpp",
			     " Cannot take 0 tmin as input. \n\n");
	}

	if (index > m_kk - 1) {
	  m_sp.resize(index+1, 0);
          m_kk = index+1;
	}
	//AssertThrow(m_sp[index] == 0, 
        //		    "Index position isn't null, duplication of assignment: " + int2str(index));

        //int nfreq = 3;
	/*
	 * Create the necessary object
	 */

	switch (type) {
	case NASA1:
	    m_sp[index] = new NasaPoly1(index, minTemp, maxTemp,
					refPressure, c);
	    break;
	case SHOMATE1:
	    m_sp[index] = new ShomatePoly(index, minTemp, maxTemp,
					  refPressure, c);
	    break;
        case CONSTANT_CP:
        case SIMPLE:
	    m_sp[index] = new ConstCpPoly(index, minTemp, maxTemp,
					  refPressure, c);
	    break;
	case MU0_INTERP:
	    m_sp[index] = new Mu0Poly(index, minTemp, maxTemp,
				      refPressure, c);
	    break;
        case SHOMATE2:
	    m_sp[index] = new ShomatePoly2(index, minTemp, maxTemp,
                                           refPressure, c);
	    break;
	case NASA2:
	    m_sp[index] = new NasaPoly2(index, minTemp, maxTemp,
					refPressure, c);
	    break;

	case STAT:
	  m_sp[index] = new StatMech(index, minTemp, maxTemp,
				     refPressure, c, name);
	  break;

#ifdef WITH_ADSORBATE
        case ADSORBATE:
	   m_sp[index] = new Adsorbate(index, minTemp, maxTemp,
                                        refPressure, c);
	    break;
#endif
	default:
	    throw UnknownSpeciesThermoModel(
                "GeneralSpeciesThermo::install",
		"unknown species type", int2str(type));
            break;
	}
        if (!m_sp[index]) {
            cout << "Null m_sp... index = " << index << endl;
            cout << "type = " << type << endl;
        }
	m_tlow_max = max(minTemp, m_tlow_max);
	m_thigh_min = min(maxTemp, m_thigh_min);
    }

  // Install a new species thermodynamic property
  // parameterization for one species.
  /*
   * @param stit_ptr Pointer to the SpeciesThermoInterpType object
   *          This will set up the thermo for one species
   */
  void GeneralSpeciesThermo::install_STIT(SpeciesThermoInterpType *stit_ptr) {
    /*
     * Resize the arrays if necessary, filling the empty
     * slots with the zero pointer.
     */
    if (!stit_ptr) {
      throw CanteraError("GeneralSpeciesThermo::install_STIT",
			 "zero pointer");
    }
    int index = stit_ptr->speciesIndex();
    if (index > m_kk - 1) {
      m_sp.resize(index+1, 0);
      m_kk = index+1;
    }
    AssertThrow(m_sp[index] == 0, 
		"Index position isn't null, duplication of assignment: " + int2str(index));
    /*
     *  Now, simply assign the position
     */
    m_sp[index] = stit_ptr;

    /*
     * Calculate max and min 
     */
    double minTemp = stit_ptr->minTemp();
    double maxTemp = stit_ptr->maxTemp();

    m_tlow_max = max(minTemp, m_tlow_max);
    m_thigh_min = min(maxTemp, m_thigh_min);
  }


  
  void GeneralSpeciesThermo::installPDSShandler(int k, PDSS *PDSS_ptr, 
						VPSSMgr *vpssmgr_ptr) {
    STITbyPDSS *stit_ptr = new STITbyPDSS(k, vpssmgr_ptr, PDSS_ptr);
    install_STIT(stit_ptr);
  }
  
  /**
   *  Update the properties for one species.
   */
  void GeneralSpeciesThermo::
  update_one(int k, doublereal t, doublereal* cp_R, 
	     doublereal* h_RT, doublereal* s_R) const {
    SpeciesThermoInterpType * sp_ptr = m_sp[k];
    if (sp_ptr) {
      sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
    }
  }

    
  /**
   *  Update the properties for all species.
   */
  void GeneralSpeciesThermo::
  update(doublereal t, doublereal* cp_R, 
	 doublereal* h_RT, doublereal* s_R) const {
    vector<SpeciesThermoInterpType *>::const_iterator _begin, _end;
    _begin  = m_sp.begin();
    _end    = m_sp.end();
    SpeciesThermoInterpType * sp_ptr = 0;
    for (; _begin != _end; ++_begin) {
      sp_ptr = *(_begin);
      if (sp_ptr) {
	//std::cout << "I am here\n\n";
	sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
      }
      // else {
      //	writelog("General::update: sp_ptr is NULL!\n");
      //}
    }
  }

  /**
   * This utility function reports the type of parameterization
   * used for the species, index.
   */
  int GeneralSpeciesThermo::reportType(int index) const {
    SpeciesThermoInterpType *sp = m_sp[index];
    if (sp) {
      return sp->reportType();
    }
    return -1;
  }

  /**
   * This utility function reports back the type of 
   * parameterization and all of the parameters for the 
   * species, index.
   *  For the NASA object, there are 15 coefficients.
   */
  void GeneralSpeciesThermo::
  reportParams(int index, int &type, doublereal * const c, 
	       doublereal &minTemp, doublereal &maxTemp, doublereal &refPressure) const {
    SpeciesThermoInterpType *sp = m_sp[index];
    int n;
    if (sp) {
      sp->reportParameters(n, type, minTemp, maxTemp, 
			   refPressure, c);      
      if (n != index) {
	throw CanteraError("GeneralSpeciesThermo::reportParams", 
			   "Internal error encountered");
      }
    } else {
      type = -1;
    }
  }

  //! Modify parameters for the standard state
  /*!
   * @param index Species index
   * @param c     Vector of coefficients used to set the
   *              parameters for the standard state.
   */
  void GeneralSpeciesThermo::
  modifyParams(int index, doublereal *c) {
    SpeciesThermoInterpType *sp = m_sp[index];
    if (sp) {
      sp->modifyParameters(c);
    }
  }


  /**
   * Return the lowest temperature at which the thermodynamic
   * parameterization is valid.  If no argument is supplied, the
   * value is the one for which all species parameterizations
   * are valid. Otherwise, if an integer argument is given, the
   * value applies only to the species with that index.
   */
  doublereal GeneralSpeciesThermo::minTemp(int k) const {
    if (k < 0)
      return m_tlow_max;
    else {
      SpeciesThermoInterpType *sp = m_sp[k];
      if (sp) {
	return sp->minTemp();
      }	      
    }
    return m_tlow_max;
  }

  doublereal GeneralSpeciesThermo::maxTemp(int k) const {
    if (k < 0) {
      return m_thigh_min;
    } else {
      SpeciesThermoInterpType *sp = m_sp[k];
      if (sp) {
	return sp->maxTemp();
      }
    }
    return m_thigh_min;
  }

  doublereal GeneralSpeciesThermo::refPressure(int k) const {
    if (k < 0) {
      return m_p0;
    } else {
      SpeciesThermoInterpType *sp = m_sp[k];
      if (sp) {
	return sp->refPressure();
      } 
    }
    return m_p0;
  }
    

  SpeciesThermoInterpType * GeneralSpeciesThermo::provideSTIT(int k) {
    return (m_sp[k]);
  }

#ifdef H298MODIFY_CAPABILITY
   
  doublereal GeneralSpeciesThermo::reportOneHf298(int k) const {
    SpeciesThermoInterpType * sp_ptr = m_sp[k];
    doublereal h = -1.0;
    if (sp_ptr) {
      h = sp_ptr->reportHf298(0);
    }
    return h;
  }

  void GeneralSpeciesThermo::modifyOneHf298(const int k, const doublereal Hf298New) {
    SpeciesThermoInterpType * sp_ptr = m_sp[k];
    if (sp_ptr) {
      sp_ptr->modifyOneHf298(k, Hf298New);
    }
  }

   
#endif


}
