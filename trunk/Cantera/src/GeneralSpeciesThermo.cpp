/**
 *  @file GeneralSpeciesThermo.cpp
 *
 */
                
// Copyright 2001-2004  California Institute of Technology
                
#include "GeneralSpeciesThermo.h"
#include "NasaPoly1.h"
#include "NasaPoly2.h"
#include "ShomatePoly.h"
#include "ConstCpPoly.h"
#include "Mu0Poly.h"
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
	m_kk(b.m_kk) {
	m_sp = b.m_sp;
    }

    const GeneralSpeciesThermo& 
    GeneralSpeciesThermo::operator=(const GeneralSpeciesThermo &b) {
	if (&b != this) {
	  m_tlow_max = b.m_tlow_max;
	  m_thigh_min = b.m_thigh_min;
	  m_kk = b.m_kk;
	  m_sp = b.m_sp;
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

    
    /**
     * Install parameterization for a species.
     * @param index    Species index
     * @param type     parameterization type
     * @param c        coefficients. The meaning of these depends on 
     *                 the parameterization.
     */
    void GeneralSpeciesThermo::install(string name,
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
	if (index > m_kk - 1) {
	  m_sp.resize(index+1, 0);
          m_kk = index+1;
	}

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
	default:
	    throw UnknownSpeciesThermoModel(
                "GeneralSpeciesThermo::install",
		"unknown species type", int2str(type));
            break;
	}
	m_tlow_max = max(minTemp, m_tlow_max);
	m_thigh_min = min(maxTemp, m_thigh_min);
    }

    /**
     *  Update the properties for one species.
     */
    void GeneralSpeciesThermo::
    update_one(int k, doublereal t, doublereal* cp_R, 
	       doublereal* h_RT, doublereal* s_R) const {
	SpeciesThermoInterpType * sp_ptr = m_sp[k];
	sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
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
	SpeciesThermoInterpType * sp_ptr;
	for (; _begin != _end; ++_begin) {
	  sp_ptr = *(_begin);
	  sp_ptr->updatePropertiesTemp(t, cp_R, h_RT, s_R);
	}
    }

    /**
     * This utility function reports the type of parameterization
     * used for the species, index.
     */
    int GeneralSpeciesThermo::reportType(int index) const {
	SpeciesThermoInterpType *sp = m_sp[index];
	return sp->reportType();
    }
    
    /**
     * This utility function reports back the type of 
     * parameterization and all of the parameters for the 
     * species, index.
     *  For the NASA object, there are 15 coefficients.
     */
    void GeneralSpeciesThermo::
    reportParams(int index, int &type, 
		 doublereal * const c, 
		 doublereal &minTemp, 
		 doublereal &maxTemp, 
		 doublereal &refPressure) {
	SpeciesThermoInterpType *sp = m_sp[index];
	int n;
	sp->reportParameters(n, type, minTemp, maxTemp, 
			     refPressure, c);      
	if (n != index) {
	  throw CanteraError("  ", "confused");
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
	  return sp->minTemp();	      
	}
    }

    doublereal GeneralSpeciesThermo::maxTemp(int k) const {
	if (k < 0) {
	  return m_thigh_min;
	} else {
	  SpeciesThermoInterpType *sp = m_sp[k];
	  return sp->maxTemp();
	}
    }

    doublereal GeneralSpeciesThermo::refPressure(int k) const {
	if (k < 0) {
	  return m_p0;
	} else {
	  SpeciesThermoInterpType *sp = m_sp[k];
	  return sp->refPressure();
	}
    }
    

}
