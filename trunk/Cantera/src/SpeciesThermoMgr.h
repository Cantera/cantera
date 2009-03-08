/**
 *  @file SpeciesThermoMgr.h
 *
 * $Author: dggoodwin $
 * $Revision: 1.7 $
 * $Date: 2005/11/22 17:59:04 $
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_SPECIESTHERMO_MGR_H
#define CT_SPECIESTHERMO_MGR_H

#include "ct_defs.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "SpeciesThermo.h"
#include <map>
using namespace std;

namespace Cantera {

    /**
     * Invokes the 'updateProperties' method of all objects in the
     * list.
     */
    template<class InputIter>
    inline void _updateAll(
        InputIter begin, 
        InputIter end,
        doublereal t,  
        vector_fp&  cp_R, 
        vector_fp& h_RT, 
        vector_fp& s_R) {
        for (; begin != end; ++begin)
            begin->updateProperties(t, cp_R, h_RT, s_R);
    }


    /**
     * Iterates through a list of objects which implement a method
     * 'minTemp()', and returns the largest 'minTemp' value.
    */
    template<class InputIter>
    doublereal _minTemp(InputIter begin, InputIter end) {
        doublereal _minT = 0.0;
        for (; begin != end; ++begin)
            _minT = fmaxx(_minT, begin->minTemp());
        return _minT;
    }
    

    /**
     * Iterates through a list of objects that implement a method
     * 'maxTemp()', and returns the smallest 'maxTemp' value.
    */
    template<class _InputIter>
    doublereal _maxTemp(_InputIter begin, _InputIter end) {
        doublereal _maxT = 1.e10;
        for (; begin != end; ++begin)
            _maxT = fminn(_maxT, begin->maxTemp());
        return _maxT;
    }


    ///////////////////////  Exceptions //////////////////////////////


    /**
     * Exception thrown if species reference pressures don't match.
     * @ingroup spthermo
     */
    class RefPressureMismatch : public CanteraError {
    public:
        RefPressureMismatch(string proc, doublereal prnew, 
            doublereal prold) : CanteraError(proc, 
                "Species reference pressure ("
                + fp2str(prnew) + ") does not match previously-defined "
                + "reference pressure (" + fp2str(prold) + ")") {}
        virtual ~RefPressureMismatch() {}
    };

    class UnknownSpeciesThermo 
        : public CanteraError {
    public:
        UnknownSpeciesThermo(string proc, int type) :
            CanteraError(proc, "Specified species "
                "parameterization type (" + int2str(type) 
                + ") does not match any known type.") {}
        virtual ~UnknownSpeciesThermo() {}
    };



    /**
     *  This species thermo manager requires that all species have one
     *  of two parameterizations.
     */
    template<class T1, class T2>
    class SpeciesThermoDuo : public SpeciesThermo {
	
    public:

	SpeciesThermoDuo() {}
        virtual ~SpeciesThermoDuo(){}
                
	virtual void install(string name, int sp, int type,
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

        virtual void update(doublereal t, doublereal* cp_R,
            doublereal* h_RT, doublereal* s_R) const {
            m_thermo1.update(t, cp_R, h_RT, s_R);
            m_thermo2.update(t, cp_R, h_RT, s_R);
        }
        
        virtual doublereal minTemp(int k = -1) const {
            doublereal tm1 = m_thermo1.minTemp();
            doublereal tm2 = m_thermo2.minTemp();
            return (tm1 < tm2 ? tm2 : tm1);
        }
        
        virtual doublereal maxTemp(int k = -1) const {
            doublereal tm1 = m_thermo1.maxTemp();
            doublereal tm2 = m_thermo2.maxTemp();
            return (tm1 < tm2 ? tm1 : tm2);
        }                        

        virtual doublereal refPressure(int k = -1) const {
            return m_p0;
        }

	virtual int reportType(int k) const {
	    map<int, int>::const_iterator p = speciesToType.find(k);
	    if (p != speciesToType.end()) {
	      const int type = p->second;
	      return type;
	    } 
	    return -1;
	}

	virtual void reportParams(int index, int &type, 
				  doublereal * const c, 
				  doublereal &minTemp, 
				  doublereal &maxTemp, 
				  doublereal &refPressure) {
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

    private:

        T1 m_thermo1;
        T2 m_thermo2;
        doublereal m_p0;
	map<int, int> speciesToType;
    };

    //#define REMOVE_FOR_V155
    //#ifndef REMOVE_FOR_V155

    /**
     *  This species thermo manager requires that all species have the
     *  same parameterization.
     */
    template<class T>
    class SpeciesThermo1 : public SpeciesThermo {
	
    public:

	SpeciesThermo1() : m_pref(0.0) {}
        virtual ~SpeciesThermo1(){}
                
	virtual void install(string name, int sp, int type, const vector_fp& c) {
	    m_thermo.push_back(T(sp, c));
	    if (m_pref) {
		if (m_thermo.begin()->refPressure() != m_pref) {
		    throw RefPressureMismatch("SpeciesThermo1:install",
                        refPressure(), m_pref);
		}
	    }
	    else  m_pref = m_thermo.begin()->refPressure();
	}

        virtual void update(doublereal t, vector_fp& cp_R, 
            vector_fp& h_RT, vector_fp& s_R) const {
            _updateAll(m_thermo.begin(),m_thermo.end(), 
                t, cp_R, h_RT, s_R);
        }

        virtual void update_one(int k, doublereal t, vector_fp& cp_R, 
            vector_fp& h_RT, vector_fp& s_R) const {
            m_thermo[k]->update(t, cp_R, h_RT, s_R);
        }
        
        virtual doublereal minTemp(int k = -1) const {
            if (k < 0)
                return _minTemp(m_thermo.begin(), m_thermo.end());
            else
                return m_thermo[k].minTemp();
        }
        
        virtual doublereal maxTemp(int k = -1) const {
            if (k < 0)
                return _maxTemp(m_thermo.begin(), m_thermo.end());
            else
                return m_thermo[k].maxTemp();
        }                        

        virtual doublereal refPressure(int k = -1) const {
            return m_pref;
        }


	virtual int reportType(int k) const {
	    return m_thermo[k]->reportType(k);
	  
	}


    private:
        vector<T> m_thermo;
        doublereal m_pref;
    };
    //#endif

}

#endif






