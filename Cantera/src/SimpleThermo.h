/**
 * $Id$
 */

#ifndef CT_SIMPLETHERMO_H
#define CT_SIMPLETHERMO_H

#include "SpeciesThermoMgr.h"

namespace Cantera {

    /**
     * A simple species thermodynamic property manager.
     */
    class SimpleThermo : public SpeciesThermo {
    
    public:

        const int ID;

        SimpleThermo() :
            ID(SIMPLE),
            m_tlow_max(0.0), 
            m_thigh_min(1.e30),
            m_p0(-1.0), m_nsp(0) {}
            

        virtual ~SimpleThermo() {}

        virtual void install(int index, int type, const doublereal* c,
            doublereal minTemp, doublereal maxTemp, doublereal refPressure) {
            m_logt0.push_back(log(c[0]));
            m_t0.push_back(c[0]);
            m_h0_R.push_back(c[1]/GasConstant);
            m_s0_R.push_back(c[2]/GasConstant);
            m_cp0_R.push_back(c[3]/GasConstant);
            m_nsp++;
            doublereal tlow  = minTemp;
            doublereal thigh = maxTemp;
            m_p0  = refPressure;
            if (tlow > m_tlow_max)    m_tlow_max = tlow;
            if (thigh < m_thigh_min)  m_thigh_min = thigh;
            m_tlow.push_back(tlow);
            m_thigh.push_back(thigh);
        }


        virtual void update(doublereal t, doublereal* cp_R, 
            doublereal* h_RT, doublereal* s_R) const {
            int k;
            doublereal logt = log(t);
            doublereal rt = 1.0/t;
            for (k = 0; k < m_nsp; k++) {
                cp_R[k] = m_cp0_R[k];
                h_RT[k] = rt*(m_h0_R[k] + (t - m_t0[k]) * m_cp0_R[k]);
                s_R[k] = m_s0_R[k] + m_cp0_R[k] * (logt - m_logt0[k]);
            }
        }

        virtual void update_one(int k, doublereal t, doublereal* cp_R, 
				doublereal* h_RT, doublereal* s_R) const {
            doublereal logt = log(t);
            doublereal rt = 1.0/t;
	    cp_R[k] = m_cp0_R[k];
	    h_RT[k] = rt*(m_h0_R[k] + (t - m_t0[k]) * m_cp0_R[k]);
	    s_R[k] = m_s0_R[k] + m_cp0_R[k] * (logt - m_logt0[k]);
        }

        virtual doublereal minTemp(int k=-1) const {
            if (k < 0)
                return m_tlow_max;
            else
                return m_tlow[k];
        }

        virtual doublereal maxTemp(int k=-1) const {
            if (k < 0)
                return m_thigh_min;
            else
                return m_thigh[k];
        }

        virtual doublereal refPressure() const {return m_p0;}

	virtual int reportType(int index) const { return SIMPLE; }

	/**
	 * This utility function reports back the type of 
	 * parameterization and all of the parameters for the 
	 * species, index.
	 *  For the SimpleThermo object, there are 4 coefficients.
	 */
	virtual void reportParams(int index, int &type, 
				  doublereal * const c, 
				  doublereal &minTemp, 
				  doublereal &maxTemp, 
				  doublereal &refPressure) {
	    type = reportType(index);
	    if (type == SIMPLE) {
	      c[0] = m_t0[index];
	      c[1] = m_h0_R[index] * GasConstant;
	      c[2] = m_s0_R[index] * GasConstant;
	      c[3] = m_cp0_R[index] * GasConstant;
	      minTemp = m_tlow[index];
	      maxTemp = m_thigh[index];
	      refPressure = m_p0;
	    }
	}


 protected:

        map<int, int>              m_index;
        doublereal                 m_tlow_max;
        doublereal                 m_thigh_min;
        vector_fp                  m_tlow;
        vector_fp                  m_thigh;
        vector_fp                  m_t0;
        vector_fp                  m_logt0;
        vector_fp                  m_h0_R;
        vector_fp                  m_s0_R;
        vector_fp                  m_cp0_R;
        doublereal                 m_p0;
        int                        m_nsp;

    };

}

#endif
