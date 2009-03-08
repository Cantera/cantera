/**
 * @file NasaThermo.h
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.15 $
 * $Date: 2006/11/06 15:20:34 $
 */


#ifndef CT_NASATHERMO_H
#define CT_NASATHERMO_H
#include <string>

#include "SpeciesThermoMgr.h"
#include "NasaPoly1.h"
#include "speciesThermoTypes.h"
#include "polyfit.h"
#include "global.h"

namespace Cantera {

    /**
     * A species thermodynamic property manager for the NASA
     * polynomial parameterization with two temperature ranges.
     *
     * This class is designed to efficiently evaluate the properties
     * of a large number of species with the NASA parameterization.
     *
     * The original NASA polynomial parameterization expressed the
     * heat capacity as a fourth-order polynomial in temperature, with
     * separate coefficients for each of two temperature ranges. (The
     * newer NASA format adds coefficients for 1/T and 1/T^2, and
     * allows multiple temperature ranges.) This class is designed for
     * use with the original parameterization, which is used, for
     * example, by the Chemkin software package.
     * 
     * In many cases, the midpoint temperature is the same for many
     * species.  To take advantage of this, class NasaThermo groups
     * species with a common midpoint temperature, so that checking
     * which range the desired temperature is in need be done only
     * once for each group.
     *
     * @note There is a special CTML element for entering the 
     * coefficients of this parameterization.
     * @see importCTML
     */
    class NasaThermo : public SpeciesThermo {
    
    public:

        const int ID;

        NasaThermo() :
            ID(NASA),  
            m_tlow_max(0.0), 
            m_thigh_min(1.e30),
            m_ngroups(0) { m_t.resize(6); }

        virtual ~NasaThermo() {}

        /**
         * Install parameterization for a species.
         * @param index    Species index
         * @param type     ignored, since only NASA type is supported
         * @param c        coefficients. These are
         * - c[0]          midpoint temperature
         * - c[1] - c[7]   coefficients for low T range
         * - c[8] - c[14]  coefficients for high T range
         */
        virtual void install(string name, int index, int type, 
			     const doublereal* c, 
			     doublereal minTemp, doublereal maxTemp,
			     doublereal refPressure) { 

            m_name[index] = name;
            int imid = int(c[0]);       // midpoint temp converted to integer
            int igrp = m_index[imid];   // has this value been seen before?
            if (igrp == 0) {            // if not, prepare new group
                vector<NasaPoly1> v;
                m_high.push_back(v);
                m_low.push_back(v);
                m_tmid.push_back(c[0]);
                m_index[imid] = igrp = static_cast<int>(m_high.size());
                m_ngroups++;
            }

	    m_group_map[index] = igrp;
	    m_posInGroup_map[index] = (int) m_low[igrp-1].size();
	    
            doublereal tlow  = minTemp;
            doublereal tmid  = c[0];
            doublereal thigh = maxTemp;
            doublereal pref  = refPressure;
            const doublereal* clow = c + 1;

            vector_fp chigh(7);
            copy(c + 8, c + 15, chigh.begin());

            m_high[igrp-1].push_back(NasaPoly1(index, tmid, thigh, 
                                         pref, &chigh[0]));
            m_low[igrp-1].push_back(NasaPoly1(index, tlow, tmid, 
                                        pref, clow));

            vector_fp clu(7), chu(7);
            clu[5] = clow[0];
            clu[6] = clow[1];
            copy(clow+2, clow+7, clu.begin());
            chu[5] = chigh[0];
            chu[6] = chigh[1];
            copy(chigh.begin()+2, chigh.begin()+7, chu.begin());

            checkContinuity(name, tmid, &clu[0], &chu[0]);

            if (tlow > m_tlow_max)    m_tlow_max = tlow;
            if (thigh < m_thigh_min)  m_thigh_min = thigh;
            m_tlow.push_back(tlow);
            m_thigh.push_back(thigh);
            m_p0 = pref;
        }

        /** 
         * update the properties for only one species.
         */
        virtual void update_one(int k, doublereal t, doublereal* cp_R, 
            doublereal* h_RT, doublereal* s_R) const {

            m_t[0] = t;
            m_t[1] = t*t;
            m_t[2] = m_t[1]*t;
            m_t[3] = m_t[2]*t;
            m_t[4] = 1.0/t;
            m_t[5] = log(t);
 
	    int grp = m_group_map[k];
	    int pos = m_posInGroup_map[k];
	    const vector<NasaPoly1> &mlg = m_low[grp-1];
	    const NasaPoly1 *nlow = &(mlg[pos]);

            doublereal tmid = nlow->maxTemp();
            if (t < tmid) {
	      nlow->updateProperties(&m_t[0], cp_R, h_RT, s_R);
	    } else {
	      const vector<NasaPoly1> &mhg = m_high[grp-1];
	      const NasaPoly1 *nhigh = &(mhg[pos]);
	      nhigh->updateProperties(&m_t[0], cp_R, h_RT, s_R);
	    }
        }


        virtual void update(doublereal t, doublereal* cp_R, 
            doublereal* h_RT, doublereal* s_R) const {
            int i;

            // load functions of temperature into m_t vector
            m_t[0] = t;
            m_t[1] = t*t;
            m_t[2] = m_t[1]*t;
            m_t[3] = m_t[2]*t;
            m_t[4] = 1.0/t;
            m_t[5] = log(t);

            // iterate over the groups
            vector<NasaPoly1>::const_iterator _begin, _end;
            for (i = 0; i != m_ngroups; i++) {
                if (t > m_tmid[i]) {
                    _begin  = m_high[i].begin();
                    _end    = m_high[i].end();
                }
                else {
                    _begin  = m_low[i].begin();
                    _end    = m_low[i].end();
                }
                for (; _begin != _end; ++_begin) 
                    _begin->updateProperties(&m_t[0], cp_R, h_RT, s_R);
            }
        }
                
        /**
         * Return the lowest temperature at which the thermodynamic
         * parameterization is valid.  If no argument is supplied, the
         * value is the one for which all species parameterizations
         * are valid. Otherwise, if an integer argument is given, the
         * value applies only to the species with that index.
         */
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

        virtual doublereal refPressure(int k = -1) const {
            return m_p0;
        }

	/**
         * This utility function reports the type of parameterization
         * used for the species, index.
         */
        virtual int reportType(int index) const { return NASA; }

	/**
	 * This utility function reports back the type of 
	 * parameterization and all of the parameters for the 
	 * species, index.
	 *  For the NASA object, there are 15 coefficients.
	 */
	virtual void reportParams(int index, int &type, 
				  doublereal * const c, 
				  doublereal &minTemp, 
				  doublereal &maxTemp, 
				  doublereal &refPressure) {
	    type = reportType(index);
	    if (type == NASA) {
                int grp = m_group_map[index];
                int pos = m_posInGroup_map[index];
                const vector<NasaPoly1> &mlg = m_low[grp-1];
                const vector<NasaPoly1> &mhg = m_high[grp-1];
                const NasaPoly1 *lowPoly  = &(mlg[pos]);
                const NasaPoly1 *highPoly = &(mhg[pos]);
                int itype = NASA;
                doublereal tmid = lowPoly->maxTemp();
	      c[0] = tmid;
	      int n;
	      double ttemp;
	      lowPoly->reportParameters(n, itype, minTemp, ttemp, refPressure,
                  c + 1);
	      if (n != index) {
                  throw CanteraError("  ", "confused");
	      }
	      if (itype != NASA1) {
                  throw CanteraError("  ", "confused");
	      }
	      highPoly->reportParameters(n, itype, ttemp, maxTemp, refPressure,
                  c + 8);
	      if (n != index) {
                  throw CanteraError("  ", "confused");
	      }
	      if (itype != NASA1) {
                  throw CanteraError("  ", "confused");
	      }
	    } else {
                throw CanteraError(" ", "confused");
	    }
	}


	/**
	 * This utility function modifies the array of coefficients.
         * The array is the same as that returned by reportParams, so
         * a call can first be made to reportParams to populate the
         * array, and then modifyParams can be called to alter
         * selected values.  For the NASA object, there are 15
         * coefficients.
	 */
	virtual void modifyParams(int index, doublereal *c) {
	    int type = reportType(index);
	    if (type == NASA) {
                int grp = m_group_map[index];
                int pos = m_posInGroup_map[index];
                vector<NasaPoly1> &mlg = m_low[grp-1];
                vector<NasaPoly1> &mhg = m_high[grp-1];
                NasaPoly1 *lowPoly  = &(mlg[pos]);
                NasaPoly1 *highPoly = &(mhg[pos]);
                doublereal tmid = lowPoly->maxTemp();
                if (c[0] != tmid) {
                    throw CanteraError(" ", "Tmid cannot be changed");
                }
                lowPoly->modifyParameters(c + 1);
                highPoly->modifyParameters(c + 8);
                checkContinuity(m_name[index], c[0], c + 1, c + 8);
	    } else {
                throw CanteraError(" ", "confused");
	    }
	}


 protected:

        vector<vector<NasaPoly1> >         m_high;
        vector<vector<NasaPoly1> >         m_low;
        map<int, int>                      m_index;
        vector_fp                          m_tmid;
        doublereal                         m_tlow_max;
        doublereal                         m_thigh_min;
        vector_fp                          m_tlow;
        vector_fp                          m_thigh;
        doublereal                         m_p0;
        int                                m_ngroups;
        mutable vector_fp                  m_t;

	/*
	 * This map takes as its index, the species index in the phase.
	 * It returns the group index, where the temperature polynomials
	 * for that species are stored. group indecises start at 1,
	 * so a decrement is always performed to access vectors.
	 */
	mutable map<int, int>              m_group_map;
	/*
	 * This map takes as its index, the species index in the phase.
	 * It returns the position index within the group, where the 
	 * temperature polynomials for that species are storred.
	 */
	mutable map<int, int>              m_posInGroup_map;
        mutable map<int, string>           m_name;

    private:

        // see SpeciesThermoFactory.cpp for the definition
        void checkContinuity(string name, double tmid, const doublereal* clow,
            doublereal* chigh);

        /// for internal use by checkContinuity
        doublereal enthalpy_RT(double t, const doublereal* c) {
            return c[0] + 0.5*c[1]*t + OneThird*c[2]*t*t 
                + 0.25*c[3]*t*t*t + 0.2*c[4]*t*t*t*t
                + c[5]/t;
        }

        /// for internal use by checkContinuity
        doublereal entropy_R(double t, const doublereal* c) {
            return c[0]*log(t) + c[1]*t + 0.5*c[2]*t*t 
                + OneThird*c[3]*t*t*t + 0.25*c[4]*t*t*t*t
                + c[6];
        }

    };

}

#endif

