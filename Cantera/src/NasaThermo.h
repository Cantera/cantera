/**
 * @file NasaThermo.h
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */


#ifndef CT_NASATHERMO_H
#define CT_NASATHERMO_H

#include "SpeciesThermoMgr.h"
#include "NasaPoly1.h"
#include "speciesThermoTypes.h"
#include "polyfit.h"

namespace Cantera {

    /**
     * A species thermodynamic property manager for the NASA
     * polynomial parameterization with two temperature ranges.
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
        virtual void install(int index, int type, const doublereal* c, 
            doublereal minTemp, doublereal maxTemp, doublereal refPressure) { 

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
            doublereal tlow  = minTemp;
            doublereal tmid  = c[0];
            doublereal thigh = maxTemp;
            doublereal pref  = refPressure;
            const doublereal* clow = c + 1;

            vector_fp chigh(7);
            copy(c + 8, c + 15, chigh.begin());

            checkContinuity(tmid, clow, chigh.begin());

            m_high[igrp-1].push_back(NasaPoly1(index, tmid, thigh, 
                                         pref, chigh.begin()));
            m_low[igrp-1].push_back(NasaPoly1(index, tlow, tmid, 
                                        pref, clow));
            if (tlow > m_tlow_max)    m_tlow_max = tlow;
            if (thigh < m_thigh_min)  m_thigh_min = thigh;
            m_tlow.push_back(tlow);
            m_thigh.push_back(thigh);
            m_p0 = pref;
            m_high_map[index] = &m_high[igrp-1].back();
            m_low_map[index] = &m_low[igrp-1].back();
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
 
            doublereal tmid = m_low_map[k]->maxTemp();
            if (t < tmid)
                m_low_map[k]->updateProperties(m_t.begin(), cp_R, h_RT, s_R);
            else
                m_high_map[k]->updateProperties(m_t.begin(), cp_R, h_RT, s_R);
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
                    _begin->updateProperties(m_t.begin(), cp_R, h_RT, s_R);
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

        virtual doublereal refPressure() const {return m_p0;}

 protected:

        mutable map<int, NasaPoly1*>       m_low_map;
        mutable map<int, NasaPoly1*>       m_high_map;
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

    private:

        // see SpeciesThermoFactory.cpp for the definition
        void checkContinuity(double tmid, const doublereal* clow,
            doublereal* chigh);

        /// for internal use by checkContinuity
        doublereal enthalpy_RT(double t, const doublereal* c) {
            return c[2] + 0.5*c[3]*t + OneThird*c[4]*t*t 
                + 0.25*c[5]*t*t*t + 0.2*c[6]*t*t*t*t
                + c[0]/t;
        }

        /// for internal use by checkContinuity
        doublereal entropy_R(double t, const doublereal* c) {
            return c[2]*log(t) + c[3]*t + 0.5*c[4]*t*t 
                + OneThird*c[5]*t*t*t + 0.25*c[6]*t*t*t*t
                + c[1];
        }

    };

}

#endif

// $Log$
// Revision 1.6  2004-07-01 23:47:45  hkmoffa
// static_cast to eliminate VC++ warnings.
//
// Revision 1.5  2004/04/24 12:53:00  dggoodwin
// *** empty log message ***
//
// Revision 1.4  2004/04/23 19:03:22  dggoodwin
// *** empty log message ***
//
// Revision 1.3  2004/04/22 21:44:36  dggoodwin
// *** empty log message ***
//
// Revision 1.2  2003/11/01 04:50:35  dggoodwin
// *** empty log message ***
//
// Revision 1.1.1.1  2003/04/14 17:57:51  dggoodwin
// Initial import.
//
// Revision 1.16  2003/01/13 10:14:32  dgg
// *** empty log message ***
//
//
// Revision 1.3  2001/12/19 03:14:27  dgg
// Added an offset to the high temperature cp constant coefficient so
// that the high and low cp fits agree precisely at Tmid. This is
// necessary to avoid spurious errors that can occur when integrators
// begin at an initial temperature precisely equal to Tmid.
//
