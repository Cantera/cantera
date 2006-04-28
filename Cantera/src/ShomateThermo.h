/**
 * @file ShomateThermo.h
 * 
 * This parameterization requires 7 coefficients A - G:
 *
 *  \f[ C_p = A + B*t + C*t2 + D*t3 + E/t^2 \f]
 *
 *  \f[ H - H_298.15= A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 - E/t + F 
 *                 - \Delta_f H_{f,298} \f]
 *
 *  \f[ S = A*ln(t) + B*t + C*t^2/2 + D*t^3/3 - E/(2*t^2) + G \f]
 *
 *    - Cp = heat capacity (J/gmol*K)
 *    - H = standard enthalpy (kJ/gmol)
 *    - \f$ \Delta_f H_298.15 \f$ = enthalpy of formation at 298.15 K (kJ/gmol)
 *    - S = standard entropy (J/gmol*K)
 *    - t = temperature (K) / 1000.
 *
 *  Note, the polynomial data (i.e., A, ... , G) is entered in dimensional form.
 *  This is in contrast to the NASA database polynomials which are entered in 
 *  nondimensional form (i.e., NASA parameterizes C_p/R, while Shomate
 *  parameterizes C_p assuming units of J/gmol*K - and kJ/gmol*K for H).
 *  Note, also that the H - H_298.15 equation has units of kJ/gmol, because of
 *  the implicit integration of (t = T 1000), which provides a 
 *  multiplier of 1000 to the Enthalpy equation.
 *
 */

#ifndef CT_SHOMATETHERMO_H
#define CT_SHOMATETHERMO_H

#include "SpeciesThermoMgr.h"
#include "ShomatePoly.h"
#include "speciesThermoTypes.h"

namespace Cantera {

    /**
     * A species thermodynamic property manager for the Shomate
     * polynomial parameterization. This is the parameterization used
     * in the NIST Chemistry WebBook (http://webbook.nist.gov/chemistry)
     */
    class ShomateThermo : public SpeciesThermo {
    
    public:

        const int ID;
 
        ShomateThermo() :
            ID(SHOMATE),
            m_tlow_max(0.0), 
            m_thigh_min(1.e30),
            m_ngroups(0) { m_t.resize(7); }

        virtual ~ShomateThermo() {}

        /**
         * Install values for a new species.
         * @param index  Species index
         * @param type   ignored, since only Shomate type is supported
         * @param c      coefficients. These are parameters A through G
         * in the same units as used in the NIST Chemistry WebBook.
         *
         */
        virtual void install(string name, int index, int type, 
			     const doublereal* c,
			     doublereal minTemp, doublereal maxTemp, 
			     doublereal refPressure) {
            int imid = int(c[0]);       // midpoint temp converted to integer
            int igrp = m_index[imid];   // has this value been seen before?
            if (igrp == 0) {            // if not, prepare new group
                vector<ShomatePoly> v;
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
            const doublereal* chigh = c + 8;
            m_high[igrp-1].push_back(ShomatePoly(index, tmid, thigh, 
                                         pref, chigh));
            m_low[igrp-1].push_back(ShomatePoly(index, tlow, tmid, 
                                        pref, clow));
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

	    doublereal tt = 1.e-3*t;
            m_t[0] = tt;
            m_t[1] = tt*tt;
            m_t[2] = m_t[1]*tt;
            m_t[3] = 1.0/m_t[1];
            m_t[4] = log(tt);
            m_t[5] = 1.0/GasConstant;
            m_t[6] = 1.0/(GasConstant * t);

	    int grp = m_group_map[k];
	    int pos = m_posInGroup_map[k];
	    const vector<ShomatePoly> &mlg = m_low[grp-1];
	    const ShomatePoly *nlow = &(mlg[pos]);

            doublereal tmid = nlow->maxTemp();
            if (t < tmid) {
	      nlow->updateProperties(&m_t[0], cp_R, h_RT, s_R);
            } else {
	      const vector<ShomatePoly> &mhg = m_high[grp-1];
	      const ShomatePoly *nhigh = &(mhg[pos]);
	      nhigh->updateProperties(&m_t[0], cp_R, h_RT, s_R);
	    }
        }

        virtual void update(doublereal t, doublereal* cp_R, 
            doublereal* h_RT, doublereal* s_R) const {
            int i;

            doublereal tt = 1.e-3*t;
            m_t[0] = tt;
            m_t[1] = tt*tt;
            m_t[2] = m_t[1]*tt;
            m_t[3] = 1.0/m_t[1];
            m_t[4] = log(tt);
            m_t[5] = 1.0/GasConstant;
            m_t[6] = 1.0/(GasConstant * t);

            vector<ShomatePoly>::const_iterator _begin, _end;
            for (i = 0; i != m_ngroups; i++) {
                if (t > m_tmid[i]) {
                    _begin  = m_high[i].begin();
                    _end    = m_high[i].end();
                }
                else {
                    _begin  = m_low[i].begin();
                    _end    = m_low[i].end();
                }
                for (; _begin != _end; ++_begin) {
                  _begin->updateProperties(&m_t[0], cp_R, h_RT, s_R);
                }
            }
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

        virtual doublereal refPressure(int k=-1) const {
            return m_p0;
        }

	virtual int reportType(int index) const { return SHOMATE; }

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
	    if (type == SHOMATE) {
	      int grp = m_group_map[index];
	      int pos = m_posInGroup_map[index];
	      int itype = SHOMATE;
	      const vector<ShomatePoly> &mlg = m_low[grp-1];
	      const vector<ShomatePoly> &mhg = m_high[grp-1];
	      const ShomatePoly *lowPoly  = &(mlg[pos]);
	      const ShomatePoly *highPoly = &(mhg[pos]);
	      doublereal tmid = lowPoly->maxTemp();
	      c[0] = tmid;
	      int n;
	      double ttemp;
	      lowPoly->reportParameters(n, itype, minTemp, ttemp, refPressure,
					c + 1);
	      if (n != index) {
		throw CanteraError("  ", "confused");
	      }
	      if (itype != SHOMATE && itype != SHOMATE1) {
		throw CanteraError("  ", "confused");
	      }
	      highPoly->reportParameters(n, itype,  ttemp, maxTemp,
					 refPressure, c + 8);
	      if (n != index) {
		throw CanteraError("  ", "confused");
	      }
	      if (itype != SHOMATE && itype != SHOMATE1) {
		throw CanteraError("  ", "confused");
	      }
	    } else {
	      throw CanteraError(" ", "confused");
	    }
	}

 protected:

        //mutable map<int, ShomatePoly*>       m_low_map;
        //mutable map<int, ShomatePoly*>       m_high_map;
        vector<vector<ShomatePoly> > m_high;
        vector<vector<ShomatePoly> > m_low;
        map<int, int>              m_index;
        vector_fp                  m_tmid;
        doublereal                 m_tlow_max;
        doublereal                 m_thigh_min;
        vector_fp                  m_tlow;
        vector_fp                  m_thigh;
        doublereal                 m_p0;
        int                        m_ngroups;
        mutable vector_fp          m_t;

	/*
	 * This map takes as its index, the species index in the phase.
	 * It returns the group index, where the temperature polynomials
	 * for that species are storred. group indecises start at 1,
	 * so a decrement is always performed to access vectors.
	 */
	mutable map<int, int>              m_group_map;

	/*
	 * This map takes as its index, the species index in the phase.
	 * It returns the position index within the group, where the 
	 * temperature polynomials for that species are storred.
	 */
	mutable map<int, int>              m_posInGroup_map;
    };

}

#endif
