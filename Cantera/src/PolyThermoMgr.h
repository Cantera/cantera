
#ifndef CT_POLYTHERMMGR_H
#define CT_POLYTHERMOMGR_H

#include "PolyThermo.h"
#include "ctexceptions.h"
#include "polyfit.h"
#include "SpeciesThermo.h"

namespace Cantera {

    /**
     *  A species thermodynamic property manager for polynomial
     *  parameterizations.
     */
    template<int N>
    class PolyThermoMgr : public SpeciesThermo {
    
    public:

        PolyThermoMgr() :
            m_nsp(0),
            m_minTemp(0.0), 
            m_maxTemp(1.e30),
            m_p0(-1.0)
            { m_t.resize(N + 1); }

        virtual ~PolyThermoMgr() {}

        virtual void install(int index, int type, const vector_fp& coeffs, 
            doublereal minTemp = 0.0, 
            doublereal maxTemp = 1.e30,
            doublereal refPressure = OneAtm) {
            m_thermo.push_back(PolyThermo<N>());
            PolyThermo<N>& th = m_thermo.back();
            if (coeffs.size() != N + 3)
                throw CanteraError("PolyThermoMgr::install",
                    string("wrong number of coefficients ")
                    +int2str(coeffs.size()));
            th.setCoefficients(coeffs);
            if (m_nsp > 0 && refPressure != m_p0) {
                throw CanteraError("PolyThermoMgr::install",
                    "reference pressure mismatch");
            }
            else if (m_nsp == 0)
                m_p0 = refPressure;

            m_index.push_back(index);
            m_nsp++;

            if (minTemp > m_minTemp) m_minTemp = minTemp;
            if (maxTemp < m_maxTemp) m_maxTemp = maxTemp;
        }

        virtual void update(doublereal t, vector_fp& cp_R, 
            vector_fp& h_RT, vector_fp& s_R) const {
            int i;

            m_t[0] = log(t);
            m_t[1] = 1.0/t;
            m_t[2] = t;
            //            int i;
            const int nmax = N + 1;
            for (i = 2; i < nmax; i++) m_t[i+1] = m_t[i]*t;
            
            doublereal* bcp = cp_R.begin();
            doublereal* bh  = h_RT.begin();
            doublereal* bs  = s_R.begin();
            doublereal* tt = m_t.begin();
            size_t k;
            for (k = 0; k < m_nsp; k++, bcp++, bh++, bs++) {
                m_thermo[k].updateProperties(tt, bcp, bh, bs);
            }
        }

        virtual doublereal minTemp(int k=-1) const {return m_minTemp;}
        virtual doublereal maxTemp(int k=-1) const {return m_maxTemp;}
        virtual doublereal refPressure() const {return m_p0;}
	virtual int reportType(int index) const { return POLYNOMIAL_4;}

    protected:

        size_t                     m_nsp;
        vector<PolyThermo<N> >     m_thermo;
        vector<int>                m_index;
        doublereal                 m_minTemp;
        doublereal                 m_maxTemp;
        doublereal                 m_p0;
        mutable vector_fp          m_t;
    };

}

#endif


