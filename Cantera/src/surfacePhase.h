#ifndef CT_SURFACE_PHASE
#define CT_SURFACE_PHASE

#include <map>
using namespace std;

#include "ctml.h"
using namespace ctml;
#include "Phase.h"

namespace Cantera {

    /**
     * Surface phases. This class is analogous to class 'Phase' for 3D phases.
     *
     * @todo Should SurfacePhase and Phase be integrated? 
     */
    class SurfacePhase : public Phase {
    public:
        SurfacePhase() : Phase(), m_s0(-1.0) {}

        virtual ~SurfacePhase() {}
        virtual void freezeSpecies() {
            Phase::freezeSpecies();
            m_work.resize(nSpecies());
        }

        /** 
         * Return the total coverage, summed over all species. 
         * Normally, this should equal 1.0, and in SurfKinetics
         * this method is used to formulate the residual equation to 
         * enforce this condition.
         */  
        doublereal totalCoverage() {
            int k;
            doublereal sum = 0.0;
            getConcentrations(m_work.begin());
            for (k = 0; k < m_kk; k++)
                sum += m_work[k]*m_size[k];
            sum /= m_s0;
            return sum;
        }

        virtual bool ready() const {
            return (Phase::ready() && m_s0 > 0.0);
        }

        // Number of surface sites per unit area.
        doublereal siteDensity() { return m_s0; }

        // Set the site density.
        void setSiteDensity(doublereal s0) { m_s0 = s0; }

        void setCoverages(const doublereal* cov) {
            int k;
            for (k = 0; k < m_kk; k++) {
                m_work[k] = cov[k]*m_s0/m_size[k];
            }
            setConcentrations(m_work.begin());
        }

        void getCoverages(doublereal* cov) const {
            int k;
            getConcentrations(m_work.begin());
            for (k = 0; k < m_kk; k++) 
                cov[k] = m_work[k]*m_size[k]/m_s0;
        }

    protected:

        doublereal m_s0;
        vector_fp m_size;
        mutable vector_fp m_work;
    };

}

#endif
