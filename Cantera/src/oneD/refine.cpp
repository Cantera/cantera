
#include <map>
#include <algorithm>
#include "Resid1D.h"

#include "refine.h"

using namespace std;

namespace Cantera {


    template<class M>
    static bool has_key(const M& m, int j) {
        if (m.find(j) != m.end()) return true;
        return false;
    }


    /**
     * Return the square root of machine precision.
     */
    static doublereal eps() {
        doublereal e = 1.0;
        while (1.0 + e != 1.0) e *= 0.5;
        return sqrt(e);
    }

    
    Refiner::Refiner(Resid1D& domain) :
        m_ratio(10.0), m_slope(0.8), m_curve(0.8), m_min_range(0.01),
        m_domain(&domain)
    {
        m_nv = m_domain->nComponents();
        m_active.resize(m_nv, true);
        m_thresh = eps();
    }
            

    int Refiner::analyze(int n, const doublereal* z, 
        const doublereal* x) {

        if (m_domain->nPoints() <= 1) return 0;
        m_nv = m_domain->nComponents();

        //m_ok = false;

        // check consistency
        if (n != m_domain->nPoints()) return -1;

        m_loc.clear();
        m_c.clear();

        /**
         * find locations where cell size ratio is too large.
         */
        int j;
        vector_fp dz(n-1, 0.0);
        dz[0] = z[1] - z[0];
        for (j = 1; j < n-1; j++) {
            dz[j] = z[j+1] - z[j];
            if (dz[j] > m_ratio*dz[j-1]) {
                m_loc[j] = 1;
                m_c["point "+int2str(j)] = 1;
            }
            if (dz[j] < dz[j-1]/m_ratio) {
                m_loc[j-1] = 1;
                m_c["point "+int2str(j-1)] = 1;                
            }
        }

        string name;
        doublereal vmin, vmax, smin, smax, aa, ss;
        doublereal dmax, r;
        vector_fp v(n), s(n-1);
        for (int i = 0; i < m_nv; i++) {
            //cout << i << "   " << m_nv << "  " << m_active[i] << endl;
            if (m_active[i]) {
                name = m_domain->componentName(i);
            
                // get component i at all points
                for (j = 0; j < n; j++) v[j] = value(x, i, j); 

                // slope of component i
                for (j = 0; j < n-1; j++)
                    s[j] = (value(x, i, j+1) - value(x, i, j))/
                           (z[j+1] - z[j]);

                // find the range of values and slopes
            
                vmin = *min_element(v.begin(), v.end());
                vmax = *max_element(v.begin(), v.end());
                smin = *min_element(s.begin(), s.end());
                smax = *max_element(s.begin(), s.end());

                // max absolute values of v and s
                aa = fmaxx(abs(vmax), abs(vmin));
                ss = fmaxx(abs(smax), abs(smin));


                // refine based on component i only if the range of v is
                // greater than a fraction 'min_range' of max |v|. This
                // eliminates components that consist of small fluctuations
                // on a constant background.

                if ((vmax - vmin) > m_min_range*aa) {

                    // maximum allowable difference in value between
                    // adjacent points.
                
                    dmax = m_slope*(vmax - vmin) + m_thresh;
                    for (j = 0; j < n-1; j++) {
                        r = abs(v[j+1] - v[j])/dmax;
                        if (r > 1.0) {
                            m_loc[j] = 1;
                            m_c[name] = 1;
                        }
                    }
                }

                        
                // refine based on the slope of component i only if the
                // range of s is greater than a fraction 'min_range' of max
                // |s|. This eliminates components that consist of small
                // fluctuations on a constant slope background.
            
                if ((smax - smin) > m_min_range*ss) {

                    // maximum allowable difference in slope between
                    // adjacent points.
                    dmax = m_curve*(smax - smin);
                    for (j = 0; j < n-2; j++) {
                        r = abs(s[j+1] - s[j]) / (dmax + m_thresh/dz[j]);
                        if (r > 1.0) {
                            m_c[name] = 1;
                            m_loc[j] = 1;
                            m_loc[j+1] = 1;
                        }
                        //cout << "at point " << j << " slope r = "
                        //     << r << " for " << name << endl
                        //     << "    threshold = " << m_thresh << endl;
                    }
                }
                //cout << name << "  " << m_curve << "  " << smax << "  " << smin << "   " << ss << "   " << m_min_range << endl;
            }
        }
        return m_loc.size();
    }

    double Refiner::value(const double* x, int i, int j) {
        return x[m_domain->index(i,j)];
    }

    void Refiner::show() {
        int nnew = m_loc.size();
        if (nnew > 0) {
            writelog("Refining grid.  "
                "New points inserted after grid points ");
            map<int, int>::const_iterator b = m_loc.begin();
            for (; b != m_loc.end(); ++b) {
                writelog(int2str(b->first)+" ");
            }
            writelog("\n");
            writelog("to resolve ");
            map<string, int>::const_iterator bb = m_c.begin();
            for (; bb != m_c.end(); ++bb) {
                writelog(string(bb->first)+" ");
            }
            writelog("\n");
        }
    }


    int Refiner::getNewGrid(int n, const doublereal* z, 
        int nn, doublereal* zn) {
        int j;
        int nnew = m_loc.size();
        if (nnew + n > nn) {
            throw CanteraError("Refine::getNewGrid",
                "array size too small.");
            return -1;
        }

        int jn = 0;
        if (m_loc.size() == 0) {
            copy(z, z + n,  zn);
            return 0;
        }

        for (j = 0; j < n - 1; j++) {
            zn[jn] = z[j];
            jn++;
            if (has_key(m_loc, j)) {
                zn[jn] = 0.5*(z[j] + z[j+1]);
                jn++;
            }
        }
        zn[jn] = z[n-1];
        return 0;
    }

//         int npts = znew.size();
//         newsoln.resize(npts*ncomp);
//         newsoln = Numeric.zeros((npts, ncomp),'d')
//         for i in range(ncomp):
//             for j in range(npts):
//                 newsoln[j,i] = interp.interp(znew[j],grid,solution[:,i])

//         return (Numeric.array(znew), Numeric.array(znew), newsoln, self.ok)
                               
    
                
            
    
}
