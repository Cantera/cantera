/**
 *  @file FalloffFactory.cpp
 */

/*  $Author: dggoodwin $
 *  $Date: 2005/11/22 17:59:04 $
 *  $Revision: 1.2 $
 */

// Copyright 2001  California Institute of Technology

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <math.h>

#include "FalloffFactory.h"

namespace Cantera {

    FalloffFactory* FalloffFactory::s_factory = 0;

    /**
     * The 3-parameter Troe falloff parameterization. 
     * This parameterization is
     * defined by
     * \f[ F = F_{cent}^{1/(1 + f_1^2)} \f]
     * where
     * \f[ F_{cent} = (1 - A)\exp(-T/T_3) + A \exp(-T/T_1) \f]
     * \f[ f_1 = (\log_{10} P_r + C) / \left(N - 0.14 
     * (\log_{10} P_r + C)\right) \f]
     * \f[ C = -0.4 - 0.67 \log_{10} F_{cent} \f]
     * \f[ N = 0.75 - 1.27 \log_{10} F_{cent} \f]
     */
    class Troe3 : public Falloff {
    public:
        
        /// Default constructor.
        Troe3() : m_a (0.0), m_rt3 (0.0), m_rt1 (0.0) {}

        // Destructor. Does nothing.
        virtual ~Troe3() {}

        /**
         * Initialize.
         * @param c Coefficient vector of length 3, 
         * with entries \f$ (A, T_3, T_1) \f$
         */
        virtual void init(const vector_fp& c) { 
            m_a  = c[0];
            m_rt3 = 1.0/c[1];
            m_rt1 = 1.0/c[2]; 
        }

        virtual void updateTemp(doublereal T, workPtr work) const {
            doublereal Fcent = (1.0 - m_a) * exp(- T * m_rt3 ) 
                               + m_a * exp(- T * m_rt1 );
            *work = log10( fmaxx( Fcent, SmallNumber ) );
        }

        virtual doublereal F(doublereal pr, const_workPtr work) const {
            doublereal lpr,f1,lgf, cc, nn;
            lpr = log10( fmaxx(pr,SmallNumber) );
            cc = -0.4 - 0.67 * (*work);
            nn = 0.75 - 1.27 * (*work);             
            f1 = ( lpr + cc )/ ( nn - 0.14 * ( lpr + cc ) );
            lgf = (*work) / ( 1.0 + f1 * f1 );
            return pow(10.0, lgf );
        }

        virtual size_t workSize() { return 1; }

    protected:

        doublereal m_a, m_rt3, m_rt1;

    private:

    };



    /**
     * The 4-parameter Troe falloff parameterization. This parameterization is
     * defined by
     *
     * \f[ F = F_{cent}^{1/(1 + f_1^2)} \f]
     * where
     * \f[ F_{cent} = (1 - A)\exp(-T/T_3) + A \exp(-T/T_1) + \exp(-T_2/T) \f]
     * \f[ f_1 = (\log_{10} P_r + C) / \left(N - 0.14 
     * (\log_{10} P_r + C)\right) \f]
     * \f[ C = -0.4 - 0.67 \log_{10} F_{cent} \f]
     * \f[ N = 0.75 - 1.27 \log_{10} F_{cent} \f]
     * 
     */

    class Troe4 : public Falloff {
    public:

        Troe4() : m_a (0.0), m_rt3 (0.0), m_rt1 (0.0),
                         m_t2 (0.0) {}
        virtual ~Troe4() {}
        
        virtual void init(const vector_fp& c) {
            m_a  = c[0];
            m_rt3 = 1.0/c[1];
            m_rt1 = 1.0/c[2]; 
            m_t2 = c[3];
        }

        virtual void updateTemp(doublereal T, workPtr work) const {
            doublereal Fcent = (1.0 - m_a) * exp(- T * m_rt3 ) 
                               + m_a * exp(- T * m_rt1 )
                               + exp(- m_t2 / T );
            *work = log10( fmaxx( Fcent, SmallNumber ) );
        }

        virtual doublereal F(doublereal pr, const_workPtr work) const {
            doublereal lpr,f1,lgf, cc, nn;
            lpr = log10( fmaxx(pr,SmallNumber) );
            cc = -0.4 - 0.67 * (*work);
            nn = 0.75 - 1.27 * (*work);             
            f1 = ( lpr + cc )/ ( nn - 0.14 * ( lpr + cc ) );
            lgf = (*work) / ( 1.0 + f1 * f1 );
            return pow(10.0, lgf );
        }

        virtual size_t workSize() { return 1; }

    protected:

        doublereal m_a, m_rt3, m_rt1;
        doublereal m_t2;

    private:
    };

    /**
     * The 3-parameter SRI falloff function.
     */
    class SRI3 : public Falloff {
      
    public:

        SRI3() {}
        virtual ~SRI3() {}
        
        virtual void init(const vector_fp& c) {
            m_a = c[0];
            m_b = c[1];
            m_c = c[2]; 
        }
  
        virtual void updateTemp(doublereal T, workPtr work) const {
            *work = m_a * exp( - m_b / T);
            if (m_c != 0.0) *work += exp( - T/m_c );
        }

        virtual doublereal F(doublereal pr, const_workPtr work) const {
            doublereal lpr = log10( fmaxx(pr,SmallNumber) );
            doublereal xx = 1.0/(1.0 + lpr*lpr);
            doublereal ff = pow( *work , xx);
            return ff;
        }

        virtual size_t workSize() { return 1; }

    protected:
        doublereal m_a, m_b, m_c;

    private:

    };


    /**
     * The 5-parameter SRI falloff function.
     */
    class SRI5 : public Falloff {
      
    public:
        SRI5() {}
        virtual ~SRI5() {}
        virtual void init(const vector_fp& c) {
            m_a = c[0];
            m_b = c[1];
            m_c = c[2]; 
            m_d  = c[3];
            m_e = c[4];
        }

        virtual void updateTemp(doublereal T, workPtr work) const {
            *work = m_a * exp( - m_b / T);
            if (m_c != 0.0) *work += exp( - T/m_c );
            work[1] = m_d * pow(T,m_e);
        }

        virtual doublereal F(doublereal pr, const_workPtr work) const {
            doublereal lpr = log10( fmaxx(pr,SmallNumber) );
            doublereal xx = 1.0/(1.0 + lpr*lpr);
            return pow( *work, xx) * work[1]; 
        }

        virtual size_t workSize() { return 2; }

    protected:

        doublereal m_a, m_b, m_c;
        doublereal m_d, m_e;

    private:

    };


    /**
     * Wang-Frenklach falloff function.  Reference: Wang, H., and
     * Frenklach, M., Chem. Phys. Lett. vol. 205, 271 (1993).
     */
    class WF93 : public Falloff {
      
    public:
        WF93() {}
        virtual ~WF93() {}

        virtual void init(const vector_fp& c) {
            m_a = c[0];
            m_rt1 = 1.0/c[1];
            m_t2 = c[2]; 
            m_rt3  = 1.0/c[3];
            m_alpha0 = c[4];
            m_alpha1 = c[5];
            m_alpha2 = c[6];
            m_sigma0 = c[7];
            m_sigma1 = c[8];
            m_sigma2 = c[9];
        }
  
        virtual void updateTemp(doublereal T, workPtr work) const {
            work[0] = m_alpha0 + (m_alpha1 + m_alpha2*T)*T; // alpha
            work[1] = m_sigma0 + (m_sigma1 + m_sigma2*T)*T; // sigma
            doublereal Fcent = (1.0 - m_a) * exp(- T * m_rt3 ) 
                               + m_a * exp(- T * m_rt1 ) + exp(-m_t2/T);
            work[2] = log10(Fcent); 
        }

        virtual doublereal F(doublereal pr, const_workPtr work) const {
            doublereal lpr = log10( fmaxx(pr, SmallNumber) );
            doublereal x = (lpr - work[0])/work[1];
            doublereal flog = work[2]/exp(x*x);
            return pow( 10.0, flog);
        }

        virtual size_t workSize() { return 3; }

    protected:

        doublereal m_alpha0, m_alpha1, m_alpha2;
        doublereal m_sigma0, m_sigma1, m_sigma2;
        doublereal m_a, m_rt1, m_t2, m_rt3;

    private:

    };


    Falloff* FalloffFactory::newFalloff(int type, const vector_fp& c) {
        Falloff* f;
        switch(type) {
        case TROE3_FALLOFF:
            f = new Troe3(); break;
        case TROE4_FALLOFF: 
            f = new Troe4(); break;
        case SRI3_FALLOFF: 
            f = new SRI3(); break;
        case SRI5_FALLOFF: 
            f = new SRI5(); break;
        case WF_FALLOFF:
            f = new WF93(); break; 
        default: return 0;
        }   
        f->init(c); 
        return f;
    }

}
