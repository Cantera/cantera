/**
 *  @file FalloffFactory.cpp
 */
// Copyright 2001  California Institute of Technology

#include "cantera/kinetics/FalloffFactory.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

namespace Cantera
{

FalloffFactory* FalloffFactory::s_factory = 0;
mutex_t FalloffFactory::falloff_mutex;

//! The 3-parameter Troe falloff parameterization.
/*!
 * The falloff function defines the value of \f$ F \f$ in the following
 * rate expression
 *
 *  \f[ k = k_{\infty} \left( \frac{P_r}{1 + P_r} \right) F \f]
 *  where
 *  \f[ P_r = \frac{k_0 [M]}{k_{\infty}} \f]
 *
 * This parameterization is  defined by
 * \f[ F = F_{cent}^{1/(1 + f_1^2)} \f]
 * where
 * \f[ F_{cent} = (1 - A)\exp(-T/T_3) + A \exp(-T/T_1) \f]
 *
 * \f[ f_1 = (\log_{10} P_r + C) / \left(N - 0.14
 *             (\log_{10} P_r + C)\right) \f]
 *
 * \f[ C = -0.4 - 0.67 \log_{10} F_{cent} \f]
 *
 * \f[ N = 0.75 - 1.27 \log_{10} F_{cent} \f]
 *
 *  - If \f$ T_3 \f$ is zero, then the corresponding term is set to zero.
 *  - If \f$ T_1 \f$ is zero, then the corresponding term is set to zero.
 *
 * @ingroup falloffGroup
 */
class Troe3 : public Falloff
{
public:
    //! Default constructor.
    Troe3() : m_a(0.0), m_rt3(0.0), m_rt1(0.0) {}

    /**
     * Initialize.
     * @param c Coefficient vector of length 3,
     * with entries \f$ (A, T_3, T_1) \f$
     */
    virtual void init(const vector_fp& c) {
        m_a  = c[0];

        if (c[1] == 0.0) {
            m_rt3 = 1000.;
        } else {
            m_rt3 = 1.0/c[1];
        }
        if (c[2] == 0.0) {
            m_rt1 = 1000.;
        } else {
            m_rt1 = 1.0/c[2];
        }
    }

    //! Update the temperature parameters in the representation
    /*!
     *   The workspace has a length of one
     *
     *   @param T         Temperature (Kelvin)
     *   @param work      Vector of working space representing
     *                    the temperature dependent part of the
     *                    parameterization.
     */
    virtual void updateTemp(doublereal T, doublereal* work) const {
        doublereal Fcent = (1.0 - m_a) * exp(- T * m_rt3)
                           + m_a * exp(- T * m_rt1);
        *work = log10(std::max(Fcent, SmallNumber));
    }

    virtual doublereal F(doublereal pr, const doublereal* work) const {
        doublereal lpr,f1,lgf, cc, nn;
        lpr = log10(std::max(pr,SmallNumber));
        cc = -0.4 - 0.67 * (*work);
        nn = 0.75 - 1.27 * (*work);
        f1 = (lpr + cc)/ (nn - 0.14 * (lpr + cc));
        lgf = (*work) / (1.0 + f1 * f1);
        return pow(10.0, lgf);
    }

    virtual size_t workSize() {
        return 1;
    }

protected:
    //! parameter a in the  4-parameter Troe falloff function. This is
    //! unitless.
    doublereal m_a;

    //! parameter 1/T_3 in the  4-parameter Troe falloff function. This has
    //! units of Kelvin-1
    doublereal m_rt3;

    //! parameter 1/T_1 in the  4-parameter Troe falloff function. This has
    //! units of Kelvin-1.
    doublereal m_rt1;
};

//! The 4-parameter Troe falloff parameterization.
/*!
 * The falloff function defines the value of \f$ F \f$ in the following
 * rate expression
 *
 *  \f[ k = k_{\infty} \left( \frac{P_r}{1 + P_r} \right) F \f]
 *  where
 *  \f[ P_r = \frac{k_0 [M]}{k_{\infty}} \f]
 *
 * This parameterization is defined by
 *
 * \f[ F = F_{cent}^{1/(1 + f_1^2)} \f]
 *    where
 * \f[ F_{cent} = (1 - A)\exp(-T/T_3) + A \exp(-T/T_1) + \exp(-T_2/T) \f]
 *
 * \f[ f_1 = (\log_{10} P_r + C) /
 *              \left(N - 0.14 (\log_{10} P_r + C)\right) \f]
 *
 * \f[ C = -0.4 - 0.67 \log_{10} F_{cent} \f]
 *
 * \f[ N = 0.75 - 1.27 \log_{10} F_{cent} \f]
 *
 *  - If \f$ T_3 \f$ is zero, then the corresponding term is set to zero.
 *  - If \f$ T_1 \f$ is zero, then the corresponding term is set to zero.
 *
 * @ingroup falloffGroup
 */
class Troe4 : public Falloff
{
public:
    //! Constructor
    Troe4() : m_a(0.0), m_rt3(0.0), m_rt1(0.0),
        m_t2(0.0) {}

    //! Initialization of the object
    /*!
     * @param c Vector of four doubles: The doubles are the parameters,
     *          a,, T_3, T_1, and T_2 of the Troe parameterization
     */
    virtual void init(const vector_fp& c) {
        m_a  = c[0];
        if (c[1] == 0.0) {
            m_rt3 = 1000.;
        } else {
            m_rt3 = 1.0/c[1];
        }
        if (c[2] == 0.0) {
            m_rt1 = 1000.;
        } else {
            m_rt1 = 1.0/c[2];
        }
        m_t2 = c[3];
    }

    //! Update the temperature parameters in the representation
    /*!
     *   The workspace has a length of one
     *
     *   @param T         Temperature (Kelvin)
     *   @param work      Vector of working space representing
     *                    the temperature dependent part of the
     *                    parameterization.
     */
    virtual void updateTemp(doublereal T, doublereal* work) const {
        doublereal Fcent = (1.0 - m_a) * exp(- T * m_rt3)
                           + m_a * exp(- T * m_rt1)
                           + exp(- m_t2 / T);
        *work = log10(std::max(Fcent, SmallNumber));
    }

    virtual doublereal F(doublereal pr, const doublereal* work) const {
        doublereal lpr,f1,lgf, cc, nn;
        lpr = log10(std::max(pr,SmallNumber));
        cc = -0.4 - 0.67 * (*work);
        nn = 0.75 - 1.27 * (*work);
        f1 = (lpr + cc)/ (nn - 0.14 * (lpr + cc));
        lgf = (*work) / (1.0 + f1 * f1);
        return pow(10.0, lgf);
    }

    virtual size_t workSize() {
        return 1;
    }

protected:
    //! parameter a in the  4-parameter Troe falloff function. This is
    //! unitless.
    doublereal m_a;

    //! parameter 1/T_3 in the  4-parameter Troe falloff function. This has
    //! units of Kelvin-1.
    doublereal m_rt3;

    //! parameter 1/T_1 in the  4-parameter Troe falloff function. This has
    //! units of Kelvin-1.
    doublereal m_rt1;

    //! parameter T_2 in the  4-parameter Troe falloff function. This has
    //! units of Kelvin.
    doublereal m_t2;
};

//! The 3-parameter SRI falloff function for <I>F</I>
/*!
 * The falloff function defines the value of \f$ F \f$ in the following
 * rate expression
 *
 *  \f[ k = k_{\infty} \left( \frac{P_r}{1 + P_r} \right) F \f]
 *  where
 *  \f[ P_r = \frac{k_0 [M]}{k_{\infty}} \f]
 *
 *  \f[ F = {\left( a \; exp(\frac{-b}{T}) + exp(\frac{-T}{c})\right)}^n \f]
 *      where
 *  \f[ n = \frac{1.0}{1.0 + {\log_{10} P_r}^2} \f]
 *
 *  \f$ c \f$ s required to greater than or equal to zero. If it is zero,
 *  then the corresponding term is set to zero.
 *
 * @ingroup falloffGroup
 */
class SRI3 : public Falloff
{
public:
    //! Constructor
    SRI3() {}

    //! Initialization of the object
    /*!
     * @param c Vector of three doubles: The doubles are the parameters,
     *          a, b, and c of the SRI parameterization
     */
    virtual void init(const vector_fp& c) {
        if (c[2] < 0.0) {
            throw CanteraError("SRI3::init()",
                               "m_c parameter is less than zero: " + fp2str(c[2]));
        }
        m_a = c[0];
        m_b = c[1];
        m_c = c[2];
    }

    //! Update the temperature parameters in the representation
    /*!
     *   The workspace has a length of one
     *
     *   @param T         Temperature (Kelvin)
     *   @param work      Vector of working space representing
     *                    the temperature dependent part of the
     *                    parameterization.
     */
    virtual void updateTemp(doublereal T, doublereal* work) const {
        *work = m_a * exp(- m_b / T);
        if (m_c != 0.0) {
            *work += exp(- T/m_c);
        }
    }

    virtual doublereal F(doublereal pr, const doublereal* work) const {
        doublereal lpr = log10(std::max(pr,SmallNumber));
        doublereal xx = 1.0/(1.0 + lpr*lpr);
        return pow(*work , xx);
    }

    virtual size_t workSize() {
        return 1;
    }

protected:
    //! parameter a in the  3-parameter SRI falloff function. This is
    //! unitless.
    doublereal m_a;

    //! parameter b in the  3-parameter SRI falloff function. This has units
    //! of Kelvin.
    doublereal m_b;

    //! parameter c in the  3-parameter SRI falloff function. This has units
    //! of Kelvin.
    doublereal m_c;
};

//! The 5-parameter SRI falloff function.
/*!
 * The falloff function defines the value of \f$ F \f$ in the following
 * rate expression
 *
 *  \f[ k = k_{\infty} \left( \frac{P_r}{1 + P_r} \right) F \f]
 *  where
 *  \f[ P_r = \frac{k_0 [M]}{k_{\infty}} \f]
 *
 *  \f[ F = {\left( a \; exp(\frac{-b}{T}) + exp(\frac{-T}{c})\right)}^n
 *              \;  d \; exp(\frac{-e}{T}) \f]
 *      where
 *  \f[ n = \frac{1.0}{1.0 + {\log_{10} P_r}^2} \f]
 *
 *  \f$ c \f$ s required to greater than or equal to zero. If it is zero, then
 *  the corresponding term is set to zero.
 *
 *  m_c is required to greater than or equal to zero. If it is zero, then the
 *  corresponding term is set to zero.
 *
 *  m_d is required to be greater than zero.
 *
 * @ingroup falloffGroup
 */
class SRI5 : public Falloff
{
public:
    //! Constructor
    SRI5() {}

    //! Initialization of the object
    /*!
     * @param c Vector of five doubles: The doubles are the parameters,
     *          a, b, c, d, and e of the SRI parameterization
     */
    virtual void init(const vector_fp& c) {
        if (c[2] < 0.0) {
            throw CanteraError("SRI5::init()",
                               "m_c parameter is less than zero: " + fp2str(c[2]));
        }
        if (c[3] < 0.0) {
            throw CanteraError("SRI5::init()",
                               "m_d parameter is less than zero: " + fp2str(c[3]));
        }
        m_a = c[0];
        m_b = c[1];
        m_c = c[2];
        m_d = c[3];
        m_e = c[4];
    }

    //! Update the temperature parameters in the representation
    /*!
     *   The workspace has a length of two
     *
     *   @param T         Temperature (Kelvin)
     *   @param work      Vector of working space representing
     *                    the temperature dependent part of the
     *                    parameterization.
     */
    virtual void updateTemp(doublereal T, doublereal* work) const {
        *work = m_a * exp(- m_b / T);
        if (m_c != 0.0) {
            *work += exp(- T/m_c);
        }
        work[1] = m_d * pow(T,m_e);
    }

    virtual doublereal F(doublereal pr, const doublereal* work) const {
        doublereal lpr = log10(std::max(pr,SmallNumber));
        doublereal xx = 1.0/(1.0 + lpr*lpr);
        return pow(*work, xx) * work[1];
    }

    virtual size_t workSize() {
        return 2;
    }

protected:
    //! parameter a in the 5-parameter SRI falloff function. This is unitless.
    doublereal m_a;

    //! parameter b in the 5-parameter SRI falloff function. This has units of
    //! Kelvin.
    doublereal m_b;

    //! parameter c in the 5-parameter SRI falloff function. This has units of
    //! Kelvin.
    doublereal m_c;

    //! parameter d in the 5-parameter SRI falloff function. This is unitless.
    doublereal m_d;

    //! parameter d in the 5-parameter SRI falloff function. This is unitless.
    doublereal m_e;
};

Falloff* FalloffFactory::newFalloff(int type, const vector_fp& c)
{
    Falloff* f;
    switch (type) {
    case SIMPLE_FALLOFF:
        f = new Falloff();
        break;
    case TROE3_FALLOFF:
        f = new Troe3();
        break;
    case TROE4_FALLOFF:
        f = new Troe4();
        break;
    case SRI3_FALLOFF:
        f = new SRI3();
        break;
    case SRI5_FALLOFF:
        f = new SRI5();
        break;
    default:
        return 0;
    }
    f->init(c);
    return f;
}

}
