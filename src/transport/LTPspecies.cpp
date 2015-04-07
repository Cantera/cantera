/**
 *  @file LTPspecies.cpp \
 *     definitions for the  LTPspecies objects and its children, which is the virtual base class
 *     for describing temperature dependence of submodels for transport parameters
 *    (see \ref tranprops and \link Cantera::LTPspecies LTPspecies \endlink) .
 */

#include "cantera/transport/LTPspecies.h"
using namespace std;
using namespace ctml;

namespace Cantera
{

//! Exception thrown if an error is encountered while reading the transport database.
class LTPError : public CanteraError
{
public:

    //! Constructor is a wrapper around CanteraError
    /*!
     *  @param msg    Informative message
     */
    explicit LTPError(const std::string& msg) :
        CanteraError("LTPspecies", "error parsing transport data: " + msg + "\n") {
    }
};

//! Parses the xml element called Arrhenius.
/*!
 * The Arrhenius expression is
 *    \f[
 *         k =  A T^(b) exp (-E_a / RT)
 *    \f]
 *
 *   @param    node       XML_Node to be read
 *   @param    A          Output pre-exponential factor. The units are variable.
 *   @param    b          output temperature power
 *   @param    E          Output activation energy in units of Kelvin
 */
static void getArrhenius(const XML_Node& node,
                         doublereal& A, doublereal& b, doublereal& E)
{
    /* parse the children for the A, b, and E components.
     */
    A = getFloat(node, "A", "toSI");
    b = getFloat(node, "b");
    E = getFloat(node, "E", "actEnergy");
    E /= GasConstant;

}

LTPspecies::LTPspecies(const XML_Node* const propNode, const std::string name,
                       TransportPropertyType tp_ind, const thermo_t* thermo) :
    m_speciesName(name),
    m_model(LTP_TD_NOTSET),
    m_property(tp_ind),
    m_thermo(thermo),
    m_mixWeight(1.0)
{
    if (propNode) {
        if (propNode->hasChild("mixtureWeighting")) {
            m_mixWeight = getFloat(*propNode, "mixtureWeighting");
        }
    }
}

LTPspecies::LTPspecies(const LTPspecies& right)
{

    *this = right;
}

LTPspecies& LTPspecies::operator=(const LTPspecies& right)
{
    if (&right != this) {
        m_speciesName = right.m_speciesName;
        m_property    = right.m_property;
        m_model       = right.m_model;
        m_coeffs      = right.m_coeffs;
        m_thermo      = right.m_thermo;
        m_mixWeight   = right.m_mixWeight;
    }
    return *this;
}

LTPspecies* LTPspecies::duplMyselfAsLTPspecies() const
{
    return new LTPspecies(*this);
}

LTPspecies::~LTPspecies()
{
}

doublereal LTPspecies::getSpeciesTransProp()
{
    return 0.0;
}

bool LTPspecies::checkPositive() const
{
    return (m_coeffs[0] > 0);

}

doublereal LTPspecies::getMixWeight() const
{
    return m_mixWeight;
}

void LTPspecies::adjustCoeffsForComposition()
{
}

LTPspecies_Const::LTPspecies_Const(const XML_Node& propNode, const std::string name,
                                   TransportPropertyType tp_ind, const thermo_t* const thermo) :
    LTPspecies(&propNode, name, tp_ind, thermo)
{
    m_model = LTP_TD_CONSTANT;
    double A_k = getFloatCurrent(propNode, "toSI");
    if (A_k > 0.0) {
        m_coeffs.push_back(A_k);
    } else {
        throw LTPError("negative or zero " + propNode.name());
    }
}

LTPspecies_Const::LTPspecies_Const(const LTPspecies_Const& right)
    : LTPspecies()
{
    *this = right; //use assignment operator to do other work
}

LTPspecies_Const& LTPspecies_Const::operator=(const LTPspecies_Const& right)
{
    if (&right != this) {
        LTPspecies::operator=(right);
    }
    return *this;
}

LTPspecies* LTPspecies_Const::duplMyselfAsLTPspecies() const
{
    return new LTPspecies_Const(*this);
}

doublereal LTPspecies_Const::getSpeciesTransProp()
{
    return m_coeffs[0];

}

LTPspecies_Arrhenius::LTPspecies_Arrhenius(const XML_Node& propNode, const std::string name,
        TransportPropertyType tp_ind,  const thermo_t* thermo) :
    LTPspecies(&propNode, name, tp_ind, thermo)
{

    m_model = LTP_TD_ARRHENIUS;
    m_temp = 0.0;
    m_prop = 0.0;

    doublereal A_k, n_k, Tact_k;
    getArrhenius(propNode, A_k, n_k, Tact_k);
    if (A_k <= 0.0) {
        throw LTPError("negative or zero " + propNode.name());
    }
    m_coeffs.push_back(A_k);
    m_coeffs.push_back(n_k);
    m_coeffs.push_back(Tact_k);
    m_coeffs.push_back(log(A_k));
}

LTPspecies_Arrhenius::LTPspecies_Arrhenius(const LTPspecies_Arrhenius& right)
    : LTPspecies()
{
    *this = right;
}

LTPspecies_Arrhenius& LTPspecies_Arrhenius::operator=(const LTPspecies_Arrhenius& right)
{
    if (&right != this) {
        LTPspecies::operator=(right);
        m_temp    = right.m_temp;
        m_logt    = right.m_logt;
        m_prop    = right.m_prop;
        m_logProp = right.m_logProp;
    }
    return *this;
}

LTPspecies* LTPspecies_Arrhenius::duplMyselfAsLTPspecies() const
{
    return new LTPspecies_Arrhenius(*this);
}

doublereal LTPspecies_Arrhenius::getSpeciesTransProp()
{
    doublereal t = m_thermo->temperature();
    //m_coeffs[0] holds A
    //m_coeffs[1] holds n
    //m_coeffs[2] holds Tact
    //m_coeffs[3] holds log(A)
    if (t != m_temp) {
        m_prop = 0;
        m_logProp = 0;
        m_temp = t;
        m_logt = log(m_temp);
        //For viscosity the sign convention on positive activation energy is swithced
        if (m_property == TP_VISCOSITY) {
            m_logProp = m_coeffs[3] + m_coeffs[1] * m_logt + m_coeffs[2] / m_temp ;
        } else {
            m_logProp = m_coeffs[3] + m_coeffs[1] * m_logt - m_coeffs[2] / m_temp ;
        }
        m_prop = exp(m_logProp);
    }
    return m_prop;
}

LTPspecies_Poly::LTPspecies_Poly(const XML_Node& propNode, const std::string name,
                                 TransportPropertyType tp_ind, const thermo_t* thermo) :
    LTPspecies(&propNode, name, tp_ind, thermo),
    m_temp(-1.0),
    m_prop(0.0)
{
    m_model = LTP_TD_POLY;
    getFloatArray(propNode, m_coeffs, "true", "toSI");
}

LTPspecies_Poly::LTPspecies_Poly(const LTPspecies_Poly& right)
    : LTPspecies()
{
    *this = right;
}

LTPspecies_Poly& LTPspecies_Poly::operator=(const LTPspecies_Poly& right)
{
    if (&right != this) {
        LTPspecies::operator=(right);
        m_temp    = right.m_temp;
        m_prop    = right.m_prop;
    }
    return *this;
}

LTPspecies* LTPspecies_Poly::duplMyselfAsLTPspecies() const
{
    return new LTPspecies_Poly(*this);
}

doublereal LTPspecies_Poly::getSpeciesTransProp()
{
    doublereal t = m_thermo->temperature();
    if (t != m_temp) {
        m_prop = 0.0;
        m_temp = t;
        double tempN = 1.0;
        for (int i = 0; i < (int) m_coeffs.size() ; i++) {
            m_prop += m_coeffs[i] * tempN;
            tempN *= m_temp;
        }
    }
    return m_prop;
}

LTPspecies_ExpT::LTPspecies_ExpT(const XML_Node& propNode, const std::string name, TransportPropertyType tp_ind,
                                 const thermo_t* thermo) :

    LTPspecies(&propNode, name, tp_ind, thermo),
    m_temp(-1.0),
    m_prop(0.0)
{
    m_model = LTP_TD_EXPT;
    getFloatArray(propNode, m_coeffs, "true", "toSI");
}

LTPspecies_ExpT::LTPspecies_ExpT(const LTPspecies_ExpT& right)
    : LTPspecies()
{
    *this = right; //use assignment operator to do other work
}

LTPspecies_ExpT& LTPspecies_ExpT::operator=(const LTPspecies_ExpT& right)
{
    if (&right != this) {
        LTPspecies::operator=(right);
        m_temp    = right.m_temp;
        m_prop    = right.m_prop;
    }
    return *this;
}

LTPspecies* LTPspecies_ExpT::duplMyselfAsLTPspecies() const
{
    return new LTPspecies_ExpT(*this);
}

doublereal LTPspecies_ExpT::getSpeciesTransProp()
{
    doublereal t = m_thermo->temperature();
    if (t != m_temp) {
        m_temp=t;
        m_prop = m_coeffs[0];
        doublereal tempN = 1.0;
        doublereal tmp = 0.0;
        for (int i = 1; i < (int) m_coeffs.size() ; i++) {
            tempN *= m_temp;
            tmp += m_coeffs[i] * tempN;
        }
        m_prop *= exp(tmp);
    }
    return m_prop;
}

}
