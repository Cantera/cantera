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
//====================================================================================================================
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
//====================================================================================================================
//! getArrhenius() parses the xml element called Arrhenius.
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
  //====================================================================================================================
  // Construct an LTPspecies object for a liquid tranport property.
  /*
   *    The species transport property is constructed from the XML node, 
   *    \verbatim <propNode>, \endverbatim that is a child of the 
   *    \verbatim <transport> \endverbatim node in the species block and specifies a type of transport
   *    property (like viscosity)
   *
   *   @param   propNode      Pointer to the XML node that contains the property information
   *   @param   name          String containing the species name
   *   @param   tp_ind        enum TransportPropertyType containing the property id that this object 
   *                          is creating a parameterization for (e.g., viscosity)
   *   @param   thermo        const pointer to the ThermoPhase object, which is used to find the temperature.
   */ 
  LTPspecies::LTPspecies(const XML_Node * const propNode, const std::string name, 
                         TransportPropertyType tp_ind, const thermo_t * thermo) :
     m_speciesName(name), 
     m_model(LTP_TD_NOTSET),
     m_property(tp_ind),
     m_thermo(thermo),
     m_mixWeight(1.0)
  {
     if (propNode) {
       if (propNode->hasChild("mixtureWeighting") ) {
         m_mixWeight = getFloat(*propNode, "mixtureWeighting");
       }
     }
  }
  //====================================================================================================================
  // Copy constructor
  LTPspecies::LTPspecies(const LTPspecies &right) 
  {

    *this = right;
}
//====================================================================================================================
// Assignment operator
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
//====================================================================================================================
// Duplication routine
/*
 *  @return  Returns a copy of this routine as a pointer to LTPspecies
 */
LTPspecies* LTPspecies::duplMyselfAsLTPspecies() const
{
    return new LTPspecies(*this);
}
//====================================================================================================================
LTPspecies::~LTPspecies()
{
}
//====================================================================================================================
// Returns the vector of pure species transport property
/*
 *  The pure species transport property (i.e. pure species viscosity)
 *  is returned.  Any temperature and composition dependence will be
 *  adjusted internally according to the information provided by the
 *  subclass object.
 */
doublereal LTPspecies::getSpeciesTransProp()
{
    return 0.0;
}
//====================================================================================================================
// Check to see if the property evaluation will be positive
bool LTPspecies::checkPositive() const
{
    return (m_coeffs[0] > 0);

  }
  //====================================================================================================================
  doublereal LTPspecies::getMixWeight() const 
  {
    return m_mixWeight; 
  }
  //====================================================================================================================
  // Internal model to adjust species-specific properties for composition.
  /*
   *  Currently just a place holder, but this method could take 
   *  the composition from the thermo object and adjust coefficients 
   *  accoding to some unspecified model.
   */
  void LTPspecies::adjustCoeffsForComposition() 
  {
  }
  //====================================================================================================================
  // Construct an LTPspecies object for a liquid tranport property 
  // expressed as a constant value.
  /* The transport property is constructed from the XML node, 
   *  \verbatim <propNode>, \endverbatim that is a child of the
   *  \verbatim <transport> \endverbatim node and specifies a type of
   *  transport property (like viscosity)
   */ 
  LTPspecies_Const::LTPspecies_Const(const XML_Node &propNode, const std::string name, 
				     TransportPropertyType tp_ind, const thermo_t * const thermo) : 
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
//====================================================================================================================
// Copy constructor
LTPspecies_Const::LTPspecies_Const(const LTPspecies_Const& right)
    : LTPspecies()
{
    *this = right; //use assignment operator to do other work
}
//====================================================================================================================
// Assignment operator
LTPspecies_Const& LTPspecies_Const::operator=(const LTPspecies_Const& right)
{
    if (&right != this) {
        LTPspecies::operator=(right);
    }
    return *this;
}
//====================================================================================================================
LTPspecies_Const::~LTPspecies_Const()
{
}
//====================================================================================================================
// Duplication routine
/*
 *  @return  Returns a copy of this routine as a pointer to LTPspecies
 */
LTPspecies* LTPspecies_Const::duplMyselfAsLTPspecies() const
{
    return new LTPspecies_Const(*this);
}
//====================================================================================================================
// Return the (constant) value for this transport property
doublereal LTPspecies_Const::getSpeciesTransProp()
{
    return m_coeffs[0];

  }
  //====================================================================================================================
  // Construct an LTPspecies object for a liquid tranport property 
  // expressed in extended Arrhenius form.
  /*
   *  The transport property is constructed from the XML node, 
   *  \verbatim <propNode>, \endverbatim that is a child of the
   *  \verbatim <transport> \endverbatim node and specifies a type of  transport property (like viscosity)
   *
   *
   *   @param   propNode      Referenc to the XML node that contains the property information.This class
   *                          is assumed to be parameterized by reading XML_Node information.
   *   @param   name          String containing the species name
   *   @param   tp_ind        enum TransportPropertyType containing the property id that this object 
   *                          is creating a parameterization for (e.g., viscosity)
   *   @param   thermo        const pointer to the ThermoPhase object, which is used to find the temperature.
   *
   */ 
  LTPspecies_Arrhenius::LTPspecies_Arrhenius(const XML_Node &propNode, const std::string name, 
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
//====================================================================================================================
// Copy constructor
LTPspecies_Arrhenius::LTPspecies_Arrhenius(const LTPspecies_Arrhenius& right)
    : LTPspecies()
{
    *this = right;
}
//====================================================================================================================
// Assignment operator
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
//====================================================================================================================
// Destructor
LTPspecies_Arrhenius::~LTPspecies_Arrhenius()
{
}
//====================================================================================================================
// Duplication routine
/*
 *  @return  Returns a copy of this routine as a pointer to LTPspecies
 */
LTPspecies* LTPspecies_Arrhenius::duplMyselfAsLTPspecies() const
{
    return new LTPspecies_Arrhenius(*this);
}
//===================================================================================================================
// Return the pure species value for this transport property evaluated
// from the Arrhenius expression
/*
 * In general the Arrhenius expression is
 *
 * \f[
 *      \mu = A T^n \exp(- E / R T).
 * \f]
 *
 * Note that for viscosity, the convention is such that
 * a positive activation energy corresponds to the typical
 * case of a positive argument to the exponential so that
 * the Arrhenius expression is
 *
 * \f[
 *      \mu = A T^n \exp(+ E / R T).
 * \f]
 *
 * Any temperature and composition dependence will be
 *  adjusted internally according to the information provided.
 */
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
  //====================================================================================================================
  // Construct an LTPspecies object for a liquid tranport property expressed as a polynomial in temperature.
  /*
   *  The transport property is constructed from the XML node, \verbatim <propNode>, \endverbatim that is a child of the
   *  \verbatim <transport> \endverbatim node and specifies a type of transport property (like viscosity).
   *
   *
   *   @param   propNode      Referenc to the XML node that contains the property information. This class
   *                          must be parameterized by reading XML_Node information.
   *   @param   name          String containing the species name
   *   @param   tp_ind        enum TransportPropertyType containing the property id that this object 
   *                          is creating a parameterization for (e.g., viscosity)
   *   @param   thermo        const pointer to the ThermoPhase object, which is used to find the temperature.
   *
   */ 
  LTPspecies_Poly::LTPspecies_Poly(const XML_Node &propNode, const std::string name, 
				   TransportPropertyType tp_ind, const thermo_t* thermo) : 
    LTPspecies(&propNode, name, tp_ind, thermo),
    m_temp(-1.0),
    m_prop(0.0)
{
    m_model = LTP_TD_POLY;
    getFloatArray(propNode, m_coeffs, "true", "toSI");
}
//====================================================================================================================
// Copy constructor
LTPspecies_Poly::LTPspecies_Poly(const LTPspecies_Poly& right)
    : LTPspecies()
{
    *this = right;
}
//====================================================================================================================
// Assignment operator
LTPspecies_Poly& LTPspecies_Poly::operator=(const LTPspecies_Poly& right)
{
    if (&right != this) {
        LTPspecies::operator=(right);
        m_temp    = right.m_temp;
        m_prop    = right.m_prop;
    }
    return *this;
}
//====================================================================================================================
LTPspecies_Poly::~LTPspecies_Poly()
{
}
//====================================================================================================================
// Duplication routine
/*
 *  @return  Returns a copy of this routine as a pointer to LTPspecies
 */
LTPspecies* LTPspecies_Poly::duplMyselfAsLTPspecies() const
{
    return new LTPspecies_Poly(*this);
}
//====================================================================================================================
// Return the value for this transport property evaluated from the polynomial expression
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
  //====================================================================================================================
  // Construct an LTPspecies object for a liquid tranport property 
  // expressed as an exponential in temperature.
  /*
   *  The transport property is constructed from the XML node, \verbatim <propNode>, \endverbatim that is a child of the
   *  \verbatim <transport> \endverbatim node and specifies a type of transport property (like viscosity).
   *
   *
   *   @param   propNode      Referenc to the XML node that contains the property information. This class
   *                          must be parameterized by reading XML_Node information.
   *   @param   name          String containing the species name
   *   @param   tp_ind        enum TransportPropertyType containing the property id that this object 
   *                          is creating a parameterization for (e.g., viscosity)
   *   @param   thermo        const pointer to the ThermoPhase object, which is used to find the temperature.
   *
   */ 
  LTPspecies_ExpT::LTPspecies_ExpT(const XML_Node &propNode, const std::string name, TransportPropertyType tp_ind, 
				   const thermo_t* thermo) : 

    LTPspecies(&propNode, name, tp_ind, thermo),
    m_temp(-1.0),
    m_prop(0.0)
{
    m_model = LTP_TD_EXPT;
    getFloatArray(propNode, m_coeffs, "true", "toSI");
}
//====================================================================================================================
// Copy constructor
LTPspecies_ExpT::LTPspecies_ExpT(const LTPspecies_ExpT& right)
    : LTPspecies()
{
    *this = right; //use assignment operator to do other work
}
//====================================================================================================================
// Assignment operator
LTPspecies_ExpT& LTPspecies_ExpT::operator=(const LTPspecies_ExpT& right)
{
    if (&right != this) {
        LTPspecies::operator=(right);
        m_temp    = right.m_temp;
        m_prop    = right.m_prop;
    }
    return *this;
}
//====================================================================================================================
LTPspecies_ExpT::~LTPspecies_ExpT()
{
}
//====================================================================================================================
// Duplication routine
/*
 *  @return  Returns a copy of this routine as a pointer to LTPspecies
 */
LTPspecies* LTPspecies_ExpT::duplMyselfAsLTPspecies() const
{
    return new LTPspecies_ExpT(*this);
}
//====================================================================================================================
// Return the value for this transport property evaluated
// from the exponential in temperature expression
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
//====================================================================================================================
}
