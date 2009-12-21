/**
 *  @file LiquidTransportData.cpp
 *  Source code for liquid transport property evaluations.
 */
/*
 *  $Author: jchewso $
 *  $Date: 2009-11-16 17:34:46 -0700 (Mon, 16 Nov 2009) $
 *  $Revision: 261 $
 *
 */

#include "LiquidTransportData.h"


namespace Cantera {



  /** 
   * Exception thrown if an error is encountered while reading the 
   * transport database.
   */
  class LTPError : public CanteraError {
  public:
    LTPError( std::string msg ) 
      : CanteraError("LTPspecies",
		     "error parsing transport data: " 
		     + msg + "\n") {}
  };


  /**
   * getArrhenius() parses the xml element called Arrhenius. 
   * The Arrhenius expression is
   * \f[        k =  A T^(b) exp (-E_a / RT). \f]
   */
  static void getArrhenius(const XML_Node& node, 
			   doublereal& A, doublereal& b, doublereal& E) {
    /* parse the children for the A, b, and E conponents.
     */
    A = getFloat(node, "A", "toSI");
    b = getFloat(node, "b");
    E = getFloat(node, "E", "actEnergy");
    E /= GasConstant;
  }                



  //! Copy constructor
  LiquidTransportData::LiquidTransportData( const LiquidTransportData &right ) 
  {
    *this = right; //use assignment operator to do other work
  }
  
  //! Assignment operator
  LiquidTransportData& LiquidTransportData::operator=(const LiquidTransportData& right ) 
  {
    if (&right != this) {
      speciesName        = right.speciesName;
      hydroRadius        = right.hydroRadius;
      viscosity          = right.viscosity;
      thermalCond        = right.thermalCond;
      electCond          = right.electCond;
      speciesDiffusivity = right.speciesDiffusivity;
    }
    return *this; 
  }




  //! Copy constructor
  LTPspecies::LTPspecies( const LTPspecies &right ) 
  {
    *this = right; //use assignment operator to do other work
  }
  
  //! Assignment operator
  LTPspecies& LTPspecies::operator=(const LTPspecies& right ) 
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


  //! Construct an LTPspecies object for a liquid tranport property 
  //! expressed as a constant value.
  /** The transport property is constructed from the XML node, 
   *  \verbatim <propNode>, \endverbatim that is a child of the
   *  \verbatim <transport> \endverbatim node and specifies a type of
   *  transport property (like viscosity)
   */ 
  LTPspecies_Const::LTPspecies_Const( const XML_Node &propNode, 
				      std::string name, 
				      TransportPropertyList tp_ind, 
				      thermo_t* thermo ) : 
    LTPspecies( propNode, name, tp_ind, thermo) 
  {
    m_model = LTR_MODEL_CONSTANT;	
    double A_k = getFloatCurrent(propNode, "toSI");
    if (A_k > 0.0) {
      m_coeffs.push_back(A_k);
    } else throw LTPError("negative or zero " + propNode.name() );
  }

  //! Copy constructor
  LTPspecies_Const::LTPspecies_Const( const LTPspecies_Const &right ) 
    : LTPspecies()
  {
    *this = right; //use assignment operator to do other work
  }
  
  //! Assignment operator
  LTPspecies_Const& LTPspecies_Const::operator=(const LTPspecies_Const& right ) 
  {
    if (&right != this) {
      //LTPspecies::operator=(right);
      m_speciesName = right.m_speciesName;
      m_property    = right.m_property;
      m_model       = right.m_model;
      m_coeffs      = right.m_coeffs;
      m_thermo      = right.m_thermo;
      m_mixWeight   = right.m_mixWeight;
    }
    return *this; 
  }

  //! Return the (constant) value for this transport property
  doublereal LTPspecies_Const::getSpeciesTransProp( ) {
    return m_coeffs[0];
  }

  ///////////////////////////////////////////////////////////////

  //! Construct an LTPspecies object for a liquid tranport property 
  //! expressed in extended Arrhenius form.
  /** The transport property is constructed from the XML node, 
   *  \verbatim <propNode>, \endverbatim that is a child of the
   *  \verbatim <transport> \endverbatim node and specifies a type of
   *  transport property (like viscosity)
   */ 
  LTPspecies_Arrhenius::LTPspecies_Arrhenius( const XML_Node &propNode, 
				      std::string name, 
				      TransportPropertyList tp_ind, 
				      thermo_t* thermo ) : 
    LTPspecies( propNode, name, tp_ind, thermo) 
  {
    m_model = LTR_MODEL_ARRHENIUS;	

    doublereal A_k, n_k, Tact_k;
    getArrhenius(propNode, A_k, n_k, Tact_k);
    if (A_k <= 0.0) {
      throw LTPError("negative or zero " + propNode.name() );
    }
    m_coeffs.push_back( A_k );
    m_coeffs.push_back( n_k );
    m_coeffs.push_back( Tact_k );
    m_coeffs.push_back( log( A_k ) );
  }
  
  //! Copy constructor
  LTPspecies_Arrhenius::LTPspecies_Arrhenius( const LTPspecies_Arrhenius &right ) 
    : LTPspecies()
  {
    *this = right; //use assignment operator to do other work
  }
  
  //! Assignment operator
  LTPspecies_Arrhenius& LTPspecies_Arrhenius::operator=(const LTPspecies_Arrhenius& right )
  {
    if (&right != this) {
      // LTPspecies::operator=(right);
      m_speciesName = right.m_speciesName;
      m_property    = right.m_property;
      m_model       = right.m_model;
      m_coeffs      = right.m_coeffs;
      m_thermo      = right.m_thermo;
      m_mixWeight   = right.m_mixWeight;

      m_temp    = right.m_temp;
      m_logt    = right.m_logt;
      m_prop    = right.m_prop;
      m_logProp = right.m_logProp;
    }
    return *this; 
  }

  //! Return the pure species value for this transport property evaluated 
  //! from the Arrhenius expression
  /**
   * In general the Arrhenius expression is 
   *
   * \f[
   *      \mu = A T^n \exp( - E / R T ).
   * \f]
   *
   * Note that for viscosity, the convention is such that 
   * a positive activation energy corresponds to the typical 
   * case of a positive argument to the exponential so that 
   * the Arrhenius expression is
   *
   * \f[
   *      \mu = A T^n \exp( + E / R T ).
   * \f]
   *
   * Any temperature and composition dependence will be 
   *  adjusted internally according to the information provided.
   */
  doublereal LTPspecies_Arrhenius::getSpeciesTransProp( ) {
    
    doublereal t = m_thermo->temperature();
    //m_coeffs[0] holds A
    //m_coeffs[1] holds n
    //m_coeffs[2] holds Tact
    //m_coeffs[3] holds log(A)
    if (t != m_temp) {    
      m_temp = t;
      m_logt = log(m_temp);
      //For viscosity the sign convention on positive activation energy is swithced
      if ( m_property == TP_VISCOSITY ) 
	m_logProp = m_coeffs[3] + m_coeffs[1] * m_logt + m_coeffs[2] / m_temp ;
      else
	m_logProp = m_coeffs[3] + m_coeffs[1] * m_logt - m_coeffs[2] / m_temp ;
      m_prop = exp( m_logProp );
    }
    return m_prop;
  }






  ///////////////////////////////////////////////////////////////

  //! Construct an LTPspecies object for a liquid tranport property 
  //! expressed as a polynomial in temperature.
  /** The transport property is constructed from the XML node, 
   *  \verbatim <propNode>, \endverbatim that is a child of the
   *  \verbatim <transport> \endverbatim node and specifies a type of
   *  transport property (like viscosity)
   */ 
  LTPspecies_Poly::LTPspecies_Poly( const XML_Node &propNode, 
				    std::string name, 
				    TransportPropertyList tp_ind, 
				    thermo_t* thermo ) : 
    LTPspecies( propNode, name, tp_ind, thermo) 
  {
    m_model = LTR_MODEL_POLY;	


    getFloatArray(propNode, m_coeffs, true); // if units labeled, convert Angstroms -> meters

    if (m_coeffs[0] <= 0.0) {
      throw LTPError("negative or zero " + propNode.name() );
    }
  }

  //! Copy constructor
  LTPspecies_Poly::LTPspecies_Poly( const LTPspecies_Poly &right ) 
    : LTPspecies()
  {
    *this = right; //use assignment operator to do other work
  }
  
  //! Assignment operator
  LTPspecies_Poly& LTPspecies_Poly::operator=(const LTPspecies_Poly& right )
  {
    if (&right != this) {
      //LTPspecies::operator=(right);
      m_speciesName = right.m_speciesName;
      m_property    = right.m_property;
      m_model       = right.m_model;
      m_coeffs      = right.m_coeffs;
      m_thermo      = right.m_thermo;
      m_mixWeight   = right.m_mixWeight;

      m_temp    = right.m_temp;
      m_prop    = right.m_prop;
    }
    return *this; 
  }

  //! Return the value for this transport property evaluated 
  //! from the polynomial expression
  doublereal LTPspecies_Poly::getSpeciesTransProp( ) {
    
    doublereal t = m_thermo->temperature();
    if (t != m_temp) {    
      double tempN = 1.0;
      for (int i = 0; i < (int) m_coeffs.size() ; i++) {
	m_prop += m_coeffs[i] * tempN;
	tempN *= m_temp;
      }
    }
    return m_prop;
  }








}
