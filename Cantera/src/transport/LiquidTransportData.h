/**
 *  @file LiquidTransportData.h
 *  Header file defining class LiquidTransportData
 */
/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *
 *
 */

#ifndef CT_LIQUIDTRANSPORTDATA_H
#define CT_LIQUIDTRANSPORTDATA_H


// STL includes
#include <vector>
#include <string>



// Cantera includes
#include "ct_defs.h"
#include "TransportBase.h"
#include "FactoryBase.h"


namespace Cantera {

  enum TransportPropertyList {
    TP_UNKNOWN = -1,
    TP_VISCOSITY = 0,
    TP_THERMALCOND,
    TP_DIFFUSIVITY,
    TP_HYDRORADIUS,
    TP_ELECTCOND
  };

  enum LiquidTR_Model {
    //! Temperature dependence type for pure (liquid) species properties
    /*!
     *  Types of temperature dependencies:
     *     0  - Independent of temperature (only one implemented so far)
     *     1  - extended arrhenius form
     *     2  - polynomial in temperature form
     */
    LTR_MODEL_NOTSET=-1,
    LTR_MODEL_CONSTANT, 
    LTR_MODEL_ARRHENIUS,
    LTR_MODEL_POLY
  };


  //! Class LTPspecies holds transport parameters for a 
  //! specific liquid-phase species.
  /** 
   * Subclasses handle different means of specifying transport properties
   * like constant, Arrhenius or polynomial fits.  In its current state, 
   * it is primarily suitable for specifying temperature dependence, but 
   * the adjustCoeffsForComposition() method can be implemented to 
   * adjust for composition dependence.  
   * Mixing rules for computing mixture transport properties are handled 
   * separately in 
   */
  class LTPspecies {

  public:

    //! Construct an LTPspecies object for a liquid tranport property.
    /** 
     *  The transport property is constructed from the 
     *  XML node, propNode, that is a child of the 
     *  <transport> node and specifies a type of
     *  transport property (like viscosity).  
     */ 
    LTPspecies( const XML_Node &propNode = 0, 
		std::string name = "-", 
		TransportPropertyList tp_ind = TP_UNKNOWN, 
		thermo_t* thermo = 0 ) :
      m_speciesName(name), 
      m_model(LTR_MODEL_NOTSET),
      m_property(tp_ind),
      m_thermo(thermo),
      m_mixWeight(1.0)
    {
      if ( propNode.hasChild("mixtureWeighting") ) 
	m_mixWeight = getFloat(propNode,"mixtureWeighting");
    }
    
    //! Copy constructor
    LTPspecies( const LTPspecies &right ); 

    //! Assignment operator
    LTPspecies&  operator=(const LTPspecies& right );

    virtual ~LTPspecies( ) { }

    //! Returns the vector of pure species tranport property
    /*!
     *  The pure species transport property (i.e. pure species viscosity)
     *  is returned.  Any temperature and composition dependence will be 
     *  adjusted internally according to the information provided by the 
     *  thermo object. 
     */
    virtual doublereal getSpeciesTransProp( ) { return 0.0; }

    virtual bool checkPositive( ) { return ( m_coeffs[0] > 0 ); }

    doublereal getMixWeight( ) {return m_mixWeight; }

  protected:
    std::string m_speciesName;
   
    //! Model type for the temperature dependence
    LiquidTR_Model m_model;

    //! enum indicating what property this is (i.e viscosity)
    TransportPropertyList m_property;

    //! Model temperature-dependence ceofficients
    vector_fp m_coeffs;

    //! pointer to thermo object to get current temperature
    thermo_t* m_thermo;

    //! Weighting used for mixing.  
    /** 
     * This weighting can be employed to allow salt transport 
     * properties to be represented by specific ions.  
     * For example, to have Li+ and Ca+ represent the mixing 
     * transport properties of LiCl and CaCl2, the weightings for
     * Li+ would be 2.0, for K+ would be 3.0 and for Cl- would be 0.0.
     * The tranport properties for Li+ would be those for LiCl and 
     * the tranport properties for Ca+ would be those for CaCl2. 
     * The transport properties for Cl- should be something innoccuous like 
     * 1.0--note that 0.0 is not innocuous if there are logarithms involved.
     */
    doublereal m_mixWeight;

    //! Internal model to adjust species-specific properties for composition.
    /** Currently just a place holder, but this method could take 
     * the composition from the thermo object and adjust coefficients 
     * accoding to some unspecified model.
     */
    virtual void adjustCoeffsForComposition() { }
  };



  //! Class LiquidTransportData holds transport parameters for a 
  //! specific liquid-phase species.
  class LiquidTransportData {

  public:

    LiquidTransportData() : 
      speciesName("-")
    {
    }
    //! copy constructor
    LiquidTransportData( const LiquidTransportData &right ) ;

    //! Assignment operator
    LiquidTransportData& operator=(const LiquidTransportData& right ); 

    std::string speciesName;   
    
    //! Model type for the hydroradius
    LTPspecies* hydroradius;

    //! Model type for the viscosity
    LTPspecies* viscosity;

    //! Model type for the thermal conductivity
    LTPspecies* thermalCond;
   
    //! Model type for the electrical conductivity
    LTPspecies* electCond;
   
    //! Model type for the speciesDiffusivity
    LTPspecies* speciesDiffusivity;
  };




  //! Class LTPspecies_Const holds transport parameters for a 
  //! specific liquid-phase species when the transport property 
  //! is just a constant value.  
  class LTPspecies_Const : public  LTPspecies{

  public:

    LTPspecies_Const( const XML_Node &propNode, 
		      std::string name, 
		      TransportPropertyList tp_ind, 
		      thermo_t* thermo ) ;
    
    //! Copy constructor
    LTPspecies_Const( const LTPspecies_Const &right ); 

    //! Assignment operator
    LTPspecies_Const&  operator=(const LTPspecies_Const& right );

    virtual ~LTPspecies_Const( ) { }

    //! Returns the pure species tranport property
    /*!
     *  The pure species transport property (i.e. pure species viscosity)
     *  is returned.  Any temperature and composition dependence will be 
     *  adjusted internally according to the information provided by the 
     *  thermo object. 
     */
    doublereal getSpeciesTransProp( );

  protected:

    //! Internal model to adjust species-specific properties for composition.
    /** Currently just a place holder, but this method could take 
     * the composition from the thermo object and adjust coefficients 
     * accoding to some unspecified model.
     */
    void adjustCoeffsForComposition( ) { }
  };



  //! Class LTPspecies_Arrhenius holds transport parameters for a 
  //! specific liquid-phase species when the transport property 
  //! is just a constant value.  
  class LTPspecies_Arrhenius : public  LTPspecies{

  public:

    LTPspecies_Arrhenius( const XML_Node &propNode, 
			  std::string name, 
			  TransportPropertyList tp_ind, 
			  thermo_t* thermo ); 
    
    //! Copy constructor
    LTPspecies_Arrhenius( const LTPspecies_Arrhenius &right ); 

    //! Assignment operator
    LTPspecies_Arrhenius&  operator=(const LTPspecies_Arrhenius& right );

    virtual ~LTPspecies_Arrhenius( ) { }

    //! Returns the pure species tranport property
    /*!
     *  The pure species transport property (i.e. pure species viscosity)
     *  is returned.  Any temperature and composition dependence will be 
     *  adjusted internally according to the information provided by the 
     *  thermo object. 
     */
    doublereal getSpeciesTransProp( );

  protected:

    //! temperature from thermo object
    doublereal m_temp;

    //! logarithm of current temperature
    doublereal m_logt;

    //! most recent evaluation of transport property
    doublereal m_prop;

    //! logarithm of most recent evaluation of transport property
    doublereal m_logProp;

    //! Internal model to adjust species-specific properties for composition.
    /** Currently just a place holder, but this method could take 
     * the composition from the thermo object and adjust coefficients 
     * accoding to some unspecified model.
     */
    void adjustCoeffsForComposition( ) { }
  };



  //! Class LTPspecies_Poly holds transport parameters for a 
  //! specific liquid-phase species when the transport property 
  //! is just a constant value.  
  class LTPspecies_Poly : public  LTPspecies{

  public:

    LTPspecies_Poly( const XML_Node &propNode, 
		     std::string name, 
		     TransportPropertyList tp_ind, 
		     thermo_t* thermo ); 
    
    //! Copy constructor
    LTPspecies_Poly( const LTPspecies_Poly &right ); 

    //! Assignment operator
    LTPspecies_Poly&  operator=(const LTPspecies_Poly& right );

    virtual ~LTPspecies_Poly( ) { }

    //! Returns the pure species tranport property
    /*!
     *  The pure species transport property (i.e. pure species viscosity)
     *  is returned.  Any temperature and composition dependence will be 
     *  adjusted internally according to the information provided by the 
     *  thermo object. 
     */
    doublereal getSpeciesTransProp( );

  protected:

    //! temperature from thermo object
    doublereal m_temp;

    //! most recent evaluation of transport property
    doublereal m_prop;

    //! Internal model to adjust species-specific properties for composition.
    /** Currently just a place holder, but this method could take 
     * the composition from the thermo object and adjust coefficients 
     * accoding to some unspecified model.
     */
    void adjustCoeffsForComposition( ){ }
  };



}
#endif
