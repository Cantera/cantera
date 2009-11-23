/**
 *  @file LiquidTransportParams.cpp
 *  Source code for liquid mixture transport property evaluations.
 */
/*
 *  $Author: jchewso $
 *  $Date: 2009-11-22 17:34:46 -0700 (Mon, 16 Nov 2009) $
 *  $Revision: 261 $
 *
 */

#include "LiquidTransportParams.h"


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




  LiquidTranInteraction::LiquidTranInteraction( 
                           const XML_Node &compModelNode, 
		           TransportPropertyList tp_ind, 
		           thermo_t* thermo ) : 
    m_model(LTI_MODEL_NOTSET),
    m_property(tp_ind),
    m_thermo(thermo)
  {
    int nsp = thermo->nSpecies();
    Aij.resize(nsp,nsp);
    Dij.resize(nsp,nsp);
    Eij.resize(nsp,nsp);
    Sij.resize(nsp,nsp);

    std::string speciesA;
    std::string speciesB;

    int num = compModelNode.nChildren();
    for (int iChild = 0; iChild < num; iChild++) {
      XML_Node &xmlChild = compModelNode.child(iChild);
      std::string nodeName = lowercase( xmlChild.name() );
      if ( nodeName != "interaction"  )	    
	throw CanteraError("TransportFactory::getLiquidInteractionsTransportData", 
			   "expected <interaction> element and got <" + nodeName + ">" );
      speciesA = xmlChild.attrib("speciesA");
      speciesB = xmlChild.attrib("speciesB");
      int iSpecies = m_thermo->speciesIndex( speciesA );
      if ( iSpecies < 0 ) 
	throw CanteraError("TransportFactory::getLiquidInteractionsTransportData", 
			   "Unknown species " + speciesA );
      int jSpecies = m_thermo->speciesIndex( speciesB );
      if ( jSpecies < 0 ) 
	throw CanteraError("TransportFactory::getLiquidInteractionsTransportData", 
			   "Unknown species " + speciesB );
      Aij(iSpecies,jSpecies) = getFloat( xmlChild, "Aij", "toSI" );
      Aij(jSpecies,iSpecies) = Aij(iSpecies,jSpecies) ;

      Eij(iSpecies,jSpecies) = getFloat( xmlChild, "Eij", "actEnergy" );
      Eij(iSpecies,jSpecies) /= GasConstant;
      Eij(jSpecies,iSpecies) = Eij(iSpecies,jSpecies) ;

      Sij(iSpecies,jSpecies) = getFloat( xmlChild, "Sij", "toSI" );
      Sij(iSpecies,jSpecies) /= GasConstant;
      Sij(jSpecies,iSpecies) = Sij(iSpecies,jSpecies) ;

      Dij(iSpecies,jSpecies) = getFloat( xmlChild, "Dij", "toSI" );
      Dij(jSpecies,iSpecies) = Dij(iSpecies,jSpecies) ;
    }
  }

  //! Copy constructor
  LiquidTranInteraction::LiquidTranInteraction( const LiquidTranInteraction &right ) {
    *this = right;  //use assignment operator to do other work
  }
   
  //! Assignment operator
  LiquidTranInteraction& LiquidTranInteraction::operator=( const LiquidTranInteraction &right )
  {
    if (&right != this) {
      m_model       = right.m_model;
      m_property    = right.m_property;
      m_thermo    = right.m_thermo;
    }
    return *this;     
  }
  
  


  LTI_MoleFracs::LTI_MoleFracs( const XML_Node &compModelNode, 
				TransportPropertyList tp_ind, 
				thermo_t* thermo ) :
    LiquidTranInteraction( compModelNode, tp_ind, thermo )
  {
    m_model = LTI_MODEL_MOLEFRACS;
  }

  
  
  
} //namespace Cantera
