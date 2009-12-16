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

  /** 
   * Exception thrown if an error is encountered while reading the 
   * transport database.
   */
  class LTPmodelError : public CanteraError {
  public:
    LTPmodelError( std::string msg ) 
      : CanteraError("LTPspecies",
		     "error parsing transport data: " 
		     + msg + "\n") {}
  };


  //!Constructor
  /**
   *  @param tp_ind          Index indicating transport property type (i.e. viscosity) 
   */
  LiquidTranInteraction::LiquidTranInteraction( TransportPropertyList tp_ind ) : 
    m_model(LTI_MODEL_NOTSET),
    m_property(tp_ind)
  {
  }




  void LiquidTranInteraction::init( const XML_Node &compModelNode, 
			       thermo_t* thermo ) 
  {
    m_thermo = thermo;

    int nsp = thermo->nSpecies();
    m_Aij.resize(nsp,nsp);
    m_Dij.resize(nsp,nsp);
    m_Eij.resize(nsp,nsp);
    m_Sij.resize(nsp,nsp);

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
      if ( xmlChild.hasChild( "Aij" ) ) {
	m_Aij(iSpecies,jSpecies) = getFloat( xmlChild, "Aij", "toSI" );
	m_Aij(jSpecies,iSpecies) = m_Aij(iSpecies,jSpecies) ;
      }

      if ( xmlChild.hasChild( "Eij" ) ) {
	m_Eij(iSpecies,jSpecies) = getFloat( xmlChild, "Eij", "actEnergy" );
	m_Eij(iSpecies,jSpecies) /= GasConstant;
	m_Eij(jSpecies,iSpecies) = m_Eij(iSpecies,jSpecies) ;
      }

      if ( xmlChild.hasChild( "Sij" ) ) {
	m_Sij(iSpecies,jSpecies) = getFloat( xmlChild, "Sij", "toSI" );
	m_Sij(iSpecies,jSpecies) /= GasConstant;
	m_Sij(jSpecies,iSpecies) = m_Sij(iSpecies,jSpecies) ;
      }

      if ( xmlChild.hasChild( "Dij" ) ) {
	m_Dij(iSpecies,jSpecies) = getFloat( xmlChild, "Dij", "toSI" );
	m_Dij(jSpecies,iSpecies) = m_Dij(iSpecies,jSpecies) ;
      }
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
      m_model     = right.m_model;
      m_property  = right.m_property;
      m_thermo    = right.m_thermo;
      //m_trParam   = right.m_trParam;
      m_Aij       = right.m_Aij;
      m_Eij       = right.m_Eij;
      m_Sij       = right.m_Sij;
      m_Dij       = right.m_Dij;
    }
    return *this;     
  }
  
  


  doublereal LTI_Solvent::getMixTransProp( doublereal *speciesValues, doublereal *speciesWeight ) {

    int nsp = m_thermo->nSpecies();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    doublereal value = 0;    
    
    //if weightings are specified, use those
    if ( speciesWeight ) {  
      for ( int k = 0; k < nsp; k++) { 
	//presume that the weighting is set to 1.0 for solvent and 0.0 for everything else.
	value += speciesValues[k] * speciesWeight[k];
      }
    }
    else {
      //This does not follow directly a solvent model 
      //although if the solvent mole fraction is dominant 
      //and the other species values are given or zero, 
      //it should work.
      for ( int k = 0; k < nsp; k++) { 
	value += speciesValues[k] * molefracs[k];
      }
    }

    for ( int i = 0; i < nsp; i++ ) 
      for ( int j = 0; j < i; j++ ) 
	value += molefracs[i] * molefracs[j] * m_Aij(i,j) ;

    return value;
  }


  doublereal LTI_Solvent::getMixTransProp( std::vector<LTPspecies*> LTPptrs ) {

    int nsp = m_thermo->nSpecies();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    doublereal value = 0;    

    for ( int k = 0; k < nsp; k++) { 
      //presume that the weighting is set to 1.0 for solvent and 0.0 for everything else.
      value += LTPptrs[k]->getSpeciesTransProp() * LTPptrs[k]->getMixWeight( ) ;
    }

    for ( int i = 0; i < nsp; i++ ) 
      for ( int j = 0; j < i; j++ ) 
	value += molefracs[i] * molefracs[j] * m_Aij(i,j) ;

    return value;
  }


  


  doublereal LTI_MoleFracs::getMixTransProp( doublereal *speciesValues, doublereal *speciesWeight ) {

    int nsp = m_thermo->nSpecies();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    doublereal value = 0;    
    
    //if weightings are specified, use those
    if ( speciesWeight ) {  
      for ( int k = 0; k < nsp; k++) { 
	value += speciesValues[k] * speciesWeight[k] * molefracs[k];
      }
    }
    else {
      for ( int k = 0; k < nsp; k++) { 
	value += speciesValues[k] * molefracs[k];
      }
    }

    for ( int i = 0; i < nsp; i++ ) 
      for ( int j = 0; j < i; j++ ) 
	value += molefracs[i] * molefracs[j] * m_Aij(i,j) ;

    return value;
  }


  doublereal LTI_MoleFracs::getMixTransProp( std::vector<LTPspecies*> LTPptrs ) {

    int nsp = m_thermo->nSpecies();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    doublereal value = 0;    

    for ( int k = 0; k < nsp; k++) { 
      value += LTPptrs[k]->getSpeciesTransProp() * LTPptrs[k]->getMixWeight( ) * molefracs[k];
    }

    for ( int i = 0; i < nsp; i++ ) 
      for ( int j = 0; j < i; j++ ) 
	value += molefracs[i] * molefracs[j] * m_Aij(i,j) ;

    return value;
  }


  


  doublereal LTI_MassFracs::getMixTransProp( doublereal *speciesValues, doublereal *speciesWeight ) {

    int nsp = m_thermo->nSpecies();
    doublereal massfracs[nsp];
    m_thermo->getMassFractions(massfracs);

    doublereal value = 0;    
    
    //if weightings are specified, use those
    if ( speciesWeight ) {  
      for ( int k = 0; k < nsp; k++) { 
	value += speciesValues[k] * speciesWeight[k] * massfracs[k];
      }
    }
    else {
      for ( int k = 0; k < nsp; k++) { 
	value += speciesValues[k] * massfracs[k];
      }
    }

    for ( int i = 0; i < nsp; i++ ) 
      for ( int j = 0; j < i; j++ ) 
	value += massfracs[i] * massfracs[j] * m_Aij(i,j) ;

    return value;
  }


  doublereal LTI_MassFracs::getMixTransProp( std::vector<LTPspecies*> LTPptrs ) {

    int nsp = m_thermo->nSpecies();
    doublereal massfracs[nsp];
    m_thermo->getMassFractions(massfracs);

    doublereal value = 0;    

    for ( int k = 0; k < nsp; k++) { 
      value += LTPptrs[k]->getSpeciesTransProp() * LTPptrs[k]->getMixWeight( ) * massfracs[k];
    }

    for ( int i = 0; i < nsp; i++ ) 
      for ( int j = 0; j < i; j++ ) 
	value += massfracs[i] * massfracs[j] * m_Aij(i,j) ;

    return value;
  }


  


  doublereal LTI_Log_MoleFracs::getMixTransProp( doublereal *speciesValues, doublereal *speciesWeight ) {

    int nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    doublereal value = 0;    
    
    //if weightings are specified, use those
    if ( speciesWeight ) {  
      for ( int k = 0; k < nsp; k++) { 
	value += log( speciesValues[k] ) * speciesWeight[k] * molefracs[k];
      }
    }
    else {
      for ( int k = 0; k < nsp; k++) { 
	value += log( speciesValues[k] ) * molefracs[k];
      }
    }

    for ( int i = 0; i < nsp; i++ ) 
      for ( int j = 0; j < i; j++ ) 
	value += molefracs[i] * molefracs[j] 
	  * ( m_Sij(i,j) + m_Eij(i,j) / temp );

    value = exp( value );
    return value;
  }


  doublereal LTI_Log_MoleFracs::getMixTransProp( std::vector<LTPspecies*> LTPptrs ) {

    int nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    doublereal value = 0;    

    for ( int k = 0; k < nsp; k++) { 
      value += log( LTPptrs[k]->getSpeciesTransProp() ) * LTPptrs[k]->getMixWeight( ) * molefracs[k];
    }

    for ( int i = 0; i < nsp; i++ ) 
      for ( int j = 0; j < i; j++ ) 
	value += molefracs[i] * molefracs[j] 
	  * ( m_Sij(i,j) + m_Eij(i,j) / temp );

    value = exp( value );

    return value;
  }


  


  void LTI_Pairwise_Interaction::setParameters( LiquidTransportParams& trParam ) {
    int nsp = m_thermo->nSpecies();
    m_diagonals.resize(nsp);

    for (int k = 0; k < nsp; k++) {
      Cantera::LiquidTransportData &ltd = trParam.LTData[k];
      if ( ltd.speciesDiffusivity )
	m_diagonals[k] = ltd.speciesDiffusivity;
    }
  }

  doublereal LTI_Pairwise_Interaction::getMixTransProp( doublereal *speciesValues, doublereal *speciesWeight ) {

    int nsp = m_thermo->nSpecies();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    doublereal value = 0;    
    
    throw LTPmodelError( "Calling LTI_Pairwise_Interaction::getMixTransProp does not make sense." );

    return value;
  }


  doublereal LTI_Pairwise_Interaction::getMixTransProp( std::vector<LTPspecies*> LTPptrs ) {

    int nsp = m_thermo->nSpecies();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    doublereal value = 0;    
    
    throw LTPmodelError( "Calling LTI_Pairwise_Interaction::getMixTransProp does not make sense." );

    return value;
  }

  DenseMatrix LTI_Pairwise_Interaction::getMatrixTransProp( doublereal *speciesValues ) {

    int nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    DenseMatrix tmp;
    tmp.resize(nsp,nsp);
    for ( int i = 0; i < nsp; i++ ) 
      for ( int j = 0; j < i; j++ ) 
	tmp(i,j) = tmp(j,i) = m_Dij(i,j) * exp( - m_Eij(i,j) / temp );
    
    for ( int i = 0; i < nsp; i++ ) 
      if ( tmp(i,i) == 0.0 && m_diagonals[i] ) 
	tmp(i,i) = m_diagonals[i]->getSpeciesTransProp() ;

    return tmp;
  }

  
  doublereal LTI_StokesEinstein::getMixTransProp( doublereal *speciesValues, doublereal *speciesWeight ) {

    int nsp = m_thermo->nSpecies();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    doublereal value = 0;    
    
    throw LTPmodelError( "Calling LTI_StokesEinstein::getMixTransProp does not make sense." );

    return value;
  }


  doublereal LTI_StokesEinstein::getMixTransProp( std::vector<LTPspecies*> LTPptrs ) {

    int nsp = m_thermo->nSpecies();
    doublereal molefracs[nsp];
    m_thermo->getMoleFractions(molefracs);

    doublereal value = 0;    
    
    throw LTPmodelError( "Calling LTI_StokesEinstein::getMixTransProp does not make sense." );

    return value;
  }



  void LTI_StokesEinstein::setParameters( LiquidTransportParams& trParam ) {
    int nsp = m_thermo->nSpecies();
    m_viscosity.resize(nsp);
    m_hydroRadius.resize(nsp);
    for (int k = 0; k < nsp; k++) {
      Cantera::LiquidTransportData &ltd = trParam.LTData[k];
      m_viscosity[k]   =  ltd.viscosity;
      m_hydroRadius[k] =  ltd.hydroRadius;
    }
  }

  DenseMatrix LTI_StokesEinstein::getMatrixTransProp( doublereal *speciesValues ) {

    int nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();

    double *viscSpec = new double(nsp);
    double *radiusSpec = new double(nsp);

    for (int k = 0; k < nsp; k++) {
      viscSpec[k] = m_viscosity[k]->getSpeciesTransProp() ;
      radiusSpec[k] = m_hydroRadius[k]->getSpeciesTransProp() ;
    } 

    DenseMatrix tmp;
    tmp.resize(nsp,nsp);
    for (int i = 0; i < nsp; i++) 
      for (int j = 0; j < nsp; j++) {
	tmp(i,j) = GasConstant * temp 
	  / ( 6.0 * Pi * radiusSpec[i] * viscSpec[j] ) ;
      }
    delete radiusSpec;
    delete viscSpec;

    return tmp;
  }



  
} //namespace Cantera
