/**
 *  @file LiquidTransportParams.h
 *  Header file defining class LiquidTransportParams
 */
/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *
 *
 */
#ifndef CT_LIQUIDTRANSPORTPARAMS_H
#define CT_LIQUIDTRANSPORTPARAMS_H

#include <vector>

#include "ct_defs.h"
#include "TransportBase.h"
#include "TransportParams.h"
#include "LiquidTransportData.h"
#include "xml.h"
#include "XML_Writer.h"

namespace Cantera {


  class NotImplemented : public CanteraError {
  public:
    NotImplemented(std::string method) : CanteraError("Transport",
						      "\n\n**** Method "+method+" not implemented. ****\n"
						      "(Did you forget to specify a transport model?)\n\n") {}
  };
  

  //! Composition dependence type for liquid mixture transport properties
  /*!
   *  Types of temperature dependencies:
   *  -   0  - Mixture calculations with this property are not allowed
   *  -   1  - Use solvent (species 0) properties
   *  -   2  - Properties weighted linearly by mole fractions
   *  -   3  - Properties weighted linearly by mass fractions
   *  -   4  - Properties weighted logarithmically by mole fractions (interaction energy weighting)
   *  -   5  - Interactions given pairwise between each possible species (i.e. D_ij)
   * 
   *   \verbatim
   *    <transport model="Liquid">
   *       <viscosity>
   *          <compositionDependence model="logMoleFractions">
   *             <interaction>
   *                <speciesA> LiCl(L) </speciesA>
   *                <speciesB> KCl(L)  </speciesB>
   *                <Eij units="J/kmol"> -1.0 </Eij>
   *                <Sij units="J/kmol/K"> 1.0E-1 </Eij>
   *             </interaction>
   *          </compositionDependence>          
   *       </viscosity>
   *       <speciesDiffusivity>
   *          <compositionDependence model="pairwiseInteraction">
   *             <interaction>
   *                <speciesA> Li+ </speciesA>
   *                <speciesB> K+  </speciesB>
   *                <Dij units="m2/s"> 1.5 </Dij>
   *             </interaction>
   *             <interaction>
   *                <speciesA> K+  </speciesA>
   *                <speciesB> Cl- </speciesB>
   *                <Dij units="m2/s"> 1.0 </Dij>
   *             </interaction>
   *             <interaction>
   *                <speciesA> Li+  </speciesA>
   *                <speciesB> Cl-  </speciesB>
   *                <Dij units="m2/s"> 1.2 </Dij>
   *             </interaction>
   *          </compositionDependence>          
   *       </speciesDiffusivity>
   *       <thermalConductivity>
   *          <compositionDependence model="massFractions"/>
   *       </thermalConductivity>
   *       <hydrodynamicRadius>
   *          <compositionDependence model="none"/>
   *       </hydrodynamicRadius>
   *    </transport>     
   *   \endverbatim
   *
   */
  enum LiquidTranMixingModel {
    LTI_MODEL_NOTSET=-1,
    LTI_MODEL_NONE, 
    LTI_MODEL_SOLVENT, 
    LTI_MODEL_MOLEFRACS,
    LTI_MODEL_MASSFRACS,
    LTI_MODEL_LOG_MOLEFRACS,
    LTI_MODEL_PAIRWISE_INTERACTION,
    LTI_MODEL_STOKES_EINSTEIN
  };
  

  //! Base class to handle transport property evaluation in a mixture.
  /** 
   * In a mixture, the mixture transport properties will generally depend on
   * the contributions of each of the pure species transport properties.
   * Many composition dependencies are possible.  This class, 
   * LiquidTranInteraction, is designed to be a base class for the
   * implementation of various models for the mixing of pure species
   * transport properties.  
   *
   * There are two very broad types of transport properties to consider.
   * First, there are properties for which a mixture value can be 
   * obtained through some mixing rule.  These are obtained using the
   * method getMixTransProp().  Viscosity is typical of this.
   * Second there are properties for which a matrix of properties may 
   * exist.  This matrix of properties is obtained from the method
   * getMatrixTransProp().  Diffusion coefficients are of this type.  
   * Subclasses should implement the appropriate one or both of 
   * these methods.
   * 
   */
  class LiquidTranInteraction {
    
  public:
    //! Constructor 
    /**
     *  @param tp_ind          Index indicating transport property type (i.e. viscosity) 
     */
    LiquidTranInteraction( TransportPropertyList tp_ind = TP_UNKNOWN );
    
    //! Copy constructor
    LiquidTranInteraction( const LiquidTranInteraction &right );
    
    //! Assignment operator
    LiquidTranInteraction& operator=( const LiquidTranInteraction &right );
    
    //! destructor
    virtual ~LiquidTranInteraction() { ; }

    //! initialize LiquidTranInteraction objects with thermo and XML node
    /**
     *  @param compModelNode   \verbatim <compositionDependence> \endverbatim XML node
     *  @param thermo          Pointer to thermo object
     */
    virtual void init( const XML_Node &compModelNode = 0, 
		       thermo_t* thermo = 0 );			  

    virtual void setParameters( LiquidTransportParams& trParam ) { ; }
    
    //! Return the mixture transport property value.
    //! (Must be implemented in subclasses.)
    virtual doublereal getMixTransProp( doublereal* speciesValues, doublereal *weightSpecies = 0 ) { 
      throw NotImplemented("LiquidTranInteraction::getMixTransProp"); 
    }

    virtual doublereal getMixTransProp( std::vector<LTPspecies*> LTPptrs ) { 
      throw NotImplemented("LiquidTranInteraction::getMixTransProp"); 
    }

    virtual void getMatrixTransProp( DenseMatrix &mat, doublereal* speciesValues = 0 ) { 
      //mat = m_Dij;
      throw NotImplemented("LiquidTranInteraction::getMixTransProp"); 
    }

  protected:
    //! Model for species interaction effects 
    //! Takes enum LiquidTranMixingModel
    LiquidTranMixingModel m_model;
    
    //! enum indicating what property this is (i.e viscosity)
    TransportPropertyList m_property;
    
    //! pointer to thermo object to get current temperature
    thermo_t* m_thermo;

    //LiquidTransportParams* m_trParam;

    //! Matrix of interactions (no temperature dependence, dimensionless)
    DenseMatrix  m_Aij; 

    //! Matrix of interactions (in energy units, 1/RT temperature dependence)
    DenseMatrix  m_Eij; 

    //! Matrix of interactions (in entropy units, divided by R)
    DenseMatrix  m_Sij;         
    
    //! Matrix of interactions 
    DenseMatrix  m_Dij;         
    
  };

  /**
   * Holds transport model parameters relevant to transport in 
   * liquids for which activated jump processes limit transport
   * (giving Arrhenius type transport properties). 
   * Used by TransportFactory.
   */
  class LiquidTransportParams : public TransportParams {
    
  public:
    
    LiquidTransportParams();
    ~LiquidTransportParams();
    LiquidTransportParams(const LiquidTransportParams &right);
    LiquidTransportParams & operator=(const LiquidTransportParams &right);
    
    //! Species transport parameters
    std::vector<Cantera::LiquidTransportData> LTData;

    LiquidTranInteraction* viscosity;
    LiquidTranInteraction* thermalCond;
    LiquidTranInteraction* speciesDiffusivity;
    LiquidTranInteraction* electCond;
    LiquidTranInteraction* hydroRadius;
    
    //! Model for species interaction effects for viscosity
    //! Takes enum LiquidTranMixingModel
    LiquidTranMixingModel model_viscosity;
    
    //! Interaction associated with linear weighting of 
    //! thermal conductivity.
    /**
     * This is used for either LTI_MODEL_MASSFRACS 
     * or LTI_MODEL_MOLEFRACS.  
     * The overall formula for the mixture viscosity is 
     *
     * \f[ \eta_{mix} = \sum_i X_i \eta_i 
     *  + \sum_i \sum_j X_i X_j A_{i,j} \f].
     */
    DenseMatrix  thermalCond_Aij; 
    
    //! Model for species interaction effects for mass diffusivity
    //! Takes enum LiquidTranMixingModel
    LiquidTranMixingModel model_speciesDiffusivity;
    
    //! Interaction associated with linear weighting of 
    //! thermal conductivity.
    /**
     * This is used for either LTI_MODEL_PAIRWISE_INTERACTION.
     * These provide species interaction coefficients associated with 
     * the Stefan-Maxwell formulation.
     */
    DenseMatrix  diff_Dij; 
    
    //! Model for species interaction effects for hydrodynamic radius
    //! Takes enum LiquidTranMixingModel
    LiquidTranMixingModel model_hydroradius;
    
    //! Interaction associated with hydrodynamic radius.
    /**
     * Not yet implemented 
     */
    DenseMatrix  radius_Aij; 
  };

  
  class LTI_Solvent;
  class LTI_MoleFracs;
  class LTI_MassFracs;
  class LTI_Log_MoleFracs;
  class LTI_Pairwise_Interaction;

  class LTI_Solvent : public LiquidTranInteraction {

  public:
    LTI_Solvent( TransportPropertyList tp_ind = TP_UNKNOWN ) :
      LiquidTranInteraction( tp_ind )
    {
      m_model = LTI_MODEL_SOLVENT;
    }
    
    //! Copy constructor
    //    LTI_Solvent( const LTI_Solvent &right );
    
    //! Assignment operator
    //    LTI_Solvent& operator=( const LTI_Solvent &right );
    
    virtual ~LTI_Solvent( ) { } 
    
    //! Return the mixture transport property value.
    /** 
     * Takes the separate species transport properties 
     * as input (this method does not know what
     * transport property it is at this point).
     */
    doublereal getMixTransProp( doublereal *valueSpecies, doublereal *weightSpecies = 0 );
    doublereal getMixTransProp( std::vector<LTPspecies*> LTPptrs ) ;

    void getMatrixTransProp( DenseMatrix &mat, doublereal* speciesValues = 0 ) { mat = m_Aij; }

  protected:    
    
  };

  //! Simple mole fraction weighting of transport properties
  /**
   * This model weights the transport property by  the mole 
   * fractions.
   * The overall formula for the mixture viscosity is 
   *
   * \f[ \eta_{mix} = \sum_i X_i \eta_i 
   *  + \sum_i \sum_j X_i X_j A_{i,j} \f].
   */
  class LTI_MoleFracs : public LiquidTranInteraction {

  public:
    LTI_MoleFracs( TransportPropertyList tp_ind = TP_UNKNOWN ) :
      LiquidTranInteraction( tp_ind )
    {
      m_model = LTI_MODEL_MOLEFRACS;
    }

    
    //! Copy constructor
    //    LTI_MoleFracs( const LTI_MoleFracs &right );
    
    //! Assignment operator
    //    LTI_MoleFracs& operator=( const LTI_MoleFracs &right );
    
    virtual ~LTI_MoleFracs( ) { } 
    
    //! Return the mixture transport property value.
    /** 
     * Takes the separate species transport properties 
     * as input (this method does not know what
     * transport property it is at this point.
     */
    doublereal getMixTransProp( doublereal *valueSpecies, doublereal *weightSpecies = 0 );
    doublereal getMixTransProp( std::vector<LTPspecies*> LTPptrs ) ;

    void getMatrixTransProp( DenseMatrix &mat, doublereal* speciesValues = 0 ) { mat = m_Aij; }

  protected:    
    
  };


  //! Simple mass fraction weighting of transport properties
  /**
   * This model weights the transport property by  the mass 
   * fractions.
   * The overall formula for the mixture viscosity is 
   *
   * \f[ \eta_{mix} = \sum_i Y_i \eta_i 
   *  + \sum_i \sum_j Y_i Y_j A_{i,j} \f].
   */
  class LTI_MassFracs : public LiquidTranInteraction {

  public:
    LTI_MassFracs( TransportPropertyList tp_ind = TP_UNKNOWN ) :
      LiquidTranInteraction( tp_ind )
    {
      m_model = LTI_MODEL_MASSFRACS;
    }

    
    //! Copy constructor
    //    LTI_MassFracs( const LTI_MassFracs &right );
    
    //! Assignment operator
    //    LTI_MassFracs& operator=( const LTI_MassFracs &right );
    
    virtual ~LTI_MassFracs( ) { } 
    
    //! Return the mixture transport property value.
    /** 
     * Takes the separate species transport properties 
     * as input (this method does not know what
     * transport property it is at this point.
     */
    doublereal getMixTransProp( doublereal *valueSpecies, doublereal *weightSpecies = 0 );
    doublereal getMixTransProp( std::vector<LTPspecies*> LTPptrs ) ;

    void getMatrixTransProp( DenseMatrix &mat, doublereal* speciesValues = 0 ) { mat = m_Aij; }

  protected:    
    
  };


  //! Mixing rule using logarithms of the mole fractions
  /** 
   * This model is based on the idea that liquid molecules are 
   * generally interacting with some energy and entropy of interaction.
   * For transport properties that depend on these energies of 
   * interaction, the mixture transport property can be written
   * in terms of its logarithm
   * 
   * \f[ \ln \eta_{mix} = \sum_i X_i \ln \eta_i 
   *  + \sum_i \sum_j X_i X_j ( S_{i,j} + E_{i,j} / T ) 
   * \f].
   * 
   * These additional interaction terms multiply the mixture property by
   *  \f[ \exp( \sum_{i} \sum_{j} X_i X_j ( S_{i,j} + E_{i,j} / T ) ) \f]
   * so that the self-interaction terms \f$ S_{i,j} \f$ and 
   * \f$ E_{i,j} \f$ should be zero.
   *
   * Note that the energies and entropies of interaction should be
   * a function of the composition themselves, but this is not yet 
   * implemented.  (We might follow the input of Margules model 
   * thermodynamic data for the purpose of implementing this.)
   * 
   *  Sample input for this method is
   *  \verbatim
   *    <transport model="Liquid">
   *       <viscosity>
   *          <compositionDependence model="logMoleFractions">
   *             <interaction speciesA="Li+" speciesB="K+">
   *	         <!-- 
   *	           interactions are from speciesA = LiCl(L)
   *		   and speciesB = KCl(L).
   *                   -->
   *                <Eij units="J/kmol"> -1.0e3 </Eij>
   *                <Sij units="J/kmol/K"> 80.0e-5 </Sij>
   *             </interaction>
   *          </compositionDependence>          
   *       </viscosity>
   *    </transport>
   *  \endverbatim
   */
  class LTI_Log_MoleFracs : public LiquidTranInteraction {

  public:
    LTI_Log_MoleFracs( TransportPropertyList tp_ind = TP_UNKNOWN ) :
      LiquidTranInteraction( tp_ind )
    {
      m_model = LTI_MODEL_LOG_MOLEFRACS;
    }

    
    //! Copy constructor
    //    LTI_Log_MoleFracs( const LTI_Log_MoleFracs &right );
    
    //! Assignment operator
    //    LTI_Log_MoleFracs& operator=( const LTI_Log_MoleFracs &right );
    
    virtual ~LTI_Log_MoleFracs( ) { } 
    
    //! Return the mixture transport property value.
    /** 
     * Takes the separate species transport properties 
     * as input (this method does not know what
     * transport property it is at this point.
     */
    doublereal getMixTransProp( doublereal *valueSpecies, doublereal *weightSpecies = 0 );
    doublereal getMixTransProp( std::vector<LTPspecies*> LTPptrs ) ;

    void getMatrixTransProp( DenseMatrix &mat, doublereal* speciesValues = 0 ) { mat = m_Eij; }

  protected:    
    
  };


  //! Transport properties that act like pairwise interactions
  //! as in binary diffusion coefficients.
  /**
   * This class holds parameters for transport properties expressed 
   * as a matrix of pairwise interaction parameters.
   * Input can be provided for constant or Arrhenius forms of the
   * separate parameters.
   *
   *  Sample input for this method is
   *  \verbatim
   *    <transport model="Liquid">
   *       <speciesDiffusivity>
   *          <compositionDependence model="pairwiseInteraction">
   *             <interaction speciesA="LiCl(L)" speciesB="KCl(L)">
   *                <Dij units="m/s"> 1.0e-8 </Dij>
   *                <Eij units="J/kmol"> 24.0e6 </Eij>
   *             </interaction>
   *          </compositionDependence>          
   *       </speciesDiffusivity>
   *    </transport>     
   *  \endverbatim
   *
   */
  class LTI_Pairwise_Interaction : public LiquidTranInteraction {

  public:
    LTI_Pairwise_Interaction( TransportPropertyList tp_ind = TP_UNKNOWN ) :
      LiquidTranInteraction( tp_ind )
    {
      m_model = LTI_MODEL_PAIRWISE_INTERACTION;
    }
    
    
    //! Copy constructor
    //    LTI_Pairwise_Interaction( const LTI_Pairwise_Interaction &right );
    
    //! Assignment operator
    //    LTI_Pairwise_Interaction& operator=( const LTI_Pairwise_Interaction &right );
    
    virtual ~LTI_Pairwise_Interaction( ) { } 
    
    void setParameters( LiquidTransportParams& trParam ) ;

    //! Return the mixture transport property value.
    /** 
     * Takes the separate species transport properties 
     * as input (this method does not know what
     * transport property it is at this point.
     */
    doublereal getMixTransProp( doublereal *valueSpecies, doublereal *weightSpecies = 0 );
    doublereal getMixTransProp( std::vector<LTPspecies*> LTPptrs ) ;

    void getMatrixTransProp( DenseMatrix &mat, doublereal* speciesValues = 0 ) ;

  protected:    

    std::vector<LTPspecies*> m_diagonals;
  };


  class LTI_StokesEinstein : public LiquidTranInteraction {

  public:
    LTI_StokesEinstein( TransportPropertyList tp_ind = TP_UNKNOWN ) :
      LiquidTranInteraction( tp_ind )
    {
      m_model = LTI_MODEL_STOKES_EINSTEIN;
    }
    
    
    //! Copy constructor
    //    LTI_StokesEinstein( const LTI_StokesEinstein &right );
    
    //! Assignment operator
    //    LTI_StokesEinstein& operator=( const LTI_StokesEinstein &right );
    
    virtual ~LTI_StokesEinstein( ) { } 

    void setParameters( LiquidTransportParams& trParam );
  
    //! Return the mixture transport property value.
    /** 
     * Takes the separate species transport properties 
     * as input (this method does not know what
     * transport property it is at this point.
     */
    doublereal getMixTransProp( doublereal *valueSpecies, doublereal *weightSpecies = 0 );
    doublereal getMixTransProp( std::vector<LTPspecies*> LTPptrs ) ;

    void getMatrixTransProp( DenseMatrix &mat, doublereal* speciesValues = 0 ) ;

  protected:    
    
    std::vector<LTPspecies*> m_viscosity;
    std::vector<LTPspecies*> m_hydroRadius;
    
  };


}

#endif
