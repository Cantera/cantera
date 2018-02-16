/**
 *  @file LiquidTranInteraction.h
 *  Header file defining the class LiquidTranInteraction and classes which
 *  derive from LiquidTranInteraction.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_LIQUIDTRANINTERACTION_H
#define CT_LIQUIDTRANINTERACTION_H

#include "TransportParams.h"
#include "LiquidTransportData.h"
#include "cantera/base/xml.h"

namespace Cantera
{
//! Composition dependence type for liquid mixture transport properties
/*!
 * @deprecated To be removed after Cantera 2.4
 *
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
 *                <Sij units="J/kmol/K"> 1.0E-1 </Sij>
 *     -or-       <Sij>
 *                  <floatArray units="J/kmol/K"> 1.0E-1, 0.001 0.01 </floatArray>
 *                </Sij>
 *     -same form for Hij,Aij,Bij-
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
 */
enum LiquidTranMixingModel {
    LTI_MODEL_NOTSET=-1,
    LTI_MODEL_SOLVENT,
    LTI_MODEL_MOLEFRACS,
    LTI_MODEL_MASSFRACS,
    LTI_MODEL_LOG_MOLEFRACS,
    LTI_MODEL_PAIRWISE_INTERACTION,
    LTI_MODEL_STEFANMAXWELL_PPN,
    LTI_MODEL_STOKES_EINSTEIN,
    LTI_MODEL_MOLEFRACS_EXPT,
    LTI_MODEL_NONE,
    LTI_MODEL_MULTIPLE
};

//! Base class to handle transport property evaluation in a mixture.
/*!
 * @deprecated To be removed after Cantera 2.4
 *
 * In a mixture, the mixture transport properties will generally depend on the
 * contributions of each of the standard state species transport properties.
 * Many composition dependencies are possible.  This class,
 * LiquidTranInteraction, is designed to be a base class for the implementation
 * of various models for the mixing of standard state species transport
 * properties.
 *
 * There are two very broad types of transport properties to consider. First,
 * there are properties for which a mixture value can be obtained through some
 * mixing rule.  These are obtained using the method getMixTransProp().
 * Viscosity is typical of this. Second, there are properties for which a matrix
 * of properties may exist.  This matrix of properties is obtained from the
 * method getMatrixTransProp().  Diffusion coefficients are of this type.
 * Subclasses should implement the appropriate one or both of these methods.
 */
class LiquidTranInteraction
{
public:
    //! Constructor
    /**
     *  @param tp_ind  Index indicating the transport property type (e.g., viscosity)
     */
    LiquidTranInteraction(TransportPropertyType tp_ind = TP_UNKNOWN);

    virtual ~LiquidTranInteraction();

    //! initialize LiquidTranInteraction objects with thermo and XML node
    /**
     *  @param compModelNode   `<compositionDependence>` XML node
     *  @param thermo          Pointer to thermo object
     */
    virtual void init(const XML_Node& compModelNode = XML_Node(),
                      thermo_t* thermo = 0);

    virtual void setParameters(LiquidTransportParams& trParam) {}

    //! Return the mixture transport property value.
    //! (Must be implemented in subclasses.)
    virtual doublereal getMixTransProp(doublereal* speciesValues, doublereal* weightSpecies = 0) {
        throw NotImplementedError("LiquidTranInteraction::getMixTransProp");
    }

    virtual doublereal getMixTransProp(std::vector<LTPspecies*> LTPptrs) {
        throw NotImplementedError("LiquidTranInteraction::getMixTransProp");
    }

    virtual void getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues = 0) {
        throw NotImplementedError("LiquidTranInteraction::getMixTransProp");
    }

protected:
    //! Model for species interaction effects. Takes enum LiquidTranMixingModel
    LiquidTranMixingModel m_model;

    //! enum indicating what property this is (i.e viscosity)
    TransportPropertyType m_property;

    //! pointer to thermo object to get current temperature
    thermo_t* m_thermo;

    //! Matrix of interaction coefficients for polynomial in molefraction*weight
    //! of speciesA (no temperature dependence, dimensionless)
    std::vector<DenseMatrix*> m_Aij;

    //! Matrix of interaction coefficients for polynomial in molefraction*weight
    //!  of speciesA (linear temperature dependence, units 1/K)
    std::vector<DenseMatrix*> m_Bij;

    //! Matrix of interactions (in energy units, 1/RT temperature dependence)
    DenseMatrix m_Eij;

    //! Matrix of interaction coefficients for polynomial in molefraction*weight
    //! of speciesA (in energy units, 1/RT temperature dependence)
    std::vector<DenseMatrix*> m_Hij;

    //! Matrix of interaction coefficients for polynomial in molefraction*weight
    //! of speciesA (in entropy units, divided by R)
    std::vector<DenseMatrix*> m_Sij;

    //! Matrix of interactions
    DenseMatrix m_Dij;
};

class LTI_Solvent : public LiquidTranInteraction
{
public:
    LTI_Solvent(TransportPropertyType tp_ind = TP_UNKNOWN);

    //! Return the mixture transport property value.
    /**
     * Takes the separate species transport properties as input (this method
     * does not know what transport property it is at this point).
     */
    doublereal getMixTransProp(doublereal* valueSpecies, doublereal* weightSpecies = 0);
    doublereal getMixTransProp(std::vector<LTPspecies*> LTPptrs);

    //! Return the matrix of binary interaction parameters.
    /**
     * Takes the proper mixing rule for the binary interaction parameters
     * and calculates them: Not implemented for this mixing rule.
     */
    void getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues = 0);
};

//! Simple mole fraction weighting of transport properties
/**
 * This model weights the transport property by the mole fractions. The
 * overall formula for the mixture viscosity is
 *
 * \f[
 *       \eta_{mix} = \sum_i X_i \eta_i  + \sum_i \sum_j X_i X_j A_{i,j}
 * \f]
 */
class LTI_MoleFracs : public LiquidTranInteraction
{
public:
    LTI_MoleFracs(TransportPropertyType tp_ind = TP_UNKNOWN) :
        LiquidTranInteraction(tp_ind) {
        m_model = LTI_MODEL_MOLEFRACS;
    }

    //! Return the mixture transport property value.
    /**
     * Takes the separate species transport properties as input (this method
     * does not know what transport property it is at this point.
     */
    doublereal getMixTransProp(doublereal* valueSpecies, doublereal* weightSpecies = 0);
    doublereal getMixTransProp(std::vector<LTPspecies*> LTPptrs);

    //! Return the matrix of binary interaction parameters.
    /**
     * Takes the proper mixing rule for the binary interaction parameters
     * and calculates them: Not Implemented for this Mixing rule;
     */
    void getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues = 0) {
        mat = (*m_Aij[0]);
    }
};

//! Simple mass fraction weighting of transport properties
/*!
 * This model weights the transport property by the mass fractions. The
 * overall formula for the mixture viscosity is
 *
 * \f[
 *   \eta_{mix} = \sum_i Y_i \eta_i
 *              + \sum_i \sum_j Y_i Y_j A_{i,j}
 * \f].
 */
class LTI_MassFracs : public LiquidTranInteraction
{
public:
    LTI_MassFracs(TransportPropertyType tp_ind = TP_UNKNOWN) :
        LiquidTranInteraction(tp_ind) {
        m_model = LTI_MODEL_MASSFRACS;
    }

    //! Return the mixture transport property value.
    /**
     * Takes the separate species transport properties as input (this method
     * does not know what transport property it is at this point.
     */
    doublereal getMixTransProp(doublereal* valueSpecies, doublereal* weightSpecies = 0);
    doublereal getMixTransProp(std::vector<LTPspecies*> LTPptrs);

    //! Return the matrix of binary interaction parameters.
    /**
     * Takes the proper mixing rule for the binary interaction parameters
     * and calculates them: Not implemented for this mixing rule.
     */
    void getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues = 0) {
        mat = (*m_Aij[0]);
    }
};

//! Mixing rule using logarithms of the mole fractions
/**
 * This model is based on the idea that liquid molecules are generally
 * interacting with some energy and entropy of interaction. For transport
 * properties that depend on these energies of interaction, the mixture
 * transport property can be written in terms of its logarithm
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
 *               <!--
 *                 interactions are from speciesA = LiCl(L)
 *                 and speciesB = KCl(L).
 *                   -->
 *                <Eij units="J/kmol"> -1.0e3 </Eij>
 *                <Sij units="J/kmol/K"> 80.0e-5 </Sij>
 *             </interaction>
 *          </compositionDependence>
 *       </viscosity>
 *    </transport>
 *  \endverbatim
 */
class LTI_Log_MoleFracs : public LiquidTranInteraction
{
public:
    LTI_Log_MoleFracs(TransportPropertyType tp_ind = TP_UNKNOWN) :
        LiquidTranInteraction(tp_ind) {
        m_model = LTI_MODEL_LOG_MOLEFRACS;
    }

    //! Return the mixture transport property value.
    /**
     * Takes the separate species transport properties as input (this method
     * does not know what transport property it is at this point.
     */
    doublereal getMixTransProp(doublereal* valueSpecies, doublereal* weightSpecies = 0);
    doublereal getMixTransProp(std::vector<LTPspecies*> LTPptrs);

    //! Return the matrix of binary interaction parameters.
    /**
     * Takes the proper mixing rule for the binary interaction parameters
     * and calculates them: Not implemented for this mixing rule.
     */
    void getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues = 0) {
        mat = m_Eij;
    }
};

//! Transport properties that act like pairwise interactions
//! as in binary diffusion coefficients.
/**
 * This class holds parameters for transport properties expressed as a matrix
 * of pairwise interaction parameters. Input can be provided for constant or
 * Arrhenius forms of the separate parameters.
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
 */
class LTI_Pairwise_Interaction : public LiquidTranInteraction
{
public:
    LTI_Pairwise_Interaction(TransportPropertyType tp_ind = TP_UNKNOWN) :
        LiquidTranInteraction(tp_ind) {
        m_model = LTI_MODEL_PAIRWISE_INTERACTION;
    }

    void setParameters(LiquidTransportParams& trParam);

    //! Return the mixture transport property value.
    /**
     * Takes the separate species transport properties as input (this method
     * does not know what transport property it is at this point.
     */
    doublereal getMixTransProp(doublereal* valueSpecies, doublereal* weightSpecies = 0);
    doublereal getMixTransProp(std::vector<LTPspecies*> LTPptrs);

    //! Return the matrix of binary interaction parameters.
    /**
     * Takes the proper mixing rule for the binary interaction parameters
     * and calculates them
     */
    void getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues = 0);
protected:

    std::vector<LTPspecies*> m_diagonals;
};


//! Stefan Maxwell Diffusion Coefficients can be solved for given
//! ion conductivity, mobility ratios, and self diffusion coeffs.
//! This class is only valid for a common anion mixture of two
//! salts with cations of equal charge.  Hence the name _PPN.
/**
 * This class requres you specify
 *
 * 1 - ion conductivity
 *
 * 2 - mobility ratio of the two cations (set all other ratios to zero)
 *
 * 3 - Self diffusion coefficients of the cations (set others to zero)
 *     is used to calculate the "mutual diffusion coefficient".  The
 *     approximation needed to do so requires the cations have equal charge.
 *
 * We than calculate the Stefan Maxwell Diffusion Coefficients by
 * \f[
 *     \frac{1}{D_{12}} = (1-\epsilon X_A)(1+\epsilon X_B)
 *     \frac{\nu_- + \nu_+}{\nu_-\nu_+^2D}
 *     + \frac{z_-z_+ F^2}{\kappa V R T}
 * \f]
 * \f[
 *     \frac{1}{D_{12}} = -\epsilon X_B(1-\epsilon X_A)
 *     \frac{\nu_- + \nu_+}{\nu_-^2\nu_+D}
 *     - \frac{z_-z_+ F^2}{\kappa V R T}
 * \f]
 * \f[
 *     \frac{1}{D_{23}} = \epsilon X_A(1+\epsilon X_B)
 *     \frac{\nu_- + \nu_+}{\nu_-^2\nu_+D}
 *     - \frac{z_-z_+ F^2}{\kappa V R T}
 * \f]
 * where F is Faraday's constant, RT is the gas constant times the
 * tempurature, and V is the molar volume (basis is moles of ions) that is
 * calculated by the ThermoPhase member. X_A and X_B are the mole fractions
 * of the salts composed of cation(1) and cation(2), respectively, that share
 * a common anion(3). \f$\nu_{+,-}\f$ are the stoichiometric coefficients in
 * the dissociation reaction of the salts to the ions with charges of
 * \f$z_{+,-}\f$.  Assuming that the cations have equal charge, the "mutual
 * diffusion coefficient" is calculated using the cation self diffusion
 * coefficients.
 * \f[
 *     \frac{1}{\nu_-\nu_+D} = \left(1+\frac{\partial \gamma_B}{\partial N_B}
 * \right)\frac{X_A}{D_2^*}+\left(1+\frac{\partial \gamma_A}{\partial N_A}
 * \right)\frac{X_B}{D_1^*}
 * \f]
 * where the self diffusion coefficients, \f$D_i^*\f$, are temperature and
 * composition parameterized inputs and the derivative of the activity
 * coefficient, \f$\frac{\partial \gamma_B}{\partial N_B}\f$, is calculated
 * by the ThermoPhase member using the excess enthalpy and entropy upon mixing.
 *
 * Finally, the deviation of the transferrence numbers from ideality,
 * \f$\epsilon\f$, is calculated from the mobility ratio of the cations.
 * \f[
 *     \epsilon = \frac{1-b_2/b_1}{X_A+X_Bb_2/b_1}
 * \f]
 * Where \f$b_i\f$ are the mobilities of the two cations.  Everywhere,
 * cation 1 corresponds with salt A and cation 2 with salt B.
 *
 *  Sample input for this method is
 *  \verbatim
 *    <transport model="Liquid">
 *       <speciesDiffusivity>
 *          <compositionDependence model="stefanMaxwell_PPN">
 *          </compositionDependence>
 *       </speciesDiffusivity>
 *    </transport>
 *  \endverbatim
 */
class LTI_StefanMaxwell_PPN : public LiquidTranInteraction
{
public:
    LTI_StefanMaxwell_PPN(TransportPropertyType tp_ind = TP_UNKNOWN) :
        LiquidTranInteraction(tp_ind) {
        m_model = LTI_MODEL_STEFANMAXWELL_PPN;
    }

    void setParameters(LiquidTransportParams& trParam);

    //! Return the mixture transport property value.
    /**
     * Takes the separate species transport properties as input (this method
     * does not know what transport property it is at this point.
     */
    doublereal getMixTransProp(doublereal* valueSpecies, doublereal* weightSpecies = 0);
    doublereal getMixTransProp(std::vector<LTPspecies*> LTPptrs);

    //! Return the matrix of binary interaction parameters.
    /**
     * Takes the proper mixing rule for the binary interaction parameters
     * and calculates them
     */
    void getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues = 0);

protected:
    doublereal m_ionCondMix;
    LiquidTranInteraction* m_ionCondMixModel;
    std::vector<LTPspecies*> m_ionCondSpecies;
    typedef std::vector<LTPspecies*> LTPvector;
    DenseMatrix m_mobRatMix;
    std::vector<LiquidTranInteraction*> m_mobRatMixModel;
    std::vector<LTPvector> m_mobRatSpecies;

    std::vector<LiquidTranInteraction*> m_selfDiffMixModel;
    vector_fp m_selfDiffMix;
    std::vector<LTPvector> m_selfDiffSpecies;
};


class LTI_StokesEinstein : public LiquidTranInteraction
{
public:
    LTI_StokesEinstein(TransportPropertyType tp_ind = TP_UNKNOWN) :
        LiquidTranInteraction(tp_ind) {
        m_model = LTI_MODEL_STOKES_EINSTEIN;
    }

    void setParameters(LiquidTransportParams& trParam);

    //! Return the mixture transport property value.
    /**
     * Takes the separate species transport properties
     * as input (this method does not know what
     * transport property it is at this point.
     */
    doublereal getMixTransProp(doublereal* valueSpecies, doublereal* weightSpecies = 0);
    doublereal getMixTransProp(std::vector<LTPspecies*> LTPptrs);

    //! Return the matrix of binary interaction parameters.
    /**
     * Takes the proper mixing rule for the binary interaction parameters
     * and calculates them
     */
    void getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues = 0);

protected:
    std::vector<LTPspecies*> m_viscosity;
    std::vector<LTPspecies*> m_hydroRadius;
};

//! Simple mole fraction weighting of transport properties
/**
 * This model weights the transport property by the mole fractions. The
 * overall formula for the mixture viscosity is
 *
 * \f[ \eta_{mix} = \sum_i X_i \eta_i
 *  + \sum_i \sum_j X_i X_j A_{i,j} \f].
 */
class LTI_MoleFracs_ExpT : public LiquidTranInteraction
{
public:
    LTI_MoleFracs_ExpT(TransportPropertyType tp_ind = TP_UNKNOWN) :
        LiquidTranInteraction(tp_ind) {
        m_model = LTI_MODEL_MOLEFRACS_EXPT;
    }

    //! Return the mixture transport property value.
    /**
     * Takes the separate species transport properties as input (this method
     * does not know what transport property it is at this point.
     */
    doublereal getMixTransProp(doublereal* valueSpecies, doublereal* weightSpecies = 0);
    doublereal getMixTransProp(std::vector<LTPspecies*> LTPptrs);

    //! Return the matrix of binary interaction parameters.
    /**
     * Takes the proper mixing rule for the binary interaction parameters
     * and calculates them: Not Implemented for this mixing rule
     */
    void getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues = 0) {
        mat = (*m_Aij[0]);
    }
};

}

#endif
