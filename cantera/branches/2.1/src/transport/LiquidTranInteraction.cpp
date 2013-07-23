/**
 *  @file LiquidTranInteraction.cpp
 *  Source code for liquid mixture transport property evaluations.
 */

#include "cantera/transport/LiquidTransportParams.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/thermo/MargulesVPSSTP.h"
#include "cantera/base/stringUtils.h"

using namespace std;
using namespace ctml;

namespace Cantera
{

/**
 * Exception thrown if an error is encountered while reading the
 * transport database.
 */
class LTPError : public CanteraError
{
public:
    explicit LTPError(const std::string& msg)
        : CanteraError("LTPspecies",
                       "error parsing transport data: "
                       + msg + "\n") {}
};

/**
 * Exception thrown if an error is encountered while reading the
 * transport database.
 */
class LTPmodelError : public CanteraError
{
public:
    explicit LTPmodelError(const std::string& msg)
        : CanteraError("LTPspecies",
                       "error parsing transport data: "
                       + msg + "\n") {}
};

LiquidTranInteraction::LiquidTranInteraction(TransportPropertyType tp_ind) :
    m_model(LTI_MODEL_NOTSET),
    m_property(tp_ind)
{
}

LiquidTranInteraction::~LiquidTranInteraction()
{
    size_t kmax = m_Aij.size();
    for (size_t k = 0; k < kmax; k++) {
        delete m_Aij[k];
    }
    kmax = m_Bij.size();
    for (size_t k = 0; k < kmax; k++) {
        delete m_Bij[k];
    }
    kmax = m_Hij.size();
    for (size_t k = 0; k < kmax; k++) {
        delete m_Hij[k];
    }
    kmax = m_Sij.size();
    for (size_t k = 0; k < kmax; k++) {
        delete m_Sij[k];
    }
}

void LiquidTranInteraction::init(const XML_Node& compModelNode,
                                 thermo_t* thermo)
{
    m_thermo = thermo;

    size_t nsp = thermo->nSpecies();
    m_Dij.resize(nsp, nsp, 0.0);
    m_Eij.resize(nsp, nsp, 0.0);
    /*
    m_Aij.resize(nsp);
    m_Bij.resize(nsp);
    m_Hij.resize(nsp);
    m_Sij.resize(nsp);
    for (int k = 0; k < nsp; k++ ){
      (*m_Aij[k]).resize(nsp, nsp, 0.0);
      (*m_Bij[k]).resize(nsp, nsp, 0.0);
      (*m_Hij[k]).resize(nsp, nsp, 0.0);
      (*m_Sij[k]).resize(nsp, nsp, 0.0);
    }
    */

    std::string speciesA;
    std::string speciesB;

    size_t num = compModelNode.nChildren();
    for (size_t iChild = 0; iChild < num; iChild++) {
        XML_Node& xmlChild = compModelNode.child(iChild);
        std::string nodeName = lowercase(xmlChild.name());
        if (nodeName != "interaction") {
            throw CanteraError("TransportFactory::getLiquidInteractionsTransportData",
                               "expected <interaction> element and got <" + nodeName + ">");
        }
        speciesA = xmlChild.attrib("speciesA");
        speciesB = xmlChild.attrib("speciesB");
        size_t iSpecies = m_thermo->speciesIndex(speciesA);
        if (iSpecies == npos) {
            throw CanteraError("TransportFactory::getLiquidInteractionsTransportData",
                               "Unknown species " + speciesA);
        }
        size_t jSpecies = m_thermo->speciesIndex(speciesB);
        if (jSpecies == npos)  {
            throw CanteraError("TransportFactory::getLiquidInteractionsTransportData",
                               "Unknown species " + speciesB);
        }
        /*      if (xmlChild.hasChild("Aij" ) ) {
        m_Aij(iSpecies,jSpecies) = getFloat(xmlChild, "Aij", "toSI" );
        m_Aij(jSpecies,iSpecies) = m_Aij(iSpecies,jSpecies) ;
        }*/

        if (xmlChild.hasChild("Eij")) {
            m_Eij(iSpecies,jSpecies) = getFloat(xmlChild, "Eij", "actEnergy");
            m_Eij(iSpecies,jSpecies) /= GasConstant;
            m_Eij(jSpecies,iSpecies) = m_Eij(iSpecies,jSpecies) ;
        }

        if (xmlChild.hasChild("Aij")) {
            vector_fp poly;
            // poly0 = getFloat(poly, xmlChild, "Aij", "toSI" );
            getFloatArray(xmlChild, poly, true, "toSI", "Aij");
            // if (!poly.size() ) poly.push_back(poly0);
            while (m_Aij.size()<poly.size()) {
                DenseMatrix* aTemp = new DenseMatrix();
                aTemp->resize(nsp, nsp, 0.0);
                m_Aij.push_back(aTemp);
            }
            for (int i = 0; i < (int)poly.size(); i++) {
                (*m_Aij[i])(iSpecies,jSpecies) = poly[i];
                //(*m_Aij[i])(jSpecies,iSpecies) = (*m_Aij[i])(iSpecies,jSpecies) ;
            }
        }

        if (xmlChild.hasChild("Bij")) {
            vector_fp poly;
            getFloatArray(xmlChild, poly, true, "toSI",  "Bij");
            //if (!poly.size() ) poly.push_back(poly0);
            while (m_Bij.size() < poly.size()) {
                DenseMatrix* bTemp = new DenseMatrix();
                bTemp->resize(nsp, nsp, 0.0);
                m_Bij.push_back(bTemp);
            }
            for (size_t i=0; i<poly.size(); i++) {
                (*m_Bij[i])(iSpecies,jSpecies) = poly[i];
                //(*m_Bij[i])(jSpecies,iSpecies) = (*m_Bij[i])(iSpecies,jSpecies) ;
            }
        }

        if (xmlChild.hasChild("Hij")) {
            vector_fp poly;
            // poly0 = getFloat(poly, xmlChild, "Hij", "actEnergy" );
            getFloatArray(xmlChild, poly, true, "actEnergy", "Hij");
            // if (!poly.size() ) poly.push_back(poly0);
            while (m_Hij.size()<poly.size()) {
                DenseMatrix* hTemp = new DenseMatrix();
                hTemp->resize(nsp, nsp, 0.0);
                m_Hij.push_back(hTemp);
            }
            for (size_t i=0; i<poly.size(); i++) {
                (*m_Hij[i])(iSpecies,jSpecies) = poly[i];
                (*m_Hij[i])(iSpecies,jSpecies) /= GasConstant;
                //(*m_Hij[i])(jSpecies,iSpecies) = (*m_Hij[i])(iSpecies,jSpecies) ;
            }
        }

        if (xmlChild.hasChild("Sij")) {
            vector_fp poly;
            //  poly0 = getFloat(poly, xmlChild, "Sij", "actEnergy" );
            getFloatArray(xmlChild, poly, true, "actEnergy", "Sij");
            // if (!poly.size() ) poly.push_back(poly0);
            while (m_Sij.size()<poly.size()) {
                DenseMatrix* sTemp = new DenseMatrix();
                sTemp->resize(nsp, nsp, 0.0);
                m_Sij.push_back(sTemp);
            }
            for (size_t i=0; i<poly.size(); i++) {
                (*m_Sij[i])(iSpecies,jSpecies) = poly[i];
                (*m_Sij[i])(iSpecies,jSpecies) /= GasConstant;
                //(*m_Sij[i])(jSpecies,iSpecies) = (*m_Sij[i])(iSpecies,jSpecies) ;
            }
        }

        /*0      if (xmlChild.hasChild("Sij" ) ) {
        m_Sij(iSpecies,jSpecies) = getFloat(xmlChild, "Sij", "toSI" );
        m_Sij(iSpecies,jSpecies) /= GasConstant;
        //m_Sij(jSpecies,iSpecies) = m_Sij(iSpecies,jSpecies) ;
        }*/

        if (xmlChild.hasChild("Dij")) {
            m_Dij(iSpecies,jSpecies) = getFloat(xmlChild, "Dij", "toSI");
            m_Dij(jSpecies,iSpecies) = m_Dij(iSpecies,jSpecies) ;
        }
    }
}

LiquidTranInteraction::LiquidTranInteraction(const LiquidTranInteraction& right)
{
    *this = right;  //use assignment operator to do other work
}

LiquidTranInteraction& LiquidTranInteraction::operator=(const LiquidTranInteraction& right)
{
    if (&right != this) {
        m_model     = right.m_model;
        m_property  = right.m_property;
        m_thermo    = right.m_thermo;
        //m_trParam   = right.m_trParam;
        m_Aij       = right.m_Aij;
        m_Bij       = right.m_Bij;
        m_Eij       = right.m_Eij;
        m_Hij       = right.m_Hij;
        m_Sij       = right.m_Sij;
        m_Dij       = right.m_Dij;
    }
    return *this;
}

LTI_Solvent::LTI_Solvent(TransportPropertyType tp_ind) :
    LiquidTranInteraction(tp_ind)
{
    m_model = LTI_MODEL_SOLVENT;
}

doublereal LTI_Solvent::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0.0;

    //if weightings are specified, use those
    if (speciesWeight) {
        for (size_t k = 0; k < nsp; k++) {
            //molefracs[k] = molefracs[k];
            // should be: molefracs[k] = molefracs[k]*speciesWeight[k]; for consistency, but weight(solvent)=1?
        }
    } else {
        throw CanteraError("LTI_Solvent::getMixTransProp","You should be specifying the speciesWeight");
        /*  //This does not follow directly a solvent model
        //although if the solvent mole fraction is dominant
        //and the other species values are given or zero,
        //it should work.
        for (int k = 0; k < nsp; k++) {
        value += speciesValues[k] * molefracs[k];
        }*/
    }

    for (size_t i = 0; i < nsp; i++) {
        //presume that the weighting is set to 1.0 for solvent and 0.0 for everything else.
        value += speciesValues[i] * speciesWeight[i];
        if (i == 0) {
            AssertTrace(speciesWeight[i] == 1.0);
        } else {
            AssertTrace(speciesWeight[i] == 0.0);
        }
        for (size_t j = 0; j < nsp; j++) {
            for (size_t k = 0; k < m_Aij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Aij[k])(i,j)*pow(molefracs[i], (int) k);
            }
            for (size_t k = 0; k < m_Bij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Bij[k])(i,j)*temp*pow(molefracs[i], (int) k);
            }
        }
    }

    return value;
}

doublereal LTI_Solvent::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0.0;

    for (size_t k = 0; k < nsp; k++) {
        //molefracs[k] = molefracs[k];
        // should be:      molefracs[k] = molefracs[k]*LTPptrs[k]->getMixWeight(); for consistency, but weight(solvent)=1?
    }

    for (size_t i = 0; i < nsp; i++) {
        //presume that the weighting is set to 1.0 for solvent and 0.0 for everything else.
        value += LTPptrs[i]->getSpeciesTransProp() * LTPptrs[i]->getMixWeight();
        for (size_t j = 0; j < nsp; j++) {
            for (size_t k = 0; k < m_Aij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Aij[k])(i,j)*pow(molefracs[i], (int) k);
            }
            for (size_t k = 0; k < m_Bij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Bij[k])(i,j)*temp*pow(molefracs[i], (int) k);
            }
        }
    }

    return value;
}

void LTI_Solvent::getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues)
{
    mat = (*m_Aij[0]);
}

doublereal LTI_MoleFracs::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0;

    //if weightings are specified, use those
    if (speciesWeight) {
        for (size_t k = 0; k < nsp; k++) {
            molefracs[k] = molefracs[k]*speciesWeight[k];
        }
    } else {
        throw CanteraError("LTI_MoleFracs::getMixTransProp","You should be specifying the speciesWeight");
    }

    for (size_t i = 0; i < nsp; i++) {
        value += speciesValues[i] * molefracs[i];
        for (size_t j = 0; j < nsp; j++) {
            for (size_t k = 0; k < m_Aij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Aij[k])(i,j)*pow(molefracs[i], (int) k);
            }
            for (size_t k = 0; k < m_Bij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Bij[k])(i,j)*temp*pow(molefracs[i], (int) k);
            }
        }
    }

    return value;
}

doublereal LTI_MoleFracs::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0;

    for (size_t k = 0; k < nsp; k++) {
        molefracs[k] = molefracs[k]*LTPptrs[k]->getMixWeight();
    }

    for (size_t i = 0; i < nsp; i++) {
        value += LTPptrs[i]->getSpeciesTransProp() * molefracs[i];
        for (size_t j = 0; j < nsp; j++) {
            for (size_t k = 0; k < m_Aij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Aij[k])(i,j)*pow(molefracs[i], (int) k);
            }
            for (size_t k = 0; k < m_Bij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Bij[k])(i,j)*temp*pow(molefracs[i], (int) k);
            }
        }
    }
    return value;
}

doublereal LTI_MassFracs::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp massfracs(nsp);
    m_thermo->getMassFractions(&massfracs[0]);

    doublereal value = 0;

    //if weightings are specified, use those
    if (speciesWeight) {
        for (size_t k = 0; k < nsp; k++) {
            massfracs[k] = massfracs[k]*speciesWeight[k];
        }
    } else {
        throw CanteraError("LTI_MassFracs::getMixTransProp","You should be specifying the speciesWeight");
    }

    for (size_t i = 0; i < nsp; i++) {
        value += speciesValues[i] * massfracs[i];
        for (size_t j = 0; j < nsp; j++) {
            for (size_t k = 0; k < m_Aij.size(); k++) {
                value += massfracs[i]*massfracs[j]*(*m_Aij[k])(i,j)*pow(massfracs[i], (int) k);
            }
            for (size_t k = 0; k < m_Bij.size(); k++) {
                value += massfracs[i]*massfracs[j]*(*m_Bij[k])(i,j)*temp*pow(massfracs[i], (int) k);
            }
        }
    }

    return value;
}

doublereal LTI_MassFracs::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp massfracs(nsp);
    m_thermo->getMassFractions(&massfracs[0]);

    doublereal value = 0;

    for (size_t k = 0; k < nsp; k++) {
        massfracs[k] = massfracs[k]*LTPptrs[k]->getMixWeight();
    }

    for (size_t i = 0; i < nsp; i++) {
        value += LTPptrs[i]->getSpeciesTransProp() * massfracs[i];
        for (size_t j = 0; j < nsp; j++) {
            for (size_t k = 0; k < m_Aij.size(); k++) {
                value += massfracs[i]*massfracs[j]*(*m_Aij[k])(i,j)*pow(massfracs[i], (int) k);
            }
            for (size_t k = 0; k < m_Bij.size(); k++) {
                value += massfracs[i]*massfracs[j]*(*m_Bij[k])(i,j)*temp*pow(massfracs[i], (int) k);
            }
        }
    }

    return value;
}

doublereal LTI_Log_MoleFracs::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);



    doublereal value = 0;

    //if weightings are specified, use those
    if (speciesWeight) {
        for (size_t k = 0; k < nsp; k++) {
            molefracs[k] = molefracs[k]*speciesWeight[k];
        }
    } else {
        throw CanteraError("LTI_Log_MoleFracs::getMixTransProp","You probably should have a speciesWeight when you call getMixTransProp to convert ion mole fractions to molecular mole fractions");
    }

    for (size_t i = 0; i < nsp; i++) {
        value += log(speciesValues[i]) * molefracs[i];
        for (size_t j = 0; j < nsp; j++) {
            for (size_t k = 0; k < m_Hij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Hij[k])(i,j)/temp*pow(molefracs[i], (int) k);
                //cout << "value = " << value << ", m_Sij = " << (*m_Sij[k])(i,j) << ", m_Hij = " << (*m_Hij[k])(i,j) << endl;
            }
            for (size_t k = 0; k < m_Sij.size(); k++) {
                value -= molefracs[i]*molefracs[j]*(*m_Sij[k])(i,j)*pow(molefracs[i], (int) k);
                //cout << "value = " << value << ", m_Sij = " << (*m_Sij[k])(i,j) << ", m_Hij = " << (*m_Hij[k])(i,j) << endl;
            }
        }
    }

    return exp(value);
}

doublereal LTI_Log_MoleFracs::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);


    doublereal value = 0;

    //if weightings are specified, use those

    for (size_t k = 0; k < nsp; k++) {
        molefracs[k] = molefracs[k]*LTPptrs[k]->getMixWeight();
    }

    for (size_t i = 0; i < nsp; i++) {
        value += log(LTPptrs[i]->getSpeciesTransProp()) * molefracs[i];
        for (size_t j = 0; j < nsp; j++) {
            for (size_t k = 0; k < m_Hij.size(); k++) {
                value +=  molefracs[i]*molefracs[j]*(*m_Hij[k])(i,j)/temp*pow(molefracs[i], (int) k);
                //cout << "1 = " << molefracs[i]+molefracs[j] << endl;
                //cout << "value = " << value << ", m_Sij = " << (*m_Sij[k])(i,j) << ", m_Hij = " << (*m_Hij[k])(i,j) << endl;
            }
            for (size_t k = 0; k < m_Sij.size(); k++) {
                value -=  molefracs[i]*molefracs[j]*(*m_Sij[k])(i,j)*pow(molefracs[i], (int) k);
                //cout << "1 = " << molefracs[i]+molefracs[j] << endl;
                //cout << "value = " << value << ", m_Sij = " << (*m_Sij[k])(i,j) << ", m_Hij = " << (*m_Hij[k])(i,j) << endl;
            }
        }
    }

    value = exp(value);
    //    cout << ", viscSpeciesA = " << LTPptrs[0]->getSpeciesTransProp() << endl;
    //cout << ", viscSpeciesB = " << LTPptrs[1]->getSpeciesTransProp() << endl;
    //cout << "value = " << value << " FINAL" << endl;
    return value;
}

void LTI_Pairwise_Interaction::setParameters(LiquidTransportParams& trParam)
{
    size_t nsp = m_thermo->nSpecies();
    m_diagonals.resize(nsp, 0);

    for (size_t k = 0; k < nsp; k++) {
        Cantera::LiquidTransportData& ltd = trParam.LTData[k];
        if (ltd.speciesDiffusivity) {
            m_diagonals[k] = ltd.speciesDiffusivity;
        }
    }
}

doublereal LTI_Pairwise_Interaction::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0;

    throw LTPmodelError("Calling LTI_Pairwise_Interaction::getMixTransProp does not make sense.");

    return value;
}

doublereal LTI_Pairwise_Interaction::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0;

    throw LTPmodelError("Calling LTI_Pairwise_Interaction::getMixTransProp does not make sense.");

    return value;
}

void LTI_Pairwise_Interaction::getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    mat.resize(nsp, nsp, 0.0);
    for (size_t i = 0; i < nsp; i++)
        for (size_t j = 0; j < i; j++) {
            mat(i,j) = mat(j,i) = exp(m_Eij(i,j) / temp) / m_Dij(i,j);
        }

    for (size_t i = 0; i < nsp; i++)
        if (mat(i,i) == 0.0 && m_diagonals[i]) {
            mat(i,i) = 1.0 / m_diagonals[i]->getSpeciesTransProp() ;
        }
}

void LTI_StefanMaxwell_PPN::setParameters(LiquidTransportParams& trParam)
{
    size_t nsp = m_thermo->nSpecies();
    size_t nsp2 = nsp*nsp;
    //vector<std

    m_ionCondMix = 0;
    m_ionCondMixModel = trParam.ionConductivity;
    //trParam.ionConductivity = 0;
    m_ionCondSpecies.resize(nsp,0);
    m_mobRatMix.resize(nsp,nsp,0.0);
    m_mobRatMixModel.resize(nsp2);
    m_mobRatSpecies.resize(nsp2);
    m_selfDiffMix.resize(nsp,0.0);
    m_selfDiffMixModel.resize(nsp);
    m_selfDiffSpecies.resize(nsp);

    for (size_t k = 0; k < nsp2; k++) {
        m_mobRatMixModel[k] = trParam.mobilityRatio[k];
        //trParam.mobilityRatio[k] = 0;
        m_mobRatSpecies[k].resize(nsp,0);
    }
    for (size_t k = 0; k < nsp; k++) {
        m_selfDiffMixModel[k] = trParam.selfDiffusion[k];
        //trParam.selfDiffusion[k] = 0;
        m_selfDiffSpecies[k].resize(nsp,0);
    }

    for (size_t k = 0; k < nsp; k++) {
        Cantera::LiquidTransportData& ltd = trParam.LTData[k];
        m_ionCondSpecies[k]   =  ltd.ionConductivity;
        //ltd.ionConductivity = 0;
        for (size_t j = 0; j < nsp2; j++) {
            m_mobRatSpecies[j][k] = ltd.mobilityRatio[j];
            //ltd.mobilityRatio[j] = 0;
        }
        for (size_t j = 0; j < nsp; j++) {
            m_selfDiffSpecies[j][k] = ltd.selfDiffusion[j];
            //ltd.selfDiffusion[j] = 0;
        }
    }
}

doublereal LTI_StefanMaxwell_PPN::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0;

    throw LTPmodelError("Calling LTI_StefanMaxwell_PPN::getMixTransProp does not make sense.");

    return value;
}

doublereal LTI_StefanMaxwell_PPN::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0;

    throw LTPmodelError("Calling LTI_StefanMaxwell_PPN::getMixTransProp does not make sense.");

    return value;
}

void LTI_StefanMaxwell_PPN::getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues)
{
    IonsFromNeutralVPSSTP* ions_thermo = dynamic_cast<IonsFromNeutralVPSSTP*>(m_thermo);
    size_t nsp = m_thermo->nSpecies();
    if (nsp != 3) {
        throw CanteraError("LTI_StefanMaxwell_PPN::getMatrixTransProp","Function may only be called with a 3-ion system");
    }
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);
    vector_fp neut_molefracs;
    ions_thermo->getNeutralMolecMoleFractions(neut_molefracs);
    vector<size_t> cation;
    vector<size_t> anion;
    ions_thermo->getCationList(cation);
    ions_thermo->getAnionList(anion);

    // Reaction Coeffs and Charges
    std::vector<double> viS(6);
    std::vector<double> charges(3);
    std::vector<size_t> neutMolIndex(3);
    ions_thermo->getDissociationCoeffs(viS,charges,neutMolIndex);

    if (anion.size() != 1) {
        throw CanteraError("LTI_StefanMaxwell_PPN::getMatrixTransProp","Must have one anion only for StefanMaxwell_PPN");
    }
    if (cation.size() != 2) {
        throw CanteraError("LTI_StefanMaxwell_PPN::getMatrixTransProp","Must have two cations of equal charge for StefanMaxwell_PPN");
    }
    if (charges[cation[0]] != charges[cation[1]]) {
        throw CanteraError("LTI_StefanMaxwell_PPN::getMatrixTransProp","Cations must be of equal charge for StefanMaxwell_PPN");
    }

    m_ionCondMix = m_ionCondMixModel->getMixTransProp(m_ionCondSpecies);

    MargulesVPSSTP* marg_thermo = dynamic_cast<MargulesVPSSTP*>(ions_thermo->neutralMoleculePhase_);
    doublereal vol = m_thermo->molarVolume();

    size_t k = 0;
    for (size_t j = 0; j < nsp; j++) {
        for (size_t i = 0; i < nsp; i++) {
            if (m_mobRatMixModel[k]) {
                m_mobRatMix(i,j) = m_mobRatMixModel[k]->getMixTransProp(m_mobRatSpecies[k]);
                if (m_mobRatMix(i,j) > 0.0) {
                    m_mobRatMix(j,i) = 1.0/m_mobRatMix(i,j);
                }
            }
            k++;
        }
    }


    for (k = 0; k < nsp; k++) {
        m_selfDiffMix[k] = m_selfDiffMixModel[k]->getMixTransProp(m_selfDiffSpecies[k]);
    }

    //! @todo Suspicious implicit conversion from double to int.
    int vP = max(viS[cation[0]],viS[cation[1]]);
    int vM = viS[anion[0]];
    int zP = charges[cation[0]];
    int zM = charges[anion[0]];
    doublereal xA, xB, eps;
    doublereal inv_vP_vM_MutualDiff;
    vector_fp dlnActCoeffdlnN_diag;
    dlnActCoeffdlnN_diag.resize(neut_molefracs.size(),0.0);
    marg_thermo->getdlnActCoeffdlnN_diag(&dlnActCoeffdlnN_diag[0]);

    xA = neut_molefracs[neutMolIndex[cation[0]]];
    xB = neut_molefracs[neutMolIndex[cation[1]]];
    eps = (1-m_mobRatMix(cation[1],cation[0]))/(xA+xB*m_mobRatMix(cation[1],cation[0]));
    inv_vP_vM_MutualDiff = (xA*(1-xB+dlnActCoeffdlnN_diag[neutMolIndex[cation[1]]])/m_selfDiffMix[cation[1]]+xB*(1-xA+dlnActCoeffdlnN_diag[neutMolIndex[cation[0]]])/m_selfDiffMix[cation[0]]);

    mat.resize(nsp, nsp, 0.0);
    mat(cation[0],cation[1]) = mat(cation[1],cation[0]) = (1+vM/vP)*(1+eps*xB)*(1-eps*xA)*inv_vP_vM_MutualDiff-zP*zP*Faraday*Faraday/GasConstant/temp/m_ionCondMix/vol;
    mat(cation[0],anion[0]) = mat(anion[0],cation[0]) = (1+vP/vM)*(-eps*xB*(1-eps*xA)*inv_vP_vM_MutualDiff)-zP*zM*Faraday*Faraday/GasConstant/temp/m_ionCondMix/vol;
    mat(cation[1],anion[0]) = mat(anion[0],cation[1]) = (1+vP/vM)*(eps*xA*(1+eps*xB)*inv_vP_vM_MutualDiff)-zP*zM*Faraday*Faraday/GasConstant/temp/m_ionCondMix/vol;
}

doublereal LTI_StokesEinstein::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0;

    throw LTPmodelError("Calling LTI_StokesEinstein::getMixTransProp does not make sense.");

    return value;
}

doublereal LTI_StokesEinstein::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0;

    throw LTPmodelError("Calling LTI_StokesEinstein::getMixTransProp does not make sense.");

    return value;
}

void LTI_StokesEinstein::setParameters(LiquidTransportParams& trParam)
{
    size_t nsp = m_thermo->nSpecies();
    m_viscosity.resize(nsp, 0);
    m_hydroRadius.resize(nsp, 0);
    for (size_t k = 0; k < nsp; k++) {
        Cantera::LiquidTransportData& ltd = trParam.LTData[k];
        m_viscosity[k]   =  ltd.viscosity;
        m_hydroRadius[k] =  ltd.hydroRadius;
    }
}

void LTI_StokesEinstein::getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();

    vector_fp viscSpec(nsp);
    vector_fp radiusSpec(nsp);

    for (size_t k = 0; k < nsp; k++) {
        viscSpec[k] = m_viscosity[k]->getSpeciesTransProp() ;
        radiusSpec[k] = m_hydroRadius[k]->getSpeciesTransProp() ;
    }

    mat.resize(nsp,nsp, 0.0);
    for (size_t i = 0; i < nsp; i++)
        for (size_t j = 0; j < nsp; j++) {
            mat(i,j) = (6.0 * Pi * radiusSpec[i] * viscSpec[j]) / GasConstant / temp;
        }
}

doublereal LTI_MoleFracs_ExpT::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0;

    //if weightings are specified, use those
    if (speciesWeight) {
        for (size_t k = 0; k < nsp; k++) {
            molefracs[k] = molefracs[k]*speciesWeight[k];
        }
    } else {
        throw CanteraError("LTI_MoleFracs_ExpT::getMixTransProp","You should be specifying the speciesWeight");
    }

    for (size_t i = 0; i < nsp; i++) {
        value += speciesValues[i] * molefracs[i];
        for (size_t j = 0; j < nsp; j++) {
            for (size_t k = 0; k < m_Aij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Aij[k])(i,j)*pow(molefracs[i], (int) k)*exp((*m_Bij[k])(i,j)*temp);
            }
        }
    }

    return value;
}

doublereal LTI_MoleFracs_ExpT::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal value = 0;

    for (size_t k = 0; k < nsp; k++) {
        molefracs[k] = molefracs[k]*LTPptrs[k]->getMixWeight();
    }

    for (size_t i = 0; i < nsp; i++) {
        value += LTPptrs[i]->getSpeciesTransProp() * molefracs[i];
        for (size_t j = 0; j < nsp; j++) {
            for (size_t k = 0; k < m_Aij.size(); k++) {
                value += molefracs[i]*molefracs[j]*(*m_Aij[k])(i,j)*pow(molefracs[i], (int) k)*exp((*m_Bij[k])(i,j)*temp);
            }
        }
    }
    return value;
}

} //namespace Cantera
