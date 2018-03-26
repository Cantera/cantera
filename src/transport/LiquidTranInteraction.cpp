/**
 *  @file LiquidTranInteraction.cpp
 *  Source code for liquid mixture transport property evaluations.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport/LiquidTransportParams.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/thermo/MargulesVPSSTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{

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
    warn_deprecated("Class LiquidTranInteraction", "To be removed after Cantera 2.4");
}

LiquidTranInteraction::~LiquidTranInteraction()
{
    for (size_t k = 0; k < m_Aij.size(); k++) {
        delete m_Aij[k];
    }
    for (size_t k = 0; k < m_Bij.size(); k++) {
        delete m_Bij[k];
    }
    for (size_t k = 0; k < m_Hij.size(); k++) {
        delete m_Hij[k];
    }
    for (size_t k = 0; k < m_Sij.size(); k++) {
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

    for (size_t iChild = 0; iChild < compModelNode.nChildren(); iChild++) {
        XML_Node& xmlChild = compModelNode.child(iChild);
        std::string nodeName = xmlChild.name();
        if (!caseInsensitiveEquals(nodeName, "interaction")) {
            throw CanteraError("TransportFactory::getLiquidInteractionsTransportData",
                               "expected <interaction> element and got <" + nodeName + ">");
        }
        string speciesA = xmlChild.attrib("speciesA");
        string speciesB = xmlChild.attrib("speciesB");
        size_t iSpecies = m_thermo->speciesIndex(speciesA);
        if (iSpecies == npos) {
            throw CanteraError("TransportFactory::getLiquidInteractionsTransportData",
                               "Unknown species " + speciesA);
        }
        size_t jSpecies = m_thermo->speciesIndex(speciesB);
        if (jSpecies == npos) {
            throw CanteraError("TransportFactory::getLiquidInteractionsTransportData",
                               "Unknown species " + speciesB);
        }

        if (xmlChild.hasChild("Eij")) {
            m_Eij(iSpecies,jSpecies) = getFloat(xmlChild, "Eij", "actEnergy");
            m_Eij(iSpecies,jSpecies) /= GasConstant;
            m_Eij(jSpecies,iSpecies) = m_Eij(iSpecies,jSpecies);
        }

        if (xmlChild.hasChild("Aij")) {
            vector_fp poly;
            getFloatArray(xmlChild, poly, true, "toSI", "Aij");
            while (m_Aij.size()<poly.size()) {
                DenseMatrix* aTemp = new DenseMatrix();
                aTemp->resize(nsp, nsp, 0.0);
                m_Aij.push_back(aTemp);
            }
            for (int i = 0; i < (int)poly.size(); i++) {
                (*m_Aij[i])(iSpecies,jSpecies) = poly[i];
            }
        }

        if (xmlChild.hasChild("Bij")) {
            vector_fp poly;
            getFloatArray(xmlChild, poly, true, "toSI", "Bij");
            while (m_Bij.size() < poly.size()) {
                DenseMatrix* bTemp = new DenseMatrix();
                bTemp->resize(nsp, nsp, 0.0);
                m_Bij.push_back(bTemp);
            }
            for (size_t i=0; i<poly.size(); i++) {
                (*m_Bij[i])(iSpecies,jSpecies) = poly[i];
            }
        }

        if (xmlChild.hasChild("Hij")) {
            vector_fp poly;
            getFloatArray(xmlChild, poly, true, "actEnergy", "Hij");
            while (m_Hij.size()<poly.size()) {
                DenseMatrix* hTemp = new DenseMatrix();
                hTemp->resize(nsp, nsp, 0.0);
                m_Hij.push_back(hTemp);
            }
            for (size_t i=0; i<poly.size(); i++) {
                (*m_Hij[i])(iSpecies,jSpecies) = poly[i];
                (*m_Hij[i])(iSpecies,jSpecies) /= GasConstant;
            }
        }

        if (xmlChild.hasChild("Sij")) {
            vector_fp poly;
            getFloatArray(xmlChild, poly, true, "actEnergy", "Sij");
            while (m_Sij.size()<poly.size()) {
                DenseMatrix* sTemp = new DenseMatrix();
                sTemp->resize(nsp, nsp, 0.0);
                m_Sij.push_back(sTemp);
            }
            for (size_t i=0; i<poly.size(); i++) {
                (*m_Sij[i])(iSpecies,jSpecies) = poly[i];
                (*m_Sij[i])(iSpecies,jSpecies) /= GasConstant;
            }
        }

        if (xmlChild.hasChild("Dij")) {
            m_Dij(iSpecies,jSpecies) = getFloat(xmlChild, "Dij", "toSI");
            m_Dij(jSpecies,iSpecies) = m_Dij(iSpecies,jSpecies);
        }
    }
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
            // should be: molefracs[k] = molefracs[k]*speciesWeight[k]; for consistency, but weight(solvent)=1?
        }
    } else {
        throw CanteraError("LTI_Solvent::getMixTransProp","You should be specifying the speciesWeight");
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
            }
            for (size_t k = 0; k < m_Sij.size(); k++) {
                value -= molefracs[i]*molefracs[j]*(*m_Sij[k])(i,j)*pow(molefracs[i], (int) k);
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
                value += molefracs[i]*molefracs[j]*(*m_Hij[k])(i,j)/temp*pow(molefracs[i], (int) k);
            }
            for (size_t k = 0; k < m_Sij.size(); k++) {
                value -= molefracs[i]*molefracs[j]*(*m_Sij[k])(i,j)*pow(molefracs[i], (int) k);
            }
        }
    }
    return exp(value);
}

void LTI_Pairwise_Interaction::setParameters(LiquidTransportParams& trParam)
{
    size_t nsp = m_thermo->nSpecies();
    m_diagonals.resize(nsp, 0);

    for (size_t k = 0; k < nsp; k++) {
        LiquidTransportData& ltd = trParam.LTData[k];
        if (ltd.speciesDiffusivity) {
            m_diagonals[k] = ltd.speciesDiffusivity;
        }
    }
}

doublereal LTI_Pairwise_Interaction::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    throw LTPmodelError("Calling LTI_Pairwise_Interaction::getMixTransProp does not make sense.");
}

doublereal LTI_Pairwise_Interaction::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    throw LTPmodelError("Calling LTI_Pairwise_Interaction::getMixTransProp does not make sense.");
}

void LTI_Pairwise_Interaction::getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    mat.resize(nsp, nsp, 0.0);
    for (size_t i = 0; i < nsp; i++) {
        for (size_t j = 0; j < i; j++) {
            mat(i,j) = mat(j,i) = exp(m_Eij(i,j) / temp) / m_Dij(i,j);
        }
    }

    for (size_t i = 0; i < nsp; i++) {
        if (mat(i,i) == 0.0 && m_diagonals[i]) {
            mat(i,i) = 1.0 / m_diagonals[i]->getSpeciesTransProp();
        }
    }
}

void LTI_StefanMaxwell_PPN::setParameters(LiquidTransportParams& trParam)
{
    size_t nsp = m_thermo->nSpecies();
    m_ionCondMix = 0;
    m_ionCondMixModel = trParam.ionConductivity;
    m_ionCondSpecies.resize(nsp,0);
    m_mobRatMix.resize(nsp,nsp,0.0);
    m_mobRatMixModel.resize(nsp*nsp);
    m_mobRatSpecies.resize(nsp*nsp);
    m_selfDiffMix.resize(nsp,0.0);
    m_selfDiffMixModel.resize(nsp);
    m_selfDiffSpecies.resize(nsp);

    for (size_t k = 0; k < nsp*nsp; k++) {
        m_mobRatMixModel[k] = trParam.mobilityRatio[k];
        m_mobRatSpecies[k].resize(nsp,0);
    }
    for (size_t k = 0; k < nsp; k++) {
        m_selfDiffMixModel[k] = trParam.selfDiffusion[k];
        m_selfDiffSpecies[k].resize(nsp,0);
    }

    for (size_t k = 0; k < nsp; k++) {
        LiquidTransportData& ltd = trParam.LTData[k];
        m_ionCondSpecies[k] = ltd.ionConductivity;
        for (size_t j = 0; j < nsp*nsp; j++) {
            m_mobRatSpecies[j][k] = ltd.mobilityRatio[j];
        }
        for (size_t j = 0; j < nsp; j++) {
            m_selfDiffSpecies[j][k] = ltd.selfDiffusion[j];
        }
    }
}

doublereal LTI_StefanMaxwell_PPN::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    throw LTPmodelError("Calling LTI_StefanMaxwell_PPN::getMixTransProp does not make sense.");
}

doublereal LTI_StefanMaxwell_PPN::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    throw LTPmodelError("Calling LTI_StefanMaxwell_PPN::getMixTransProp does not make sense.");
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
    vector_fp viS(6);
    vector_fp charges(3);
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
    MargulesVPSSTP* marg_thermo = dynamic_cast<MargulesVPSSTP*>(ions_thermo->getNeutralMoleculePhase().get());
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

    double vP = max(viS[cation[0]],viS[cation[1]]);
    double vM = viS[anion[0]];
    double zP = charges[cation[0]];
    double zM = charges[anion[0]];
    vector_fp dlnActCoeffdlnN_diag(neut_molefracs.size(),0.0);
    marg_thermo->getdlnActCoeffdlnN_diag(&dlnActCoeffdlnN_diag[0]);

    double xA = neut_molefracs[neutMolIndex[cation[0]]];
    double xB = neut_molefracs[neutMolIndex[cation[1]]];
    double eps = (1-m_mobRatMix(cation[1],cation[0]))/(xA+xB*m_mobRatMix(cation[1],cation[0]));
    double inv_vP_vM_MutualDiff = (xA*(1-xB+dlnActCoeffdlnN_diag[neutMolIndex[cation[1]]])/m_selfDiffMix[cation[1]]+xB*(1-xA+dlnActCoeffdlnN_diag[neutMolIndex[cation[0]]])/m_selfDiffMix[cation[0]]);

    mat.resize(nsp, nsp, 0.0);
    mat(cation[0],cation[1]) = mat(cation[1],cation[0]) = (1+vM/vP)*(1+eps*xB)*(1-eps*xA)*inv_vP_vM_MutualDiff-zP*zP*Faraday*Faraday/GasConstant/temp/m_ionCondMix/vol;
    mat(cation[0],anion[0]) = mat(anion[0],cation[0]) = (1+vP/vM)*(-eps*xB*(1-eps*xA)*inv_vP_vM_MutualDiff)-zP*zM*Faraday*Faraday/GasConstant/temp/m_ionCondMix/vol;
    mat(cation[1],anion[0]) = mat(anion[0],cation[1]) = (1+vP/vM)*(eps*xA*(1+eps*xB)*inv_vP_vM_MutualDiff)-zP*zM*Faraday*Faraday/GasConstant/temp/m_ionCondMix/vol;
}

doublereal LTI_StokesEinstein::getMixTransProp(doublereal* speciesValues, doublereal* speciesWeight)
{
    throw LTPmodelError("Calling LTI_StokesEinstein::getMixTransProp does not make sense.");
}

doublereal LTI_StokesEinstein::getMixTransProp(std::vector<LTPspecies*> LTPptrs)
{
    throw LTPmodelError("Calling LTI_StokesEinstein::getMixTransProp does not make sense.");
}

void LTI_StokesEinstein::setParameters(LiquidTransportParams& trParam)
{
    size_t nsp = m_thermo->nSpecies();
    m_viscosity.resize(nsp, 0);
    m_hydroRadius.resize(nsp, 0);
    for (size_t k = 0; k < nsp; k++) {
        LiquidTransportData& ltd = trParam.LTData[k];
        m_viscosity[k] = ltd.viscosity;
        m_hydroRadius[k] = ltd.hydroRadius;
    }
}

void LTI_StokesEinstein::getMatrixTransProp(DenseMatrix& mat, doublereal* speciesValues)
{
    size_t nsp = m_thermo->nSpecies();
    doublereal temp = m_thermo->temperature();
    vector_fp viscSpec(nsp);
    vector_fp radiusSpec(nsp);

    for (size_t k = 0; k < nsp; k++) {
        viscSpec[k] = m_viscosity[k]->getSpeciesTransProp();
        radiusSpec[k] = m_hydroRadius[k]->getSpeciesTransProp();
    }

    mat.resize(nsp,nsp, 0.0);
    for (size_t i = 0; i < nsp; i++) {
        for (size_t j = 0; j < nsp; j++) {
            mat(i,j) = (6.0 * Pi * radiusSpec[i] * viscSpec[j]) / GasConstant / temp;
        }
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
