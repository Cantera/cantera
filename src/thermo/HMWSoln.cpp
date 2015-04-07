/**
 *  @file HMWSoln.cpp
 *    Definitions for the %HMWSoln ThermoPhase object, which
 *    models concentrated electrolyte solutions
 *    (see \ref thermoprops and \link Cantera::HMWSoln HMWSoln \endlink) .
 *
 * Class %HMWSoln represents a concentrated liquid electrolyte phase which
 * obeys the Pitzer formulation for nonideality using molality-based
 * standard states.
 *
 * This version of the code was modified to have the binary Beta2 Pitzer
 * parameter consistent with the temperature expansions used for Beta0,
 * Beta1, and Cphi.(CFJC, SNL)
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/HMWSoln.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/WaterProps.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/base/stringUtils.h"
#include <cstdio>

namespace Cantera
{

HMWSoln::HMWSoln() :
    MolalityVPSSTP(),
    m_formPitzer(PITZERFORM_BASE),
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_WATER),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0),
    m_molalitiesAreCropped(false),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.2),
    IMS_gamma_o_min_(1.0E-5),
    IMS_gamma_k_min_(10.0),
    IMS_cCut_(0.05),
    IMS_slopefCut_(0.6),
    IMS_slopegCut_(0.0),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0),
    MC_X_o_cutoff_(0.0),
    MC_X_o_min_(0.0),
    MC_slopepCut_(0.0),
    MC_dpCut_(0.0),
    MC_epCut_(0.0),
    MC_apCut_(0.0),
    MC_bpCut_(0.0),
    MC_cpCut_(0.0),
    CROP_ln_gamma_o_min(-6.0),
    CROP_ln_gamma_o_max(3.0),
    CROP_ln_gamma_k_min(-5.0),
    CROP_ln_gamma_k_max(15.0),
    m_last_is(-1.0),
    m_debugCalc(0)
{
    for (size_t i = 0; i < 17; i++) {
        elambda[i] = 0.0;
        elambda1[i] = 0.0;
    }
}

HMWSoln::HMWSoln(const std::string& inputFile, const std::string& id_) :
    MolalityVPSSTP(),
    m_formPitzer(PITZERFORM_BASE),
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_WATER),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0),
    m_molalitiesAreCropped(false),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.2),
    IMS_gamma_o_min_(1.0E-5),
    IMS_gamma_k_min_(10.0),
    IMS_cCut_(0.05),
    IMS_slopefCut_(0.6),
    IMS_slopegCut_(0.0),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0),
    MC_X_o_cutoff_(0.0),
    MC_X_o_min_(0.0),
    MC_slopepCut_(0.0),
    MC_dpCut_(0.0),
    MC_epCut_(0.0),
    MC_apCut_(0.0),
    MC_bpCut_(0.0),
    MC_cpCut_(0.0),
    CROP_ln_gamma_o_min(-6.0),
    CROP_ln_gamma_o_max(3.0),
    CROP_ln_gamma_k_min(-5.0),
    CROP_ln_gamma_k_max(15.0),
    m_last_is(-1.0),
    m_debugCalc(0)
{
    for (int i = 0; i < 17; i++) {
        elambda[i] = 0.0;
        elambda1[i] = 0.0;
    }
    initThermoFile(inputFile, id_);
}

HMWSoln::HMWSoln(XML_Node& phaseRoot, const std::string& id_) :
    MolalityVPSSTP(),
    m_formPitzer(PITZERFORM_BASE),
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_WATER),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0),
    m_molalitiesAreCropped(false),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.2),
    IMS_gamma_o_min_(1.0E-5),
    IMS_gamma_k_min_(10.0),
    IMS_cCut_(0.05),
    IMS_slopefCut_(0.6),
    IMS_slopegCut_(0.0),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0),
    MC_X_o_cutoff_(0.0),
    MC_X_o_min_(0.0),
    MC_slopepCut_(0.0),
    MC_dpCut_(0.0),
    MC_epCut_(0.0),
    MC_apCut_(0.0),
    MC_bpCut_(0.0),
    MC_cpCut_(0.0),
    CROP_ln_gamma_o_min(-6.0),
    CROP_ln_gamma_o_max(3.0),
    CROP_ln_gamma_k_min(-5.0),
    CROP_ln_gamma_k_max(15.0),
    m_last_is(-1.0),
    m_debugCalc(0)
{
    for (int i = 0; i < 17; i++) {
        elambda[i] = 0.0;
        elambda1[i] = 0.0;
    }
    importPhase(*findXMLPhase(&phaseRoot, id_), this);
}

HMWSoln::HMWSoln(const HMWSoln& b) :
    MolalityVPSSTP(),
    m_formPitzer(PITZERFORM_BASE),
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_WATER),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0),
    m_molalitiesAreCropped(false),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.2),
    IMS_gamma_o_min_(1.0E-5),
    IMS_gamma_k_min_(10.0),
    IMS_cCut_(0.05),
    IMS_slopefCut_(0.6),
    IMS_slopegCut_(0.0),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0),
    MC_X_o_cutoff_(0.0),
    MC_X_o_min_(0.0),
    MC_slopepCut_(0.0),
    MC_dpCut_(0.0),
    MC_epCut_(0.0),
    MC_apCut_(0.0),
    MC_bpCut_(0.0),
    MC_cpCut_(0.0),
    CROP_ln_gamma_o_min(-6.0),
    CROP_ln_gamma_o_max(3.0),
    CROP_ln_gamma_k_min(-5.0),
    CROP_ln_gamma_k_max(15.0),
    m_last_is(-1.0),
    m_debugCalc(0)
{
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = b;
}

HMWSoln& HMWSoln::operator=(const HMWSoln& b)
{
    if (&b != this) {
        MolalityVPSSTP::operator=(b);
        m_formPitzer          = b.m_formPitzer;
        m_formPitzerTemp      = b.m_formPitzerTemp;
        m_formGC              = b.m_formGC;
        m_Aionic              = b.m_Aionic;
        m_IionicMolality      = b.m_IionicMolality;
        m_maxIionicStrength   = b.m_maxIionicStrength;
        m_TempPitzerRef       = b.m_TempPitzerRef;
        m_IionicMolalityStoich= b.m_IionicMolalityStoich;
        m_form_A_Debye        = b.m_form_A_Debye;
        m_A_Debye             = b.m_A_Debye;

        // This is an internal shallow copy of the PDSS_Water pointer
        m_waterSS = providePDSS(0);
        if (!m_waterSS) {
            throw CanteraError("HMWSoln::operator=()", "Dynamic cast to PDSS_Water failed");
        }

        m_densWaterSS         = b.m_densWaterSS;

        delete m_waterProps;
        m_waterProps = 0;
        if (b.m_waterProps) {
            m_waterProps = new WaterProps(dynamic_cast<PDSS_Water*>(m_waterSS));
        }

        m_pp                  = b.m_pp;
        m_tmpV                = b.m_tmpV;
        m_speciesCharge_Stoich= b.m_speciesCharge_Stoich;
        m_Beta0MX_ij          = b.m_Beta0MX_ij;
        m_Beta0MX_ij_L        = b.m_Beta0MX_ij_L;
        m_Beta0MX_ij_LL       = b.m_Beta0MX_ij_LL;
        m_Beta0MX_ij_P        = b.m_Beta0MX_ij_P;
        m_Beta0MX_ij_coeff    = b.m_Beta0MX_ij_coeff;
        m_Beta1MX_ij          = b.m_Beta1MX_ij;
        m_Beta1MX_ij_L        = b.m_Beta1MX_ij_L;
        m_Beta1MX_ij_LL       = b.m_Beta1MX_ij_LL;
        m_Beta1MX_ij_P        = b.m_Beta1MX_ij_P;
        m_Beta1MX_ij_coeff    = b.m_Beta1MX_ij_coeff;
        m_Beta2MX_ij          = b.m_Beta2MX_ij;
        m_Beta2MX_ij_L        = b.m_Beta2MX_ij_L;
        m_Beta2MX_ij_LL       = b.m_Beta2MX_ij_LL;
        m_Beta2MX_ij_P        = b.m_Beta2MX_ij_P;
        m_Beta2MX_ij_coeff    = b.m_Beta2MX_ij_coeff;
        m_Alpha1MX_ij         = b.m_Alpha1MX_ij;
        m_Alpha2MX_ij         = b.m_Alpha2MX_ij;
        m_CphiMX_ij           = b.m_CphiMX_ij;
        m_CphiMX_ij_L         = b.m_CphiMX_ij_L;
        m_CphiMX_ij_LL        = b.m_CphiMX_ij_LL;
        m_CphiMX_ij_P         = b.m_CphiMX_ij_P;
        m_CphiMX_ij_coeff     = b.m_CphiMX_ij_coeff;
        m_Theta_ij            = b.m_Theta_ij;
        m_Theta_ij_L          = b.m_Theta_ij_L;
        m_Theta_ij_LL         = b.m_Theta_ij_LL;
        m_Theta_ij_P          = b.m_Theta_ij_P;
        m_Theta_ij_coeff      = b.m_Theta_ij_coeff;
        m_Psi_ijk             = b.m_Psi_ijk;
        m_Psi_ijk_L           = b.m_Psi_ijk_L;
        m_Psi_ijk_LL          = b.m_Psi_ijk_LL;
        m_Psi_ijk_P           = b.m_Psi_ijk_P;
        m_Psi_ijk_coeff       = b.m_Psi_ijk_coeff;
        m_Lambda_nj           = b.m_Lambda_nj;
        m_Lambda_nj_L         = b.m_Lambda_nj_L;
        m_Lambda_nj_LL        = b.m_Lambda_nj_LL;
        m_Lambda_nj_P         = b.m_Lambda_nj_P;
        m_Lambda_nj_coeff     = b.m_Lambda_nj_coeff;

        m_Mu_nnn              = b.m_Mu_nnn;
        m_Mu_nnn_L            = b.m_Mu_nnn_L;
        m_Mu_nnn_LL           = b.m_Mu_nnn_LL;
        m_Mu_nnn_P            = b.m_Mu_nnn_P;
        m_Mu_nnn_coeff        = b.m_Mu_nnn_coeff;

        m_lnActCoeffMolal_Scaled       = b.m_lnActCoeffMolal_Scaled;
        m_lnActCoeffMolal_Unscaled     = b.m_lnActCoeffMolal_Unscaled;

        m_dlnActCoeffMolaldT_Scaled    = b.m_dlnActCoeffMolaldT_Scaled;
        m_dlnActCoeffMolaldT_Unscaled  = b.m_dlnActCoeffMolaldT_Unscaled;

        m_d2lnActCoeffMolaldT2_Scaled  = b.m_d2lnActCoeffMolaldT2_Scaled;
        m_d2lnActCoeffMolaldT2_Unscaled= b.m_d2lnActCoeffMolaldT2_Unscaled;

        m_dlnActCoeffMolaldP_Scaled    = b.m_dlnActCoeffMolaldP_Scaled;
        m_dlnActCoeffMolaldP_Unscaled  = b.m_dlnActCoeffMolaldP_Unscaled;

        m_molalitiesCropped    = b.m_molalitiesCropped;
        m_molalitiesAreCropped = b.m_molalitiesAreCropped;
        m_CounterIJ            = b.m_CounterIJ;

        m_gfunc_IJ            = b.m_gfunc_IJ;
        m_g2func_IJ           = b.m_g2func_IJ;
        m_hfunc_IJ            = b.m_hfunc_IJ;
        m_h2func_IJ           = b.m_h2func_IJ;
        m_BMX_IJ              = b.m_BMX_IJ;
        m_BMX_IJ_L            = b.m_BMX_IJ_L;
        m_BMX_IJ_LL           = b.m_BMX_IJ_LL;
        m_BMX_IJ_P            = b.m_BMX_IJ_P;
        m_BprimeMX_IJ         = b.m_BprimeMX_IJ;
        m_BprimeMX_IJ_L       = b.m_BprimeMX_IJ_L;
        m_BprimeMX_IJ_LL      = b.m_BprimeMX_IJ_LL;
        m_BprimeMX_IJ_P       = b.m_BprimeMX_IJ_P;
        m_BphiMX_IJ           = b.m_BphiMX_IJ;
        m_BphiMX_IJ_L         = b.m_BphiMX_IJ_L;
        m_BphiMX_IJ_LL        = b.m_BphiMX_IJ_LL;
        m_BphiMX_IJ_P         = b.m_BphiMX_IJ_P;
        m_Phi_IJ              = b.m_Phi_IJ;
        m_Phi_IJ_L            = b.m_Phi_IJ_L;
        m_Phi_IJ_LL           = b.m_Phi_IJ_LL;
        m_Phi_IJ_P            = b.m_Phi_IJ_P;
        m_Phiprime_IJ         = b.m_Phiprime_IJ;
        m_PhiPhi_IJ           = b.m_PhiPhi_IJ;
        m_PhiPhi_IJ_L         = b.m_PhiPhi_IJ_L;
        m_PhiPhi_IJ_LL        = b.m_PhiPhi_IJ_LL;
        m_PhiPhi_IJ_P         = b.m_PhiPhi_IJ_P;
        m_CMX_IJ              = b.m_CMX_IJ;
        m_CMX_IJ_L            = b.m_CMX_IJ_L;
        m_CMX_IJ_LL           = b.m_CMX_IJ_LL;
        m_CMX_IJ_P            = b.m_CMX_IJ_P;
        m_gamma_tmp           = b.m_gamma_tmp;

        IMS_lnActCoeffMolal_  = b.IMS_lnActCoeffMolal_;
        IMS_typeCutoff_       = b.IMS_typeCutoff_;
        IMS_X_o_cutoff_       = b.IMS_X_o_cutoff_;
        IMS_gamma_o_min_      = b.IMS_gamma_o_min_;
        IMS_gamma_k_min_      = b.IMS_gamma_k_min_;
        IMS_cCut_             = b.IMS_cCut_;
        IMS_slopefCut_        = b.IMS_slopefCut_;
        IMS_dfCut_            = b.IMS_dfCut_;
        IMS_efCut_            = b.IMS_efCut_;
        IMS_afCut_            = b.IMS_afCut_;
        IMS_bfCut_            = b.IMS_bfCut_;
        IMS_slopegCut_        = b.IMS_slopegCut_;
        IMS_dgCut_            = b.IMS_dgCut_;
        IMS_egCut_            = b.IMS_egCut_;
        IMS_agCut_            = b.IMS_agCut_;
        IMS_bgCut_            = b.IMS_bgCut_;
        MC_X_o_cutoff_        = b.MC_X_o_cutoff_;
        MC_X_o_min_           = b.MC_X_o_min_;
        MC_slopepCut_         = b.MC_slopepCut_;
        MC_dpCut_             = b.MC_dpCut_;
        MC_epCut_             = b.MC_epCut_;
        MC_apCut_             = b.MC_apCut_;
        MC_bpCut_             = b.MC_bpCut_;
        MC_cpCut_             = b.MC_cpCut_;
        CROP_ln_gamma_o_min   = b.CROP_ln_gamma_o_min;
        CROP_ln_gamma_o_max   = b.CROP_ln_gamma_o_max;
        CROP_ln_gamma_k_min   = b.CROP_ln_gamma_k_min;
        CROP_ln_gamma_k_max   = b.CROP_ln_gamma_k_max;
        CROP_speciesCropped_  = b.CROP_speciesCropped_;

        m_debugCalc           = b.m_debugCalc;
    }
    return *this;
}

HMWSoln::HMWSoln(int testProb) :
    MolalityVPSSTP(),
    m_formPitzer(PITZERFORM_BASE),
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_WATER),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0),
    m_molalitiesAreCropped(false),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.2),
    IMS_gamma_o_min_(1.0E-5),
    IMS_gamma_k_min_(10.0),
    IMS_cCut_(0.05),
    IMS_slopefCut_(0.6),
    IMS_slopegCut_(0.0),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0),
    MC_X_o_cutoff_(0.0),
    MC_X_o_min_(0.0),
    MC_slopepCut_(0.0),
    MC_dpCut_(0.0),
    MC_epCut_(0.0),
    MC_apCut_(0.0),
    MC_bpCut_(0.0),
    MC_cpCut_(0.0),
    CROP_ln_gamma_o_min(-6.0),
    CROP_ln_gamma_o_max(3.0),
    CROP_ln_gamma_k_min(-5.0),
    CROP_ln_gamma_k_max(15.0),
    m_last_is(-1.0),
    m_debugCalc(0)
{
    if (testProb != 1) {
        throw CanteraError("HMWSoln::HMWSoln", "unknown test problem");
    }

    initThermoFile("HMW_NaCl.xml", "");

    size_t i = speciesIndex("Cl-");
    size_t j = speciesIndex("H+");
    size_t n = i * m_kk + j;
    size_t ct = m_CounterIJ[n];
    m_Beta0MX_ij[ct] = 0.1775;
    m_Beta1MX_ij[ct] = 0.2945;
    m_CphiMX_ij[ct]  = 0.0008;
    m_Alpha1MX_ij[ct]= 2.000;


    i = speciesIndex("Cl-");
    j = speciesIndex("Na+");
    n = i * m_kk + j;
    ct = m_CounterIJ[n];
    m_Beta0MX_ij[ct] = 0.0765;
    m_Beta1MX_ij[ct] = 0.2664;
    m_CphiMX_ij[ct]  = 0.00127;
    m_Alpha1MX_ij[ct]= 2.000;


    i = speciesIndex("Cl-");
    j = speciesIndex("OH-");
    n = i * m_kk + j;
    ct = m_CounterIJ[n];
    m_Theta_ij[ct] = -0.05;

    i = speciesIndex("H+");
    j = speciesIndex("Na+");
    n = i * m_kk + j;
    ct = m_CounterIJ[n];
    m_Theta_ij[ct] = 0.036;

    i = speciesIndex("Na+");
    j = speciesIndex("OH-");
    n = i * m_kk + j;
    ct = m_CounterIJ[n];
    m_Beta0MX_ij[ct] = 0.0864;
    m_Beta1MX_ij[ct] = 0.253;
    m_CphiMX_ij[ct]  = 0.0044;
    m_Alpha1MX_ij[ct]= 2.000;

    i = speciesIndex("Cl-");
    j = speciesIndex("H+");
    size_t k = speciesIndex("Na+");
    double param = -0.004;
    n = i * m_kk *m_kk + j * m_kk + k ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;
    n = i * m_kk *m_kk + k * m_kk + j ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;
    n = j * m_kk *m_kk + i * m_kk + k ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;
    n = j * m_kk *m_kk + k * m_kk + i ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;
    n = k * m_kk *m_kk + j * m_kk + i ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;
    n = k * m_kk *m_kk + i * m_kk + j ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;

    i = speciesIndex("Cl-");
    j = speciesIndex("Na+");
    k = speciesIndex("OH-");
    param = -0.006;
    n = i * m_kk *m_kk + j * m_kk + k ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;
    n = i * m_kk *m_kk + k * m_kk + j ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;
    n = j * m_kk *m_kk + i * m_kk + k ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;
    n = j * m_kk *m_kk + k * m_kk + i ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;
    n = k * m_kk *m_kk + j * m_kk + i ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;
    n = k * m_kk *m_kk + i * m_kk + j ;
    m_Psi_ijk[n] = param;
    m_Psi_ijk_coeff(0,n) = param;

    printCoeffs();
}

HMWSoln::~HMWSoln()
{
    delete m_waterProps;
    m_waterProps = 0;
}

ThermoPhase* HMWSoln::duplMyselfAsThermoPhase() const
{
    return new HMWSoln(*this);
}

int HMWSoln::eosType() const
{
    int res;
    switch (m_formGC) {
    case 0:
        res = cHMWSoln0;
        break;
    case 1:
        res = cHMWSoln1;
        break;
    case 2:
        res = cHMWSoln2;
        break;
    default:
        throw CanteraError("eosType", "Unknown type");
    }
    return res;
}

//
// -------- Molar Thermodynamic Properties of the Solution ---------------
//
doublereal HMWSoln::enthalpy_mole() const
{
    getPartialMolarEnthalpies(DATA_PTR(m_tmpV));
    getMoleFractions(DATA_PTR(m_pp));
    return mean_X(DATA_PTR(m_tmpV));
}

doublereal HMWSoln::relative_enthalpy() const
{
    getPartialMolarEnthalpies(DATA_PTR(m_tmpV));
    double hbar = mean_X(DATA_PTR(m_tmpV));
    getEnthalpy_RT(DATA_PTR(m_gamma_tmp));
    double RT = GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        m_gamma_tmp[k] *= RT;
    }
    double h0bar = mean_X(DATA_PTR(m_gamma_tmp));
    return hbar - h0bar;
}

doublereal HMWSoln::relative_molal_enthalpy() const
{
    double L = relative_enthalpy();
    getMoleFractions(DATA_PTR(m_tmpV));
    double xanion = 0.0;
    size_t kcation = npos;
    double xcation = 0.0;
    size_t kanion = npos;
    for (size_t k = 0; k < m_kk; k++) {
        if (charge(k) > 0.0) {
            if (m_tmpV[k] > xanion) {
                xanion = m_tmpV[k];
                kanion = k;
            }
        } else if (charge(k) < 0.0) {
            if (m_tmpV[k] > xcation) {
                xcation = m_tmpV[k];
                kcation = k;
            }
        }
    }
    if (kcation == npos || kanion == npos) {
        return L;
    }
    double xuse = xcation;
    double factor = 1;
    if (xanion < xcation) {
        xuse = xanion;
        if (charge(kcation) != 1.0) {
            factor = charge(kcation);
        }
    } else {
        if (charge(kanion) != 1.0) {
            factor = charge(kanion);
        }
    }
    xuse = xuse / factor;
    return L / xuse;
}

doublereal HMWSoln::entropy_mole() const
{
    getPartialMolarEntropies(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
}

doublereal HMWSoln::gibbs_mole() const
{
    getChemPotentials(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
}

doublereal HMWSoln::cp_mole() const
{
    getPartialMolarCp(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
}

doublereal HMWSoln::cv_mole() const
{
    double kappa_t = isothermalCompressibility();
    double beta = thermalExpansionCoeff();
    double cp = cp_mole();
    double tt = temperature();
    double molarV = molarVolume();
    return cp - beta * beta * tt * molarV / kappa_t;
}

//
// ------- Mechanical Equation of State Properties ------------------------
//

doublereal HMWSoln::pressure() const
{
    return m_Pcurrent;
}

void HMWSoln::setPressure(doublereal p)
{
    setState_TP(temperature(), p);
}

void HMWSoln::calcDensity()
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if(cached.validate(temperature(), pressure(), stateMFNumber())) {
        return;
    }

    double* vbar = &m_pp[0];
    getPartialMolarVolumes(vbar);
    double* x = &m_tmpV[0];
    getMoleFractions(x);
    doublereal vtotal = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        vtotal += vbar[i] * x[i];
    }
    doublereal dd = meanMolecularWeight() / vtotal;
    Phase::setDensity(dd);
}

double HMWSoln::density() const
{
    return Phase::density();
}

void HMWSoln::setDensity(const doublereal rho)
{
    double dens_old = density();

    if (rho != dens_old) {
        throw CanteraError("HMWSoln::setDensity",
                           "Density is not an independent variable");
    }
}

void HMWSoln::setMolarDensity(const doublereal rho)
{
    throw CanteraError("HMWSoln::setMolarDensity",
                       "Density is not an independent variable");
}

void HMWSoln::setTemperature(const doublereal temp)
{
    setState_TP(temp, m_Pcurrent);
}

void HMWSoln::setState_TP(doublereal temp, doublereal pres)
{
    Phase::setTemperature(temp);
    /*
     * Store the current pressure
     */
    m_Pcurrent = pres;

    /*
     * update the standard state thermo
     * -> This involves calling the water function and setting the pressure
     */
    updateStandardStateThermo();

    /*
     * Store the internal density of the water SS.
     * Note, we would have to do this for all other
     * species if they had pressure dependent properties.
     */
    m_densWaterSS = m_waterSS->density();
    /*
     * Calculate all of the other standard volumes
     * -> note these are constant for now
     */
    calcDensity();
}

//
// ------- Activities and Activity Concentrations
//

void HMWSoln::getActivityConcentrations(doublereal* c) const
{
    double cs_solvent = standardConcentration();
    getActivities(c);
    c[0] *= cs_solvent;
    if (m_kk > 1) {
        double cs_solute = standardConcentration(1);
        for (size_t k = 1; k < m_kk; k++) {
            c[k] *= cs_solute;
        }
    }
}

doublereal HMWSoln::standardConcentration(size_t k) const
{
    getStandardVolumes(DATA_PTR(m_tmpV));
    double mvSolvent = m_tmpV[m_indexSolvent];
    if (k > 0) {
        return m_Mnaught / mvSolvent;
    }
    return 1.0 / mvSolvent;
}

void HMWSoln::getUnitsStandardConc(double* uA, int k, int sizeUA) const
{
    for (int i = 0; i < sizeUA; i++) {
        if (i == 0) {
            uA[0] = 1.0;
        }
        if (i == 1) {
            uA[1] = -int(nDim());
        }
        if (i == 2) {
            uA[2] = 0.0;
        }
        if (i == 3) {
            uA[3] = 0.0;
        }
        if (i == 4) {
            uA[4] = 0.0;
        }
        if (i == 5) {
            uA[5] = 0.0;
        }
    }
}

void HMWSoln::getActivities(doublereal* ac) const
{
    updateStandardStateThermo();
    /*
     * Update the molality array, m_molalities()
     *   This requires an update due to mole fractions
     */
    s_update_lnMolalityActCoeff();
    /*
     * Now calculate the array of activities.
     */
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_indexSolvent) {
            ac[k] = m_molalities[k] * exp(m_lnActCoeffMolal_Scaled[k]);
        }
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    ac[m_indexSolvent] =
        exp(m_lnActCoeffMolal_Scaled[m_indexSolvent]) * xmolSolvent;
}

void HMWSoln::getUnscaledMolalityActivityCoefficients(doublereal* acMolality) const
{
    updateStandardStateThermo();
    A_Debye_TP(-1.0, -1.0);
    s_update_lnMolalityActCoeff();
    std::copy(m_lnActCoeffMolal_Unscaled.begin(), m_lnActCoeffMolal_Unscaled.end(), acMolality);
    for (size_t k = 0; k < m_kk; k++) {
        acMolality[k] = exp(acMolality[k]);
    }
}

//
// ------ Partial Molar Properties of the Solution -----------------
//

void HMWSoln::getChemPotentials(doublereal* mu) const
{
    double xx;
    /*
     * First get the standard chemical potentials in
     * molar form.
     *  -> this requires updates of standard state as a function
     *     of T and P
     */
    getStandardChemPotentials(mu);
    /*
     * Update the activity coefficients
     * This also updates the internal molality array.
     */
    s_update_lnMolalityActCoeff();
    doublereal RT = GasConstant * temperature();
    double xmolSolvent = moleFraction(m_indexSolvent);
    for (size_t k = 0; k < m_kk; k++) {
        if (m_indexSolvent != k) {
            xx = std::max(m_molalities[k], SmallNumber);
            mu[k] += RT * (log(xx) + m_lnActCoeffMolal_Scaled[k]);
        }
    }
    xx = std::max(xmolSolvent, SmallNumber);
    mu[m_indexSolvent] +=
        RT * (log(xx) + m_lnActCoeffMolal_Scaled[m_indexSolvent]);
}

void HMWSoln::getPartialMolarEnthalpies(doublereal* hbar) const
{
    /*
     * Get the nondimensional standard state enthalpies
     */
    getEnthalpy_RT(hbar);
    /*
     * dimensionalize it.
     */
    double T = temperature();
    double RT = GasConstant * T;
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnMolalityActCoeff();
    s_update_dlnMolalityActCoeff_dT();
    double RTT = RT * T;
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= RTT * m_dlnActCoeffMolaldT_Scaled[k];
    }
}

void HMWSoln::getPartialMolarEntropies(doublereal* sbar) const
{
    /*
     * Get the standard state entropies at the temperature
     * and pressure of the solution.
     */
    getEntropy_R(sbar);
    /*
     * Dimensionalize the entropies
     */
    doublereal R = GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= R;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnMolalityActCoeff();
    /*
     * First we will add in the obvious dependence on the T
     * term out front of the log activity term
     */
    doublereal mm;
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_indexSolvent) {
            mm = std::max(SmallNumber, m_molalities[k]);
            sbar[k] -= R * (log(mm) + m_lnActCoeffMolal_Scaled[k]);
        }
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    mm = std::max(SmallNumber, xmolSolvent);
    sbar[m_indexSolvent] -= R *(log(mm) + m_lnActCoeffMolal_Scaled[m_indexSolvent]);
    /*
     * Check to see whether activity coefficients are temperature
     * dependent. If they are, then calculate the their temperature
     * derivatives and add them into the result.
     */
    s_update_dlnMolalityActCoeff_dT();
    double RT = R * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -= RT * m_dlnActCoeffMolaldT_Scaled[k];
    }
}

void HMWSoln::getPartialMolarVolumes(doublereal* vbar) const
{
    /*
     * Get the standard state values in m^3 kmol-1
     */
    getStandardVolumes(vbar);
    /*
     * Update the derivatives wrt the activity coefficients.
     */
    s_update_lnMolalityActCoeff();
    s_update_dlnMolalityActCoeff_dP();
    double T = temperature();
    double RT = GasConstant * T;
    for (size_t k = 0; k < m_kk; k++) {
        vbar[k] += RT * m_dlnActCoeffMolaldP_Scaled[k];
    }
}

void HMWSoln::getPartialMolarCp(doublereal* cpbar) const
{
    /*
     * Get the nondimensional gibbs standard state of the
     * species at the T and P of the solution.
     */
    getCp_R(cpbar);

    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnMolalityActCoeff();
    s_update_dlnMolalityActCoeff_dT();
    s_update_d2lnMolalityActCoeff_dT2();
    double T = temperature();
    double RT = GasConstant * T;
    double RTT = RT * T;
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] -= (2.0 * RT * m_dlnActCoeffMolaldT_Scaled[k] +
                     RTT * m_d2lnActCoeffMolaldT2_Scaled[k]);
    }
}

/*
 *  -------------- Utilities -------------------------------
 */

doublereal HMWSoln::satPressure(doublereal t) {
    double p_old = pressure();
    double t_old = temperature();
    double pres = m_waterSS->satPressure(t);
    /*
     * Set the underlying object back to its original state.
     */
    m_waterSS->setState_TP(t_old, p_old);
    return pres;
}

double HMWSoln::A_Debye_TP(double tempArg, double presArg) const
{
    double T = temperature();
    double A;
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }

    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if(cached.validate(T, P)) {
        return m_A_Debye;
    }

    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        A = m_A_Debye;
        break;
    case A_DEBYE_WATER:
        A = m_waterProps->ADebye(T, P, 0);
        m_A_Debye = A;
        break;
    default:
        throw CanteraError("HMWSoln::A_Debye_TP", "shouldn't be here");
    }
    return A;
}

double HMWSoln::dA_DebyedT_TP(double tempArg, double presArg) const
{
    doublereal T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    doublereal P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }
    doublereal dAdT;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        dAdT = 0.0;
        break;
    case A_DEBYE_WATER:
        dAdT = m_waterProps->ADebye(T, P, 1);
        break;
    default:
        throw CanteraError("HMWSoln::dA_DebyedT_TP", "shouldn't be here");
    }
    return dAdT;
}

double HMWSoln::dA_DebyedP_TP(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }

    double dAdP;
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        dAdP = 0.0;
        break;
    case A_DEBYE_WATER:
        if(cached.validate(T, P)) {
            dAdP = cached.value;
        } else {
            dAdP = m_waterProps->ADebye(T, P, 3);
            cached.value  = dAdP;
        }
        break;
    default:
        throw CanteraError("HMWSoln::dA_DebyedP_TP", "shouldn't be here");
    }
    return dAdP;
}

double HMWSoln::ADebye_L(double tempArg, double presArg) const
{
    double dAdT = dA_DebyedT_TP();
    double dAphidT = dAdT /3.0;
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    return dAphidT * (4.0 * GasConstant * T * T);
}

double HMWSoln::ADebye_V(double tempArg, double presArg) const
{
    double dAdP = dA_DebyedP_TP();
    double dAphidP = dAdP /3.0;
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    return - dAphidP * (4.0 * GasConstant * T);
}

double HMWSoln::ADebye_J(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double A_L = ADebye_L(T, presArg);
    double d2 = d2A_DebyedT2_TP(T, presArg);
    double d2Aphi = d2 / 3.0;
    return 2.0 * A_L / T + 4.0 * GasConstant * T * T *d2Aphi;
}

double HMWSoln::d2A_DebyedT2_TP(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }
    double d2AdT2;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        d2AdT2 = 0.0;
        break;
    case A_DEBYE_WATER:
        d2AdT2 = m_waterProps->ADebye(T, P, 2);
        break;
    default:
        throw CanteraError("HMWSoln::d2A_DebyedT2_TP", "shouldn't be here");
    }
    return d2AdT2;
}

/*
 * ---------- Other Property Functions
 */
double HMWSoln::AionicRadius(int k) const
{
    return m_Aionic[k];
}

/*
 * ------------ Private and Restricted Functions ------------------
 */

void HMWSoln::initLengths()
{
    m_kk = nSpecies();

    /*
     * Resize lengths equal to the number of species in
     * the phase.
     */
    m_electrolyteSpeciesType.resize(m_kk, cEST_polarNeutral);
    m_speciesSize.resize(m_kk);
    m_speciesCharge_Stoich.resize(m_kk, 0.0);
    m_Aionic.resize(m_kk, 0.0);

    m_pp.resize(m_kk, 0.0);
    m_tmpV.resize(m_kk, 0.0);
    m_molalitiesCropped.resize(m_kk, 0.0);

    size_t maxCounterIJlen = 1 + (m_kk-1) * (m_kk-2) / 2;

    /*
     * Figure out the size of the temperature coefficient
     * arrays
     */
    int TCoeffLength = 1;
    if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
        TCoeffLength = 2;
    } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        TCoeffLength = 5;
    }

    m_Beta0MX_ij.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_L.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_P.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    m_Beta1MX_ij.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_L.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_P.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    m_Beta2MX_ij.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_L.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_P.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    m_CphiMX_ij.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_L.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_P.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    m_Alpha1MX_ij.resize(maxCounterIJlen, 2.0);
    m_Alpha2MX_ij.resize(maxCounterIJlen, 12.0);
    m_Theta_ij.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_L.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_P.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    size_t n = m_kk*m_kk*m_kk;
    m_Psi_ijk.resize(n, 0.0);
    m_Psi_ijk_L.resize(n, 0.0);
    m_Psi_ijk_LL.resize(n, 0.0);
    m_Psi_ijk_P.resize(n, 0.0);
    m_Psi_ijk_coeff.resize(TCoeffLength, n, 0.0);

    m_Lambda_nj.resize(m_kk, m_kk, 0.0);
    m_Lambda_nj_L.resize(m_kk, m_kk, 0.0);
    m_Lambda_nj_LL.resize(m_kk, m_kk, 0.0);
    m_Lambda_nj_P.resize(m_kk, m_kk, 0.0);
    m_Lambda_nj_coeff.resize(TCoeffLength, m_kk * m_kk, 0.0);

    m_Mu_nnn.resize(m_kk, 0.0);
    m_Mu_nnn_L.resize(m_kk, 0.0);
    m_Mu_nnn_LL.resize(m_kk, 0.0);
    m_Mu_nnn_P.resize(m_kk, 0.0);
    m_Mu_nnn_coeff.resize(TCoeffLength, m_kk, 0.0);

    m_lnActCoeffMolal_Scaled.resize(m_kk, 0.0);
    m_dlnActCoeffMolaldT_Scaled.resize(m_kk, 0.0);
    m_d2lnActCoeffMolaldT2_Scaled.resize(m_kk, 0.0);
    m_dlnActCoeffMolaldP_Scaled.resize(m_kk, 0.0);

    m_lnActCoeffMolal_Unscaled.resize(m_kk, 0.0);
    m_dlnActCoeffMolaldT_Unscaled.resize(m_kk, 0.0);
    m_d2lnActCoeffMolaldT2_Unscaled.resize(m_kk, 0.0);
    m_dlnActCoeffMolaldP_Unscaled.resize(m_kk, 0.0);

    m_CounterIJ.resize(m_kk*m_kk, 0);

    m_gfunc_IJ.resize(maxCounterIJlen, 0.0);
    m_g2func_IJ.resize(maxCounterIJlen, 0.0);
    m_hfunc_IJ.resize(maxCounterIJlen, 0.0);
    m_h2func_IJ.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ_P.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ_P.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ_P.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ_L.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ_P.resize(maxCounterIJlen, 0.0);
    m_Phiprime_IJ.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ_L.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ_P.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ_P.resize(maxCounterIJlen, 0.0);

    m_gamma_tmp.resize(m_kk, 0.0);

    IMS_lnActCoeffMolal_.resize(m_kk, 0.0);
    CROP_speciesCropped_.resize(m_kk, 0);

    counterIJ_setup();
}

void HMWSoln::s_update_lnMolalityActCoeff() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if( cached.validate(temperature(), pressure(), stateMFNumber()) ) {
        return;
    }

    /*
     * Calculate the molalities. Currently, the molalities
     * may not be current with respect to the contents of the
     * State objects' data.
     */
    calcMolalities();
    /*
     *  Calculate a cropped set of molalities that will be used
     *  in all activity coefficient calculations.
     */
    calcMolalitiesCropped();

    /*
     * Calculate the stoichiometric ionic charge. This isn't used in the
     * Pitzer formulation.
     */
    m_IionicMolalityStoich = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        double z_k = charge(k);
        double zs_k1 =  m_speciesCharge_Stoich[k];
        if (z_k == zs_k1) {
            m_IionicMolalityStoich += m_molalities[k] * z_k * z_k;
        } else {
            double zs_k2 = z_k - zs_k1;
            m_IionicMolalityStoich
            += m_molalities[k] * (zs_k1 * zs_k1 + zs_k2 * zs_k2);
        }
    }

    /*
     * Update the temperature dependence of the pitzer coefficients
     * and their derivatives
     */
    s_updatePitzer_CoeffWRTemp();

    /*
     * Calculate the IMS cutoff factors
     */
    s_updateIMS_lnMolalityActCoeff();

    /*
     * Now do the main calculation.
     */
    s_updatePitzer_lnMolalityActCoeff();

    double xmolSolvent = moleFraction(m_indexSolvent);
    double xx = std::max(m_xmolSolventMIN, xmolSolvent);
    double lnActCoeffMolal0 = - log(xx) + (xx - 1.0)/xx;
    double lnxs = log(xx);

    for (size_t k = 1; k < m_kk; k++) {
        CROP_speciesCropped_[k] = 0;
        m_lnActCoeffMolal_Unscaled[k] += IMS_lnActCoeffMolal_[k];
        if (m_lnActCoeffMolal_Unscaled[k] > (CROP_ln_gamma_k_max- 2.5 *lnxs)) {
            CROP_speciesCropped_[k] = 2;
            m_lnActCoeffMolal_Unscaled[k] = CROP_ln_gamma_k_max - 2.5 * lnxs;
        }
        if (m_lnActCoeffMolal_Unscaled[k] < (CROP_ln_gamma_k_min - 2.5 *lnxs)) {
            // -1.0 and -1.5 caused multiple solutions
            CROP_speciesCropped_[k] = 2;
            m_lnActCoeffMolal_Unscaled[k] = CROP_ln_gamma_k_min - 2.5 * lnxs;
        }
    }
    CROP_speciesCropped_[0] = 0;
    m_lnActCoeffMolal_Unscaled[0] += (IMS_lnActCoeffMolal_[0] - lnActCoeffMolal0);
    if (m_lnActCoeffMolal_Unscaled[0] < CROP_ln_gamma_o_min) {
        CROP_speciesCropped_[0] = 2;
        m_lnActCoeffMolal_Unscaled[0] = CROP_ln_gamma_o_min;
    }
    if (m_lnActCoeffMolal_Unscaled[0] > CROP_ln_gamma_o_max) {
        CROP_speciesCropped_[0] = 2;
        // -0.5 caused multiple solutions
        m_lnActCoeffMolal_Unscaled[0] = CROP_ln_gamma_o_max;
    }
    if (m_lnActCoeffMolal_Unscaled[0] > CROP_ln_gamma_o_max - 0.5 * lnxs) {
        CROP_speciesCropped_[0] = 2;
        m_lnActCoeffMolal_Unscaled[0] = CROP_ln_gamma_o_max - 0.5 * lnxs;
    }

    /*
     * Now do the pH Scaling
     */
    s_updateScaling_pHScaling();
}

void HMWSoln::calcMolalitiesCropped() const
{
    doublereal Imax = 0.0, Itmp;
    doublereal Iac_max;
    m_molalitiesAreCropped = false;

    for (size_t k = 0; k < m_kk; k++) {
        m_molalitiesCropped[k] = m_molalities[k];
        Imax = std::max(m_molalities[k] * charge(k) * charge(k), Imax);
    }

    int cropMethod = 1;


    if (cropMethod == 0) {

        /*
         * Quick return
         */
        if (Imax < m_maxIionicStrength) {
            return;
        }

        m_molalitiesAreCropped = true;

        for (size_t i = 1; i < (m_kk - 1); i++) {
            double charge_i = charge(i);
            double abs_charge_i = fabs(charge_i);
            if (charge_i == 0.0) {
                continue;
            }
            for (size_t j = (i+1); j < m_kk; j++) {
                double charge_j = charge(j);
                double abs_charge_j = fabs(charge_j);
                /*
                 * Only loop over oppositely charge species
                 */
                if (charge_i * charge_j < 0) {
                    Iac_max = m_maxIionicStrength;

                    if (m_molalitiesCropped[i] > m_molalitiesCropped[j]) {
                        Imax = m_molalitiesCropped[i] * abs_charge_i * abs_charge_i;
                        if (Imax > Iac_max) {
                            m_molalitiesCropped[i] = Iac_max / (abs_charge_i * abs_charge_i);
                        }
                        Imax = m_molalitiesCropped[j] * fabs(abs_charge_j * abs_charge_i);
                        if (Imax > Iac_max) {
                            m_molalitiesCropped[j] = Iac_max / (abs_charge_j * abs_charge_i);
                        }
                    } else {
                        Imax = m_molalitiesCropped[j] * abs_charge_j * abs_charge_j;
                        if (Imax > Iac_max) {
                            m_molalitiesCropped[j] = Iac_max / (abs_charge_j * abs_charge_j);
                        }
                        Imax = m_molalitiesCropped[i] * abs_charge_j * abs_charge_i;
                        if (Imax > Iac_max) {
                            m_molalitiesCropped[i] = Iac_max / (abs_charge_j * abs_charge_i);
                        }
                    }
                }
            }
        }

        /*
         * Do this loop 10 times until we have achieved charge neutrality
         * in the cropped molalities
         */
        for (int times = 0; times< 10; times++) {
            double anion_charge = 0.0;
            double cation_charge = 0.0;
            size_t anion_contrib_max_i = npos;
            double anion_contrib_max = -1.0;
            size_t cation_contrib_max_i = npos;
            double cation_contrib_max = -1.0;
            for (size_t i = 0; i < m_kk; i++) {
                double charge_i = charge(i);
                if (charge_i < 0.0) {
                    double anion_contrib =  - m_molalitiesCropped[i] * charge_i;
                    anion_charge += anion_contrib ;
                    if (anion_contrib > anion_contrib_max) {
                        anion_contrib_max = anion_contrib;
                        anion_contrib_max_i = i;
                    }
                } else if (charge_i > 0.0) {
                    double cation_contrib = m_molalitiesCropped[i] * charge_i;
                    cation_charge += cation_contrib ;
                    if (cation_contrib > cation_contrib_max) {
                        cation_contrib_max = cation_contrib;
                        cation_contrib_max_i = i;
                    }
                }
            }
            double total_charge = cation_charge - anion_charge;
            if (total_charge > 1.0E-8) {
                double desiredCrop = total_charge/charge(cation_contrib_max_i);
                double maxCrop =  0.66 * m_molalitiesCropped[cation_contrib_max_i];
                if (desiredCrop < maxCrop) {
                    m_molalitiesCropped[cation_contrib_max_i] -= desiredCrop;
                    break;
                } else {
                    m_molalitiesCropped[cation_contrib_max_i] -= maxCrop;
                }
            } else if (total_charge < -1.0E-8) {
                double desiredCrop = total_charge/charge(anion_contrib_max_i);
                double maxCrop =  0.66 * m_molalitiesCropped[anion_contrib_max_i];
                if (desiredCrop < maxCrop) {
                    m_molalitiesCropped[anion_contrib_max_i] -= desiredCrop;
                    break;
                } else {
                    m_molalitiesCropped[anion_contrib_max_i] -= maxCrop;
                }
            } else {
                break;
            }
        }
    }

    if (cropMethod == 1) {
        double* molF = DATA_PTR(m_gamma_tmp);
        getMoleFractions(molF);
        double xmolSolvent = molF[m_indexSolvent];
        if (xmolSolvent >= MC_X_o_cutoff_) {
            return;
        }

        m_molalitiesAreCropped = true;

        double poly = MC_apCut_ + MC_bpCut_ * xmolSolvent + MC_dpCut_* xmolSolvent * xmolSolvent;
        double  p =  xmolSolvent + MC_epCut_ + exp(- xmolSolvent/ MC_cpCut_) * poly;
        double denomInv = 1.0/ (m_Mnaught * p);

        for (size_t k = 0; k < m_kk; k++) {
            m_molalitiesCropped[k] = molF[k] * denomInv ;
        }

        // Do a further check to see if the Ionic strength is below a max value
        // Reduce the molalities to enforce this. Note, this algorithm preserves
        // the charge neutrality of the solution after cropping.
        Itmp = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            Itmp += m_molalitiesCropped[k] * charge(k) * charge(k);
        }
        if (Itmp > m_maxIionicStrength) {
            double ratio = Itmp / m_maxIionicStrength;
            for (size_t k = 0; k < m_kk; k++) {
                if (charge(k) != 0.0) {
                    m_molalitiesCropped[k] *= ratio;
                }
            }
        }
    }

}

void HMWSoln::counterIJ_setup(void) const
{
    size_t n, nc, i, j;
    m_CounterIJ.resize(m_kk * m_kk);
    int counter = 0;
    for (i = 0; i < m_kk; i++) {
        n = i;
        nc = m_kk * i;
        m_CounterIJ[n] = 0;
        m_CounterIJ[nc] = 0;
    }
    for (i = 1; i < (m_kk - 1); i++) {
        n = m_kk * i + i;
        m_CounterIJ[n] = 0;
        for (j = (i+1); j < m_kk; j++) {
            n = m_kk * j + i;
            nc = m_kk * i + j;
            counter++;
            m_CounterIJ[n] = counter;
            m_CounterIJ[nc] = counter;
        }
    }
}

void HMWSoln::s_updatePitzer_CoeffWRTemp(int doDerivs) const
{

    size_t i, j, n, counterIJ;
    const double* beta0MX_coeff;
    const double* beta1MX_coeff;
    const double* beta2MX_coeff;
    const double* CphiMX_coeff;
    const double* Theta_coeff;
    double T = temperature();
    const double twoT = 2.0 * T;
    const double invT = 1.0 / T;
    const double invT2 = invT * invT;
    const double twoinvT3 = 2.0 * invT * invT2;
    double Tr = m_TempPitzerRef;
    double tinv = 0.0, tln = 0.0, tlin = 0.0, tquad = 0.0;
    if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
        tlin = T - Tr;
    } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        tlin = T - Tr;
        tquad = T * T - Tr * Tr;
        tln = log(T/ Tr);
        tinv = 1.0/T - 1.0/Tr;
    }

    for (i = 1; i < (m_kk - 1); i++) {
        for (j = (i+1); j < m_kk; j++) {

            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];

            beta0MX_coeff = m_Beta0MX_ij_coeff.ptrColumn(counterIJ);
            beta1MX_coeff = m_Beta1MX_ij_coeff.ptrColumn(counterIJ);
            beta2MX_coeff = m_Beta2MX_ij_coeff.ptrColumn(counterIJ);
            CphiMX_coeff = m_CphiMX_ij_coeff.ptrColumn(counterIJ);
            Theta_coeff = m_Theta_ij_coeff.ptrColumn(counterIJ);

            switch (m_formPitzerTemp) {
            case PITZER_TEMP_CONSTANT:
                break;
            case PITZER_TEMP_LINEAR:

                m_Beta0MX_ij[counterIJ] = beta0MX_coeff[0]
                                          + beta0MX_coeff[1]*tlin;
                m_Beta0MX_ij_L[counterIJ] = beta0MX_coeff[1];
                m_Beta0MX_ij_LL[counterIJ] = 0.0;

                m_Beta1MX_ij[counterIJ]   = beta1MX_coeff[0]
                                            + beta1MX_coeff[1]*tlin;
                m_Beta1MX_ij_L[counterIJ] = beta1MX_coeff[1];
                m_Beta1MX_ij_LL[counterIJ] = 0.0;

                m_Beta2MX_ij[counterIJ]    = beta2MX_coeff[0]
                                             + beta2MX_coeff[1]*tlin;
                m_Beta2MX_ij_L[counterIJ]  = beta2MX_coeff[1];
                m_Beta2MX_ij_LL[counterIJ] = 0.0;

                m_CphiMX_ij[counterIJ]     = CphiMX_coeff[0]
                                             + CphiMX_coeff[1]*tlin;
                m_CphiMX_ij_L[counterIJ]   = CphiMX_coeff[1];
                m_CphiMX_ij_LL[counterIJ]  = 0.0;

                m_Theta_ij[counterIJ]      = Theta_coeff[0] + Theta_coeff[1]*tlin;
                m_Theta_ij_L[counterIJ]    = Theta_coeff[1];
                m_Theta_ij_LL[counterIJ]   = 0.0;

                break;

            case PITZER_TEMP_COMPLEX1:
                m_Beta0MX_ij[counterIJ] = beta0MX_coeff[0]
                                          + beta0MX_coeff[1]*tlin
                                          + beta0MX_coeff[2]*tquad
                                          + beta0MX_coeff[3]*tinv
                                          + beta0MX_coeff[4]*tln;

                m_Beta1MX_ij[counterIJ] = beta1MX_coeff[0]
                                          + beta1MX_coeff[1]*tlin
                                          + beta1MX_coeff[2]*tquad
                                          + beta1MX_coeff[3]*tinv
                                          + beta1MX_coeff[4]*tln;

                m_Beta2MX_ij[counterIJ] = beta2MX_coeff[0]
                                          + beta2MX_coeff[1]*tlin
                                          + beta2MX_coeff[2]*tquad
                                          + beta2MX_coeff[3]*tinv
                                          + beta2MX_coeff[4]*tln;

                m_CphiMX_ij[counterIJ] = CphiMX_coeff[0]
                                         + CphiMX_coeff[1]*tlin
                                         + CphiMX_coeff[2]*tquad
                                         + CphiMX_coeff[3]*tinv
                                         + CphiMX_coeff[4]*tln;

                m_Theta_ij[counterIJ] = Theta_coeff[0]
                                        + Theta_coeff[1]*tlin
                                        + Theta_coeff[2]*tquad
                                        + Theta_coeff[3]*tinv
                                        + Theta_coeff[4]*tln;

                m_Beta0MX_ij_L[counterIJ] =  beta0MX_coeff[1]
                                             + beta0MX_coeff[2]*twoT
                                             - beta0MX_coeff[3]*invT2
                                             + beta0MX_coeff[4]*invT;

                m_Beta1MX_ij_L[counterIJ] =  beta1MX_coeff[1]
                                             + beta1MX_coeff[2]*twoT
                                             - beta1MX_coeff[3]*invT2
                                             + beta1MX_coeff[4]*invT;

                m_Beta2MX_ij_L[counterIJ] =  beta2MX_coeff[1]
                                             + beta2MX_coeff[2]*twoT
                                             - beta2MX_coeff[3]*invT2
                                             + beta2MX_coeff[4]*invT;

                m_CphiMX_ij_L[counterIJ] =  CphiMX_coeff[1]
                                            + CphiMX_coeff[2]*twoT
                                            - CphiMX_coeff[3]*invT2
                                            + CphiMX_coeff[4]*invT;

                m_Theta_ij_L[counterIJ] =   Theta_coeff[1]
                                            + Theta_coeff[2]*twoT
                                            - Theta_coeff[3]*invT2
                                            + Theta_coeff[4]*invT;

                doDerivs = 2;
                if (doDerivs > 1) {
                    m_Beta0MX_ij_LL[counterIJ] =
                        + beta0MX_coeff[2]*2.0
                        + beta0MX_coeff[3]*twoinvT3
                        - beta0MX_coeff[4]*invT2;

                    m_Beta1MX_ij_LL[counterIJ] =
                        + beta1MX_coeff[2]*2.0
                        + beta1MX_coeff[3]*twoinvT3
                        - beta1MX_coeff[4]*invT2;

                    m_Beta2MX_ij_LL[counterIJ] =
                        + beta2MX_coeff[2]*2.0
                        + beta2MX_coeff[3]*twoinvT3
                        - beta2MX_coeff[4]*invT2;

                    m_CphiMX_ij_LL[counterIJ] =
                        + CphiMX_coeff[2]*2.0
                        + CphiMX_coeff[3]*twoinvT3
                        - CphiMX_coeff[4]*invT2;

                    m_Theta_ij_LL[counterIJ] =
                        + Theta_coeff[2]*2.0
                        + Theta_coeff[3]*twoinvT3
                        - Theta_coeff[4]*invT2;
                }
                break;
            }
        }
    }

    // Lambda interactions and Mu_nnn
    // i must be neutral for this term to be nonzero. We take advantage of this
    // here to lower the operation count.
    for (i = 1; i < m_kk; i++) {
        if (charge(i) == 0.0) {
            for (j = 1; j < m_kk; j++) {
                n = i * m_kk + j;
                const double* Lambda_coeff = m_Lambda_nj_coeff.ptrColumn(n);
                switch (m_formPitzerTemp) {
                case PITZER_TEMP_CONSTANT:
                    m_Lambda_nj(i,j) = Lambda_coeff[0];
                    break;
                case PITZER_TEMP_LINEAR:
                    m_Lambda_nj(i,j)      = Lambda_coeff[0] + Lambda_coeff[1]*tlin;
                    m_Lambda_nj_L(i,j)    = Lambda_coeff[1];
                    m_Lambda_nj_LL(i,j)   = 0.0;
                    break;
                case PITZER_TEMP_COMPLEX1:
                    m_Lambda_nj(i,j) = Lambda_coeff[0]
                                       + Lambda_coeff[1]*tlin
                                       + Lambda_coeff[2]*tquad
                                       + Lambda_coeff[3]*tinv
                                       + Lambda_coeff[4]*tln;

                    m_Lambda_nj_L(i,j) = Lambda_coeff[1]
                                         + Lambda_coeff[2]*twoT
                                         - Lambda_coeff[3]*invT2
                                         + Lambda_coeff[4]*invT;

                    m_Lambda_nj_LL(i,j) =
                        Lambda_coeff[2]*2.0
                        + Lambda_coeff[3]*twoinvT3
                        - Lambda_coeff[4]*invT2;
                }

                if (j == i) {
                    const double* Mu_coeff = m_Mu_nnn_coeff.ptrColumn(i);
                    switch (m_formPitzerTemp) {
                    case PITZER_TEMP_CONSTANT:
                        m_Mu_nnn[i] = Mu_coeff[0];
                        break;
                    case PITZER_TEMP_LINEAR:
                        m_Mu_nnn[i]      = Mu_coeff[0] + Mu_coeff[1]*tlin;
                        m_Mu_nnn_L[i]    = Mu_coeff[1];
                        m_Mu_nnn_LL[i]   = 0.0;
                        break;
                    case PITZER_TEMP_COMPLEX1:
                        m_Mu_nnn[i] = Mu_coeff[0]
                                      + Mu_coeff[1]*tlin
                                      + Mu_coeff[2]*tquad
                                      + Mu_coeff[3]*tinv
                                      + Mu_coeff[4]*tln;

                        m_Mu_nnn_L[i] = Mu_coeff[1]
                                        + Mu_coeff[2]*twoT
                                        - Mu_coeff[3]*invT2
                                        + Mu_coeff[4]*invT;

                        m_Mu_nnn_LL[i] =
                            Mu_coeff[2]*2.0
                            + Mu_coeff[3]*twoinvT3
                            - Mu_coeff[4]*invT2;
                    }
                }
            }
        }
    }

    switch(m_formPitzerTemp) {
    case PITZER_TEMP_CONSTANT:
      for (i = 1; i < m_kk; i++) {
          for (j = 1; j < m_kk; j++) {
              for (size_t k = 1; k < m_kk; k++) {
                  n = i * m_kk *m_kk + j * m_kk + k ;
                  const double* Psi_coeff = m_Psi_ijk_coeff.ptrColumn(n);
                  m_Psi_ijk[n] = Psi_coeff[0];
              }
          }
      }
      break;
    case PITZER_TEMP_LINEAR:
      for (i = 1; i < m_kk; i++) {
          for (j = 1; j < m_kk; j++) {
              for (size_t k = 1; k < m_kk; k++) {
                  n = i * m_kk *m_kk + j * m_kk + k ;
                  const double* Psi_coeff = m_Psi_ijk_coeff.ptrColumn(n);
                  m_Psi_ijk[n]      = Psi_coeff[0] + Psi_coeff[1]*tlin;
                  m_Psi_ijk_L[n]    = Psi_coeff[1];
                  m_Psi_ijk_LL[n]   = 0.0;
              }
          }
      }
      break;
    case PITZER_TEMP_COMPLEX1:
      for (i = 1; i < m_kk; i++) {
          for (j = 1; j < m_kk; j++) {
              for (size_t k = 1; k < m_kk; k++) {
                  n = i * m_kk *m_kk + j * m_kk + k ;
                  const double* Psi_coeff = m_Psi_ijk_coeff.ptrColumn(n);
                  m_Psi_ijk[n] = Psi_coeff[0]
                                 + Psi_coeff[1]*tlin
                                 + Psi_coeff[2]*tquad
                                 + Psi_coeff[3]*tinv
                                 + Psi_coeff[4]*tln;

                  m_Psi_ijk_L[n] = Psi_coeff[1]
                                   + Psi_coeff[2]*twoT
                                   - Psi_coeff[3]*invT2
                                   + Psi_coeff[4]*invT;

                  m_Psi_ijk_LL[n] =
                      Psi_coeff[2]*2.0
                      + Psi_coeff[3]*twoinvT3
                      - Psi_coeff[4]*invT2;
              }
          }
      }
      break;
    }

}

void HMWSoln::s_updatePitzer_lnMolalityActCoeff() const
{
    /*
     * HKM -> Assumption is made that the solvent is
     *        species 0.
     */
    if (m_indexSolvent != 0) {
        throw CanteraError("HMWSoln::s_updatePitzer_lnMolalityActCoeff",
                           "Wrong index solvent value!");
    }

    std::string sni,  snj, snk;

    /*
     * Use the CROPPED molality of the species in solution.
     */
    const double* molality = DATA_PTR(m_molalitiesCropped);

    /*
     * These are data inputs about the Pitzer correlation. They come
     * from the input file for the Pitzer model.
     */
    const double* beta0MX =  DATA_PTR(m_Beta0MX_ij);
    const double* beta1MX =  DATA_PTR(m_Beta1MX_ij);
    const double* beta2MX =  DATA_PTR(m_Beta2MX_ij);
    const double* CphiMX  =  DATA_PTR(m_CphiMX_ij);
    const double* thetaij =  DATA_PTR(m_Theta_ij);
    const double* alpha1MX =  DATA_PTR(m_Alpha1MX_ij);
    const double* alpha2MX =  DATA_PTR(m_Alpha2MX_ij);

    const double* psi_ijk =  DATA_PTR(m_Psi_ijk);
    //n = k + j * m_kk + i * m_kk * m_kk;


    double* gamma_Unscaled = DATA_PTR(m_gamma_tmp);
    /*
     * Local variables defined by Coltrin
     */
    double etheta[5][5], etheta_prime[5][5], sqrtIs;
    /*
     * Molality based ionic strength of the solution
     */
    double Is = 0.0;
    /*
     * Molarcharge of the solution: In Pitzer's notation,
     * this is his variable called "Z".
     */
    double molarcharge = 0.0;
    /*
     * molalitysum is the sum of the molalities over all solutes,
     * even those with zero charge.
     */
    double molalitysumUncropped = 0.0;

    double* gfunc    =  DATA_PTR(m_gfunc_IJ);
    double* g2func   =  DATA_PTR(m_g2func_IJ);
    double* hfunc    =  DATA_PTR(m_hfunc_IJ);
    double* h2func   =  DATA_PTR(m_h2func_IJ);
    double* BMX      =  DATA_PTR(m_BMX_IJ);
    double* BprimeMX =  DATA_PTR(m_BprimeMX_IJ);
    double* BphiMX   =  DATA_PTR(m_BphiMX_IJ);
    double* Phi      =  DATA_PTR(m_Phi_IJ);
    double* Phiprime =  DATA_PTR(m_Phiprime_IJ);
    double* Phiphi   =  DATA_PTR(m_PhiPhi_IJ);
    double* CMX      =  DATA_PTR(m_CMX_IJ);


    double x1, x2;
    double Aphi, F, zsqF;
    double sum1, sum2, sum3, sum4, sum5, term1;
    double sum_m_phi_minus_1, osmotic_coef, lnwateract;

    int z1, z2;
    size_t n, i, j, m, counterIJ,  counterIJ2;

    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf("\n Debugging information from hmw_act \n");
    }
    /*
     * Make sure the counter variables are setup
     */
    counterIJ_setup();

    /*
     * ---------- Calculate common sums over solutes ---------------------
     */
    for (n = 1; n < m_kk; n++) {
        //      ionic strength
        Is += charge(n) * charge(n) * molality[n];
        //      total molar charge
        molarcharge +=  fabs(charge(n)) * molality[n];
        molalitysumUncropped += m_molalities[n];
    }
    Is *= 0.5;

    /*
     * Store the ionic molality in the object for reference.
     */
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 1: \n");
        printf(" ionic strenth      = %14.7le \n total molar "
               "charge = %14.7le \n", Is, molarcharge);
    }

    /*
     * The following call to calc_lambdas() calculates all 16 elements
     * of the elambda and elambda1 arrays, given the value of the
     * ionic strength (Is)
     */
    calc_lambdas(Is);

    /*
     * ----- Step 2:  Find the coefficients E-theta and -------------------
     *                E-thetaprime for all combinations of positive
     *                unlike charges up to 4
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 2: \n");
    }
    for (z1 = 1; z1 <=4; z1++) {
        for (z2 =1; z2 <=4; z2++) {
            calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf(" z1=%3d z2=%3d E-theta(I) = %f, E-thetaprime(I) = %f\n",
                       z1, z2, etheta[z1][z2], etheta_prime[z1][z2]);
            }
        }
    }

    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 3: \n");
        printf(" Species          Species            g(x) "
               " hfunc(x)   \n");
    }

    /*
     *  calculate g(x) and hfunc(x) for each cation-anion pair MX
     *   In the original literature, hfunc, was called gprime. However,
     *   it's not the derivative of g(x), so I renamed it.
     */
    for (i = 1; i < (m_kk - 1); i++) {
        for (j = (i+1); j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];

            /*
             * Only loop over oppositely charge species
             */
            if (charge(i)*charge(j) < 0) {
                /*
                 * x is a reduced function variable
                 */
                x1 = sqrtIs * alpha1MX[counterIJ];
                if (x1 > 1.0E-100) {
                    gfunc[counterIJ] =  2.0*(1.0-(1.0 + x1) * exp(-x1)) / (x1 * x1);
                    hfunc[counterIJ] = -2.0 *
                                       (1.0-(1.0 + x1 + 0.5 * x1 * x1) * exp(-x1)) / (x1 * x1);
                } else {
                    gfunc[counterIJ] = 0.0;
                    hfunc[counterIJ] = 0.0;
                }

                if (beta2MX[counterIJ] != 0.0) {
                    x2 = sqrtIs * alpha2MX[counterIJ];
                    if (x2 > 1.0E-100) {
                        g2func[counterIJ] =  2.0*(1.0-(1.0 + x2) * exp(-x2)) / (x2 * x2);
                        h2func[counterIJ] = -2.0 *
                                            (1.0-(1.0 + x2 + 0.5 * x2 * x2) * exp(-x2)) / (x2 * x2);
                    } else {
                        g2func[counterIJ] = 0.0;
                        h2func[counterIJ] = 0.0;
                    }
                }
            } else {
                gfunc[counterIJ] = 0.0;
                hfunc[counterIJ] = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %9.5f %9.5f \n", sni.c_str(), snj.c_str(),
                       gfunc[counterIJ], hfunc[counterIJ]);
            }
        }
    }

    /*
     * --------- SUBSECTION TO CALCULATE BMX, BprimeMX, BphiMX ----------
     * --------- Agrees with Pitzer, Eq. (49), (51), (55)
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 4: \n");
        printf(" Species          Species            BMX    "
               "BprimeMX    BphiMX   \n");
    }

    for (i = 1; i < m_kk - 1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];

            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0.0) {
                BMX[counterIJ]  = beta0MX[counterIJ]
                                  + beta1MX[counterIJ] * gfunc[counterIJ]
                                  + beta2MX[counterIJ] * g2func[counterIJ];

                if (DEBUG_MODE_ENABLED && m_debugCalc) {
                    printf("%d %g: %g %g %g %g\n",
                           (int) counterIJ,  BMX[counterIJ], beta0MX[counterIJ],
                           beta1MX[counterIJ], beta2MX[counterIJ], gfunc[counterIJ]);
                }
                if (Is > 1.0E-150) {
                    BprimeMX[counterIJ] = (beta1MX[counterIJ] * hfunc[counterIJ]/Is +
                                           beta2MX[counterIJ] * h2func[counterIJ]/Is);
                } else {
                    BprimeMX[counterIJ] = 0.0;
                }
                BphiMX[counterIJ]   = BMX[counterIJ] + Is*BprimeMX[counterIJ];
            } else {
                BMX[counterIJ]      = 0.0;
                BprimeMX[counterIJ] = 0.0;
                BphiMX[counterIJ]   = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %11.7f %11.7f %11.7f \n",
                       sni.c_str(), snj.c_str(),
                       BMX[counterIJ], BprimeMX[counterIJ], BphiMX[counterIJ]);
            }
        }
    }

    /*
     * --------- SUBSECTION TO CALCULATE CMX ----------
     * --------- Agrees with Pitzer, Eq. (53).
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 5: \n");
        printf(" Species          Species            CMX \n");
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0.0) {
                CMX[counterIJ] = CphiMX[counterIJ]/
                                 (2.0* sqrt(fabs(charge(i)*charge(j))));
            } else {
                CMX[counterIJ] = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %11.7f \n", sni.c_str(), snj.c_str(),
                       CMX[counterIJ]);
            }
        }
    }

    /*
     * ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
     * --------- Agrees with Pitzer, Eq. 72, 73, 74
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 6: \n");
        printf(" Species          Species            Phi_ij "
               " Phiprime_ij  Phi^phi_ij \n");
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) > 0) {
                z1 = (int) fabs(charge(i));
                z2 = (int) fabs(charge(j));
                Phi[counterIJ] = thetaij[counterIJ] + etheta[z1][z2];
                Phiprime[counterIJ] = etheta_prime[z1][z2];
                Phiphi[counterIJ] = Phi[counterIJ] + Is * Phiprime[counterIJ];
            } else {
                Phi[counterIJ]      = 0.0;
                Phiprime[counterIJ] = 0.0;
                Phiphi[counterIJ]   = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %10.6f %10.6f %10.6f \n",
                       sni.c_str(), snj.c_str(),
                       Phi[counterIJ], Phiprime[counterIJ], Phiphi[counterIJ]);
            }
        }
    }

    /*
     * ------------- SUBSECTION FOR CALCULATION OF F ----------------------
     * ------------ Agrees with Pitzer Eqn. (65) --------------------------
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 7: \n");
    }
    Aphi = A_Debye_TP() / 3.0;
    F = -Aphi * (sqrt(Is) / (1.0 + 1.2*sqrt(Is))
                 + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" initial value of F = %10.6f \n", F);
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0) {
                F = F + molality[i]*molality[j] * BprimeMX[counterIJ];
            }
            /*
             * Both species have a non-zero charge, and they
             * have the same sign
             */
            if (charge(i)*charge(j) > 0) {
                F = F + molality[i]*molality[j] * Phiprime[counterIJ];
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf(" F = %10.6f \n", F);
            }
        }
    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 8: Summing in All Contributions to Activity Coefficients \n");
    }

    for (i = 1; i < m_kk; i++) {

        /*
         * -------- SUBSECTION FOR CALCULATING THE ACTCOEFF FOR CATIONS -----
         * -------- -> equations agree with my notes, Eqn. (118).
         *          -> Equations agree with Pitzer, eqn.(63)
         */
        if (charge(i) > 0.0) {
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf("  Contributions to ln(ActCoeff_%s):\n", sni.c_str());
            }
            // species i is the cation (positive) to calc the actcoeff
            zsqF = charge(i)*charge(i)*F;
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf("      Unary term:                                      z*z*F = %10.5f\n", zsqF);
            }
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            sum4 = 0.0;
            sum5 = 0.0;
            for (j = 1; j < m_kk; j++) {
                /*
                 * Find the counterIJ for the symmetric binary interaction
                 */
                n = m_kk*i + j;
                counterIJ = m_CounterIJ[n];

                if (charge(j) < 0.0) {
                    // sum over all anions
                    sum1 = sum1 + molality[j]*
                           (2.0*BMX[counterIJ] + molarcharge*CMX[counterIJ]);
                    if (DEBUG_MODE_ENABLED && m_debugCalc) {
                        snj = speciesName(j) + ":";
                        printf("      Bin term with %-13s                  2 m_j BMX = %10.5f\n", snj.c_str(),
                               molality[j]*2.0*BMX[counterIJ]);
                        printf("                                                   m_j Z CMX = %10.5f\n",
                               molality[j]* molarcharge*CMX[counterIJ]);
                    }
                    if (j < m_kk-1) {
                        /*
                         * This term is the ternary interaction involving the
                         * non-duplicate sum over double anions, j, k, with
                         * respect to the cation, i.
                         */
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all anions
                            if (charge(k) < 0.0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 = sum3 + molality[j]*molality[k]*psi_ijk[n];
                                if (DEBUG_MODE_ENABLED && m_debugCalc) {
                                    if (psi_ijk[n] != 0.0) {
                                        snj = speciesName(j) + "," + speciesName(k) + ":";
                                        printf("      Psi term on %-16s           m_j m_k psi_ijk = %10.5f\n", snj.c_str(),
                                               molality[j]*molality[k]*psi_ijk[n]);
                                    }
                                }
                            }
                        }
                    }
                }


                if (charge(j) > 0.0) {
                    // sum over all cations
                    if (j != i) {
                        sum2 = sum2 + molality[j]*(2.0*Phi[counterIJ]);
                        if (DEBUG_MODE_ENABLED && m_debugCalc) {
                            if ((molality[j] * Phi[counterIJ])!= 0.0) {
                                snj = speciesName(j) + ":";
                                printf("      Phi term with %-12s                2 m_j Phi_cc = %10.5f\n", snj.c_str(),
                                       molality[j]*(2.0*Phi[counterIJ]));
                            }
                        }
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            // two inner sums over anions

                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 = sum2 + molality[j]*molality[k]*psi_ijk[n];
                            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                                if (psi_ijk[n] != 0.0) {
                                    snj = speciesName(j) + "," + speciesName(k) + ":";
                                    printf("      Psi term on %-16s           m_j m_k psi_ijk = %10.5f\n", snj.c_str(),
                                           molality[j]*molality[k]*psi_ijk[n]);
                                }
                            }
                            /*
                             * Find the counterIJ for the j,k interaction
                             */
                            n = m_kk*j + k;
                            counterIJ2 = m_CounterIJ[n];
                            sum4 = sum4 + (fabs(charge(i))*
                                           molality[j]*molality[k]*CMX[counterIJ2]);
                            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                                if ((molality[j]*molality[k]*CMX[counterIJ2]) != 0.0) {
                                    snj = speciesName(j) + "," + speciesName(k) + ":";
                                    printf("      Tern CMX term on %-16s abs(z_i) m_j m_k CMX = %10.5f\n", snj.c_str(),
                                           fabs(charge(i))* molality[j]*molality[k]*CMX[counterIJ2]);
                                }
                            }
                        }
                    }
                }

                /*
                 * Handle neutral j species
                 */
                if (charge(j) == 0) {
                    sum5 = sum5 + molality[j]*2.0*m_Lambda_nj(j,i);
                    if (DEBUG_MODE_ENABLED && m_debugCalc) {
                        if ((molality[j]*2.0*m_Lambda_nj(j,i)) != 0.0) {
                            snj = speciesName(j) + ":";
                            printf("      Lambda term with %-12s                 2 m_j lam_ji = %10.5f\n", snj.c_str(),
                                   molality[j]*2.0*m_Lambda_nj(j,i));
                        }
                    }
                    /*
                     * Zeta interaction term
                     */
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            size_t izeta = j;
                            size_t jzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + k;
                            double zeta = psi_ijk[n];
                            if (zeta != 0.0) {
                                sum5 = sum5 + molality[j]*molality[k]*zeta;
                                if (DEBUG_MODE_ENABLED && m_debugCalc) {
                                    snj = speciesName(j) + "," + speciesName(k) + ":";
                                    printf("      Zeta term on %-16s         m_n m_a zeta_nMa = %10.5f\n", snj.c_str(),
                                           molality[j]*molality[k]*psi_ijk[n]);
                                }
                            }
                        }
                    }
                }
            }
            /*
             * Add all of the contributions up to yield the log of the
             * solute activity coefficients (molality scale)
             */
            m_lnActCoeffMolal_Unscaled[i] = zsqF + sum1 + sum2 + sum3 + sum4 + sum5;
            gamma_Unscaled[i] = exp(m_lnActCoeffMolal_Unscaled[i]);
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf("      Net %-16s                        lngamma[i] =  %9.5f         gamma[i]=%10.6f \n",
                       sni.c_str(), m_lnActCoeffMolal_Unscaled[i], gamma_Unscaled[i]);
            }
        }

        /*
         * -------- SUBSECTION FOR CALCULATING THE ACTCOEFF FOR ANIONS ------
         * -------- -> equations agree with my notes, Eqn. (119).
         *          -> Equations agree with Pitzer, eqn.(64)
         */
        if (charge(i) < 0) {
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf("  Contributions to ln(ActCoeff_%s):\n", sni.c_str());
            }
            //          species i is an anion (negative)
            zsqF = charge(i)*charge(i)*F;
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf("      Unary term:                                      z*z*F = %10.5f\n", zsqF);
            }
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            sum4 = 0.0;
            sum5 = 0.0;
            for (j = 1; j < m_kk; j++) {
                /*
                 * Find the counterIJ for the symmetric binary interaction
                 */
                n = m_kk*i + j;
                counterIJ = m_CounterIJ[n];

                /*
                 * For Anions, do the cation interactions.
                 */
                if (charge(j) > 0) {
                    sum1 = sum1 + molality[j]*
                           (2.0*BMX[counterIJ]+molarcharge*CMX[counterIJ]);
                    if (DEBUG_MODE_ENABLED && m_debugCalc) {
                        snj = speciesName(j) + ":";
                        printf("      Bin term with %-13s                  2 m_j BMX = %10.5f\n", snj.c_str(),
                               molality[j]*2.0*BMX[counterIJ]);
                        printf("                                                   m_j Z CMX = %10.5f\n",
                               molality[j]* molarcharge*CMX[counterIJ]);
                    }
                    if (j < m_kk-1) {
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all cations
                            if (charge(k) > 0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 = sum3 + molality[j]*molality[k]*psi_ijk[n];
                                if (DEBUG_MODE_ENABLED && m_debugCalc) {
                                    if (psi_ijk[n] != 0.0) {
                                        snj = speciesName(j) + "," + speciesName(k) + ":";
                                        printf("      Psi term on %-16s           m_j m_k psi_ijk = %10.5f\n", snj.c_str(),
                                               molality[j]*molality[k]*psi_ijk[n]);
                                    }
                                }
                            }
                        }
                    }
                }

                /*
                 * For Anions, do the other anion interactions.
                 */
                if (charge(j) < 0.0) {
                    //  sum over all anions
                    if (j != i) {
                        sum2 = sum2 + molality[j]*(2.0*Phi[counterIJ]);
                        if (DEBUG_MODE_ENABLED && m_debugCalc) {
                            if ((molality[j] * Phi[counterIJ])!= 0.0) {
                                snj = speciesName(j) + ":";
                                printf("      Phi term with %-12s                2 m_j Phi_aa = %10.5f\n", snj.c_str(),
                                       molality[j]*(2.0*Phi[counterIJ]));
                            }
                        }
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            // two inner sums over cations
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 = sum2 + molality[j]*molality[k]*psi_ijk[n];
                            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                                if (psi_ijk[n] != 0.0) {
                                    snj = speciesName(j) + "," + speciesName(k) + ":";
                                    printf("      Psi term on %-16s           m_j m_k psi_ijk = %10.5f\n", snj.c_str(),
                                           molality[j]*molality[k]*psi_ijk[n]);
                                }
                            }
                            /*
                             * Find the counterIJ for the symmetric binary interaction
                             */
                            n = m_kk*j + k;
                            counterIJ2 = m_CounterIJ[n];
                            sum4 = sum4 +
                                   (fabs(charge(i))*
                                    molality[j]*molality[k]*CMX[counterIJ2]);
                            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                                if ((molality[j]*molality[k]*CMX[counterIJ2]) != 0.0) {
                                    snj = speciesName(j) + "," + speciesName(k) + ":";
                                    printf("      Tern CMX term on %-16s abs(z_i) m_j m_k CMX = %10.5f\n", snj.c_str(),
                                           fabs(charge(i))* molality[j]*molality[k]*CMX[counterIJ2]);
                                }
                            }
                        }
                    }
                }

                /*
                 * for Anions, do the neutral species interaction
                 */
                if (charge(j) == 0.0) {
                    sum5 = sum5 + molality[j]*2.0*m_Lambda_nj(j,i);
                    if (DEBUG_MODE_ENABLED && m_debugCalc) {
                        if ((molality[j]*2.0*m_Lambda_nj(j,i)) != 0.0) {
                            snj = speciesName(j) + ":";
                            printf("      Lambda term with %-12s                 2 m_j lam_ji = %10.5f\n", snj.c_str(),
                                   molality[j]*2.0*m_Lambda_nj(j,i));
                        }
                    }
                    /*
                     * Zeta interaction term
                     */
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            size_t izeta = j;
                            size_t jzeta = k;
                            size_t kzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + kzeta;
                            double zeta = psi_ijk[n];
                            if (zeta != 0.0) {
                                sum5 = sum5 + molality[j]*molality[k]*zeta;
                                if (DEBUG_MODE_ENABLED && m_debugCalc) {
                                    snj = speciesName(j) + "," + speciesName(k) + ":";
                                    printf("      Zeta term on %-16s         m_n m_c zeta_ncX = %10.5f\n", snj.c_str(),
                                           molality[j]*molality[k]*psi_ijk[n]);
                                }
                            }
                        }
                    }
                }
            }
            m_lnActCoeffMolal_Unscaled[i] = zsqF + sum1 + sum2 + sum3 + sum4 + sum5;
            gamma_Unscaled[i] = exp(m_lnActCoeffMolal_Unscaled[i]);
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf("      Net %-16s                        lngamma[i] =  %9.5f             gamma[i]=%10.6f\n",
                       sni.c_str(), m_lnActCoeffMolal_Unscaled[i], gamma_Unscaled[i]);
            }
        }
        /*
         * ------ SUBSECTION FOR CALCULATING NEUTRAL SOLUTE ACT COEFF -------
         * ------ -> equations agree with my notes,
         *        -> Equations agree with Pitzer,
         */
        if (charge(i) == 0.0) {
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf("  Contributions to ln(ActCoeff_%s):\n", sni.c_str());
            }
            sum1 = 0.0;
            sum3 = 0.0;
            for (j = 1; j < m_kk; j++) {
                sum1 = sum1 + molality[j]*2.0*m_Lambda_nj(i,j);
                if (DEBUG_MODE_ENABLED && m_debugCalc) {
                    if (m_Lambda_nj(i,j) != 0.0) {
                        snj = speciesName(j) + ":";
                        printf("      Lambda_n term on %-16s     2 m_j lambda_n_j = %10.5f\n", snj.c_str(),
                               molality[j]*2.0*m_Lambda_nj(i,j));
                    }
                }
                /*
                 * Zeta term -> we piggyback on the psi term
                 */
                if (charge(j) > 0.0) {
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum3 = sum3 + molality[j]*molality[k]*psi_ijk[n];
                            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                                if (psi_ijk[n] != 0.0) {
                                    snj = speciesName(j) + "," + speciesName(k) + ":";
                                    printf("      Zeta term on %-16s           m_j m_k psi_ijk = %10.5f\n", snj.c_str(),
                                           molality[j]*molality[k]*psi_ijk[n]);
                                }
                            }
                        }
                    }
                }
            }
            sum2 = 3.0 * molality[i]* molality[i] * m_Mu_nnn[i];
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                if (m_Mu_nnn[i] != 0.0) {
                    printf("      Mu_nnn term              3 m_n m_n Mu_n_n = %10.5f\n",
                           3.0 * molality[i]* molality[i] * m_Mu_nnn[i]);
                }
            }
            m_lnActCoeffMolal_Unscaled[i] = sum1 + sum2 + sum3;
            gamma_Unscaled[i] = exp(m_lnActCoeffMolal_Unscaled[i]);
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf("      Net %-16s                        lngamma[i] =  %9.5f             gamma[i]=%10.6f\n",
                       sni.c_str(), m_lnActCoeffMolal_Unscaled[i], gamma_Unscaled[i]);
            }
        }

    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 9: \n");
    }
    /*
     * -------- SUBSECTION FOR CALCULATING THE OSMOTIC COEFF ---------
     * -------- -> equations agree with my notes, Eqn. (117).
     *          -> Equations agree with Pitzer, eqn.(62)
     */
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    double sum6 = 0.0;
    double sum7 = 0.0;
    /*
     * term1 is the DH term in the osmotic coefficient expression
     * b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer
     *                          implementations.
     * Is = Ionic strength on the molality scale (units of (gmol/kg))
     * Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
     */
    term1 = -Aphi * pow(Is,1.5) / (1.0 + 1.2 * sqrt(Is));

    for (j = 1; j < m_kk; j++) {
        /*
         * Loop Over Cations
         */
        if (charge(j) > 0.0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];

                    sum1 = sum1 + molality[j]*molality[k]*
                           (BphiMX[counterIJ] + molarcharge*CMX[counterIJ]);
                }
            }

            for (size_t k = j+1; k < m_kk; k++) {
                if (j == (m_kk-1)) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_lnMolalityActCoeff",
                                       "logic error 1 in Step 9 of hmw_act");
                }
                if (charge(k) > 0.0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     * between 2 cations.
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];
                    sum2 = sum2 + molality[j]*molality[k]*Phiphi[counterIJ];
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) < 0.0) {
                            // species m is an anion
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum2 = sum2 +
                                   molality[j]*molality[k]*molality[m]*psi_ijk[n];
                        }
                    }
                }
            }
        }

        /*
         * Loop Over Anions
         */
        if (charge(j) < 0) {
            for (size_t k = j+1; k < m_kk; k++) {
                if (j == m_kk-1) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_lnMolalityActCoeff",
                                       "logic error 2 in Step 9 of hmw_act");
                }
                if (charge(k) < 0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     * between two anions
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];

                    sum3 = sum3 + molality[j]*molality[k]*Phiphi[counterIJ];
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum3 = sum3 +
                                   molality[j]*molality[k]*molality[m]*psi_ijk[n];
                        }
                    }
                }
            }
        }

        /*
         * Loop Over Neutral Species
         */
        if (charge(j) == 0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    sum4 = sum4 + molality[j]*molality[k]*m_Lambda_nj(j,k);
                }
                if (charge(k) > 0.0) {
                    sum5 = sum5 + molality[j]*molality[k]*m_Lambda_nj(j,k);
                }
                if (charge(k) == 0.0) {
                    if (k > j) {
                        sum6 = sum6 + molality[j]*molality[k]*m_Lambda_nj(j,k);
                    } else if (k == j) {
                        sum6 = sum6 + 0.5 * molality[j]*molality[k]*m_Lambda_nj(j,k);
                    }
                }
                if (charge(k) < 0.0) {
                    size_t izeta = j;
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            size_t jzeta = m;
                            n = k + jzeta * m_kk + izeta * m_kk * m_kk;
                            double zeta = psi_ijk[n];
                            if (zeta != 0.0) {
                                sum7 += molality[izeta]*molality[jzeta]*molality[k]*zeta;
                            }
                        }
                    }
                }
            }
            sum7 += molality[j]*molality[j]*molality[j]*m_Mu_nnn[j];
        }
    }
    sum_m_phi_minus_1 = 2.0 *
                        (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7);
    /*
     * Calculate the osmotic coefficient from
     *       osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
     */
    if (molalitysumUncropped > 1.0E-150) {
        osmotic_coef = 1.0 + (sum_m_phi_minus_1 / molalitysumUncropped);
    } else {
        osmotic_coef = 1.0;
    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" term1=%10.6f sum1=%10.6f sum2=%10.6f "
               "sum3=%10.6f sum4=%10.6f sum5=%10.6f\n",
               term1, sum1, sum2, sum3, sum4, sum5);
        printf("     sum_m_phi_minus_1=%10.6f        osmotic_coef=%10.6f\n",
               sum_m_phi_minus_1, osmotic_coef);
        printf(" Step 10: \n");
    }
    lnwateract = -(m_weightSolvent/1000.0) * molalitysumUncropped * osmotic_coef;

    /*
     * In Cantera, we define the activity coefficient of the solvent as
     *
     *     act_0 = actcoeff_0 * Xmol_0
     *
     * We have just computed act_0. However, this routine returns
     *  ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
     */
    double xmolSolvent = moleFraction(m_indexSolvent);
    double xx = std::max(m_xmolSolventMIN, xmolSolvent);
    m_lnActCoeffMolal_Unscaled[0] = lnwateract - log(xx);
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        double wateract = exp(lnwateract);
        printf(" Weight of Solvent = %16.7g\n", m_weightSolvent);
        printf(" molalitySumUncropped = %16.7g\n", molalitysumUncropped);
        printf(" ln_a_water=%10.6f a_water=%10.6f\n\n",
               lnwateract, wateract);
    }
}

void HMWSoln::s_update_dlnMolalityActCoeff_dT() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if( cached.validate(temperature(), pressure(), stateMFNumber()) ) {
        return;
    }

    /*
     *  Zero the unscaled 2nd derivatives
     */
    m_dlnActCoeffMolaldT_Unscaled.assign(m_kk, 0.0);
    /*
     *  Do the actual calculation of the unscaled temperature derivatives
     */
    s_updatePitzer_dlnMolalityActCoeff_dT();

    for (size_t k = 1; k < m_kk; k++) {
        if (CROP_speciesCropped_[k] == 2) {
            m_dlnActCoeffMolaldT_Unscaled[k] = 0.0;
        }
    }

    if (CROP_speciesCropped_[0]) {
        m_dlnActCoeffMolaldT_Unscaled[0] = 0.0;
    }

    /*
     *  Do the pH scaling to the derivatives
     */
    s_updateScaling_pHScaling_dT();


}

void HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dT() const
{
    /*
     * It may be assumed that the Pitzer activity coefficient routine is
     * called immediately preceding the calling of this routine. Therefore,
     * some quantities do not need to be recalculated in this routine.
     */

    /*
     * HKM -> Assumption is made that the solvent is
     *        species 0.
     */
#ifdef DEBUG_MODE
    m_debugCalc = 0;
#endif
    if (m_indexSolvent != 0) {
        throw CanteraError("HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dT",
                           "Wrong index solvent value!");
    }

    std::string sni, snj, snk;

    const double* molality  =  DATA_PTR(m_molalitiesCropped);
    const double* beta0MX_L =  DATA_PTR(m_Beta0MX_ij_L);
    const double* beta1MX_L =  DATA_PTR(m_Beta1MX_ij_L);
    const double* beta2MX_L =  DATA_PTR(m_Beta2MX_ij_L);
    const double* CphiMX_L  =  DATA_PTR(m_CphiMX_ij_L);
    const double* thetaij_L =  DATA_PTR(m_Theta_ij_L);
    const double* alpha1MX   =  DATA_PTR(m_Alpha1MX_ij);
    const double* alpha2MX   =  DATA_PTR(m_Alpha2MX_ij);
    const double* psi_ijk_L =  DATA_PTR(m_Psi_ijk_L);
    double* d_gamma_dT_Unscaled   =  DATA_PTR(m_gamma_tmp);
    /*
     * Local variables defined by Coltrin
     */
    double etheta[5][5], etheta_prime[5][5], sqrtIs;
    /*
     * Molality based ionic strength of the solution
     */
    double Is = 0.0;
    /*
     * Molarcharge of the solution: In Pitzer's notation,
     * this is his variable called "Z".
     */
    double molarcharge = 0.0;
    /*
     * molalitysum is the sum of the molalities over all solutes,
     * even those with zero charge.
     */
    double molalitysum = 0.0;

    double* gfunc    =  DATA_PTR(m_gfunc_IJ);
    double* g2func   =  DATA_PTR(m_g2func_IJ);
    double* hfunc    =  DATA_PTR(m_hfunc_IJ);
    double* h2func   =  DATA_PTR(m_h2func_IJ);
    double* BMX_L    =  DATA_PTR(m_BMX_IJ_L);
    double* BprimeMX_L= DATA_PTR(m_BprimeMX_IJ_L);
    double* BphiMX_L =  DATA_PTR(m_BphiMX_IJ_L);
    double* Phi_L    =  DATA_PTR(m_Phi_IJ_L);
    double* Phiprime =  DATA_PTR(m_Phiprime_IJ);
    double* Phiphi_L =  DATA_PTR(m_PhiPhi_IJ_L);
    double* CMX_L    =  DATA_PTR(m_CMX_IJ_L);

    double x1, x2;
    double dFdT, zsqdFdT;
    double sum1, sum2, sum3, sum4, sum5, term1;
    double sum_m_phi_minus_1, d_osmotic_coef_dT, d_lnwateract_dT;

    int z1, z2;
    size_t n, i, j, m, counterIJ,  counterIJ2;

    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf("\n Debugging information from "
               "s_Pitzer_dlnMolalityActCoeff_dT()\n");
    }
    /*
     * Make sure the counter variables are setup
     */
    counterIJ_setup();

    /*
     * ---------- Calculate common sums over solutes ---------------------
     */
    for (n = 1; n < m_kk; n++) {
        //      ionic strength
        Is += charge(n) * charge(n) * molality[n];
        //      total molar charge
        molarcharge +=  fabs(charge(n)) * molality[n];
        molalitysum += molality[n];
    }
    Is *= 0.5;

    /*
     * Store the ionic molality in the object for reference.
     */
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 1: \n");
        printf(" ionic strenth      = %14.7le \n total molar "
               "charge = %14.7le \n", Is, molarcharge);
    }

    /*
     * The following call to calc_lambdas() calculates all 16 elements
     * of the elambda and elambda1 arrays, given the value of the
     * ionic strength (Is)
     */
    calc_lambdas(Is);

    /*
     * ----- Step 2:  Find the coefficients E-theta and -------------------
     *                E-thetaprime for all combinations of positive
     *                unlike charges up to 4
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 2: \n");
    }
    for (z1 = 1; z1 <=4; z1++) {
        for (z2 =1; z2 <=4; z2++) {
            calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf(" z1=%3d z2=%3d E-theta(I) = %f, E-thetaprime(I) = %f\n",
                       z1, z2, etheta[z1][z2], etheta_prime[z1][z2]);
            }
        }
    }

    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 3: \n");
        printf(" Species          Species            g(x) "
               " hfunc(x)   \n");
    }

    /*
     *  calculate g(x) and hfunc(x) for each cation-anion pair MX
     *   In the original literature, hfunc, was called gprime. However,
     *   it's not the derivative of g(x), so I renamed it.
     */
    for (i = 1; i < (m_kk - 1); i++) {
        for (j = (i+1); j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * Only loop over oppositely charge species
             */
            if (charge(i)*charge(j) < 0) {
                /*
                 * x is a reduced function variable
                 */
                x1 = sqrtIs * alpha1MX[counterIJ];
                if (x1 > 1.0E-100) {
                    gfunc[counterIJ]     =  2.0*(1.0-(1.0 + x1) * exp(-x1)) / (x1 * x1);
                    hfunc[counterIJ] = -2.0 *
                                       (1.0-(1.0 + x1 + 0.5 * x1 *x1) * exp(-x1)) / (x1 * x1);
                } else {
                    gfunc[counterIJ] = 0.0;
                    hfunc[counterIJ] = 0.0;
                }

                if (beta2MX_L[counterIJ] != 0.0) {
                    x2 = sqrtIs * alpha2MX[counterIJ];
                    if (x2 > 1.0E-100) {
                        g2func[counterIJ] =  2.0*(1.0-(1.0 + x2) * exp(-x2)) / (x2 * x2);
                        h2func[counterIJ] = -2.0 *
                                            (1.0-(1.0 + x2 + 0.5 * x2 * x2) * exp(-x2)) / (x2 * x2);
                    } else {
                        g2func[counterIJ] = 0.0;
                        h2func[counterIJ] = 0.0;
                    }
                }
            } else {
                gfunc[counterIJ] = 0.0;
                hfunc[counterIJ] = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %9.5f %9.5f \n", sni.c_str(), snj.c_str(),
                       gfunc[counterIJ], hfunc[counterIJ]);
            }
        }
    }

    /*
     * ------- SUBSECTION TO CALCULATE BMX_L, BprimeMX_L, BphiMX_L ----------
     * ------- These are now temperature derivatives of the
     *         previously calculated quantities.
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 4: \n");
        printf(" Species          Species            BMX    "
               "BprimeMX    BphiMX   \n");
    }

    for (i = 1; i < m_kk - 1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0.0) {
                BMX_L[counterIJ]  = beta0MX_L[counterIJ]
                                    + beta1MX_L[counterIJ] * gfunc[counterIJ]
                                    + beta2MX_L[counterIJ] * gfunc[counterIJ];
                if (DEBUG_MODE_ENABLED && m_debugCalc) {
                    printf("%d %g: %g %g %g %g\n",
                           (int) counterIJ,  BMX_L[counterIJ], beta0MX_L[counterIJ],
                           beta1MX_L[counterIJ],  beta2MX_L[counterIJ], gfunc[counterIJ]);
                }
                if (Is > 1.0E-150) {
                    BprimeMX_L[counterIJ] = (beta1MX_L[counterIJ] * hfunc[counterIJ]/Is +
                                             beta2MX_L[counterIJ] * h2func[counterIJ]/Is);
                } else {
                    BprimeMX_L[counterIJ] = 0.0;
                }
                BphiMX_L[counterIJ] = BMX_L[counterIJ] + Is*BprimeMX_L[counterIJ];
            } else {
                BMX_L[counterIJ]      = 0.0;
                BprimeMX_L[counterIJ] = 0.0;
                BphiMX_L[counterIJ]     = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %11.7f %11.7f %11.7f \n",
                       sni.c_str(), snj.c_str(),
                       BMX_L[counterIJ], BprimeMX_L[counterIJ], BphiMX_L[counterIJ]);
            }
        }
    }

    /*
     * --------- SUBSECTION TO CALCULATE CMX_L ----------
     * ---------
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 5: \n");
        printf(" Species          Species            CMX \n");
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0.0) {
                CMX_L[counterIJ] = CphiMX_L[counterIJ]/
                                   (2.0* sqrt(fabs(charge(i)*charge(j))));
            } else {
                CMX_L[counterIJ] = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %11.7f \n", sni.c_str(), snj.c_str(),
                       CMX_L[counterIJ]);
            }
        }
    }

    /*
     * ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
     * --------
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 6: \n");
        printf(" Species          Species            Phi_ij "
               " Phiprime_ij  Phi^phi_ij \n");
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) > 0) {
                z1 = (int) fabs(charge(i));
                z2 = (int) fabs(charge(j));
                //Phi[counterIJ] = thetaij_L[counterIJ] + etheta[z1][z2];
                Phi_L[counterIJ] = thetaij_L[counterIJ];
                //Phiprime[counterIJ] = etheta_prime[z1][z2];
                Phiprime[counterIJ] = 0.0;
                Phiphi_L[counterIJ] = Phi_L[counterIJ] + Is * Phiprime[counterIJ];
            } else {
                Phi_L[counterIJ]      = 0.0;
                Phiprime[counterIJ] = 0.0;
                Phiphi_L[counterIJ]   = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %10.6f %10.6f %10.6f \n",
                       sni.c_str(), snj.c_str(),
                       Phi_L[counterIJ], Phiprime[counterIJ], Phiphi_L[counterIJ]);
            }
        }
    }

    /*
     * ----------- SUBSECTION FOR CALCULATION OF dFdT ---------------------
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 7: \n");
    }
    double dA_DebyedT = dA_DebyedT_TP();
    double dAphidT = dA_DebyedT /3.0;
    dFdT = -dAphidT * (sqrt(Is) / (1.0 + 1.2*sqrt(Is))
                       + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" initial value of dFdT = %10.6f \n", dFdT);
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0) {
                dFdT = dFdT + molality[i]*molality[j] * BprimeMX_L[counterIJ];
            }
            /*
             * Both species have a non-zero charge, and they
             * have the same sign, e.g., both positive or both negative.
             */
            if (charge(i)*charge(j) > 0) {
                dFdT = dFdT + molality[i]*molality[j] * Phiprime[counterIJ];
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf(" dFdT = %10.6f \n", dFdT);
            }
        }
    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 8: \n");
    }

    for (i = 1; i < m_kk; i++) {

        /*
         * -------- SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR CATIONS -----
         * --
         */
        if (charge(i) > 0) {
            // species i is the cation (positive) to calc the actcoeff
            zsqdFdT = charge(i)*charge(i)*dFdT;
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            sum4 = 0.0;
            sum5 = 0.0;
            for (j = 1; j < m_kk; j++) {
                /*
                 * Find the counterIJ for the symmetric binary interaction
                 */
                n = m_kk*i + j;
                counterIJ = m_CounterIJ[n];

                if (charge(j) < 0.0) {
                    // sum over all anions
                    sum1 = sum1 + molality[j]*
                           (2.0*BMX_L[counterIJ] + molarcharge*CMX_L[counterIJ]);
                    if (j < m_kk-1) {
                        /*
                         * This term is the ternary interaction involving the
                         * non-duplicate sum over double anions, j, k, with
                         * respect to the cation, i.
                         */
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all anions
                            if (charge(k) < 0.0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 = sum3 + molality[j]*molality[k]*psi_ijk_L[n];
                            }
                        }
                    }
                }


                if (charge(j) > 0.0) {
                    // sum over all cations
                    if (j != i) {
                        sum2 = sum2 + molality[j]*(2.0*Phi_L[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            // two inner sums over anions

                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 = sum2 + molality[j]*molality[k]*psi_ijk_L[n];
                            /*
                             * Find the counterIJ for the j,k interaction
                             */
                            n = m_kk*j + k;
                            counterIJ2 = m_CounterIJ[n];
                            sum4 = sum4 + (fabs(charge(i))*
                                           molality[j]*molality[k]*CMX_L[counterIJ2]);
                        }
                    }
                }

                /*
                 * Handle neutral j species
                 */
                if (charge(j) == 0) {
                    sum5 = sum5 + molality[j]*2.0*m_Lambda_nj_L(j,i);
                }
                /*
                 * Zeta interaction term
                 */
                for (size_t k = 1; k < m_kk; k++) {
                    if (charge(k) < 0.0) {
                        size_t izeta = j;
                        size_t jzeta = i;
                        n = izeta * m_kk * m_kk + jzeta * m_kk + k;
                        double zeta_L = psi_ijk_L[n];
                        if (zeta_L != 0.0) {
                            sum5 = sum5 + molality[j]*molality[k]*zeta_L;
                        }
                    }
                }
            }
            /*
             * Add all of the contributions up to yield the log of the
             * solute activity coefficients (molality scale)
             */
            m_dlnActCoeffMolaldT_Unscaled[i] =
                zsqdFdT + sum1 + sum2 + sum3 + sum4 + sum5;
            d_gamma_dT_Unscaled[i] = exp(m_dlnActCoeffMolaldT_Unscaled[i]);
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf(" %-16s lngamma[i]=%10.6f gamma[i]=%10.6f \n",
                       sni.c_str(), m_dlnActCoeffMolaldT_Unscaled[i], d_gamma_dT_Unscaled[i]);
                printf("                   %12g %12g %12g %12g %12g %12g\n",
                       zsqdFdT, sum1, sum2, sum3, sum4, sum5);
            }
        }

        /*
         * ------ SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR ANIONS ------
         *
         */
        if (charge(i) < 0) {
            //          species i is an anion (negative)
            zsqdFdT = charge(i)*charge(i)*dFdT;
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            sum4 = 0.0;
            sum5 = 0.0;
            for (j = 1; j < m_kk; j++) {
                /*
                 * Find the counterIJ for the symmetric binary interaction
                 */
                n = m_kk*i + j;
                counterIJ = m_CounterIJ[n];

                /*
                 * For Anions, do the cation interactions.
                 */
                if (charge(j) > 0) {
                    sum1 = sum1 + molality[j]*
                           (2.0*BMX_L[counterIJ] + molarcharge*CMX_L[counterIJ]);
                    if (j < m_kk-1) {
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all cations
                            if (charge(k) > 0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 = sum3 + molality[j]*molality[k]*psi_ijk_L[n];
                            }
                        }
                    }
                }

                /*
                 * For Anions, do the other anion interactions.
                 */
                if (charge(j) < 0.0) {
                    //  sum over all anions
                    if (j != i) {
                        sum2 = sum2 + molality[j]*(2.0*Phi_L[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            // two inner sums over cations
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 = sum2 + molality[j]*molality[k]*psi_ijk_L[n];
                            /*
                             * Find the counterIJ for the symmetric binary interaction
                             */
                            n = m_kk*j + k;
                            counterIJ2 = m_CounterIJ[n];
                            sum4 = sum4 +
                                   (fabs(charge(i))*
                                    molality[j]*molality[k]*CMX_L[counterIJ2]);
                        }
                    }
                }

                /*
                 * for Anions, do the neutral species interaction
                 */
                if (charge(j) == 0.0) {
                    sum5 = sum5 + molality[j]*2.0*m_Lambda_nj_L(j,i);
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            size_t izeta = j;
                            size_t jzeta = k;
                            size_t kzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + kzeta;
                            double zeta_L = psi_ijk_L[n];
                            if (zeta_L != 0.0) {
                                sum5 = sum5 + molality[j]*molality[k]*zeta_L;
                            }
                        }
                    }
                }
            }
            m_dlnActCoeffMolaldT_Unscaled[i] =
                zsqdFdT + sum1 + sum2 + sum3 + sum4 + sum5;
            d_gamma_dT_Unscaled[i] = exp(m_dlnActCoeffMolaldT_Unscaled[i]);
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf(" %-16s lngamma[i]=%10.6f gamma[i]=%10.6f\n",
                       sni.c_str(), m_dlnActCoeffMolaldT_Unscaled[i], d_gamma_dT_Unscaled[i]);
                printf("                   %12g %12g %12g %12g %12g %12g\n",
                       zsqdFdT, sum1, sum2, sum3, sum4, sum5);
            }
        }
        /*
         * ------ SUBSECTION FOR CALCULATING NEUTRAL SOLUTE ACT COEFF -------
         * ------ -> equations agree with my notes,
         *        -> Equations agree with Pitzer,
         */
        if (charge(i) == 0.0) {
            sum1 = 0.0;
            sum3 = 0.0;
            for (j = 1; j < m_kk; j++) {
                sum1 = sum1 + molality[j]*2.0*m_Lambda_nj_L(i,j);
                /*
                 * Zeta term -> we piggyback on the psi term
                 */
                if (charge(j) > 0.0) {
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum3 = sum3 + molality[j]*molality[k]*psi_ijk_L[n];
                        }
                    }
                }
            }
            sum2 = 3.0 * molality[i] * molality[i] * m_Mu_nnn_L[i];
            m_dlnActCoeffMolaldT_Unscaled[i] = sum1 + sum2 + sum3;
            d_gamma_dT_Unscaled[i] = exp(m_dlnActCoeffMolaldT_Unscaled[i]);
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf(" %-16s lngamma[i]=%10.6f gamma[i]=%10.6f \n",
                       sni.c_str(), m_dlnActCoeffMolaldT_Unscaled[i], d_gamma_dT_Unscaled[i]);
            }
        }

    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 9: \n");
    }
    /*
     * ------ SUBSECTION FOR CALCULATING THE d OSMOTIC COEFF dT ---------
     *
     */
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    double sum6 = 0.0;
    double sum7 = 0.0;
    /*
     * term1 is the temperature derivative of the
     * DH term in the osmotic coefficient expression
     * b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer
     *                          implementations.
     * Is = Ionic strength on the molality scale (units of (gmol/kg))
     * Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
     */
    term1 = -dAphidT * Is * sqrt(Is) / (1.0 + 1.2 * sqrt(Is));

    for (j = 1; j < m_kk; j++) {
        /*
         * Loop Over Cations
         */
        if (charge(j) > 0.0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];

                    sum1 = sum1 + molality[j]*molality[k]*
                           (BphiMX_L[counterIJ] + molarcharge*CMX_L[counterIJ]);
                }
            }

            for (size_t k = j+1; k < m_kk; k++) {
                if (j == (m_kk-1)) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dT",
                                       "logic error 1 in Step 9 of hmw_act");
                }
                if (charge(k) > 0.0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     * between 2 cations.
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];
                    sum2 = sum2 + molality[j]*molality[k]*Phiphi_L[counterIJ];
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) < 0.0) {
                            // species m is an anion
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum2 = sum2 +
                                   molality[j]*molality[k]*molality[m]*psi_ijk_L[n];
                        }
                    }
                }
            }
        }

        /*
         * Loop Over Anions
         */
        if (charge(j) < 0) {
            for (size_t k = j+1; k < m_kk; k++) {
                if (j == m_kk-1) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dT",
                                       "logic error 2 in Step 9 of hmw_act");
                }
                if (charge(k) < 0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     * between two anions
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];

                    sum3 = sum3 + molality[j]*molality[k]*Phiphi_L[counterIJ];
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum3 = sum3 +
                                   molality[j]*molality[k]*molality[m]*psi_ijk_L[n];
                        }
                    }
                }
            }
        }

        /*
         * Loop Over Neutral Species
         */
        if (charge(j) == 0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    sum4 = sum4 + molality[j]*molality[k]*m_Lambda_nj_L(j,k);
                }
                if (charge(k) > 0.0) {
                    sum5 = sum5 + molality[j]*molality[k]*m_Lambda_nj_L(j,k);
                }
                if (charge(k) == 0.0) {
                    if (k > j) {
                        sum6 = sum6 + molality[j]*molality[k]*m_Lambda_nj_L(j,k);
                    } else if (k == j) {
                        sum6 = sum6 + 0.5 * molality[j]*molality[k]*m_Lambda_nj_L(j,k);
                    }
                }
                if (charge(k) < 0.0) {
                    size_t izeta = j;
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            size_t jzeta = m;
                            n = k + jzeta * m_kk + izeta * m_kk * m_kk;
                            double zeta_L = psi_ijk_L[n];
                            if (zeta_L != 0.0) {
                                sum7 += molality[izeta]*molality[jzeta]*molality[k]*zeta_L;
                            }
                        }
                    }
                }
            }
            sum7 += molality[j]*molality[j]*molality[j]*m_Mu_nnn_L[j];
        }
    }
    sum_m_phi_minus_1 = 2.0 *
                        (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7);
    /*
     * Calculate the osmotic coefficient from
     *       osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
     */
    if (molalitysum > 1.0E-150) {
        d_osmotic_coef_dT = 0.0 + (sum_m_phi_minus_1 / molalitysum);
    } else {
        d_osmotic_coef_dT = 0.0;
    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" term1=%10.6f sum1=%10.6f sum2=%10.6f "
               "sum3=%10.6f sum4=%10.6f sum5=%10.6f\n",
               term1, sum1, sum2, sum3, sum4, sum5);
        printf("     sum_m_phi_minus_1=%10.6f        d_osmotic_coef_dT =%10.6f\n",
               sum_m_phi_minus_1, d_osmotic_coef_dT);
        printf(" Step 10: \n");
    }
    d_lnwateract_dT = -(m_weightSolvent/1000.0) * molalitysum * d_osmotic_coef_dT;

    /*
     * In Cantera, we define the activity coefficient of the solvent as
     *
     *     act_0 = actcoeff_0 * Xmol_0
     *
     * We have just computed act_0. However, this routine returns
     *  ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
     */
    m_dlnActCoeffMolaldT_Unscaled[0] = d_lnwateract_dT;
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        double d_wateract_dT = exp(d_lnwateract_dT);
        printf(" d_ln_a_water_dT = %10.6f d_a_water_dT=%10.6f\n\n",
               d_lnwateract_dT, d_wateract_dT);
    }
}

void HMWSoln::s_update_d2lnMolalityActCoeff_dT2() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if( cached.validate(temperature(), pressure(), stateMFNumber()) ) {
        return;
    }

    /*
     *  Zero the unscaled 2nd derivatives
     */
    m_d2lnActCoeffMolaldT2_Unscaled.assign(m_kk, 0.0);
    /*
     * Calculate the unscaled 2nd derivatives
     */
    s_updatePitzer_d2lnMolalityActCoeff_dT2();

    for (size_t k = 1; k < m_kk; k++) {
        if (CROP_speciesCropped_[k] == 2) {
            m_d2lnActCoeffMolaldT2_Unscaled[k] = 0.0;
        }
    }

    if (CROP_speciesCropped_[0]) {
        m_d2lnActCoeffMolaldT2_Unscaled[0] = 0.0;
    }

    /*
     * Scale the 2nd derivatives
     */
    s_updateScaling_pHScaling_dT2();
}

void HMWSoln::s_updatePitzer_d2lnMolalityActCoeff_dT2() const
{
    /*
     * HKM -> Assumption is made that the solvent is species 0.
     */
#ifdef DEBUG_MODE
    m_debugCalc = 0;
#endif
    if (m_indexSolvent != 0) {
        throw CanteraError("HMWSoln::s_updatePitzer_d2lnMolalityActCoeff_dT2",
                           "Wrong index solvent value!");
    }

    std::string sni, snj, snk;

    const double* molality  =  DATA_PTR(m_molalitiesCropped);
    const double* beta0MX_LL=  DATA_PTR(m_Beta0MX_ij_LL);
    const double* beta1MX_LL=  DATA_PTR(m_Beta1MX_ij_LL);
    const double* beta2MX_LL=  DATA_PTR(m_Beta2MX_ij_LL);
    const double* CphiMX_LL =  DATA_PTR(m_CphiMX_ij_LL);
    const double* thetaij_LL=  DATA_PTR(m_Theta_ij_LL);
    const double* alpha1MX  =  DATA_PTR(m_Alpha1MX_ij);
    const double* alpha2MX  =  DATA_PTR(m_Alpha2MX_ij);
    const double* psi_ijk_LL=  DATA_PTR(m_Psi_ijk_LL);

    /*
     * Local variables defined by Coltrin
     */
    double etheta[5][5], etheta_prime[5][5], sqrtIs;
    /*
     * Molality based ionic strength of the solution
     */
    double Is = 0.0;
    /*
     * Molarcharge of the solution: In Pitzer's notation,
     * this is his variable called "Z".
     */
    double molarcharge = 0.0;
    /*
     * molalitysum is the sum of the molalities over all solutes,
     * even those with zero charge.
     */
    double molalitysum = 0.0;

    double* gfunc    =  DATA_PTR(m_gfunc_IJ);
    double* g2func   =  DATA_PTR(m_g2func_IJ);
    double* hfunc    =  DATA_PTR(m_hfunc_IJ);
    double* h2func   =  DATA_PTR(m_h2func_IJ);
    double* BMX_LL   =  DATA_PTR(m_BMX_IJ_LL);
    double* BprimeMX_LL=DATA_PTR(m_BprimeMX_IJ_LL);
    double* BphiMX_LL=  DATA_PTR(m_BphiMX_IJ_LL);
    double* Phi_LL   =  DATA_PTR(m_Phi_IJ_LL);
    double* Phiprime =  DATA_PTR(m_Phiprime_IJ);
    double* Phiphi_LL=  DATA_PTR(m_PhiPhi_IJ_LL);
    double* CMX_LL   =  DATA_PTR(m_CMX_IJ_LL);


    double x1, x2;
    double d2FdT2, zsqd2FdT2;
    double sum1, sum2, sum3, sum4, sum5, term1;
    double sum_m_phi_minus_1, d2_osmotic_coef_dT2, d2_lnwateract_dT2;

    int z1, z2;
    size_t n, i, j, m, counterIJ,  counterIJ2;

    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf("\n Debugging information from "
               "s_Pitzer_d2lnMolalityActCoeff_dT2()\n");
    }
    /*
     * Make sure the counter variables are setup
     */
    counterIJ_setup();


    /*
     * ---------- Calculate common sums over solutes ---------------------
     */
    for (n = 1; n < m_kk; n++) {
        //      ionic strength
        Is += charge(n) * charge(n) * molality[n];
        //      total molar charge
        molarcharge +=  fabs(charge(n)) * molality[n];
        molalitysum += molality[n];
    }
    Is *= 0.5;

    /*
     * Store the ionic molality in the object for reference.
     */
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 1: \n");
        printf(" ionic strenth      = %14.7le \n total molar "
               "charge = %14.7le \n", Is, molarcharge);
    }

    /*
     * The following call to calc_lambdas() calculates all 16 elements
     * of the elambda and elambda1 arrays, given the value of the
     * ionic strength (Is)
     */
    calc_lambdas(Is);

    /*
     * ----- Step 2:  Find the coefficients E-theta and -------------------
     *                E-thetaprime for all combinations of positive
     *                unlike charges up to 4
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 2: \n");
    }
    for (z1 = 1; z1 <=4; z1++) {
        for (z2 =1; z2 <=4; z2++) {
            calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf(" z1=%3d z2=%3d E-theta(I) = %f, E-thetaprime(I) = %f\n",
                       z1, z2, etheta[z1][z2], etheta_prime[z1][z2]);
            }
        }
    }

    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 3: \n");
        printf(" Species          Species            g(x) "
               " hfunc(x)   \n");
    }

    /*
     *
     *  calculate gfunc(x) and hfunc(x) for each cation-anion pair MX
     *   In the original literature, hfunc, was called gprime. However,
     *   it's not the derivative of gfunc(x), so I renamed it.
     */
    for (i = 1; i < (m_kk - 1); i++) {
        for (j = (i+1); j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * Only loop over oppositely charge species
             */
            if (charge(i)*charge(j) < 0) {
                /*
                 * x is a reduced function variable
                 */
                x1 = sqrtIs * alpha1MX[counterIJ];
                if (x1 > 1.0E-100) {
                    gfunc[counterIJ] =  2.0*(1.0-(1.0 + x1) * exp(-x1)) / (x1 *x1);
                    hfunc[counterIJ] = -2.0*
                                       (1.0-(1.0 + x1 + 0.5*x1 * x1) * exp(-x1)) / (x1 * x1);
                } else {
                    gfunc[counterIJ] = 0.0;
                    hfunc[counterIJ] = 0.0;
                }

                if (beta2MX_LL[counterIJ] != 0.0) {
                    x2 = sqrtIs * alpha2MX[counterIJ];
                    if (x2 > 1.0E-100) {
                        g2func[counterIJ] =  2.0*(1.0-(1.0 + x2) * exp(-x2)) / (x2 * x2);
                        h2func[counterIJ] = -2.0 *
                                            (1.0-(1.0 + x2 + 0.5 * x2 * x2) * exp(-x2)) / (x2 * x2);
                    } else {
                        g2func[counterIJ] = 0.0;
                        h2func[counterIJ] = 0.0;
                    }
                }
            } else {
                gfunc[counterIJ] = 0.0;
                hfunc[counterIJ] = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %9.5f %9.5f \n", sni.c_str(), snj.c_str(),
                       gfunc[counterIJ], hfunc[counterIJ]);
            }
        }
    }
    /*
     * ------- SUBSECTION TO CALCULATE BMX_L, BprimeMX_LL, BphiMX_L ----------
     * ------- These are now temperature derivatives of the
     *         previously calculated quantities.
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 4: \n");
        printf(" Species          Species            BMX    "
               "BprimeMX    BphiMX   \n");
    }

    for (i = 1; i < m_kk - 1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0.0) {
                BMX_LL[counterIJ]  = beta0MX_LL[counterIJ]
                                     + beta1MX_LL[counterIJ] * gfunc[counterIJ]
                                     + beta2MX_LL[counterIJ] * g2func[counterIJ];
                if (DEBUG_MODE_ENABLED && m_debugCalc) {
                    printf("%d %g: %g %g %g %g\n",
                           (int) counterIJ,  BMX_LL[counterIJ], beta0MX_LL[counterIJ],
                           beta1MX_LL[counterIJ], beta2MX_LL[counterIJ], gfunc[counterIJ]);
                }
                if (Is > 1.0E-150) {
                    BprimeMX_LL[counterIJ] = (beta1MX_LL[counterIJ] * hfunc[counterIJ]/Is +
                                              beta2MX_LL[counterIJ] * h2func[counterIJ]/Is);
                } else {
                    BprimeMX_LL[counterIJ] = 0.0;
                }
                BphiMX_LL[counterIJ] = BMX_LL[counterIJ] + Is*BprimeMX_LL[counterIJ];
            } else {
                BMX_LL[counterIJ]      = 0.0;
                BprimeMX_LL[counterIJ] = 0.0;
                BphiMX_LL[counterIJ]     = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %11.7f %11.7f %11.7f \n",
                       sni.c_str(), snj.c_str(),
                       BMX_LL[counterIJ], BprimeMX_LL[counterIJ], BphiMX_LL[counterIJ]);
            }
        }
    }

    /*
     * --------- SUBSECTION TO CALCULATE CMX_LL ----------
     * ---------
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 5: \n");
        printf(" Species          Species            CMX \n");
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0.0) {
                CMX_LL[counterIJ] = CphiMX_LL[counterIJ]/
                                    (2.0* sqrt(fabs(charge(i)*charge(j))));
            } else {
                CMX_LL[counterIJ] = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %11.7f \n", sni.c_str(), snj.c_str(),
                       CMX_LL[counterIJ]);
            }
        }
    }

    /*
     * ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
     * --------
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 6: \n");
        printf(" Species          Species            Phi_ij "
               " Phiprime_ij  Phi^phi_ij \n");
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) > 0) {
                z1 = (int) fabs(charge(i));
                z2 = (int) fabs(charge(j));
                Phi_LL[counterIJ] = thetaij_LL[counterIJ];
                Phiprime[counterIJ] = 0.0;
                Phiphi_LL[counterIJ] = Phi_LL[counterIJ];
            } else {
                Phi_LL[counterIJ]      = 0.0;
                Phiprime[counterIJ] = 0.0;
                Phiphi_LL[counterIJ]   = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %10.6f %10.6f %10.6f \n",
                       sni.c_str(), snj.c_str(),
                       Phi_LL[counterIJ], Phiprime[counterIJ], Phiphi_LL[counterIJ]);
            }
        }
    }

    /*
     * ----------- SUBSECTION FOR CALCULATION OF d2FdT2 ---------------------
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 7: \n");
    }
    double d2AphidT2 = d2A_DebyedT2_TP() / 3.0;
    d2FdT2 = -d2AphidT2 * (sqrt(Is) / (1.0 + 1.2*sqrt(Is))
                           + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" initial value of d2FdT2 = %10.6f \n", d2FdT2);
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0) {
                d2FdT2 = d2FdT2 + molality[i]*molality[j] * BprimeMX_LL[counterIJ];
            }
            /*
             * Both species have a non-zero charge, and they
             * have the same sign, e.g., both positive or both negative.
             */
            if (charge(i)*charge(j) > 0) {
                d2FdT2 = d2FdT2 + molality[i]*molality[j] * Phiprime[counterIJ];
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf(" d2FdT2 = %10.6f \n", d2FdT2);
            }
        }
    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 8: \n");
    }

    for (i = 1; i < m_kk; i++) {

        /*
         * -------- SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR CATIONS -----
         * --
         */
        if (charge(i) > 0) {
            // species i is the cation (positive) to calc the actcoeff
            zsqd2FdT2 = charge(i)*charge(i)*d2FdT2;
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            sum4 = 0.0;
            sum5 = 0.0;
            for (j = 1; j < m_kk; j++) {
                /*
                 * Find the counterIJ for the symmetric binary interaction
                 */
                n = m_kk*i + j;
                counterIJ = m_CounterIJ[n];

                if (charge(j) < 0.0) {
                    // sum over all anions
                    sum1 = sum1 + molality[j]*
                           (2.0*BMX_LL[counterIJ] + molarcharge*CMX_LL[counterIJ]);
                    if (j < m_kk-1) {
                        /*
                         * This term is the ternary interaction involving the
                         * non-duplicate sum over double anions, j, k, with
                         * respect to the cation, i.
                         */
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all anions
                            if (charge(k) < 0.0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 = sum3 + molality[j]*molality[k]*psi_ijk_LL[n];
                            }
                        }
                    }
                }


                if (charge(j) > 0.0) {
                    // sum over all cations
                    if (j != i) {
                        sum2 = sum2 + molality[j]*(2.0*Phi_LL[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            // two inner sums over anions

                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 = sum2 + molality[j]*molality[k]*psi_ijk_LL[n];
                            /*
                             * Find the counterIJ for the j,k interaction
                             */
                            n = m_kk*j + k;
                            counterIJ2 = m_CounterIJ[n];
                            sum4 = sum4 + (fabs(charge(i))*
                                           molality[j]*molality[k]*CMX_LL[counterIJ2]);
                        }
                    }
                }

                /*
                 * Handle neutral j species
                 */
                if (charge(j) == 0) {
                    sum5 = sum5 + molality[j]*2.0*m_Lambda_nj_LL(j,i);
                    /*
                     * Zeta interaction term
                     */
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            size_t izeta = j;
                            size_t jzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + k;
                            double zeta_LL = psi_ijk_LL[n];
                            if (zeta_LL != 0.0) {
                                sum5 = sum5 + molality[j]*molality[k]*zeta_LL;
                            }
                        }
                    }
                }
            }
            /*
             * Add all of the contributions up to yield the log of the
             * solute activity coefficients (molality scale)
             */
            m_d2lnActCoeffMolaldT2_Unscaled[i] =
                zsqd2FdT2 + sum1 + sum2 + sum3 + sum4 + sum5;
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf(" %-16s d2lngammadT2[i]=%10.6f \n",
                       sni.c_str(), m_d2lnActCoeffMolaldT2_Unscaled[i]);
                printf("                   %12g %12g %12g %12g %12g %12g\n",
                       zsqd2FdT2, sum1, sum2, sum3, sum4, sum5);
            }
        }


        /*
         * ------ SUBSECTION FOR CALCULATING THE d2ACTCOEFFdT2 FOR ANIONS ------
         *
         */
        if (charge(i) < 0) {
            //          species i is an anion (negative)
            zsqd2FdT2 = charge(i)*charge(i)*d2FdT2;
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            sum4 = 0.0;
            sum5 = 0.0;
            for (j = 1; j < m_kk; j++) {
                /*
                 * Find the counterIJ for the symmetric binary interaction
                 */
                n = m_kk*i + j;
                counterIJ = m_CounterIJ[n];

                /*
                 * For Anions, do the cation interactions.
                 */
                if (charge(j) > 0) {
                    sum1 = sum1 + molality[j]*
                           (2.0*BMX_LL[counterIJ] + molarcharge*CMX_LL[counterIJ]);
                    if (j < m_kk-1) {
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all cations
                            if (charge(k) > 0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 = sum3 + molality[j]*molality[k]*psi_ijk_LL[n];
                            }
                        }
                    }
                }

                /*
                 * For Anions, do the other anion interactions.
                 */
                if (charge(j) < 0.0) {
                    //  sum over all anions
                    if (j != i) {
                        sum2 = sum2 + molality[j]*(2.0*Phi_LL[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            // two inner sums over cations
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 = sum2 + molality[j]*molality[k]*psi_ijk_LL[n];
                            /*
                             * Find the counterIJ for the symmetric binary interaction
                             */
                            n = m_kk*j + k;
                            counterIJ2 = m_CounterIJ[n];
                            sum4 = sum4 +
                                   (fabs(charge(i))*
                                    molality[j]*molality[k]*CMX_LL[counterIJ2]);
                        }
                    }
                }

                /*
                 * for Anions, do the neutral species interaction
                 */
                if (charge(j) == 0.0) {
                    sum5 = sum5 + molality[j]*2.0*m_Lambda_nj_LL(j,i);
                    /*
                     * Zeta interaction term
                     */
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            size_t izeta = j;
                            size_t jzeta = k;
                            size_t kzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + kzeta;
                            double zeta_LL = psi_ijk_LL[n];
                            if (zeta_LL != 0.0) {
                                sum5 = sum5 + molality[j]*molality[k]*zeta_LL;
                            }
                        }
                    }
                }
            }
            m_d2lnActCoeffMolaldT2_Unscaled[i] =
                zsqd2FdT2 + sum1 + sum2 + sum3 + sum4 + sum5;
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf(" %-16s d2lngammadT2[i]=%10.6f\n",
                       sni.c_str(), m_d2lnActCoeffMolaldT2_Unscaled[i]);
                printf("                   %12g %12g %12g %12g %12g %12g\n",
                       zsqd2FdT2, sum1, sum2, sum3, sum4, sum5);
            }
        }
        /*
         * ------ SUBSECTION FOR CALCULATING NEUTRAL SOLUTE ACT COEFF -------
         * ------ -> equations agree with my notes,
         *        -> Equations agree with Pitzer,
         */
        if (charge(i) == 0.0) {
            sum1 = 0.0;
            sum3 = 0.0;
            for (j = 1; j < m_kk; j++) {
                sum1 = sum1 + molality[j]*2.0*m_Lambda_nj_LL(i,j);
                /*
                 * Zeta term -> we piggyback on the psi term
                 */
                if (charge(j) > 0.0) {
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum3 = sum3 + molality[j]*molality[k]*psi_ijk_LL[n];
                        }
                    }
                }
            }
            sum2 = 3.0 * molality[i] * molality[i] * m_Mu_nnn_LL[i];
            m_d2lnActCoeffMolaldT2_Unscaled[i] = sum1 + sum2 + sum3;
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf(" %-16s d2lngammadT2[i]=%10.6f \n",
                       sni.c_str(), m_d2lnActCoeffMolaldT2_Unscaled[i]);
            }
        }

    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 9: \n");
    }

    /*
     * ------ SUBSECTION FOR CALCULATING THE d2 OSMOTIC COEFF dT2 ---------
     */
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    double sum6 = 0.0;
    double sum7 = 0.0;
    /*
     * term1 is the temperature derivative of the
     * DH term in the osmotic coefficient expression
     * b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer
     *                          implementations.
     * Is = Ionic strength on the molality scale (units of (gmol/kg))
     * Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
     */
    term1 = -d2AphidT2 * Is * sqrt(Is) / (1.0 + 1.2 * sqrt(Is));

    for (j = 1; j < m_kk; j++) {
        /*
         * Loop Over Cations
         */
        if (charge(j) > 0.0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];

                    sum1 = sum1 + molality[j]*molality[k]*
                           (BphiMX_LL[counterIJ] + molarcharge*CMX_LL[counterIJ]);
                }
            }

            for (size_t k = j+1; k < m_kk; k++) {
                if (j == (m_kk-1)) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_d2lnMolalityActCoeff_dT2",
                                       "logic error 1 in Step 9 of hmw_act");
                }
                if (charge(k) > 0.0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     * between 2 cations.
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];
                    sum2 = sum2 + molality[j]*molality[k]*Phiphi_LL[counterIJ];
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) < 0.0) {
                            // species m is an anion
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum2 = sum2 +
                                   molality[j]*molality[k]*molality[m]*psi_ijk_LL[n];
                        }
                    }
                }
            }
        }

        /*
         * Loop Over Anions
         */
        if (charge(j) < 0) {
            for (size_t k = j+1; k < m_kk; k++) {
                if (j == m_kk-1) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_d2lnMolalityActCoeff_dT2",
                                       "logic error 2 in Step 9 of hmw_act");
                }
                if (charge(k) < 0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     * between two anions
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];

                    sum3 = sum3 + molality[j]*molality[k]*Phiphi_LL[counterIJ];
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum3 = sum3 +
                                   molality[j]*molality[k]*molality[m]*psi_ijk_LL[n];
                        }
                    }
                }
            }
        }

        /*
         * Loop Over Neutral Species
         */
        if (charge(j) == 0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    sum4 = sum4 + molality[j]*molality[k]*m_Lambda_nj_LL(j,k);
                }
                if (charge(k) > 0.0) {
                    sum5 = sum5 + molality[j]*molality[k]*m_Lambda_nj_LL(j,k);
                }
                if (charge(k) == 0.0) {
                    if (k > j) {
                        sum6 = sum6 + molality[j]*molality[k]*m_Lambda_nj_LL(j,k);
                    } else if (k == j) {
                        sum6 = sum6 + 0.5 * molality[j]*molality[k]*m_Lambda_nj_LL(j,k);
                    }
                }
                if (charge(k) < 0.0) {
                    size_t izeta = j;
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            size_t jzeta = m;
                            n = k + jzeta * m_kk + izeta * m_kk * m_kk;
                            double zeta_LL = psi_ijk_LL[n];
                            if (zeta_LL != 0.0) {
                                sum7 += molality[izeta]*molality[jzeta]*molality[k]*zeta_LL;
                            }
                        }
                    }
                }
            }

            sum7 += molality[j] * molality[j] * molality[j] * m_Mu_nnn_LL[j];
        }
    }
    sum_m_phi_minus_1 = 2.0 *
                        (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7);
    /*
     * Calculate the osmotic coefficient from
     *       osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
     */
    if (molalitysum > 1.0E-150) {
        d2_osmotic_coef_dT2 = 0.0 + (sum_m_phi_minus_1 / molalitysum);
    } else {
        d2_osmotic_coef_dT2 = 0.0;
    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" term1=%10.6f sum1=%10.6f sum2=%10.6f "
               "sum3=%10.6f sum4=%10.6f sum5=%10.6f\n",
               term1, sum1, sum2, sum3, sum4, sum5);
        printf("     sum_m_phi_minus_1=%10.6f        d2_osmotic_coef_dT2=%10.6f\n",
               sum_m_phi_minus_1, d2_osmotic_coef_dT2);
        printf(" Step 10: \n");
    }
    d2_lnwateract_dT2 = -(m_weightSolvent/1000.0) * molalitysum * d2_osmotic_coef_dT2;

    /*
     * In Cantera, we define the activity coefficient of the solvent as
     *
     *     act_0 = actcoeff_0 * Xmol_0
     *
     * We have just computed act_0. However, this routine returns
     *  ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
     */
    m_d2lnActCoeffMolaldT2_Unscaled[0] = d2_lnwateract_dT2;

    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        double d2_wateract_dT2 = exp(d2_lnwateract_dT2);
        printf(" d2_ln_a_water_dT2 = %10.6f d2_a_water_dT2=%10.6f\n\n",
               d2_lnwateract_dT2, d2_wateract_dT2);
    }
}

void HMWSoln::s_update_dlnMolalityActCoeff_dP() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    if( cached.validate(temperature(), pressure(), stateMFNumber()) ) {
        return;
    }

    m_dlnActCoeffMolaldP_Unscaled.assign(m_kk, 0.0);
    s_updatePitzer_dlnMolalityActCoeff_dP();

    for (size_t k = 1; k < m_kk; k++) {
        if (CROP_speciesCropped_[k] == 2) {
            m_dlnActCoeffMolaldP_Unscaled[k] = 0.0;
        }
    }

    if (CROP_speciesCropped_[0]) {
        m_dlnActCoeffMolaldP_Unscaled[0] = 0.0;
    }

    s_updateScaling_pHScaling_dP();
}

void HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dP() const
{
    /*
     * HKM -> Assumption is made that the solvent is species 0.
     */
#ifdef DEBUG_MODE
    m_debugCalc = 0;
#endif
    if (m_indexSolvent != 0) {
        throw CanteraError("HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dP",
                           "Wrong index solvent value!");
    }

    std::string sni, snj, snk;

    const double* molality  =  DATA_PTR(m_molalitiesCropped);
    const double* beta0MX_P =  DATA_PTR(m_Beta0MX_ij_P);
    const double* beta1MX_P =  DATA_PTR(m_Beta1MX_ij_P);
    const double* beta2MX_P =  DATA_PTR(m_Beta2MX_ij_P);
    const double* CphiMX_P  =  DATA_PTR(m_CphiMX_ij_P);
    const double* thetaij_P =  DATA_PTR(m_Theta_ij_P);
    const double* alpha1MX  =  DATA_PTR(m_Alpha1MX_ij);
    const double* alpha2MX  =  DATA_PTR(m_Alpha2MX_ij);
    const double* psi_ijk_P =  DATA_PTR(m_Psi_ijk_P);

    /*
     * Local variables defined by Coltrin
     */
    double etheta[5][5], etheta_prime[5][5], sqrtIs;
    /*
     * Molality based ionic strength of the solution
     */
    double Is = 0.0;
    /*
     * Molarcharge of the solution: In Pitzer's notation,
     * this is his variable called "Z".
     */
    double molarcharge = 0.0;
    /*
     * molalitysum is the sum of the molalities over all solutes,
     * even those with zero charge.
     */
    double molalitysum = 0.0;

    double* gfunc    =  DATA_PTR(m_gfunc_IJ);
    double* g2func   =  DATA_PTR(m_g2func_IJ);
    double* hfunc    =  DATA_PTR(m_hfunc_IJ);
    double* h2func   =  DATA_PTR(m_h2func_IJ);
    double* BMX_P    =  DATA_PTR(m_BMX_IJ_P);
    double* BprimeMX_P= DATA_PTR(m_BprimeMX_IJ_P);
    double* BphiMX_P =  DATA_PTR(m_BphiMX_IJ_P);
    double* Phi_P    =  DATA_PTR(m_Phi_IJ_P);
    double* Phiprime =  DATA_PTR(m_Phiprime_IJ);
    double* Phiphi_P =  DATA_PTR(m_PhiPhi_IJ_P);
    double* CMX_P    =  DATA_PTR(m_CMX_IJ_P);

    double x1, x2;
    double dFdP, zsqdFdP;
    double sum1, sum2, sum3, sum4, sum5, term1;
    double sum_m_phi_minus_1, d_osmotic_coef_dP, d_lnwateract_dP;

    int z1, z2;
    size_t n, i, j, m, counterIJ,  counterIJ2;

    double currTemp = temperature();
    double currPres = pressure();

    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf("\n Debugging information from "
               "s_Pitzer_dlnMolalityActCoeff_dP()\n");
    }
    /*
     * Make sure the counter variables are setup
     */
    counterIJ_setup();

    /*
     * ---------- Calculate common sums over solutes ---------------------
     */
    for (n = 1; n < m_kk; n++) {
        //      ionic strength
        Is += charge(n) * charge(n) * molality[n];
        //      total molar charge
        molarcharge +=  fabs(charge(n)) * molality[n];
        molalitysum += molality[n];
    }
    Is *= 0.5;

    /*
     * Store the ionic molality in the object for reference.
     */
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 1: \n");
        printf(" ionic strenth      = %14.7le \n total molar "
               "charge = %14.7le \n", Is, molarcharge);
    }

    /*
     * The following call to calc_lambdas() calculates all 16 elements
     * of the elambda and elambda1 arrays, given the value of the
     * ionic strength (Is)
     */
    calc_lambdas(Is);


    /*
     * ----- Step 2:  Find the coefficients E-theta and -------------------
     *                E-thetaprime for all combinations of positive
     *                unlike charges up to 4
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 2: \n");
    }
    for (z1 = 1; z1 <=4; z1++) {
        for (z2 =1; z2 <=4; z2++) {
            calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf(" z1=%3d z2=%3d E-theta(I) = %f, E-thetaprime(I) = %f\n",
                       z1, z2, etheta[z1][z2], etheta_prime[z1][z2]);
            }
        }
    }

    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 3: \n");
        printf(" Species          Species            g(x) "
               " hfunc(x)   \n");
    }

    /*
     *
     *  calculate g(x) and hfunc(x) for each cation-anion pair MX
     *   In the original literature, hfunc, was called gprime. However,
     *   it's not the derivative of g(x), so I renamed it.
     */
    for (i = 1; i < (m_kk - 1); i++) {
        for (j = (i+1); j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * Only loop over oppositely charge species
             */
            if (charge(i)*charge(j) < 0) {
                /*
                 * x is a reduced function variable
                 */
                x1 = sqrtIs * alpha1MX[counterIJ];
                if (x1 > 1.0E-100) {
                    gfunc[counterIJ] =  2.0*(1.0-(1.0 + x1) * exp(-x1)) / (x1 * x1);
                    hfunc[counterIJ] = -2.0*
                                       (1.0-(1.0 + x1 + 0.5 * x1 * x1) * exp(-x1)) / (x1 * x1);
                } else {
                    gfunc[counterIJ] = 0.0;
                    hfunc[counterIJ] = 0.0;
                }

                if (beta2MX_P[counterIJ] != 0.0) {
                    x2 = sqrtIs * alpha2MX[counterIJ];
                    if (x2 > 1.0E-100) {
                        g2func[counterIJ] =  2.0*(1.0-(1.0 + x2) * exp(-x2)) / (x2 * x2);
                        h2func[counterIJ] = -2.0 *
                                            (1.0-(1.0 + x2 + 0.5 * x2 * x2) * exp(-x2)) / (x2 * x2);
                    } else {
                        g2func[counterIJ] = 0.0;
                        h2func[counterIJ] = 0.0;
                    }
                }
            } else {
                gfunc[counterIJ] = 0.0;
                hfunc[counterIJ] = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %9.5f %9.5f \n", sni.c_str(), snj.c_str(),
                       gfunc[counterIJ], hfunc[counterIJ]);
            }
        }
    }

    /*
     * ------- SUBSECTION TO CALCULATE BMX_P, BprimeMX_P, BphiMX_P ----------
     * ------- These are now temperature derivatives of the
     *         previously calculated quantities.
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 4: \n");
        printf(" Species          Species            BMX    "
               "BprimeMX    BphiMX   \n");
    }

    for (i = 1; i < m_kk - 1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0.0) {
                BMX_P[counterIJ]  = beta0MX_P[counterIJ]
                                    + beta1MX_P[counterIJ] * gfunc[counterIJ]
                                    + beta2MX_P[counterIJ] * g2func[counterIJ];
                if (DEBUG_MODE_ENABLED && m_debugCalc) {
                    printf("%d %g: %g %g %g %g\n",
                           (int) counterIJ,  BMX_P[counterIJ], beta0MX_P[counterIJ],
                           beta1MX_P[counterIJ], beta2MX_P[counterIJ], gfunc[counterIJ]);
                }
                if (Is > 1.0E-150) {
                    BprimeMX_P[counterIJ] = (beta1MX_P[counterIJ] * hfunc[counterIJ]/Is +
                                             beta2MX_P[counterIJ] * h2func[counterIJ]/Is);
                } else {
                    BprimeMX_P[counterIJ] = 0.0;
                }
                BphiMX_P[counterIJ] = BMX_P[counterIJ] + Is*BprimeMX_P[counterIJ];
            } else {
                BMX_P[counterIJ]      = 0.0;
                BprimeMX_P[counterIJ] = 0.0;
                BphiMX_P[counterIJ]     = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %11.7f %11.7f %11.7f \n",
                       sni.c_str(), snj.c_str(),
                       BMX_P[counterIJ], BprimeMX_P[counterIJ], BphiMX_P[counterIJ]);
            }
        }
    }

    /*
     * --------- SUBSECTION TO CALCULATE CMX_P ----------
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 5: \n");
        printf(" Species          Species            CMX \n");
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0.0) {
                CMX_P[counterIJ] = CphiMX_P[counterIJ]/
                                   (2.0* sqrt(fabs(charge(i)*charge(j))));
            } else {
                CMX_P[counterIJ] = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %11.7f \n", sni.c_str(), snj.c_str(),
                       CMX_P[counterIJ]);
            }
        }
    }

    /*
     * ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
     * --------
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 6: \n");
        printf(" Species          Species            Phi_ij "
               " Phiprime_ij  Phi^phi_ij \n");
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) > 0) {
                z1 = (int) fabs(charge(i));
                z2 = (int) fabs(charge(j));
                Phi_P[counterIJ] = thetaij_P[counterIJ];
                Phiprime[counterIJ] = 0.0;
                Phiphi_P[counterIJ] = Phi_P[counterIJ] + Is * Phiprime[counterIJ];
            } else {
                Phi_P[counterIJ]      = 0.0;
                Phiprime[counterIJ] = 0.0;
                Phiphi_P[counterIJ]   = 0.0;
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                snj = speciesName(j);
                printf(" %-16s %-16s %10.6f %10.6f %10.6f \n",
                       sni.c_str(), snj.c_str(),
                       Phi_P[counterIJ], Phiprime[counterIJ], Phiphi_P[counterIJ]);
            }
        }
    }

    /*
     * ----------- SUBSECTION FOR CALCULATION OF dFdT ---------------------
     */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 7: \n");
    }
    double dA_DebyedP = dA_DebyedP_TP(currTemp, currPres);
    double dAphidP = dA_DebyedP /3.0;
    dFdP = -dAphidP * (sqrt(Is) / (1.0 + 1.2*sqrt(Is))
                       + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" initial value of dFdP = %10.6f \n", dFdP);
    }
    for (i = 1; i < m_kk-1; i++) {
        for (j = i+1; j < m_kk; j++) {
            /*
             * Find the counterIJ for the symmetric binary interaction
             */
            n = m_kk*i + j;
            counterIJ = m_CounterIJ[n];
            /*
             * both species have a non-zero charge, and one is positive
             * and the other is negative
             */
            if (charge(i)*charge(j) < 0) {
                dFdP = dFdP + molality[i]*molality[j] * BprimeMX_P[counterIJ];
            }
            /*
             * Both species have a non-zero charge, and they
             * have the same sign, e.g., both positive or both negative.
             */
            if (charge(i)*charge(j) > 0) {
                dFdP = dFdP + molality[i]*molality[j] * Phiprime[counterIJ];
            }
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf(" dFdP = %10.6f \n", dFdP);
            }
        }
    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 8: \n");
    }

    for (i = 1; i < m_kk; i++) {

        /*
         * -------- SUBSECTION FOR CALCULATING THE dACTCOEFFdP FOR CATIONS -----
         */
        if (charge(i) > 0) {
            // species i is the cation (positive) to calc the actcoeff
            zsqdFdP = charge(i)*charge(i)*dFdP;
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            sum4 = 0.0;
            sum5 = 0.0;
            for (j = 1; j < m_kk; j++) {
                /*
                 * Find the counterIJ for the symmetric binary interaction
                 */
                n = m_kk*i + j;
                counterIJ = m_CounterIJ[n];

                if (charge(j) < 0.0) {
                    // sum over all anions
                    sum1 = sum1 + molality[j]*
                           (2.0*BMX_P[counterIJ] + molarcharge*CMX_P[counterIJ]);
                    if (j < m_kk-1) {
                        /*
                         * This term is the ternary interaction involving the
                         * non-duplicate sum over double anions, j, k, with
                         * respect to the cation, i.
                         */
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all anions
                            if (charge(k) < 0.0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 = sum3 + molality[j]*molality[k]*psi_ijk_P[n];
                            }
                        }
                    }
                }


                if (charge(j) > 0.0) {
                    // sum over all cations
                    if (j != i) {
                        sum2 = sum2 + molality[j]*(2.0*Phi_P[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            // two inner sums over anions

                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 = sum2 + molality[j]*molality[k]*psi_ijk_P[n];
                            /*
                             * Find the counterIJ for the j,k interaction
                             */
                            n = m_kk*j + k;
                            counterIJ2 = m_CounterIJ[n];
                            sum4 = sum4 + (fabs(charge(i))*
                                           molality[j]*molality[k]*CMX_P[counterIJ2]);
                        }
                    }
                }

                /*
                 * for Anions, do the neutral species interaction
                 */
                if (charge(j) == 0) {
                    sum5 = sum5 + molality[j]*2.0*m_Lambda_nj_P(j,i);
                    /*
                     * Zeta interaction term
                     */
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            size_t izeta = j;
                            size_t jzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + k;
                            double zeta_P = psi_ijk_P[n];
                            if (zeta_P != 0.0) {
                                sum5 = sum5 + molality[j]*molality[k]*zeta_P;
                            }
                        }
                    }
                }
            }

            /*
             * Add all of the contributions up to yield the log of the
             * solute activity coefficients (molality scale)
             */
            m_dlnActCoeffMolaldP_Unscaled[i] =
                zsqdFdP + sum1 + sum2 + sum3 + sum4 + sum5;

            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf(" %-16s lngamma[i]=%10.6f \n",
                       sni.c_str(), m_dlnActCoeffMolaldP_Unscaled[i]);
                printf("                   %12g %12g %12g %12g %12g %12g\n",
                       zsqdFdP, sum1, sum2, sum3, sum4, sum5);
            }
        }


        /*
         * ------ SUBSECTION FOR CALCULATING THE dACTCOEFFdP FOR ANIONS ------
         */
        if (charge(i) < 0) {
            //          species i is an anion (negative)
            zsqdFdP = charge(i)*charge(i)*dFdP;
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            sum4 = 0.0;
            sum5 = 0.0;
            for (j = 1; j < m_kk; j++) {
                /*
                 * Find the counterIJ for the symmetric binary interaction
                 */
                n = m_kk*i + j;
                counterIJ = m_CounterIJ[n];

                /*
                 * For Anions, do the cation interactions.
                 */
                if (charge(j) > 0) {
                    sum1 = sum1 + molality[j]*
                           (2.0*BMX_P[counterIJ] + molarcharge*CMX_P[counterIJ]);
                    if (j < m_kk-1) {
                        for (size_t k = j+1; k < m_kk; k++) {
                            // an inner sum over all cations
                            if (charge(k) > 0) {
                                n = k + j * m_kk + i * m_kk * m_kk;
                                sum3 = sum3 + molality[j]*molality[k]*psi_ijk_P[n];
                            }
                        }
                    }
                }

                /*
                 * For Anions, do the other anion interactions.
                 */
                if (charge(j) < 0.0) {
                    //  sum over all anions
                    if (j != i) {
                        sum2 = sum2 + molality[j]*(2.0*Phi_P[counterIJ]);
                    }
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            // two inner sums over cations
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum2 = sum2 + molality[j]*molality[k]*psi_ijk_P[n];
                            /*
                             * Find the counterIJ for the symmetric binary interaction
                             */
                            n = m_kk*j + k;
                            counterIJ2 = m_CounterIJ[n];
                            sum4 = sum4 +
                                   (fabs(charge(i))*
                                    molality[j]*molality[k]*CMX_P[counterIJ2]);
                        }
                    }
                }

                /*
                 * for Anions, do the neutral species interaction
                 */
                if (charge(j) == 0.0) {
                    sum5 = sum5 + molality[j]*2.0*m_Lambda_nj_P(j,i);
                    /*
                     * Zeta interaction term
                     */
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) > 0.0) {
                            size_t izeta = j;
                            size_t jzeta = k;
                            size_t kzeta = i;
                            n = izeta * m_kk * m_kk + jzeta * m_kk + kzeta;
                            double zeta_P = psi_ijk_P[n];
                            if (zeta_P != 0.0) {
                                sum5 = sum5 + molality[j]*molality[k]*zeta_P;
                            }
                        }
                    }
                }
            }
            m_dlnActCoeffMolaldP_Unscaled[i] =
                zsqdFdP + sum1 + sum2 + sum3 + sum4 + sum5;
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf(" %-16s lndactcoeffmolaldP[i]=%10.6f \n",
                       sni.c_str(), m_dlnActCoeffMolaldP_Unscaled[i]);
                printf("                   %12g %12g %12g %12g %12g %12g\n",
                       zsqdFdP, sum1, sum2, sum3, sum4, sum5);
            }
        }

        /*
         * ------ SUBSECTION FOR CALCULATING d NEUTRAL SOLUTE ACT COEFF dP -------
         */
        if (charge(i) == 0.0) {
            sum1 = 0.0;
            sum3 = 0.0;
            for (j = 1; j < m_kk; j++) {
                sum1 +=  molality[j]*2.0*m_Lambda_nj_P(i,j);
                /*
                 * Zeta term -> we piggyback on the psi term
                 */
                if (charge(j) > 0.0) {
                    for (size_t k = 1; k < m_kk; k++) {
                        if (charge(k) < 0.0) {
                            n = k + j * m_kk + i * m_kk * m_kk;
                            sum3 = sum3 + molality[j]*molality[k]*psi_ijk_P[n];
                        }
                    }
                }
            }
            sum2 = 3.0 * molality[i] * molality[i] * m_Mu_nnn_P[i];
            m_dlnActCoeffMolaldP_Unscaled[i] = sum1 + sum2 + sum3;
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                sni = speciesName(i);
                printf(" %-16s dlnActCoeffMolaldP[i]=%10.6f \n",
                       sni.c_str(), m_dlnActCoeffMolaldP_Unscaled[i]);
            }
        }

    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Step 9: \n");
    }

    /*
     * ------ SUBSECTION FOR CALCULATING THE d OSMOTIC COEFF dP ---------
     */
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    double sum6 = 0.0;
    double sum7 = 0.0;
    /*
     * term1 is the temperature derivative of the
     * DH term in the osmotic coefficient expression
     * b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer
     *                          implementations.
     * Is = Ionic strength on the molality scale (units of (gmol/kg))
     * Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
     */
    term1 = -dAphidP * Is * sqrt(Is) / (1.0 + 1.2 * sqrt(Is));

    for (j = 1; j < m_kk; j++) {
        /*
         * Loop Over Cations
         */
        if (charge(j) > 0.0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];

                    sum1 = sum1 + molality[j]*molality[k]*
                           (BphiMX_P[counterIJ] + molarcharge*CMX_P[counterIJ]);
                }
            }

            for (size_t k = j+1; k < m_kk; k++) {
                if (j == (m_kk-1)) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dP",
                                       "logic error 1 in Step 9 of hmw_act");
                }
                if (charge(k) > 0.0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     * between 2 cations.
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];
                    sum2 = sum2 + molality[j]*molality[k]*Phiphi_P[counterIJ];
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) < 0.0) {
                            // species m is an anion
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum2 = sum2 +
                                   molality[j]*molality[k]*molality[m]*psi_ijk_P[n];
                        }
                    }
                }
            }
        }


        /*
         * Loop Over Anions
         */
        if (charge(j) < 0) {
            for (size_t k = j+1; k < m_kk; k++) {
                if (j == m_kk-1) {
                    // we should never reach this step
                    throw CanteraError("HMWSoln::s_updatePitzer_dlnMolalityActCoeff_dP",
                                       "logic error 2 in Step 9 of hmw_act");
                }
                if (charge(k) < 0) {
                    /*
                     * Find the counterIJ for the symmetric j,k binary interaction
                     * between two anions
                     */
                    n = m_kk*j + k;
                    counterIJ = m_CounterIJ[n];

                    sum3 = sum3 + molality[j]*molality[k]*Phiphi_P[counterIJ];
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            n = m + k * m_kk + j * m_kk * m_kk;
                            sum3 = sum3 +
                                   molality[j]*molality[k]*molality[m]*psi_ijk_P[n];
                        }
                    }
                }
            }
        }

        /*
         * Loop Over Neutral Species
         */
        if (charge(j) == 0) {
            for (size_t k = 1; k < m_kk; k++) {
                if (charge(k) < 0.0) {
                    sum4 = sum4 + molality[j]*molality[k]*m_Lambda_nj_P(j,k);
                }
                if (charge(k) > 0.0) {
                    sum5 = sum5 + molality[j]*molality[k]*m_Lambda_nj_P(j,k);
                }
                if (charge(k) == 0.0) {
                    if (k > j) {
                        sum6 = sum6 + molality[j]*molality[k]*m_Lambda_nj_P(j,k);
                    } else if (k == j) {
                        sum6 = sum6 + 0.5 * molality[j]*molality[k]*m_Lambda_nj_P(j,k);
                    }
                }
                if (charge(k) < 0.0) {
                    size_t izeta = j;
                    for (m = 1; m < m_kk; m++) {
                        if (charge(m) > 0.0) {
                            size_t jzeta = m;
                            n = k + jzeta * m_kk + izeta * m_kk * m_kk;
                            double zeta_P = psi_ijk_P[n];
                            if (zeta_P != 0.0) {
                                sum7 += molality[izeta]*molality[jzeta]*molality[k]*zeta_P;
                            }
                        }
                    }
                }
            }

            sum7 += molality[j] * molality[j] * molality[j] * m_Mu_nnn_P[j];
        }
    }
    sum_m_phi_minus_1 = 2.0 *
                        (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7);

    /*
     * Calculate the osmotic coefficient from
     *       osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
     */
    if (molalitysum > 1.0E-150) {
        d_osmotic_coef_dP = 0.0 + (sum_m_phi_minus_1 / molalitysum);
    } else {
        d_osmotic_coef_dP = 0.0;
    }
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" term1=%10.6f sum1=%10.6f sum2=%10.6f "
               "sum3=%10.6f sum4=%10.6f sum5=%10.6f\n",
               term1, sum1, sum2, sum3, sum4, sum5);
        printf("     sum_m_phi_minus_1=%10.6f        d_osmotic_coef_dP =%10.6f\n",
               sum_m_phi_minus_1, d_osmotic_coef_dP);
        printf(" Step 10: \n");
    }
    d_lnwateract_dP = -(m_weightSolvent/1000.0) * molalitysum * d_osmotic_coef_dP;


    /*
     * In Cantera, we define the activity coefficient of the solvent as
     *
     *     act_0 = actcoeff_0 * Xmol_0
     *
     * We have just computed act_0. However, this routine returns
     *  ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
     */
    m_dlnActCoeffMolaldP_Unscaled[0] = d_lnwateract_dP;
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        double d_wateract_dP = exp(d_lnwateract_dP);
        printf(" d_ln_a_water_dP = %10.6f d_a_water_dP=%10.6f\n\n",
               d_lnwateract_dP, d_wateract_dP);
    }
}

void HMWSoln::calc_lambdas(double is) const
{
    if( m_last_is == is ) {
      return;
    }
    m_last_is = is;

    double aphi, dj, jfunc, jprime, t, x, zprod;
    int i, ij, j;
    /*
     * Coefficients c1-c4 are used to approximate
     * the integral function "J";
     * aphi is the Debye-Huckel constant at 25 C
     */

    double c1 = 4.581, c2 = 0.7237, c3 = 0.0120, c4 = 0.528;

    aphi = 0.392;   /* Value at 25 C */
    if (DEBUG_MODE_ENABLED && m_debugCalc) {
        printf(" Is = %g\n", is);
    }
    if (is < 1.0E-150) {
        for (i = 0; i < 17; i++) {
            elambda[i] = 0.0;
            elambda1[i] = 0.0;
        }
        return;
    }
    /*
     * Calculate E-lambda terms for charge combinations of like sign,
     * using method of Pitzer (1975). Charges up to 4 are calculated.
     */

    for (i=1; i<=4; i++) {
        for (j=i; j<=4; j++) {
            ij = i*j;
            /*
             * calculate the product of the charges
             */
            zprod = (double)ij;
            /*
             * calculate Xmn (A1) from Harvie, Weare (1980).
             */
            x = 6.0* zprod * aphi * sqrt(is);                      /* eqn 23 */

            jfunc = x / (4.0 + c1*pow(x,-c2)*exp(-c3*pow(x,c4)));  /* eqn 47 */

            t = c3 * c4 * pow(x,c4);
            dj = c1* pow(x,(-c2-1.0)) * (c2+t) * exp(-c3*pow(x,c4));
            jprime = (jfunc/x)*(1.0 + jfunc*dj);

            elambda[ij] = zprod*jfunc / (4.0*is);                  /* eqn 14 */
            elambda1[ij] = (3.0*zprod*zprod*aphi*jprime/(4.0*sqrt(is))
                            - elambda[ij])/is;
            if (DEBUG_MODE_ENABLED && m_debugCalc) {
                printf(" ij = %d, elambda = %g, elambda1 = %g\n",
                       ij, elambda[ij], elambda1[ij]);
            }
        }
    }
}

void HMWSoln::calc_thetas(int z1, int z2,
                          double* etheta, double* etheta_prime) const
{
    int i, j;
    double f1, f2;

    /*
     * Calculate E-theta(i) and E-theta'(I) using method of
     * Pitzer (1987)
     */
    i = abs(z1);
    j = abs(z2);

    AssertThrowMsg(i <= 4 && j <= 4, "HMWSoln::calc_thetas",
                   "we shouldn't be here");
    AssertThrowMsg(i != 0 && j != 0, "HMWSoln::calc_thetas",
                   "called with one species being neutral");

    /*
     *  Check to see if the charges are of opposite sign. If they are of
     *  opposite sign then their etheta interaction is zero.
     */
    if (z1*z2 < 0) {
        *etheta = 0.0;
        *etheta_prime = 0.0;
    }
    /*
     * Actually calculate the interaction.
     */
    else {
        f1 = (double)i / (2.0  * j);
        f2 = (double)j / (2.0  * i);
        *etheta = elambda[i*j] - f1*elambda[j*j] - f2*elambda[i*i];
        *etheta_prime = elambda1[i*j] - f1*elambda1[j*j] - f2*elambda1[i*i];
    }
}

void  HMWSoln::s_updateIMS_lnMolalityActCoeff() const
{
    double tmp;
    /*
     * Calculate the molalities. Currently, the molalities
     * may not be current with respect to the contents of the
     * State objects' data.
     */
    calcMolalities();

    double xmolSolvent = moleFraction(m_indexSolvent);
    double xx = std::max(m_xmolSolventMIN, xmolSolvent);
    if (IMS_typeCutoff_ == 0) {
        for (size_t k = 1; k < m_kk; k++) {
            IMS_lnActCoeffMolal_[k]= 0.0;
        }
        IMS_lnActCoeffMolal_[m_indexSolvent] = - log(xx) + (xx - 1.0)/xx;
        return;
    } else if (IMS_typeCutoff_ == 1) {
        if (xmolSolvent > 3.0 * IMS_X_o_cutoff_/2.0) {
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= 0.0;
            }
            IMS_lnActCoeffMolal_[m_indexSolvent] = - log(xx) + (xx - 1.0)/xx;
            return;
        } else if (xmolSolvent < IMS_X_o_cutoff_/2.0) {
            tmp = log(xx * IMS_gamma_k_min_);
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= tmp;
            }
            IMS_lnActCoeffMolal_[m_indexSolvent] = log(IMS_gamma_o_min_);
            return;
        } else {
            /*
                 * If we are in the middle region, calculate the connecting polynomials
                 */
            double xminus  = xmolSolvent - IMS_X_o_cutoff_/2.0;
            double xminus2 = xminus * xminus;
            double xminus3 = xminus2 * xminus;
            double x_o_cut2 = IMS_X_o_cutoff_ * IMS_X_o_cutoff_;
            double x_o_cut3 =  x_o_cut2 * IMS_X_o_cutoff_;

            double h2 = 3.5 * xminus2 /  IMS_X_o_cutoff_ - 2.0 * xminus3 / x_o_cut2;
            double h2_prime = 7.0 * xminus /  IMS_X_o_cutoff_ - 6.0 * xminus2 /  x_o_cut2;

            double h1 = (1.0 - 3.0 * xminus2 /  x_o_cut2 + 2.0 *  xminus3/ x_o_cut3);
            double h1_prime = (- 6.0 * xminus /  x_o_cut2 + 6.0 *  xminus2/ x_o_cut3);

            double h1_g = h1 / IMS_gamma_o_min_;
            double h1_g_prime  = h1_prime / IMS_gamma_o_min_;

            double alpha = 1.0 / (exp(1.0) * IMS_gamma_k_min_);
            double h1_f = h1 * alpha;
            double h1_f_prime  = h1_prime * alpha;

            double f = h2 + h1_f;
            double f_prime = h2_prime + h1_f_prime;

            double g = h2 + h1_g;
            double g_prime = h2_prime + h1_g_prime;

            tmp = (xmolSolvent/ g * g_prime + (1.0-xmolSolvent) / f * f_prime);
            double lngammak = -1.0 - log(f) + tmp * xmolSolvent;
            double lngammao =-log(g) - tmp * (1.0-xmolSolvent);

            tmp = log(xmolSolvent) + lngammak;
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= tmp;
            }
            IMS_lnActCoeffMolal_[m_indexSolvent] = lngammao;
        }
    }
    // Exponentials - trial 2
    else if (IMS_typeCutoff_ == 2) {
        if (xmolSolvent > IMS_X_o_cutoff_) {
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= 0.0;
            }
            IMS_lnActCoeffMolal_[m_indexSolvent] = - log(xx) + (xx - 1.0)/xx;
            return;
        } else {

            double xoverc = xmolSolvent/IMS_cCut_;
            double eterm = std::exp(-xoverc);

            double fptmp = IMS_bfCut_  - IMS_afCut_ / IMS_cCut_ - IMS_bfCut_*xoverc
                           + 2.0*IMS_dfCut_*xmolSolvent - IMS_dfCut_*xmolSolvent*xoverc;
            double f_prime = 1.0 + eterm*fptmp;
            double f = xmolSolvent + IMS_efCut_
                       + eterm * (IMS_afCut_ + xmolSolvent * (IMS_bfCut_ + IMS_dfCut_*xmolSolvent));

            double gptmp = IMS_bgCut_  - IMS_agCut_ / IMS_cCut_ - IMS_bgCut_*xoverc
                           + 2.0*IMS_dgCut_*xmolSolvent - IMS_dgCut_*xmolSolvent*xoverc;
            double g_prime = 1.0 + eterm*gptmp;
            double g = xmolSolvent + IMS_egCut_
                       + eterm * (IMS_agCut_ + xmolSolvent * (IMS_bgCut_ + IMS_dgCut_*xmolSolvent));

            tmp = (xmolSolvent / g * g_prime + (1.0 - xmolSolvent) / f * f_prime);
            double lngammak = -1.0 - log(f) + tmp * xmolSolvent;
            double lngammao =-log(g) - tmp * (1.0-xmolSolvent);

            tmp = log(xx) + lngammak;
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= tmp;
            }
            IMS_lnActCoeffMolal_[m_indexSolvent] = lngammao;
        }
    }
    return;
}

void HMWSoln::printCoeffs() const
{
    size_t i, j, k;
    std::string sni, snj;
    calcMolalities();
    double* molality = DATA_PTR(m_molalitiesCropped);
    double* moleF = DATA_PTR(m_tmpV);
    /*
     * Update the coefficients wrt Temperature
     * Calculate the derivatives as well
     */
    s_updatePitzer_CoeffWRTemp(2);
    getMoleFractions(moleF);

    printf("Index  Name                  MoleF   MolalityCropped  Charge\n");
    for (k = 0; k < m_kk; k++) {
        sni = speciesName(k);
        printf("%2s     %-16s %14.7le %14.7le %5.1f \n",
               int2str(k).c_str(), sni.c_str(), moleF[k], molality[k], charge(k));
    }

    printf("\n Species          Species            beta0MX  "
           "beta1MX   beta2MX   CphiMX    alphaMX thetaij    \n");
    for (i = 1; i < m_kk - 1; i++) {
        sni = speciesName(i);
        for (j = i+1; j < m_kk; j++) {
            snj = speciesName(j);
            size_t n  = i * m_kk + j;
            size_t ct = m_CounterIJ[n];
            printf(" %-16s %-16s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f \n",
                   sni.c_str(), snj.c_str(),
                   m_Beta0MX_ij[ct], m_Beta1MX_ij[ct],
                   m_Beta2MX_ij[ct], m_CphiMX_ij[ct],
                   m_Alpha1MX_ij[ct], m_Theta_ij[ct]);


        }
    }

    printf("\n Species          Species          Species       "
           "psi   \n");
    for (i = 1; i < m_kk; i++) {
        sni = speciesName(i);
        for (j = 1; j < m_kk; j++) {
            snj = speciesName(j);
            for (k = 1; k < m_kk; k++) {
                std::string snk = speciesName(k);
                size_t n = k + j * m_kk + i * m_kk * m_kk;
                if (m_Psi_ijk[n] != 0.0) {
                    printf(" %-16s %-16s %-16s %9.5f \n",
                           sni.c_str(), snj.c_str(),
                           snk.c_str(), m_Psi_ijk[n]);
                }
            }
        }
    }
}

void HMWSoln::applyphScale(doublereal* acMolality) const
{
    if (m_pHScalingType == PHSCALE_PITZER) {
        return;
    }
    AssertTrace(m_pHScalingType == PHSCALE_NBS);
    doublereal lnGammaClMs2 = s_NBS_CLM_lnMolalityActCoeff();
    doublereal lnGammaCLMs1 = m_lnActCoeffMolal_Unscaled[m_indexCLM];
    doublereal afac = -1.0 *(lnGammaClMs2 - lnGammaCLMs1);
    for (size_t k = 0; k < m_kk; k++) {
        acMolality[k] *= exp(charge(k) * afac);
    }
}

void HMWSoln::s_updateScaling_pHScaling() const
{
    if (m_pHScalingType == PHSCALE_PITZER) {
        m_lnActCoeffMolal_Scaled = m_lnActCoeffMolal_Unscaled;
        return;
    }
    AssertTrace(m_pHScalingType == PHSCALE_NBS);
    doublereal lnGammaClMs2 = s_NBS_CLM_lnMolalityActCoeff();
    doublereal lnGammaCLMs1 = m_lnActCoeffMolal_Unscaled[m_indexCLM];
    doublereal afac = -1.0 *(lnGammaClMs2 - lnGammaCLMs1);
    for (size_t k = 0; k < m_kk; k++) {
        m_lnActCoeffMolal_Scaled[k] = m_lnActCoeffMolal_Unscaled[k] + charge(k) * afac;
    }
}

void HMWSoln::s_updateScaling_pHScaling_dT() const
{
    if (m_pHScalingType == PHSCALE_PITZER) {
        m_dlnActCoeffMolaldT_Scaled = m_dlnActCoeffMolaldT_Unscaled;
        return;
    }
    AssertTrace(m_pHScalingType == PHSCALE_NBS);
    doublereal dlnGammaClM_dT_s2 = s_NBS_CLM_dlnMolalityActCoeff_dT();
    doublereal dlnGammaCLM_dT_s1 = m_dlnActCoeffMolaldT_Unscaled[m_indexCLM];
    doublereal afac = -1.0 *(dlnGammaClM_dT_s2 - dlnGammaCLM_dT_s1);
    for (size_t k = 0; k < m_kk; k++) {
        m_dlnActCoeffMolaldT_Scaled[k] = m_dlnActCoeffMolaldT_Unscaled[k] + charge(k) * afac;
    }
}

void HMWSoln::s_updateScaling_pHScaling_dT2() const
{
    if (m_pHScalingType == PHSCALE_PITZER) {
        m_d2lnActCoeffMolaldT2_Scaled = m_d2lnActCoeffMolaldT2_Unscaled;
        return;
    }
    AssertTrace(m_pHScalingType == PHSCALE_NBS);
    doublereal d2lnGammaClM_dT2_s2 = s_NBS_CLM_d2lnMolalityActCoeff_dT2();
    doublereal d2lnGammaCLM_dT2_s1 = m_d2lnActCoeffMolaldT2_Unscaled[m_indexCLM];
    doublereal afac = -1.0 *(d2lnGammaClM_dT2_s2 - d2lnGammaCLM_dT2_s1);
    for (size_t k = 0; k < m_kk; k++) {
        m_d2lnActCoeffMolaldT2_Scaled[k] = m_d2lnActCoeffMolaldT2_Unscaled[k] + charge(k) * afac;
    }
}

void HMWSoln::s_updateScaling_pHScaling_dP() const
{
    if (m_pHScalingType == PHSCALE_PITZER) {
        m_dlnActCoeffMolaldP_Scaled = m_dlnActCoeffMolaldP_Unscaled;
        return;
    }
    AssertTrace(m_pHScalingType == PHSCALE_NBS);
    doublereal dlnGammaClM_dP_s2 = s_NBS_CLM_dlnMolalityActCoeff_dP();
    doublereal dlnGammaCLM_dP_s1 = m_dlnActCoeffMolaldP_Unscaled[m_indexCLM];
    doublereal afac = -1.0 *(dlnGammaClM_dP_s2 - dlnGammaCLM_dP_s1);
    for (size_t k = 0; k < m_kk; k++) {
        m_dlnActCoeffMolaldP_Scaled[k] = m_dlnActCoeffMolaldP_Unscaled[k] + charge(k) * afac;
    }
}

doublereal HMWSoln::s_NBS_CLM_lnMolalityActCoeff() const
{
    doublereal sqrtIs = sqrt(m_IionicMolality);
    doublereal A = A_Debye_TP();
    doublereal lnGammaClMs2 = - A * sqrtIs /(1.0 + 1.5 * sqrtIs);
    return  lnGammaClMs2;
}

doublereal HMWSoln::s_NBS_CLM_dlnMolalityActCoeff_dT() const
{
    doublereal sqrtIs = sqrt(m_IionicMolality);
    doublereal dAdT = dA_DebyedT_TP();
    return - dAdT * sqrtIs /(1.0 + 1.5 * sqrtIs);
}

doublereal HMWSoln::s_NBS_CLM_d2lnMolalityActCoeff_dT2() const
{
    doublereal sqrtIs = sqrt(m_IionicMolality);
    doublereal d2AdT2 = d2A_DebyedT2_TP();
    return - d2AdT2 * sqrtIs /(1.0 + 1.5 * sqrtIs);
}

doublereal HMWSoln::s_NBS_CLM_dlnMolalityActCoeff_dP() const
{
    doublereal sqrtIs = sqrt(m_IionicMolality);
    doublereal dAdP = dA_DebyedP_TP();
    return - dAdP * sqrtIs /(1.0 + 1.5 * sqrtIs);
}

int HMWSoln::debugPrinting()
{
#ifdef DEBUG_MODE
    return m_debugCalc;
#else
    return 0;
#endif
}

}
