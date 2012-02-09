/**
 * @file electrolyteThermo.h
 *
 * Support for thermo property calculation from C++ application programs.
 * This header file includes several headers from the Cantera kernel needed
 * to evaluate thermo properties.
 */

#ifndef CT_ELECTROLYTETHERMO_INCL
#define CT_ELECTROLYTETHERMO_INCL

#include "thermo.h"

#include "kernel/electrolytes.h"
#include "kernel/MolalityVPSSTP.h"
#include "kernel/VPStandardStateTP.h"
#include "kernel/IdealMolalSoln.h"
#include "kernel/WaterPropsIAPWS.h"
#include "kernel/WaterProps.h"
#include "kernel/PDSS.h"
#include "kernel/PDSS_Water.h"
#include "kernel/PDSS_HKFT.h"
#include "kernel/HMWSoln.h"
#include "kernel/DebyeHuckel.h"
#include "kernel/WaterSSTP.h"
#include "kernel/VPSSMgr_Water_HKFT.h"
#include "kernel/VPSSMgr_Water_ConstVol.h"
#endif
