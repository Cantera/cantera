/**
 * @file electrolyteThermo.h
 *
 * Support for thermo property calculation from C++ application programs.
 * This header file includes several headers from the Cantera kernel needed
 * to evaluate thermo properties.
 * @deprecated To be removed after Cantera 2.3. Include relevant headers directly.
 */

#ifndef CT_ELECTROLYTETHERMO_INCL
#define CT_ELECTROLYTETHERMO_INCL

#pragma message "Deprecated. electrolyteThermo.h will be removed after Cantera 2.3. Include relevant headers directly."

#include "thermo/electrolytes.h"
#include "thermo/MolalityVPSSTP.h"
#include "thermo/VPStandardStateTP.h"
#include "thermo/IdealMolalSoln.h"
#include "thermo/WaterPropsIAPWS.h"
#include "thermo/WaterProps.h"
#include "thermo/PDSS.h"
#include "thermo/PDSS_Water.h"
#include "thermo/PDSS_HKFT.h"
#include "thermo/HMWSoln.h"
#include "thermo/DebyeHuckel.h"
#include "thermo/WaterSSTP.h"
#include "thermo/VPSSMgr_Water_HKFT.h"
#include "thermo/VPSSMgr_Water_ConstVol.h"

#endif
