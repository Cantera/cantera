#ifndef CT_MIX_DEFS_H
#define CT_MIX_DEFS_H

namespace Cantera {

  /**
   *  This generic id is used as the default in virtual base
   *  classes that employ id's. It is used to indicate the lack
   *  of an inherited class that would define the id.
   */
  const int cNone = 0;

  // species thermo types
  const int cNASA = 1;
  const int cShomate = 2;
  const int cNASA96 = 3;
  const int cHarmonicOsc = 4;

  /**
   * Equation of state types:
   *
   *  These types are used in the member function eosType() of
   *  the virtual base class ThermoPhase. They are used to
   *  distinguish different types of equation of states. Also, they
   *  may be used for upcasting from the ThermoPhase class.  Their
   *  id's should be distinct.
   *
   *  Users who wish to define their own equation of states which
   *  derive from ThermoPhase should define a unique id which
   *  doesn't conflict with those listed below. The Cantera Kernel
   *  however, will not be know about the class and will therefore
   *  not be able to initialize the class within its "factory"
   *  routines. 
   */
  const int cIdealGas = 1;       //  IdealGasPhase in IdealGasPhase.h
  const int cIncompressible = 2; //  ConstDensityThermo in ConstDensityThermo.h
  /// A surface phase. Used by class SurfPhase.
  const int cSurf = 3;           

  /// A metal phase. 
  const int cMetal = 4;          //  MetalPhase in MetalPhase.h
  //    const int cSolidCompound = 5;  //  SolidCompound in SolidCompound.h
  const int cStoichSubstance = 5; // StoichSubstance.h
  const int cSemiconductor = 7;

  const int cMineralEQ3 = 8; // MineralEQ3 in MineralEQ3.h
  const int cMetalSHEelectrons = 9; // SHE electrode electrons

  const int cLatticeSolid = 20; // LatticeSolidPhase.h
  const int cLattice = 21; 

  // pure fluids with liquid/vapor eqs of state
  const int cPureFluid = 10;

  /// An edge between two 2D surfaces    
  const int cEdge = 6;

  /// Constant partial molar volume solution IdealSolidSolnPhase.h
  const int cIdealSolidSolnPhase = 5009;

  //! HMW - Strong electrolyte using the Pitzer formulation
  const int cHMW = 40;

  //! DebyeHuckel - Weak electrolyte using various Debye-Huckel formulations
  const int cDebyeHuckel = 50;

  //! IdealMolalSoln - molality based solution with molality-based act coeffs of 1
  const int cIdealMolalSoln = 60;

  const int cIdealSolnGasVPSS = 500;
  const int cIdealSolnGasVPSS_iscv = 501;

  const int cMargulesVPSSTP = 301;

  // perfect gas -- companion of PecosTransport for high speed flows
  const int cPerfectGas = 666;

  const int cIonsFromNeutral = 2000;

  //! Variable Pressure Standard State ThermoPhase objects
  const int cVPSS_IdealGas     = 1001;
  const int cVPSS_ConstVol     = 1002;
  const int cVPSS_PureFluid    = 1010;
  const int cVPSS_HMW          = 1040;
  const int cVPSS_DebyeHuckel = 1050;
  const int cVPSS_MolalSoln   = 1060;

  //! Types of general formulations for the specification of the standard state volume
  enum SSVolume_Model_enumType {
    //! This approximation is for a constant volume
    cSSVOLUME_CONSTANT = 0,
    //! This approximation is for a species with a quadratic polynomial in temperature
    /*!
     *       V^ss_i = ai + bi T + ci T2
     */
    cSSVOLUME_TPOLY,
    //! This approximation is for a species where the density is expressed as a
    //! quadratic polynomial in temperature
    /*!
     *       V^ss_i = M_i / (ai + bi T + ci T2)
     */
    cSSVOLUME_DENSITY_TPOLY
  };

  //! Types of PDSS's
  enum PDSS_enumType {
    cPDSS_UNDEF = 100,
    cPDSS_IDEALGAS,
    cPDSS_CONSTVOL,
    cPDSS_SSVOL,
    cPDSS_MOLAL_CONSTVOL,
    cPDSS_WATER,
    cPDSS_MOLAL_HKFT,
    cPDSS_IONSFROMNEUTRAL
  };


  //! enum for VPSSMgr types
  enum VPSSMgr_enumType {
    cVPSSMGR_UNDEF = 1000,
    cVPSSMGR_IDEALGAS,
    cVPSSMGR_CONSTVOL , 
    cVPSSMGR_PUREFLUID,
    cVPSSMGR_WATER_CONSTVOL, 
    cVPSSMGR_WATER_HKFT,
    cVPSSMGR_GENERAL
  };
     

  // kinetic manager types
  const int cGasKinetics = 2;
  const int cGRI30 = 3;
  const int cInterfaceKinetics = 4;
  const int cLineKinetics = 5;
  const int cEdgeKinetics = 6;
  const int cSolidKinetics = 7;
  const int cAqueousKinetics = 8;
}

#endif
 
