/**
 * @file MMCollisionInt.h
 *  Monk and Monchick collision integrals
 */
/*
 * $Author: hkmoffa $
 * $Revision: 1.5 $
 * $Date: 2007/07/26 21:55:32 $
 */

// Copyright 2001  California Institute of Technology 


#ifndef CT_MMCOLLISIONINT_H
#define CT_MMCOLLISIONINT_H

#include <vector>


#include "ct_defs.h"

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

namespace Cantera {

  class XML_Writer;

  class MMCollisionIntError {
  public:
    MMCollisionIntError(std::ostream& logfile, std::string msg) {
      logfile << "#### ERROR ####" << std::endl;
      logfile << "MMCollisionInt: " << msg << std::endl;
      std::cerr << "Error in fitting collision integrals. "
	   << "Execution terminated." << std::endl 
	   << "See transport log file for more information." << std::endl; 
    }
  };


  /**
   * Collision integrals. This class provides functions that
   * interpolate the tabulated collision integrals in Monchick and
   * Mason, "Transport Properties of Polar Gases," J. Chem. Phys. (1961)
   *
   * @ingroup transportgroup
   */
  class MMCollisionInt {

  public:

    MMCollisionInt(){}
    virtual ~MMCollisionInt();
    void init(XML_Writer* xml, doublereal tsmin, 
	      doublereal tsmax, int loglevel = 0);

    doublereal omega22(double ts, double deltastar);
    doublereal astar(double ts, double deltastar);
    doublereal bstar(double ts, double deltastar);
    doublereal cstar(double ts, double deltastar);

    void fit(std::ostream& logfile, int degree, doublereal deltastar,
	     doublereal* astar, doublereal* bstar, doublereal* cstar);

    void fit_omega22(std::ostream& logfile, int degree, doublereal deltastar, doublereal* om22);
    doublereal omega11(double ts, double deltastar) {
      return omega22(ts, deltastar)/astar(ts, deltastar); 
    }    

  private:

    doublereal fitDelta(int table, int ntstar, 
			int degree, doublereal* c);

    std::vector<vector_fp>  m_o22poly;
    std::vector<vector_fp>  m_apoly;
    std::vector<vector_fp>  m_bpoly;
    std::vector<vector_fp>  m_cpoly;

    static doublereal delta[8];
    static doublereal tstar22[37];
    static doublereal omega22_table[37*8];
    static doublereal tstar[39];
    static doublereal astar_table[39*8];
    static doublereal bstar_table[39*8];
    static doublereal cstar_table[39*8];

    vector_fp  m_logTemp;
    int m_nmin, m_nmax;
    XML_Writer* m_xml;
    int m_loglevel;
  };
}
#endif


