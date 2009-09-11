
#include "WaterPropsIAPWS.h"
#include "TransportBase.h"
#include "DenseMatrix.h"
#include "LiquidTransportParams.h"


#include "WaterTransport.h"

#include <iostream>
using namespace std;


namespace Cantera {

  //! default constructor
  WaterTransport::WaterTransport(thermo_t* thermo, int ndim) :
    Transport(thermo, ndim)
  {

  }

  //  Copy Constructor for the %WaterThermo object.
  /* 
   *    @param right  ThermoPhase to be copied
   */
  WaterTransport::WaterTransport(const WaterTransport &right) :
    Transport(right.m_thermo, right.m_nDim)
  {
    *this = right;

  }

  // Assignment operator
  /*
   *
   * @param right    Reference to %WaterTransport object to be copied into the
   *                 current one.
   */
  WaterTransport&  WaterTransport::operator=(const  WaterTransport& right)
  {
    if (&right != this) {
      return *this;
    }
    Transport::operator=(right);
 
    return *this;
  }

  // Duplication routine for objects which inherit from %Transport
  /*
   *  This virtual routine can be used to duplicate %Transport objects
   *  inherited from %Transport even if the application only has
   *  a pointer to %Transport to work with.
   *
   *  These routines are basically wrappers around the derived copy
   *  constructor.
   */
  Transport *  WaterTransport::duplMyselfAsTransport() const {
    WaterTransport* tr = new WaterTransport(*this);
    return dynamic_cast<Transport *> (tr);
  }


  // virtual destructor
  WaterTransport::~WaterTransport() {
  }





  const double H[4] = {1.,
		       0.978197,
		       0.579829,
		       -0.202354};



  const double Hij[6][7] =
    {
      { 0.5132047, 0.2151778, -0.2818107,  0.1778064, -0.04176610,          0.,           0.},
      { 0.3205656, 0.7317883, -1.070786 ,  0.4605040,          0., -0.01578386,           0.},
      { 0.,        1.241044 , -1.263184 ,  0.2340379,          0.,          0.,           0.},
      { 0.,        1.476783 ,         0., -0.4924179,   0.1600435,          0., -0.003629481},
      {-0.7782567,      0.0 ,         0.,  0.       ,          0.,          0.,           0.},
      { 0.1885447,      0.0 ,         0.,  0.       ,          0.,          0.,           0.},
    };
  const double TStar = 647.27;       // Kelvin
  const double rhoStar = 317.763;    // kg / m3
  const double presStar = 22.115E6;  // Pa
  const double muStar = 55.071E-6;   //Pa s




  double WaterTransport::viscosity() {

    double temp = m_thermo->temperature();
    double dens = m_thermo->density();

    WaterPropsIAPWS *waterP = new WaterPropsIAPWS();
    waterP->setState_TR(temp, dens);

 
    double pressure = waterP->pressure();
    printf("pressure = %g\n", pressure);

    //dens = 18.02 * pressure / (GasConstant * temp);
    //printf ("mod dens = %g\n", dens);

    double rhobar = dens/rhoStar;

    double tbar = temp / TStar;
    // double pbar = pressure / presStar;

    

    double tbar2 = tbar * tbar;
    double tbar3 = tbar2 * tbar;

    double mu0bar = sqrt(tbar) / (H[0] + H[1]/tbar + H[2]/tbar2 + H[3]/tbar3);

    printf("mu0bar = %g\n", mu0bar);
    printf("mu0 = %g\n", mu0bar * muStar);
    //double tfac0 = 1.0;
    double tfac1 = 1.0 / tbar - 1.0;
    double tfac2 = tfac1 * tfac1;
    double tfac3 = tfac2 * tfac1;
    double tfac4 = tfac3 * tfac1;
    double tfac5 = tfac4 * tfac1;

    //double rfac0 = 1.0;
    double rfac1 = rhobar - 1.0;
    double rfac2 = rfac1 * rfac1;
    double rfac3 = rfac2 * rfac1;
    double rfac4 = rfac3 * rfac1;
    double rfac5 = rfac4 * rfac1;
    double rfac6 = rfac5 * rfac1;

    double sum = (Hij[0][0]       + Hij[1][0]*tfac1       + Hij[4][0]*tfac4       + Hij[5][0]*tfac5 +
		  Hij[0][1]*rfac1 + Hij[1][1]*tfac1*rfac1 + Hij[2][1]*tfac2*rfac1 + Hij[3][1]*tfac3*rfac1 +
		  Hij[0][2]*rfac2 + Hij[1][2]*tfac1*rfac2 + Hij[2][2]*tfac2*rfac2 +
		  Hij[0][3]*rfac3 + Hij[1][3]*tfac1*rfac3 + Hij[2][3]*tfac2*rfac3 + Hij[3][3]*tfac3*rfac3 +
		  Hij[0][4]*rfac4 + Hij[3][4]*tfac3*rfac4 + 
		  Hij[1][5]*tfac1*rfac5 + Hij[3][6]*tfac3*rfac6 
		  );
    double mu1bar = exp(rhobar * sum);
    double mu2bar = 1.0;


    double mubar = mu0bar * mu1bar * mu2bar;




    return mubar * muStar;

  }

}
