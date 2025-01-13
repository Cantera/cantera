/**
 * @file CustomKinetics.h
 *
 * @ingroup chemkinetics
 */

// Author: Q. Cazeres, A. Felden, P. Pepiot

#include "GasKinetics.h"

namespace Cantera
{
  /**
   *  Kinetics manager implementing reaction mechanism GRI-Mech 3.0
   *  @deprecated
   */
  class CustomKinetics : public GasKinetics
  {
  public:
    /// Default constructor.
    CustomKinetics();

    // virtual int type() const {
    //     return cCustomKinetics;
    // }

    void getNetProductionRates(double* net) {
          get_wdot_reduced(net);
    }

    string kineticsType() const override {
       return "custom";
    }

    void closeDynamicLib() {
      close_dl();
    }

    void* handle;

private:
    typedef void (*ck_t)(double*,double*,const double*, double*);
    // New yarc2 format
    //typedef void (*ck_t)(const double*,const double*,const double*,const double*,double*);
    ck_t ck;
    void get_wdot_reduced(double* wdot);
    void close_dl();
  };
}
