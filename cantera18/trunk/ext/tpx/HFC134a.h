#ifndef TPX_HFC134_H
#define TPX_HFC134_H

#include "Sub.h"

namespace tpx {
    class HFC134a : public Substance {
    public:
	HFC134a(){
          m_name = "HFC-134a";
          m_formula = "C2F4H2";
        }
        ~HFC134a(){}

	double MolWt();
	double Tcrit();
        double Pcrit();
	double Vcrit();
	double Tmin();
	double Tmax();
	char * name();
	char * formula();

	double Pp();
	double fp();
	double up();
	double sp() 
            { return (up() - fp())/T + m_entropy_offset; }
	double Psat();
        //    double dPsatdT();
    private:
	double ldens();
    };
}
#endif // ! HFC134_H

