#ifndef TPX_HFC134_H
#define TPX_HFC134_H

#include "Sub.h"

namespace tpx {
    class HFC134a : public Substance {
    public:
	HFC134a(){}
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

