#ifndef CT_FLAMELET_H
#define CT_FLAMELET_H

// #include "Domain1D.h"
#include "cantera/base/Array.h"
#include "cantera/oneD/Flow1D.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/ThermoPhase.h"



#include<iostream>
/**
 * A class for counter-flow diffusion flames in mixture fraction space
 * @ingroup onedim
 * L. Esclapez - 01/2018
 */

namespace Cantera
{

class Flamelet : public Flow1D
{
public:
    Flamelet(IdealGasPhase* ph = 0, size_t nsp = 1, size_t nsoot = 0, size_t neq = c_offset_Yflamelet, size_t points = 1);

    //! Delegating constructor
    Flamelet(shared_ptr<ThermoPhase> th, size_t nsp = 1, size_t nsoot = 0, size_t neq = c_offset_Yflamelet, size_t points = 1);

    //! Create a new flow domain.
    //! @param sol  Solution object used to evaluate all thermodynamic, kinetic, and
    //!     transport properties
    //! @param id  name of flow domain
    //! @param points  initial number of grid points
    Flamelet(shared_ptr<Solution> sol, const string& id="", size_t nsoot = 0, size_t neq = c_offset_Yflamelet, size_t points=1);

    virtual ~Flamelet();

    virtual void eval(size_t j, double* x, double* r, integer* mask, double rdt) override;

    string domainType() const override{
        return "flamelet-flow";
    }

    virtual void resize(size_t components, size_t points) override;

    virtual string componentName(size_t n) const override;

    virtual size_t componentIndex(const string& name) const override;

    vector<double> dCpdT(const Array2D sol);

    vector<double> d2CpdT2(const Array2D sol);

    virtual void setGas(const double* x, size_t j) override;

    virtual void setGasAtMidpoint(const double* x, size_t j) override;

    virtual void resetBadValues(double* xg) override;

    virtual void setChiSt(double chiSt) {
	    m_chiSt = chiSt;
        m_updateChi = true;
    }

    virtual void setzSt(double zSt) {
	    m_zSt = zSt;
        m_updateChi = true;
    }

    virtual double chiSt() const { 
        return m_chiSt; 
    }

    virtual double zSt() const { 
        return m_zSt;
    }

    void setChi();

    void setUnityLewisNumber(bool doLewisNumber) {
	    m_do_unityLewisNumber = doLewisNumber;
    }

    virtual void updateTransport(double* x, size_t j0, size_t j1) override;

    virtual void _getInitialSoln(double* x) override{
	     for (size_t j = 0; j < m_points; j++) {
	        T(x,j) = m_thermo->temperature();
    	    m_thermo->getMassFractions(&Y(x, 0, j));
       }
    }

    void showChiProfile() {
	    std::cout << "Chi profile " << std::endl;
	    for (size_t i = 0; i < m_points; i++) {
	        std::cout << z(i) << "  " << m_chi[i] << std::endl;
	    }
    }

    double chi(size_t j) const {
	    return m_chi[j];
    }

    size_t indexSt() const
    {
	    auto low = std::lower_bound(m_z.begin(), m_z.end(), m_zSt);
	    return (low - m_z.begin());
    }

    virtual double T(const double* x, size_t j) const override {
        return x[index(c_offset_Tflamelet, j)];
    }
    virtual double& T(double* x, size_t j) override {
        return x[index(c_offset_Tflamelet, j)];
    }
    virtual double T_prev(size_t j) const override {
        return prevSoln(c_offset_Tflamelet, j);
    }

    virtual double Y(const double* x, size_t k, size_t j) const override {
        return x[index(c_offset_Yflamelet + k, j)];
    }

    virtual double& Y(double* x, size_t k, size_t j) override {
        return x[index(c_offset_Yflamelet + k, j)];
    }

    virtual double Y_prev(size_t k, size_t j) const override {
        return prevSoln(c_offset_Yflamelet + k, j);
    }

    double d2Ydz2(const double* x, size_t k, size_t j) const {
	double c1 = (Y(x,k,j) - Y(x,k,j-1));
	double c2 = (Y(x,k,j+1) - Y(x,k,j));
	return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    double d2Tdz2(const double* x, size_t j) const {
	double c1 = (T(x,j) - T(x,j-1));
	double c2 = (T(x,j+1) - T(x,j));
	return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    double Lek(size_t k, size_t j) {
       return m_Lek[j*m_nsp+k];
    }

    double exactChi(double zz);

    double erfinv(double p1, double q1);

    // Lewis number stuff
    bool m_do_unityLewisNumber = true;
    vector<double> m_Lek;

    // Thermo stuff
    vector<double> h_RT;
    vector<double> cp_R;


    // Scalar dissipation rate stuff
    vector<double> m_chi;
    vector<double> m_exactchi;
    double m_chiSt = 0.0;
    double m_zSt = 0.0;
    bool m_updateChi = true;

// private:
    // vector<double> m_ybar;
};
    
}

#endif