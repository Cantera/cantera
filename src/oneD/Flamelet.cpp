// ! Flamelet definitions
//

#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"


#include "cantera/base/SolutionArray.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/refine.h"
#include "cantera/transport/Transport.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/numerics/funcs.h"
#include "cantera/base/global.h"
#include "cantera/base/Parser.h"

#include <iostream>

#include "cantera/oneD/Flamelet.h"
// import erf function
using namespace std;

namespace Cantera
{

Flamelet::Flamelet(IdealGasPhase* ph, size_t nsp, size_t nsoot, size_t neq, size_t points) :
	StFlow(ph, nsp, nsoot, neq, points) {
	m_dovisc = true;
	m_updateChi = true;
	m_do_unityLewisNumber = true;
    }

Flamelet::Flamelet(shared_ptr<ThermoPhase> th, size_t nsp, size_t nsoot, size_t neq, size_t points)
    : StFlow(th.get(), nsp, nsoot, neq, points) {
    m_solution = Solution::create();
    m_solution->setThermo(th);
    }

Flamelet::Flamelet(shared_ptr<Solution> sol, const string& id, size_t nsoot, size_t neq, size_t points)
    : StFlow(sol->thermo().get(), sol->thermo()->nSpecies(), nsoot, neq, points) {
    m_solution = sol;
    m_id = id;
    m_kin = m_solution->kinetics().get();
    m_trans = m_solution->transport().get();
        
    m_solution->registerChangedCallback(this, [this]() {
        setKinetics(m_solution->kinetics());
        setTransport(m_solution->transport());
    });

    if (m_trans->transportModel() == "none") {
        // @deprecated
        warn_deprecated("StFlow",
            "An appropriate transport model\nshould be set when instantiating the "
            "Solution ('gas') object.\nImplicit setting of the transport model "
            "is deprecated and\nwill be removed after Cantera 3.0.");
        setTransportModel("mixture-averaged");
        }
    }
    
Flamelet::~Flamelet() {}

void Flamelet::resize(size_t ncomponents, size_t points)
{
    Domain1D::resize(ncomponents, points);
    m_rho.resize(m_points, 0.0);
    m_wtm.resize(m_points, 0.0);
    m_cp.resize(m_points, 0.0);
    m_hr.resize(m_points, 0.0);
    m_h.resize(m_points, 0.0);
    m_tcon.resize(m_points, 0.0);
    m_visc.resize(m_points, 0.0);
    h_RT.resize(m_nsp, 0.0);
    cp_R.resize(m_nsp, 0.0);
    m_hk.resize(m_nsp, m_points, 0.0);
    m_Lek.resize(m_nsp*m_points);

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }

    m_wdot.resize(m_nsp,m_points, 0.0);
    m_do_energy.resize(m_points,false);
    m_qdotRadiation.resize(m_points, 0.0);
    m_fixedtemp.resize(m_points);

    m_dz.resize(m_points-1);
    m_z.resize(m_points);
    m_chi.resize(m_points);
    m_updateChi = true;
}

void Flamelet::eval(size_t jg, double* xg,
                        double* rg, integer* diagg, double rdt)
{
    // if evaluating a Jacobian, and the global point is outside
    // the domain of influence for this domain, then skip
    // evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // if evaluating a Jacobian, compute the steady-state residual
    if (jg != npos) {
        rdt = 0.0;
    }

    // start of local part of global arrays
    double* x = xg + loc();
    double* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;

    if (jg == npos) {      // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else {          // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    // properties are computed for grid points from j0 to j1
    size_t j0 = std::max<size_t>(jmin, 1) - 1;
    size_t j1 = std::min(jmax+1,m_points-1);

    //-----------------------------------------------------
    //              update properties
    //-----------------------------------------------------
    // update thermodynamic properties
    updateThermo(x, j0, j1);

    // update transport properties
    if (jg == npos) {
      updateTransport(x, j0, j1);
    }

    //---------------------------------------------------
    // evaluate the chi profile
    //---------------------------------------------------
    if ( m_updateChi ) setChi();

    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {

            // these may be modified by a boundary object
            rsd[index(c_offset_Tflamelet, j)] =  T(x,j);

            // The default boundary condition for species is zero
            // flux. However, the boundary object may modify this.
            double sum = 0.0;
            for ( size_t k = 0; k < m_nsp; k++) {
                sum += Y(x,k,0);
                rsd[index(c_offset_Yflamelet + k, j)] = Y(x,k,j);
            }
            rsd[index(c_offset_Yflamelet + leftExcessSpecies(), 0)] = 1.0 - sum;
        }

        //----------------------------------------------
        //         right boundary
        //----------------------------------------------

        else if (j == m_points - 1) {

            // the boundary object connected to the right of this
            // one may modify or replace these equations. The
            // default boundary conditions are zero u, V, and T,
            // and zero diffusive flux for all species.
            rsd[index(c_offset_Tflamelet,j)] =  T(x,j);

            double sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
        	sum += Y(x,k,j);
                rsd[index(c_offset_Yflamelet + k ,j)] = Y(x,k,j);
            }
    	    rsd[index(c_offset_Yflamelet + rightExcessSpecies(), j)] = 1.0 - sum;
    	    diag[index(c_offset_Yflamelet + rightExcessSpecies(), j)] = 0;
        }

        //------------------------------------------
        //     interior points
        //------------------------------------------

        else {

            //-------------------------------------------------
            //    Species equations
            //
            //   rho * dY/dt = 1/2 * 1/Le_k *  \rho * \chi * d^2Y/dz^2
	        //         + 1/4 * (1/Le_k - 1 ) * ( d\rho\chi/dz + \rho * \chi * Cp * 1/\lambda * d(\lambda * 1/Cp)/dz )  
            //         + M_k * \omega_k
            //
            //    Temperature equation
            //
            //   rho * dT/dt = 1/2 * \rho * \chi * d^2T/dz^2 + 1/2 * \rho * \chi 1/Cp * dCp/dz * dT/dz + HR
            //                 + 1/2 * \rho * \chi 1/Cp * dT/dz * ∑_k Cp_k dYk/dz
            //-------------------------------------------------
            getWdot(x,j);

	        double drhoChidZ = (m_rho[j+1]*chi(j+1) - m_rho[j-1]*chi(j-1))/(z(j+1) - z(j-1));
	        double dlambdaoverCpdZ = (m_tcon[j+1]/m_cp[j+1] - m_tcon[j-1]/m_cp[j-1])/(z(j+1) - z(j-1));

            for (size_t k = 0; k < m_nsp; k++) {
                rsd[index(c_offset_Yflamelet + k, j)] =
                    m_wt[k]*wdot(k,j)/m_rho[j]
		            + 0.5*chi(j)/Lek(k,j)*d2Ydz2(x,k,j)
		            + 0.25 / m_rho[j] * (1.0/Lek(k,j) - 1.0)
		            * ( drhoChidZ + m_rho[j] * chi(j) * m_cp[j] / m_tcon[j] * dlambdaoverCpdZ )
		            * cdif1(x,c_offset_Yflamelet+k,j)
		            - rdt*(Y(x,k,j) - Y_prev(k,j));
                diag[index(c_offset_Yflamelet + k, j)] = 1;
            }

            // heat release term
            setGas(x,j);
            
            m_thermo->getEnthalpy_RT_ref(&h_RT[0]);
            m_thermo->getCp_R_ref(&cp_R[0]);

            double dCpdZ;
            double sum = 0.0;
            double sum2 = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
               sum  += wdot(k,j) * h_RT[k];
               sum2 += cp_R[k] / m_wt[k] * GasConstant / Lek(k,j)
		       * cdif1(x,c_offset_Yflamelet+k,j);
            }
            sum *= -T(x,j) * GasConstant / m_cp[j];      // consistent units
            m_hr[j] = sum*m_cp[j];
            dCpdZ = (m_cp[j+1] - m_cp[j-1])/(z(j+1) - z(j-1));

            rsd[index(c_offset_Tflamelet, j)] =
                sum/m_rho[j]
	            + 0.5*chi(j)*d2Tdz2(x,j)
                + 0.5*chi(j)/m_cp[j]*dCpdZ*cdif1(x,c_offset_Tflamelet,j)
                + 0.5*chi(j)/m_cp[j]*sum2*cdif1(x,c_offset_Tflamelet,j)
                - rdt*(T(x,j) - T_prev(j));
            diag[index(c_offset_Tflamelet, j)] = 1;
        }
    }
}

string Flamelet::componentName(size_t n) const
{
    switch (n) {
    case c_offset_Tflamelet:
        return "T";
    default:
        if (n >= c_offset_Yflamelet && n < (c_offset_Yflamelet + m_nsp)) {
            return m_thermo->speciesName(n - c_offset_Yflamelet);
        } else {
            return "<unknown>";
        }
    }
}

size_t Flamelet::componentIndex(const string &name) const
{
	if (name=="T") {
		return c_offset_Tflamelet;
	} else {
		for (size_t n=m_neq; n<m_nsp+m_neq; n++) {
			if (componentName(n)==name) {
				return n;
			}
		}
	}
	return npos;
}

/**
 * Set the gas object state to be consistent with the solution at
 * point j.
 */
void Flamelet::setGas(const double* x, size_t j)
{   
	m_thermo->setTemperature(T(x,j));
	const double* yy = x + m_nv*j + c_offset_Yflamelet;
	m_thermo->setMassFractions_NoNorm(yy);
	m_thermo->setPressure(m_press);
}

void Flamelet::resetBadValues(double* xg)
{
    double* x = xg + loc();
    for (size_t j = 0; j < m_points; j++) {
        double* Y = x + m_nv*j + c_offset_Yflamelet;
        m_thermo->setMassFractions(Y);
        m_thermo->getMassFractions(Y);
    }
}

vector<double> Flamelet::dCpdT(const Array2D sol)
{
    vector<double> dCpdT;
    double delta_T = 1.0;
    double cp_l, cp_u;
    double * yspec = new double[m_nsp];
	m_thermo->setPressure(m_press);
    for (size_t n = 0; n < m_points; n++) {
      for (size_t k = 0; k < m_nsp; k++) {
        yspec[k] = sol(k+1,n);
      }
      m_thermo->setMassFractions_NoNorm(yspec);
	  m_thermo->setTemperature(sol(0,n)-delta_T);
      cp_l = m_thermo->cp_mass();
	  m_thermo->setTemperature(sol(0,n)+delta_T);
      cp_u = m_thermo->cp_mass();
      dCpdT.push_back( ( cp_u - cp_l ) / ( 2.0 * delta_T ));
    }
    return dCpdT;
    delete[] yspec;
}

vector<double> Flamelet::d2CpdT2(const Array2D sol)
{
    vector<double> d2CpdT2;
    double delta_T = 1.0;
    double cp_l, cp_c, cp_u;
    double * yspec = new double[m_nsp];
	m_thermo->setPressure(m_press);
    for (size_t n = 0; n < m_points; n++) {
      for (size_t k = 0; k < m_nsp; k++) {
        yspec[k] = sol(k+1,n);
      }
      m_thermo->setMassFractions_NoNorm(yspec);
	  m_thermo->setTemperature(sol(0,n));
      cp_c = m_thermo->cp_mass();
	  m_thermo->setTemperature(sol(0,n)-2.0*delta_T);
      cp_l = m_thermo->cp_mass();
	  m_thermo->setTemperature(sol(0,n)+2.0*delta_T);
      cp_u = m_thermo->cp_mass();
      d2CpdT2.push_back( ( cp_u + cp_l - 2.0*cp_c ) / ( 4.0 * std::pow(delta_T,2.0) ));
    }
    return d2CpdT2;
    delete[] yspec;
}

/**
 * Set the gas state to be consistent with the solution at the
 * midpoint between j and j + 1.
 */
void Flamelet::setGasAtMidpoint(const double* x, size_t j)
{
	m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
	const double* yyj = x + m_nv*j + c_offset_Yflamelet;
	const double* yyjp = x + m_nv*(j+1) + c_offset_Yflamelet;
	for (size_t k = 0; k < m_nsp; k++) {
		m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
	}

	m_thermo->setMassFractions_NoNorm(m_ybar.data());
	m_thermo->setPressure(m_press);
}

void Flamelet::updateTransport(double* x, size_t j0, size_t j1)
{
   for (size_t j = j0; j < j1; j++) {
      setGasAtMidpoint(x,j);
      m_visc[j] = m_trans->viscosity();
      m_tcon[j] = m_trans->thermalConductivity();
      m_trans->getMixDiffCoeffs(&m_diff[j*m_nsp]);
      for (size_t k = 0; k < m_nsp; k++) {
        if ( m_do_unityLewisNumber ) {
          m_Lek[m_nsp*j+k] = 1.0;
        } else {
          m_Lek[m_nsp*j+k] = m_tcon[j] / ( m_rho[j] * m_cp[j] * m_diff[j*m_nsp+k] );
        }
      }
   }
}

void Flamelet::setChi()
{
	m_chi[0]=0.0;
	m_chi[m_points-1]=0.0;
	for (size_t j = 1; j < m_points-1; j++) {
		m_chi[j] = m_chiSt * exactChi(z(j)) / exactChi(m_zSt);
	}
	m_updateChi = false;
}

double Flamelet::exactChi(double zz)
{
	double p1 = fabs(2.0*zz - 1.0);
	double q1 = 1.0 - p1;
	return exp(-2.0 * erfinv(p1,q1) * erfinv(p1,q1));
}

double Flamelet::erfinv(double p1,double q1)
{
	double X,S,T,LNQ,F;
	T = 0.0;
	S = 0.0;
	const double C  = 0.5625;
	const double C1 = 0.87890625;
	const double C2 =-0.2302585092994045684017991454684364e3;
	const double R  = 0.8862269254527580136490837416705726;
	const double A[7] = {0.841467547194693616e-01,0.160499904248262200e1,0.809451641478547505e1,0.164273396973002581e2,0.154297507839223692e2,0.669584134660994039e1,0.108455979679682472e1};
	const double B[7] = {.352281538790042405e-2,0.293409069065309557,0.326709873508963100e1,0.123611641257633210e2,0.207984023857547070e2,0.170791197367677668e2,0.669253523595376683e1};
	const double A1[7] = {0.552755110179178015e2,0.657347545992519152e3,0.124276851197202733e4,0.818859792456464820e3,0.234425632359410093e3,0.299942187305427917e2,0.140496035731853946e1};
	const double B1[6] = {0.179209835890172156e3,0.991315839349539886e3,0.138271033653003487e4,0.764020340925985926e3,0.194354053300991923e3,0.228139510050586581e2};

	const double A2[7] ={0.500926197430588206e1,0.111349802614499199e3,0.353872732756132161e3,0.356000407341490731e3,0.143264457509959760e3,0.240823237485307567e2,0.140496035273226366e1};
	const double B2[6]= {0.209004294324106981e2,0.198607335199741185e3,0.439311287748524270e3,0.355415991280861051e3,0.123303672628828521e3,0.186060775181898848e2};
	const double A3[11]={0.237121026548776092e4,0.732899958728969905e6,0.182063754893444775e7,0.269191299062422172e7,0.304817224671614253e7,0.130643103351072345e7,0.296799076241952125e6,0.457006532030955554e5,0.373449801680687213e4,0.118062255483596543e3,0.100000329157954960e1};
	const double B3[10]={0.851911109952055378e6,0.194746720192729966e7,0.373640079258593694e7,0.397271370110424145e7,0.339457682064283712e7,0.136888294898155938e7,0.303357770911491406e6,0.459721480357533823e5,0.373762573565814355e4,0.118064334590001264e3};
	const double A4[9]={0.154269429680540807e12,0.430207405012067454e12,0.182623446525965017e12,0.248740194409838713e11,0.133506080294978121e10,0.302446226073105850e08,0.285909602878724425e6,0.101789226017835707e04,0.100000004821118676e1};
	const double B4[9]={0.220533001293836387e12,0.347822938010402687e12,0.468373326975152250e12,0.185251723580351631e12,0.249464490520921771e11,0.133587491840784926e10,0.302480682561295591e08,0.285913799407861384e6,0.101789250893050230e04};

	if (p1 <= 0.0 || q1 <= 0.0) {
		return 0.0;
	}
	if (p1 <= 0.75) {
		X = C - p1*p1;
		S =  (((((A[0]*X + A[1])*X + A[2])*X + A[3])*X + A[4])*X + A[5])*X + A[6];
		T = ((((((B[0]*X + B[1])*X + B[2])*X + B[3])*X + B[4])*X + B[5])*X + B[6])*X+1.0;
		X=p1*(S/T);
		F=erf(X)-p1;
		return X-R*exp(X*X)*F;
	}
	//  p = [0,0.75]
	else if (p1 >0.75 && p1 <= 0.9375) {
		X = C1 - p1*p1;
		if (X < 0.1) {
			S = ((((((A1[0]*X + A1[1])*X + A1[2])*X + A1[3])*X+ A1[4])*X + A1[5])*X + A1[6]);
			T = ((((((B1[0]*X + B1[1])*X + B1[2])*X + B1[3])*X+ B1[4])*X + B1[5])*X +1.0);
		}
		else if (X>0.1) {
			S = ((((((A2[0]*X + A2[1])*X + A2[2])*X + A2[3])*X+ A2[4])*X + A2[5])*X + A2[6]);
			T = ((((((B2[0]*X + B2[1])*X + B2[2])*X + B2[3])*X+ B2[4])*X + B2[5])*X + 1.0);
		}
		X = p1*(S/T);
		T = exp(X*X)*erfc(X)-exp(X*X)*q1;
		return X + R*T;
	}
	else if (p1 > 0.9375) {
		LNQ = log(q1);
		X= 1.0 / pow(-LNQ,0.5);
		if (LNQ > C2) {
			S = (((((((((A3[0]*X + A3[1])*X + A3[2])*X + A3[3])*X + A3[4])*X+ A3[5])*X + A3[6])*X + A3[7])*X + A3[8])*X+ A3[9])*X + A3[10];
			T = (((((((((B3[0]*X + B3[1])*X + B3[2])*X + B3[3])*X + B3[4])*X+ B3[5])*X + B3[6])*X + B3[7])*X + B3[8])*X + B3[9])*X + 1.0;
		}
		else if (LNQ < C2) {
			S = (((((((A4[0]*X + A4[1])*X + A4[2])*X + A4[3])*X + A4[4])*X+ A4[5])*X + A4[6])*X + A4[7])*X + A4[8];
			T = ((((((((B4[0]*X + B4[1])*X + B4[2])*X + B4[3])*X + B4[4])*X+ B4[5])*X + B4[6])*X + B4[7])*X + B4[8])*X + 1.0;
		}
		X=S/(X*T);
		T = exp(X*X)*erfc(X);
		F = log(T) - LNQ - X*X;
		return X + R*T*F;
	}
	else {
		return -2.0;
	}

}

}