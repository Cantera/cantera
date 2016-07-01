/**
 *  @file MixTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */
// copyright 2001 California Institute of Technology

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/MixTransport.h"
#include "cantera/numerics/polyfit.h"
#include "cantera/base/utilities.h"
#include "cantera/transport/TransportParams.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"
#include "TransportCharged.h"
#include <cmath>
#include "MMCollisionInt.h"
#include "MMCollisionIntCharged.h"

using namespace std;

namespace Cantera
{
MixTransport::MixTransport() :
    m_condcoeffs(0),
    m_cond(0),
    m_omega12(0),
    m_omega13(0),
    m_omega14(0),
    m_omega15(0),
    m_omega23(0),
    m_omega24(0),
    alpha(0),
    beta(0),
    V(0,0),
    H(0,0),
    DIJ(0,0),
    G(0,0),
    m_lambda(0.0),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_debug(false)
{

}

MixTransport::MixTransport(const MixTransport& right) :
    GasTransport(right),
    m_condcoeffs(0),
    m_cond(0),
    m_omega12(0),
    m_omega13(0),
    m_omega14(0),
    m_omega15(0),
    m_omega23(0),
    m_omega24(0),
    alpha(0),
    beta(0),
    V(0,0),
    H(0,0),
    DIJ(0,0),
    G(0,0),
    m_lambda(0.0),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_debug(false)
{
    *this = right;
}

MixTransport&  MixTransport::operator=(const MixTransport& right)
{
    if (&right == this) {
        return *this;
    }
    GasTransport::operator=(right);

    m_condcoeffs = right.m_condcoeffs;
    m_cond = right.m_cond;
    m_omega12 = right.m_omega12;
    m_omega13 = right.m_omega13;
    m_omega14 = right.m_omega14;
    m_omega15 = right.m_omega15;
    m_omega23 = right.m_omega23;
    m_omega24 = right.m_omega24;
    alpha = right.alpha;
    beta = right.beta;
    V=right.V;
    H=right.H;
    DIJ=right.DIJ;
    G=right.G;
    m_lambda = right.m_lambda;
    m_spcond_ok = right.m_spcond_ok;
    m_condmix_ok = right.m_condmix_ok;
    m_debug = right.m_debug;

    return *this;
}

Transport* MixTransport::duplMyselfAsTransport() const
{
    return new MixTransport(*this);
}

bool MixTransport::initGas(GasTransportParams& tr)
{

    GasTransport::initGas(tr);
    m_eps = tr.eps;
    m_sigma = tr.sigma;
    m_alpha = tr.alpha;
    m_dipole = tr.dipole;
    m_zrot = tr.zrot;
    m_crot = tr.crot;

    // copy polynomials and parameters into local storage
    m_condcoeffs = tr.condcoeffs;
    m_cond.resize(m_nsp);
    m_omega12.resize(m_nsp);
    m_omega13.resize(m_nsp);
    m_omega14.resize(m_nsp);
    m_omega15.resize(m_nsp);
    m_omega23.resize(m_nsp);
    m_omega24.resize(m_nsp);

    // set flags all false
    m_spcond_ok = false;
    m_condmix_ok = false;

    return true;
}

void MixTransport::getMobilities(doublereal* const mobil)
{
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}

doublereal MixTransport::thermalConductivity()
{
    update_T();
    update_C();

    doublereal trThermCond = 0.0;
    // call the new function for the updated translational thermal conductivity
    trThermCond = translationalThermalConductivity();

    return trThermCond;

}



doublereal MixTransport::ElectronTranslationalThermalConductivity()
{

	update_T();
   	update_C();

        MMCollisionIntCharged integrals;
        double t = m_thermo->temperature();
        double p = m_thermo->pressure();

        const int a = m_nsp*(m_nsp+1)/2;
        const int b = m_nsp;
        double cp[m_nsp];
	doublereal m_lambdaEL = 0.0;
        double sum;
        int ij = 0;
        double om11;
        double om22;
        double om12;
        double om13;
        double om14;
        double om15;
        double om23;
        double om24;
        double mol_fract[m_nsp];
        doublereal fac3 = 0.0;

        doublereal D_omega00EE = 0.0;
        doublereal D_omega10EE = 0.0;
        doublereal D_omega20EE = 0.0;
        doublereal D_omega11EE = 0.0;
        doublereal D_omega12EE = 0.0;
        doublereal D_omega22EE = 0.0;

        int numberSpecies;
        numberSpecies = m_nsp;

        size_t ic2 = 0;

        // CK mode not used
        if (m_mode == CK_Mode) {
        for (size_t i = 0; i < m_nsp; i++) {

                om12 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega12Coeff[i]));
                om13 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega13Coeff[i]));
                om14 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega14Coeff[i]));
                om15 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega15Coeff[i]));
                om23 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega23Coeff[i]));
                om24 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega24Coeff[i]));


            for (size_t j = i; j < m_nsp; j++) {


                        fac3 = 8*m_molefracs[i];

                        om22 = exp(dot4(m_polytempvec, m_omega22Coeff[ic2]));
                        om11 = exp(dot4(m_polytempvec, m_omega11Coeff[ic2]));
                        ic2++;

                        if ( j ==  m_thermo->speciesIndex("E") )
                        {

                                D_omega11EE = D_omega11EE + fac3 / 4 * (25 * om11 - 60*om12 + 48*om13);
                                D_omega12EE = D_omega12EE + fac3 / 16 * (175 * om11 - 630*om12 + 912*om13 - 480*om14);
                                D_omega22EE = D_omega22EE + fac3/ 64 * (1225*om11 - 5885*om12 + 12768*om13 - 13440*om14 + 5760*om15);

                        }
                        }
        }

        fac3 = 8*pow(2, 0.5)*m_molefracs[m_thermo->speciesIndex("E")];

	// N.B. omegas at this point are in between el-el
        D_omega11EE = D_omega11EE + fac3*om22;
        D_omega12EE = D_omega12EE + fac3/4 * (7*om22 - 8*om23);
        D_omega22EE = D_omega22EE + fac3/16 * (77*om22 - 112*om23 + 80*om24);


    } else {
      for (size_t i = 0; i < m_nsp; i++) {

                mol_fract[i] = m_molefracs[i];

                if (m_thermo->charge(i) != 0)
                {
                        om12 = integrals.omega12_charged(m_thermo->speciesName(i), "E", m_thermo->charge(i), -1, t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                        om13 = integrals.omega13_charged(m_thermo->speciesName(i), "E", m_thermo->charge(i), -1, t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                        om14 = integrals.omega14_charged(m_thermo->speciesName(i), "E", m_thermo->charge(i), -1, t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                        om15 = integrals.omega15_charged(m_thermo->speciesName(i), "E", m_thermo->charge(i), -1, t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                        om23 = integrals.omega23_charged(m_thermo->speciesName(i), "E", m_thermo->charge(i), -1, t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                        om24 = integrals.omega24_charged(m_thermo->speciesName(i), "E", m_thermo->charge(i), -1, t, t, m_molefracs[m_thermo->speciesIndex("E")],p);

                }

		else
		{

                        om12 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega12Coeff[i]));
                        om13 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega13Coeff[i]));
                        om14 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega14Coeff[i]));
                        om15 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega15Coeff[i]));
                        om23 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega23Coeff[i]));
                        om24 = m_temp * m_temp * exp(dot4(m_polytempvec, m_omega24Coeff[i]));

		}

                if ( mol_fract[i] <= 1e-20)
                { 
			mol_fract[i] = 0.0;
		}

                fac3 = 8*mol_fract[i];

                for (size_t j = i; j < m_nsp; j++) {

			if ( (m_thermo->charge(i) != 0) and (m_thermo->charge(j) != 0) )
                                {

					om22 = integrals.omega22_charged(m_thermo->speciesName(i), m_thermo->speciesName(j), m_thermo->charge(i), m_thermo->charge(j), t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                        		om11 = integrals.omega11_charged(m_thermo->speciesName(i), m_thermo->speciesName(j), m_thermo->charge(i), m_thermo->charge(j), t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                                }

                        else
                                {
                                        om22 = dot5(m_polytempvec, m_omega22Coeff[ic2]);
                                        om11 = dot5(m_polytempvec, m_omega11Coeff[ic2]);
                                }

                        ic2++;



                        if ( (i != m_thermo->speciesIndex("E")) and (j == m_thermo->speciesIndex("E")) )
                        {

                                D_omega11EE = D_omega11EE + fac3/4 * (25*om11 - 60*om12 + 48*om13);
                                D_omega12EE = D_omega12EE + fac3/16 * (175*om11 - 630*om12 + 912*om13 - 480*om14);
                                D_omega22EE = D_omega22EE + fac3/64 * (1225*om11 - 5880*om12 + 12768*om13 - 13440*om14 + 5760*om15);

                        }
                }
        }


	// at this point: om22, om23, om24 are in between electron-electron
        fac3 = 8*pow(2, 0.5)*m_molefracs[m_thermo->speciesIndex("E")];

        D_omega11EE = D_omega11EE + fac3*om22;
        D_omega12EE = D_omega12EE + fac3/4 * (7*om22 - 8*om23);
        D_omega22EE = D_omega22EE + fac3/16 * (77*om22 - 112*om23 + 80*om24);

    }



        const double pre = 9.375;	// 75/8;
        double A;
        double num;

        A = 1/(D_omega11EE - D_omega12EE*D_omega12EE/D_omega22EE);
        num = 2*Pi*GasConstant/1000*m_temp/(0.00000055);
        m_lambdaEL = pre * 1.3806488*pow(10, -23) * m_molefracs[m_thermo->speciesIndex("E")]*sqrt(num)*A;

    return m_lambdaEL;


}


doublereal MixTransport::translationalThermalConductivity()
{

	update_T();
    	update_C();

        double rp = 1/(Boltzmann*m_thermo->temperature());
        const int a = m_nsp*(m_nsp+1)/2;
        const int b = m_nsp;
        double GIJ[a];
        double BIJ[a];

        double mIJ = 0;

        double om22;
        double om11;


        MMCollisionIntCharged integrals;
        double t = m_thermo->temperature();
        double p = m_thermo->pressure();

        DenseMatrix m_astar(b,b);
        DenseMatrix m_bstar(b,b);

        // obtain the Astar and Bstar to compute thermal conductivity with Gupta-Yos mixture rule
        size_t ic2 = 0;
        if (m_mode == CK_Mode) {
        for (size_t i = 0; i < m_nsp; i++) {
            for (size_t j = i; j < m_nsp; j++) {


		if ( ( m_thermo->charge(i) != 0) and ( m_thermo->charge(j) != 0 ) )
                {
                  om22 = integrals.omega22_charged(m_thermo->speciesName(i), m_thermo->speciesName(j), m_thermo->charge(i), m_thermo->charge(j), t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                  om11 = integrals.omega11_charged(m_thermo->speciesName(i), m_thermo->speciesName(j), m_thermo->charge(i), m_thermo->charge(j), t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
		  m_bstar(i,j) = integrals.bstar_charged(m_thermo->speciesName(i), m_thermo->speciesName(j), m_thermo->charge(i), m_thermo->charge(j), t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                  m_bstar(j,i) = m_bstar(i,j);
                  m_astar(i,j) = om22/om11;
                  m_astar(j,i) = m_astar(i,j);
                }

                else
                {
                  m_astar(i,j) = exp(dot4(m_polytempvec, m_omega22Coeff[ic2])) / exp(dot4(m_polytempvec, m_omega11Coeff[ic2]));
                  m_astar(j,i) = m_astar(i,j);
		  m_bstar(i,j) = exp(dot4(m_polytempvec, m_bstarCoeff[ic2]));
                  m_bstar(j,i) = m_bstar(i,j);
                  om22 = exp(dot4(m_polytempvec, m_omega22Coeff[ic2]));
                  om11 = exp(dot4(m_polytempvec, m_omega11Coeff[ic2]));
                }

		ic2++;
            }
        }
    } else {
        for (size_t i = 0; i < m_nsp; i++) {
            for (size_t j = i; j < m_nsp; j++) {

                if ( ( m_thermo->charge(i) != 0) and ( m_thermo->charge(j) != 0 ) )
                {
                  om22 = integrals.omega22_charged(m_thermo->speciesName(i), m_thermo->speciesName(j), m_thermo->charge(i), m_thermo->charge(j), t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                  om11 = integrals.omega11_charged(m_thermo->speciesName(i), m_thermo->speciesName(j), m_thermo->charge(i), m_thermo->charge(j), t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                  m_bstar(i,j) = integrals.bstar_charged(m_thermo->speciesName(i), m_thermo->speciesName(j), m_thermo->charge(i), m_thermo->charge(j), t, t, m_molefracs[m_thermo->speciesIndex("E")],p);
                  m_bstar(j,i) = m_bstar(i,j);
                  m_astar(i,j) = om22/om11;
                  m_astar(j,i) = m_astar(i,j);
                }

                else
                {
                  m_astar(i,j) = dot5(m_polytempvec, m_omega22Coeff[ic2])  / dot5(m_polytempvec, m_omega11Coeff[ic2]);
                  m_astar(j,i) = m_astar(i,j);
                  m_bstar(i,j) = exp(dot4(m_polytempvec, m_bstarCoeff[ic2]));
                  m_bstar(j,i) = m_bstar(i,j);
                  om22 = dot5(m_polytempvec, m_omega22Coeff[ic2])/(m_temp);
                  om11 = dot5(m_polytempvec, m_omega11Coeff[ic2])/(m_temp);
                }

                ic2++;
            }
        }
    }


        int numberSpecies;
        numberSpecies = m_nsp;
	if (m_molefracs[m_thermo->speciesIndex("E")] >= 1e-20)
	{m_nsp = m_nsp-1;}
	else
	{m_nsp = m_nsp;}


        int ij =0;
        int ic = 0;
        for (size_t i = 0; i < m_nsp; i++) {
                for (size_t j = i; j < m_nsp; j++) {

                ij = ((i+1-1)*(2*m_nsp-i-1)+2*(j+1))/2;
		ij = ij - 1;

                ic = ij;
                mIJ = (m_mw[i]*m_mw[j]) / ( (m_mw[i] + m_mw[j]) * (m_mw[i] + m_mw[j]) );
		GIJ[ic] = mIJ / (75 * Boltzmann * rp*m_bdiff(i,j) ) * (165 - 36*m_bstar(i,j) - 48*m_astar(i,j) ) ;

        	}
       }


        DenseMatrix B(m_nsp,m_nsp);
        ij = 0;
        for (int i = 0; i < m_nsp; i++) {

                for (int j = 0; j < i; j++) {
                        ij = ((j+1-1)*(2*m_nsp-j-1)+2*(i+1))/2;
                	ij = ij - 1;

                        mIJ = (m_mw[i]*m_mw[j]) / ( (m_mw[i] + m_mw[j]) * (m_mw[i] + m_mw[j]) );
                        B(i,j) = mIJ / (rp*m_bdiff(i,j)) * (9.6 * m_astar(i,j) + (m_mw[i] - m_mw[j]) * (9/m_mw[j] + (-7.5 + 3.6*m_bstar(i,j))/m_mw[i]));

                }

                for (int j = i; j < m_nsp; j++) {
                        ij = ((i+1-1)*(2*m_nsp-i-1)+2*(j+1))/2;
                	ij = ij - 1;

                        mIJ = (m_mw[i]*m_mw[j]) / ( (m_mw[i] + m_mw[j]) * (m_mw[i] + m_mw[j]) );
                        B(i,j) = mIJ / (rp*m_bdiff(i,j)) * (9.6 * m_astar(i,j) + (m_mw[i] - m_mw[j]) * (9/m_mw[j] + (-7.5 + 3.6*m_bstar(i,j))/m_mw[i]));

                }
        }


        const double FAC = 2/(15*Boltzmann);
        double GII[m_nsp];

        for (size_t i = 0; i < m_nsp; i++) {
                GII[i] = 0;

                for (size_t j = 0; j < m_nsp; j++) {

                GII[i] = GII[i] + B(i,j)*m_molefracs[j]*FAC;

        	}
        }


        double fac = 0;
        double sum1 = 0;
        double sum2 = 0;
        ic = 0;

        for (size_t i = 0; i < m_nsp-1; i++) {                          
                for (size_t j = i+1; j < m_nsp; j++) {                 

                ij = ((i+1-1)*(2*m_nsp-i-1)+2*(j+1))/2;
                ij = ij - 1;

                ic = ij;

                fac = 2 * m_molefracs[i]*m_molefracs[j]*( 1/GII[i] - 1/GII[j])*( 1/GII[i] - 1/GII[j]);
                sum1 = sum1 + fac;
                sum2 = sum2 + fac*GIJ[ic];

                }
        }


        double gav = 0;
        gav = sum2/sum1;
        double sum;
        sum =0;
        for (size_t i = 0; i < m_nsp; i++) {
                sum = sum + m_molefracs[i]/(GII[i] + gav);
        }

        m_nsp = numberSpecies;

        return sum / (1- gav*sum);


}




void MixTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    for (size_t k = 0; k < m_nsp; k++) {
        dt[k] = 0.0;
    }
}

void MixTransport::getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                    size_t ldx, const doublereal* const grad_X,
                                    size_t ldf, doublereal* const fluxes)
{

    update_T();
    update_C();


    getMixDiffCoeffs(DATA_PTR(m_spwork));

    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* y  = m_thermo->massFractions();
    doublereal rhon = m_thermo->molarDensity();

    const int numS= m_thermo->spPosIndex.size();	// number of ionized species
    double gamma_s[numS];              // species molar concentration
    double denom[numS];
    double num[numS];
	for (size_t k = 0; k < numS; k++) {
		gamma_s[k] = 0.0;
        	denom[k] = 0.0;
        	num[k] = 0.0;
	}
        double numerator = 0.0;
        double denominator = 0.0;
	int c=0;



// correction for Ambipolar Diffusion
vector_fp sum(ndim,0.0);
    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {

		if ( m_thermo->charge(k) > 0 )
        	{
        		m_spwork[k] = m_spwork[k]*2;
        		gamma_s[c] = m_thermo->concentration(k)/m_thermo->density();
        		denom[c] = gamma_s[c]*m_thermo->molecularWeight(k)/Avogadro;
        		num[c] = m_spwork[k] * gamma_s[c];
        		c = c + 1;
        	}
	}
    }

    for (size_t j = 0; j < numS; j++)
    {
    	numerator = numerator + num[j];
    	denominator = denominator + denom[j];
    }

    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
		
		if (m_thermo->charge(k) > 0 )
             	{
              		fluxes[n*ldf + k] = -rhon * mw[k] * m_spwork[k] * grad_X[n*ldx + k];
                	sum[n] += fluxes[n*ldf + k];
                }

		else if  ( k == m_thermo->speciesIndex("E") )
		{
			m_spwork[k] = (m_thermo->molecularWeight(k)/Avogadro) * numerator / denominator;
                                fluxes[n*ldf + k] = -rhon * mw[k] * m_spwork[k] * grad_X[n*ldx + k];
                                sum[n] += fluxes[n*ldf + k];
		}

		else 
		{
            		fluxes[n*ldf + k] = -rhon * mw[k] * m_spwork[k] * grad_X[n*ldx + k];
           		sum[n] += fluxes[n*ldf + k];
		}

        }
    }
    // add correction flux to enforce sum to zero
    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] -= y[k]*sum[n];

        }
    }

}


void MixTransport::getSpeciesFluxesSM(size_t ndim, const doublereal* const grad_T,
                                    size_t ldx, const doublereal* const grad_X,
                                    size_t ldf, doublereal* const fluxes)
{


    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

	double nd =0;
	double t = m_thermo->temperature();
	double rp = 1/(Boltzmann*t);
	TransportCharged trChCh;
	nd = trChCh.getNumberDensity(t, t, m_molefracs[m_thermo->speciesIndex("E")], m_thermo->pressure());
	double atomicW[m_nsp];
	m_thermo->getMolecularWeights(atomicW);		// kg/kmol
	double MM = 0;
        doublereal mol_fract[m_nsp];
	double charge[m_nsp];
	doublereal TOL = 1e-15;
	doublereal SUM = 0.0;


		for (size_t k = 0; k < m_nsp; k++) {

			mol_fract[k] = m_molefracs[k] + TOL;
			SUM = SUM + mol_fract[k];
			charge[k] = m_thermo->charge(k);
		}


		for (size_t k = 0; k < m_nsp; k++) {

			mol_fract[k] = mol_fract[k]/SUM;

		}



	MM = m_thermo->meanMolecularWeight();

	doublereal QT = 0.0;
	doublereal M[m_nsp];
	doublereal R[m_nsp];
        DIJ.resize(m_nsp,m_nsp);
        G.resize(m_nsp+1,m_nsp+1);
	beta.resize(m_nsp+1);
	doublereal init = 0.0;

// Initialize to zero the new vectors and matrices
for (int count = 0; count < m_nsp; count++)
{
	M[count] = init;
	R[count] = init;
	
	for (int count2 = 0; count2 < m_nsp; count2++)
	{
		DIJ(count,count2) = init;
	}
}

for (int count = 0; count < m_nsp+1; count++)
{
        beta[count] = init;

        for (int count2 = 0; count2 < m_nsp+1; count2++)
        {
                G(count,count2) = init;
        }
}


// System for V_i
    for (size_t n = 0; n < ndim; n++) {
	for (size_t k = 0; k < m_nsp; k++) {

		beta[k] = - grad_X[n*ldx + k];
		M[k] = mol_fract[k]*atomicW[k]/MM;		
		R[k] = 1.0;
		QT = QT + mol_fract[k]*charge[k];
	
        }
		beta[m_nsp] = 0.0;
    }	

    doublereal SF = 0.0;
    doublereal KA[m_nsp];
    const doublereal E = 1.60219099e-19 ;  // electron charge

	for (size_t i = 0; i < m_nsp; i++) {
		
		KA[i] =  mol_fract[i]* (charge[i] - atomicW[i]/MM * QT) * E / ( Boltzmann * t);
		SF = SF + KA[i]*KA[i];
		DIJ(i,i) = 1.0;
		G(i,i) = 0.0;
	
		for (size_t j = i+1; j < m_nsp; j++) {

			DIJ(i,j) = rp*m_bdiff(i,j)/nd;
			DIJ(j,i) = DIJ(i,j);	
			G(i,j) = 0.0;
			G(j,i) = 0.0;

		}

	}

	SF = sqrt(SF);


	// Electron Contribution
	size_t i = m_thermo->speciesIndex("E");
	for (size_t j = 0; j < m_nsp - 1; j++) {	

		G(i,j) = - mol_fract[i]*mol_fract[j]/DIJ(i,j);
		G(i,i) = G(i,i) - G(i,j);
		G(j,i) = G(i,j);
		G(j,j) = G(j,j) - G(i,j);


	}
	size_t j = m_nsp;				
		G(i,j) = -KA[i]/(SF);


	// Heavy particle contribution
	for ( i = 0; i < m_nsp - 1; i++) {
		for ( j = i+1; j < m_nsp - 1; j++) {

			G(i,j) = - mol_fract[i]*mol_fract[j]/DIJ(i,j);
			G(j,i) = G(i,j);
			G(i,i) = G(i,i) - G(i,j);
			G(j,j) = G(j,j) - G(i,j);
		}

        j = m_nsp;                                              
                G(i,j) = -KA[i]/(SF);

	}

	// Ambipolar constraint
	i = m_nsp;					
	for ( j = 0; j < m_nsp; j++) {
		G(i,j) = -KA[j]/(SF);

	}

	G(m_nsp,m_nsp) = 0.0;				


	// GMRES iterative solution
	doublereal P[m_nsp+1];
	for ( i = 0; i < m_nsp; i++) {
		P[i] = G(i,i);
	}
	P[m_nsp] = 1.0;				


	const int iterMax =50;
	const double tolRes =0.000001;

	int N = m_nsp + 1;

	alpha.resize(N);
	doublereal norm0 = 0.0;

	V.resize(N,iterMax+1);
	H.resize(iterMax+1,iterMax);


	// Initialize to zero the new vectors and matrices
	for ( int count=0; count < N; count++)
	{
		alpha[count] = 0.0;

		for ( int count2=0; count2 < iterMax+1; count2++)
		{
			V(count,count2) =init;
		}
	}

	for ( int count=0; count < iterMax+1; count++)
        {
                for ( int count2=0; count2 < iterMax; count2++)
                {	
                        H(count,count2) =init;
                }
        }


	// Initial guess is zero
	for ( j = 0; j < N; j++) {

		alpha[j] = init;
		V(j,0) = beta[j];
		norm0 = norm0 + V(j,0)*V(j,0);

	}
	norm0 = sqrt(norm0);


	// First term of Hessenberg system
	if (norm0 != 0)
	{
		for (j = 0; j < N; j++) {
			V(j,0) = V(j,0) / norm0;
		}
	
	}

	for ( j = 0; j < iterMax+1; j++) {
		for (size_t k = 0; k < iterMax; k++) {
			H(j,k) =init;
		}
	}

	// Iterative loop
	doublereal norm = 0.0;
	norm = norm0;
	int I = -1;
	doublereal R2[iterMax+1];
	R2[0] = norm;
	int I1 = 0;
	doublereal Z[N];
	j = 0;
	doublereal C[iterMax];
	doublereal S[iterMax];
	doublereal lowN =init;
	doublereal GAM =init;
        doublereal T = 0.0;
	int k1 = 0;

	for (int l=0; l<iterMax; l++)
	{
		C[l] = lowN;
		S[l] = lowN;

	}


	while ((I<iterMax-1) and (norm > norm0*tolRes)) {

		I = I+1;
		I1 = I+1;

		for ( j = 0; j < N; j++) {
			Z[j] = init;
                	Z[j] = V(j,I)/P[j];

        	}


		for ( j = 0; j < N; j++) {

                        V(j,I1) =init;
			for (size_t k = 0; k < N; k++) {
				V(j,I1) = V(j,I1) + G(j,k)*Z[k];
			}
		}


		// Modified Graham-Schmidt and Arnoldi step
		for ( j = 0; j <= I; j++) {		
			T = 0.0;

			for (size_t k = 0; k < N; k++) {
				T = T + V(k,j)*V(k,I1);

			}


			H(j,I) = T;


			for (size_t k = 0; k < N; k++) {

                                V(k, I1) = V(k, I1) -T*V(k,j);
                        
			}

		}

		T = 0.0;
		for ( j = 0; j < N; j++) {
			T = T + V(j, I1)*V(j, I1);

		}

		T = sqrt(T);
		H(I1,I) = T;

		if ( T!= 0.0)
		{
			for ( j = 0; j < N; j++) {
                        	V(j, I1) = V(j, I1)/T;
                	}
		}		



		// Factorization of H
		// Perform previous transformation on i-th column of H
		for (size_t k = 1; k <= I; k++) {	
				k1 = k-1;
				T = H(k1,I);
				H(k1,I) = C[k1]*T + S[k1]*H(k,I);
				H(k,I) = -S[k1]*T + C[k1]*H(k,I);
                        }
			GAM = sqrt(H(I,I)*H(I,I) + H(I1,I)*H(I1,I));

			if (GAM == 0.0)
			{
				GAM = 1e-32;
			}

			// Next plane rotation
			C[I] = H(I,I)/GAM;
			S[I] = H(I1,I)/GAM;
			R2[I1] = -S[I]*R2[I];
			R2[I] = C[I]*R2[I];


			// Residual norm
			H(I,I) = C[I]*H(I,I) + S[I]*H(I1,I);
			norm = abs(R2[I1]);


		}


	// New solution: solve for upper triangular system
	for ( int k = I; k>=0; k--)			
	{
		T = 0;
		for ( int l = k+1; l<=I; l++)
        	{
			T = T + H(k,l)*R2[l];


		}
		for ( int l = I+1; l<=k; l++)
                {
                        T = T + H(k,l)*R2[l];


                }
		R2[k] = (R2[k] - T) / H(k,k);


	}

	// Form linear combination of V(*,I)'s to get solution
	for ( j = 0; j < N; j++) {

		T = 0.0;
		for ( size_t k = 0; k<=I; k++)
        	{			
			T = T + V(j,k)*R2[k];
		}
		
		Z[j] = T / P[j];


	}

	for ( j = 0; j < N; j++) {
		alpha[j] = alpha[j] + Z[j];
	}

	
	// Projection
	doublereal RM = 0.0;
	doublereal MX = 0.0;

	for (i=0; i < m_nsp; i++)
	{		
		RM = RM + M[i]*R[i];
		MX = MX + M[i]*alpha[i];
	}

	doublereal PROJ =0;
	PROJ = MX/RM;


	for (i=0; i < m_nsp; i++)
        {				
		fluxes[i] = alpha[i] - PROJ*R[i];
	}

	// Mass diffusion fluxes
        for (i=0; i < m_nsp; i++)
        {
		fluxes[i] = mol_fract[i]*fluxes[i]*(atomicW[i]/1000)/(6.0221367e23) * nd;
	}



	// Ambipolar field (scaling factor correction)
	doublereal EAMB=0;
      	EAMB = alpha[N-1]/SF;
}



void MixTransport::getSpeciesFluxesNeutSM(size_t ndim, const doublereal* const grad_T,
                                    size_t ldx, const doublereal* const grad_X,
                                    size_t ldf, doublereal* const fluxes)
{

    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

        double nd = 0;
        double t = m_thermo->temperature();
        double rp = 1/(Boltzmann*t);

	TransportCharged trChCh;
        nd = trChCh.getNumberDensity(t, t, m_molefracs[m_thermo->speciesIndex("E")], m_thermo->pressure());
        double atomicW[m_nsp];
        m_thermo->getMolecularWeights(atomicW);         // kg/kmol
        double MM = 0;
        doublereal mol_fract[m_nsp];
        doublereal P[m_nsp];
        doublereal V[m_nsp];
        doublereal TOL = 1e-15;
        doublereal SUM = 0.0;
	doublereal beta[m_nsp];
	doublereal G[(m_nsp+1)*m_nsp/2];
	int II = 0;
	int IJ = 0;
	int JJ = 0;
	doublereal FAC = 0.0;

                for (size_t k = 0; k < m_nsp; k++) {

                        mol_fract[k] = m_molefracs[k] + TOL;
                        SUM = SUM + mol_fract[k];

                }


                for (size_t k = 0; k < m_nsp; k++) {

                        mol_fract[k] = mol_fract[k]/SUM;

                }



        MM = m_thermo->meanMolecularWeight();

        doublereal M[m_nsp];
        doublereal R[m_nsp];
        doublereal init = 0.0;


// Initialize to zero the new vectors and matrices
for (int count = 0; count < m_nsp; count++)
{
        M[count] = init;
        R[count] = init;
	beta[count] = init;
}


// System for V_i
    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
                beta[k] = - grad_X[n*ldx + k];
                M[k] = mol_fract[k]*atomicW[k]/MM;
                R[k] = 1.0;
		II = ((k+1-1)*(2*m_nsp-k-1)+2*(k+1))/2;
		II = II - 1;

 		G[II] = 0.0;
        }
                beta[m_nsp] = beta[m_nsp];
    }


        for (size_t i = 0; i < m_nsp; i++) {

                for (size_t j = i+1; j < m_nsp; j++) {

			IJ = ((i+1-1)*(2*m_nsp-i-1)+2*(j+1))/2;
			II = ((i+1-1)*(2*m_nsp-i-1)+2*(i+1))/2;
			JJ = ((j+1-1)*(2*m_nsp-j-1)+2*(j+1))/2;
			IJ = IJ - 1;
			II = II - 1;
			JJ = JJ - 1;

			FAC = mol_fract[i]*mol_fract[j]*nd/(rp*m_bdiff(i,j));
			G[IJ] = - FAC;
			G[II] = G[II] + FAC;
			G[JJ] = G[JJ] + FAC;
                }
        }

	
	// Krylov iterative solution
	for (size_t i = 0; i < m_nsp; i++) {
		
			II = ((i+1-1)*(2*m_nsp-i-1)+2*(i+1))/2;
			II = II - 1;
			P[i] = G[II];
	}

	int iterMax = 100;
	double tolRes = 0.0001;


	// Conjugate Gradient with preconditioning
	doublereal tol = tolRes*tolRes;
	int iter = 0;
	doublereal ZN[m_nsp];
        doublereal DMI[m_nsp];
        doublereal arnold[m_nsp];
        doublereal TEMP[m_nsp];
	bool loopIN;
	doublereal res = 1.0;
	doublereal error = 1.0;
	

        for (size_t i = 0; i < m_nsp; i++) {
	
		V[i] = 0.0;
		ZN[i] = 0.0;
		DMI[i] = 1/P[i];
	}

	doublereal betaN = 0.0;
	doublereal AAA = 0.0;
        doublereal BBB = 0.0;
        doublereal CCC = 0.0;
	doublereal dif = 0.0;


	for (size_t i = 0; i < m_nsp; i++) {

		AAA = AAA + DMI[i]*beta[i]*beta[i];
	}

	doublereal RESINI = 0.0;
	RESINI = AAA;

	if (RESINI < 1e-32)
		{ loopIN = false;}
	else
        	{ loopIN = true;}


	while ((iter < iterMax) and (loopIN == true) and (res > tol)) {

		iter = iter +1;

		for (size_t i = 0; i < m_nsp; i++) {
		
			arnold[i] = V[i];

		}

		for (size_t i = 0; i < m_nsp; i++) {

                        ZN[i] = DMI[i]*beta[i] + betaN*ZN[i];

                }

		// Matrix vector multiplication
		int III = 0;
		for (size_t i = 0; i < m_nsp; i++) {
			
			III = m_nsp*(i+1-1) - ((i+1)*(i+1-1))/2 + (i+1);
			III = III - 1;

			TEMP[i] = G[III]*ZN[i];
		}

		for (size_t i = 0; i < m_nsp; i++) {

			for (size_t j = i+1; j < m_nsp; j++) {

				III = m_nsp*(i+1-1) - ((i+1)*(i+1-1))/2 + (j+1);
				III = III - 1;

				TEMP[i] = TEMP[i] + G[III]*ZN[j];
				TEMP[j] = TEMP[j] + G[III]*ZN[i];

			}
                }


		// Scalar product
		doublereal product = 0.0;
		for (size_t i = 0; i < m_nsp; i++) {
		
			product = product + ZN[i]*TEMP[i];

		}

		BBB = product;

		if ( BBB <= 0.0)
		{
			iter = iter -1;
			loopIN = false;
		}
		
		else
		{
			for (size_t i = 0; i < m_nsp; i++) {
			
				V[i] = V[i] + AAA/BBB*ZN[i];
				beta[i] = beta[i] - AAA/BBB*TEMP[i];
			}

		}

		CCC = 0.0;
		for (size_t i = 0; i < m_nsp; i++) {

			CCC = CCC + DMI[i]*beta[i]*beta[i];
		}
	
		res = CCC/RESINI;

		if ( res > tol)
		{
			betaN = CCC/AAA;

		}

		AAA = CCC;
		error = 0.0;

		for (size_t i = 0; i < m_nsp; i++) {
		
			dif = V[i] - arnold[i];
			error = error + dif*dif; 

		}


	} // end while


	// Constraint 
	doublereal RM = 0.0;
        doublereal MX = 0.0;
        doublereal PROJ = 0.0;

	for (size_t i = 0; i < m_nsp; i++) {

                RM = RM + M[i]*R[i];
                MX = MX + M[i]*V[i];


	}

	PROJ = MX/RM;

	for (size_t i = 0; i < m_nsp; i++) {

		V[i] = V[i] - PROJ*R[i];

	}


	// Mass diffusion fluxes
	for (size_t i = 0; i < m_nsp; i++) {

		fluxes[i] = V[i] * mol_fract[i]*(atomicW[i]/1000)/(6.0221367e23) * nd;
	}

}


void MixTransport::getSpeciesFluxesSMEamb(size_t ndim, const doublereal* const grad_T,
                                    size_t ldx, const doublereal* const grad_X,
                                    size_t ldf, doublereal Eamb, doublereal* const fluxes)
{

    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

        double nd = 0;
        double t = m_thermo->temperature();
        double rp = 1/(Boltzmann*t);

	TransportCharged trChCh;
        nd = trChCh.getNumberDensity(t, t, m_molefracs[m_thermo->speciesIndex("E")], m_thermo->pressure());
        double atomicW[m_nsp];
        m_thermo->getMolecularWeights(atomicW);         // kg/kmol
        double MM = 0;
        MM = m_thermo->meanMolecularWeight();
	const doublereal Q = 1.60219099e-19 ;  // electron charge
        doublereal mol_fract[m_nsp];
        doublereal mass_fract[m_nsp];
        doublereal P[m_nsp];
        doublereal V[m_nsp];
        doublereal TOL = 1e-15;
        doublereal SUM = 0.0;
	doublereal beta[m_nsp];
	doublereal G[(m_nsp+1)*m_nsp/2];
	int II = 0;
	int IJ = 0;
	int JJ = 0;
	doublereal FAC = 0.0;
	doublereal k_amb[m_nsp];
        doublereal charge[m_nsp];
	doublereal tot_charge = 0.0;
	doublereal Y_i[m_nsp];
	doublereal N_ns[m_nsp];		// vector null space
	doublereal normY = 0.0;		// normalization for N_ns
	DenseMatrix NN(m_nsp,m_nsp);
	doublereal k_proj[m_nsp];

                for (size_t k = 0; k < m_nsp; k++) {

                        mol_fract[k] = m_molefracs[k] + TOL;
                        SUM = SUM + mol_fract[k];

                }


                for (size_t k = 0; k < m_nsp; k++) {

                        mol_fract[k] = mol_fract[k]/SUM;
			mass_fract[k] = mol_fract[k] * atomicW[k]/MM;
			charge[k] = m_thermo->charge(k);
			tot_charge = tot_charge + mol_fract[k]*charge[k];
                }

		for (size_t k = 0; k < m_nsp; k++) {

                        k_amb[k] = (mol_fract[k]*charge[k] - mass_fract[k]*tot_charge)*Q*rp;
			Y_i[k] = mass_fract[k];
                        normY = normY + Y_i[k]*Y_i[k];

                }

		normY = sqrt(normY);

		for (size_t k = 0; k < m_nsp; k++) {

			N_ns[k] = Y_i[k] / normY;	

		}

		for (size_t i = 0; i < m_nsp; i++) {
            		for (size_t j = i; j < m_nsp; j++) {
			
				NN(i,j) = N_ns[i]*N_ns[j];
				NN(j,i) = NN(i,j);
			}

                }

		multiply(NN, k_amb, k_proj);

        doublereal M[m_nsp];
        doublereal R[m_nsp];
        doublereal init = 0.0;

// Initialize to zero the new vectors and matrices
for (int count = 0; count < m_nsp; count++)
{
        M[count] = init;
        R[count] = init;
	beta[count] = init;
}



// System for V_i
    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {

		beta[k] = - grad_X[n*ldx + k] + (k_amb[k]*Eamb -  k_proj[k]*Eamb);
                M[k] = mol_fract[k]*atomicW[k]/MM;
                R[k] = 1.0;
		II = ((k+1-1)*(2*m_nsp-k-1)+2*(k+1))/2;
		II = II -1;
 		G[II] = 0.0;

        }
                beta[m_nsp-1] = beta[m_nsp-1];
    }


        for (size_t i = 0; i < m_nsp; i++) {
        	for (size_t j = i+1; j < m_nsp; j++) {

			IJ = ((i+1-1)*(2*m_nsp-i-1)+2*(j+1))/2;
			II = ((i+1-1)*(2*m_nsp-i-1)+2*(i+1))/2;
			JJ = ((j+1-1)*(2*m_nsp-j-1)+2*(j+1))/2;
			IJ = IJ -1;
			II = II -1;
			JJ = JJ -1;
			FAC = mol_fract[i]*mol_fract[j]*nd/(rp*m_bdiff(i,j));

			G[IJ] = - FAC;
			G[II] = G[II] + FAC;
			G[JJ] = G[JJ] + FAC;
                }
        }


size_t countG = 0;
        for (size_t i = 0; i < m_nsp; i++) {
                for (size_t j = i; j < m_nsp; j++) {
			G[countG] = G[countG] +  NN(i,j);
			countG++;
		}
	}

	
	// Krylov iterative solution
	for (size_t i = 0; i < m_nsp; i++) {
		
			II = ((i+1-1)*(2*m_nsp-i-1)+2*(i+1))/2;
			II = II -1;

			P[i] = G[II];
	}

	int iterMax = 50;
	doublereal tolRes =  0.000001;


	// Conjugate Gradient with preconditioning
	doublereal tol = tolRes*tolRes;
	int iter = 0;
	doublereal ZN[m_nsp];
        doublereal DMI[m_nsp];
        doublereal arnold[m_nsp];
        doublereal TEMP[m_nsp];
	bool loopIN;
	doublereal res = 1.0;
	doublereal error = 1.0;

        for (size_t i = 0; i < m_nsp; i++) {
	
		V[i] = 0.0;
		ZN[i] = 0.0;
		DMI[i] = 1.0/P[i];
	}

	doublereal betaN = 0.0;
	doublereal AAA = 0.0;
        doublereal BBB = 0.0;
        doublereal CCC = 0.0;
	doublereal dif = 0.0;


	for (size_t i = 0; i < m_nsp; i++) {

		AAA = AAA + DMI[i]*beta[i]*beta[i];
	}

	doublereal RESINI = 0.0;
	RESINI = AAA;

	if (RESINI < 1e-32)
		{ loopIN = false;}
	else
        	{ loopIN = true;}


	while ((iter < iterMax) and (loopIN == true) and (res > tol)) {

		iter = iter +1;

		for (size_t i = 0; i < m_nsp; i++) {
		
			arnold[i] = V[i];

		}

		for (size_t i = 0; i < m_nsp; i++) {

                        ZN[i] = DMI[i]*beta[i] + betaN*ZN[i];

                }

		// Matrix vector multiplication
		int III = 0;
		for (size_t i = 0; i < m_nsp; i++) {
			
			III = m_nsp*(i+1-1) - ((i+1)*(i+1-1))/2 + (i+1);
			III = III - 1;

			TEMP[i] = G[III]*ZN[i];
		}

		for (size_t i = 0; i < m_nsp; i++) {

			for (size_t j = i+1; j < m_nsp; j++) {

				III = m_nsp*(i+1-1) - ((i+1)*(i+1-1))/2 + (j+1);
				III = III - 1;

				TEMP[i] = TEMP[i] + G[III]*ZN[j];
				TEMP[j] = TEMP[j] + G[III]*ZN[i];

			}
                }


		// Scalar product
		doublereal product = 0.0;
		for (size_t i = 0; i < m_nsp; i++) {
		
			product = product + ZN[i]*TEMP[i];

		}

		BBB = product;

		if ( BBB <= 0.0)
		{
			iter = iter -1;
			loopIN = false;
		}
		
		else
		{
			for (size_t i = 0; i < m_nsp; i++) {
				V[i] = V[i] + AAA/BBB*ZN[i];
				beta[i] = beta[i] - AAA/BBB*TEMP[i];
			}
		}

		CCC = 0.0;
		for (size_t i = 0; i < m_nsp; i++) {

			CCC = CCC + DMI[i]*beta[i]*beta[i];
		}
		res = CCC/RESINI;

		if ( res > tol)
		{
			betaN = CCC/AAA;
		}

		AAA = CCC;
		error = 0.0;

		for (size_t i = 0; i < m_nsp; i++) {
		
			dif = V[i] - arnold[i];
			error = error + dif*dif; 
		}


	} // end while


	// Constraint 
	doublereal RM = 0.0;
        doublereal MX = 0.0;
        doublereal PROJ = 0.0;

	for (size_t i = 0; i < m_nsp; i++) {

                RM = RM + M[i]*R[i];
                MX = MX + M[i]*V[i];

	}

	PROJ = MX/RM;

	for (size_t i = 0; i < m_nsp; i++) {

		V[i] = V[i] - PROJ*R[i];

	}

	// Mass diffusion fluxes
	for (size_t i = 0; i < m_nsp; i++) {

		fluxes[i] = V[i] * mol_fract[i]*(atomicW[i]/1000)/(6.0221367e23) * nd;
	}
}


void MixTransport::update_T()
{
    doublereal t = m_thermo->temperature();
    if (t == m_temp) {
        return;
    }
    if (t < 0.0) {
        throw CanteraError("MixTransport::update_T",
                           "negative temperature "+fp2str(t));
    }
    GasTransport::update_T();
    // temperature has changed, so polynomial fits will need to be redone.
    m_spcond_ok = false;
    m_bindiff_ok = false;
    m_condmix_ok = false;
}




// TO BE CHANGED FOR m_cond because the interactions charged-charged are not initialized in TransportFactory
void MixTransport::update_C()
{
    // signal that concentration-dependent quantities will need to
    // be recomputed before use, and update the local mole
    // fractions.

    m_visc_ok = false;
    m_condmix_ok = false;

    m_thermo->getMoleFractions(DATA_PTR(m_molefracs));

    // add an offset to avoid a pure species condition
    for (size_t k = 0; k < m_nsp; k++) {
        m_molefracs[k] = std::max(Tiny, m_molefracs[k]);
    }
}

void MixTransport::updateCond_T()
{
    if (m_mode == CK_Mode) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = exp(dot4(m_polytempvec, m_condcoeffs[k]));



	m_omega12[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega12Coeff[k]) );
        m_omega13[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega13Coeff[k]) );
        m_omega14[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega14Coeff[k]) );
        m_omega15[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega15Coeff[k]) );
        m_omega23[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega23Coeff[k]) );
        m_omega24[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega24Coeff[k]) );


        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {

	if ( m_thermo->charge(k) != 0 )
            {

		m_cond[k] = m_temp * m_temp * dot5(m_polytempvec, m_condcoeffs[k]);

        	m_omega12[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega12Coeff[k]) );
        	m_omega13[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega13Coeff[k]) );
        	m_omega14[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega14Coeff[k]) );
        	m_omega15[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega15Coeff[k]) );
        	m_omega23[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega23Coeff[k]) );
        	m_omega24[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega24Coeff[k]) );

	    }

	else
	    {
            	m_cond[k] = m_sqrt_t * dot5(m_polytempvec, m_condcoeffs[k]);

        	m_omega12[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega12Coeff[k]) );
        	m_omega13[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega13Coeff[k]) );
        	m_omega14[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega14Coeff[k]) );
        	m_omega15[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega15Coeff[k]) );
        	m_omega23[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega23Coeff[k]) );
        	m_omega24[k] = m_temp * m_temp * exp( dot5(m_polytempvec,m_omega24Coeff[k]) );

	    }

        }

    }
    m_spcond_ok = true;
    m_condmix_ok = false;
}

}
