/**
 *
 *  @file TransportFactory.cpp
 *
 *  Implementation file for class TransportFactory.
 *
 *
 */

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif


// known transport models
#include "MultiTransport.h"
#include "MixTransport.h"
#include "SolidTransport.h"

#include "TransportFactory.h"

#include "polyfit.h"
#include "MMCollisionInt.h"
#include <stdio.h>
#include "xml.h"
#include "TransportParams.h"
#include "global.h"
#include "IdealGasPhase.h"
#include "ctml.h"

#include <iostream>


/**
 * polynomial degree used for fitting collision integrals
 * except in CK mode, where the degree is 6.
 */
#define COLL_INT_POLY_DEGREE 8


namespace Cantera {

    TransportFactory* TransportFactory::__factory = 0;

    
    ////////////////////////// exceptions /////////////////////////

    /** 
     * Exception thrown if an error is encountered while reading the 
     * transport database.
     */
    class TransportDBError : public CanteraError {
    public:
        TransportDBError(int linenum, string msg) 
            : CanteraError("getTransportData",
                "error reading transport data: " 
                + msg + "\n") {}
    };



    /////////////////////////// constants //////////////////////////

    const doublereal ThreeSixteenths = 3.0/16.0;
    const doublereal TwoOverPi       = 2.0/Pi;
    const doublereal FiveThirds      = 5.0/3.0;



    //////////////////// class TransportFactory methods //////////////


    /**
     * Calculate second-order corrections to binary diffusion
     * coefficient pair (dkj, djk). At first order, the binary
     * diffusion coefficients are independent of composition, and
     * d(k,j) = d(j,k). But at second order, there is a weak
     * dependence on composition, with the result that d(k,j) !=
     * d(j,k). This method computes the multiplier by which the
     * first-order binary diffusion coefficient should be multiplied
     * to produce the value correct to second order. The expressions
     * here are taken from Marerro and Mason,
     * J. Phys. Chem. Ref. Data, vol. 1, p. 3 (1972). 
     *
     * @param t   Temperature (K)
     * @param tr  Transport parameters
     * @param k   index of first species
     * @param j   index of second species
     * @param xmk mole fraction of species k
     * @param xmj mole fraction of species j
     * @param fkj multiplier for d(k,j)
     * @param fjk multiplier for d(j,k) 
     */
    void TransportFactory::getBinDiffCorrection(doublereal t, 
        const TransportParams& tr, int k, int j, doublereal xk, doublereal xj, 
        doublereal& fkj, doublereal& fjk) {

        doublereal w1, w2, wsum, sig1, sig2, sig12, sigratio, sigratio2,
            sigratio3, tstar1, tstar2, tstar12,
            om22_1, om22_2, om22_12, om11_12, astar_12, bstar_12, cstar_12,
            cnst, wmwp, sqw12, p1, p2, p12, q1, q2, q12;

        w1 = tr.mw[k];
        w2 = tr.mw[j];
        wsum = w1 + w2;
        wmwp = (w1 - w2)/wsum;
        sqw12 = sqrt(w1*w2);
        
        sig1 = tr.sigma[k];
        sig2 = tr.sigma[j];
        sig12 = 0.5*(tr.sigma[k] + tr.sigma[j]);
        sigratio = sig1*sig1/(sig2*sig2);
        sigratio2 = sig1*sig1/(sig12*sig12);
        sigratio3 = sig2*sig2/(sig12*sig12);

        tstar1 = Boltzmann * t / tr.eps[k];
        tstar2 = Boltzmann * t / tr.eps[j];
        tstar12 = Boltzmann * t / sqrt(tr.eps[k] * tr.eps[j]);

        om22_1 = m_integrals->omega22(tstar1, tr.delta(k,k));
        om22_2 = m_integrals->omega22(tstar2, tr.delta(j,j));
        om22_12 = m_integrals->omega22(tstar12, tr.delta(k,j));
        om11_12 = m_integrals->omega11(tstar12, tr.delta(k,j));
        astar_12 = m_integrals->astar(tstar12, tr.delta(k,j));
        bstar_12 = m_integrals->bstar(tstar12, tr.delta(k,j));
        cstar_12 = m_integrals->cstar(tstar12, tr.delta(k,j));

        cnst = sigratio * sqrt(2.0*w2/wsum) * 2.0 * 
                 w1*w1/(wsum * w2);
        p1 = cnst * om22_1 / om11_12;

        cnst = (1.0/sigratio) * sqrt(2.0*w1/wsum) * 2.0*w2*w2/(wsum*w1);
        p2 = cnst * om22_2 / om11_12;

        p12 = 15.0 * wmwp*wmwp + 8.0*w1*w2*astar_12/(wsum*wsum);

        cnst = (2.0/(w2*wsum))*sqrt(2.0*w2/wsum)*sigratio2;
        q1 = cnst*((2.5 - 1.2*bstar_12)*w1*w1 + 3.0*w2*w2
            + 1.6*w1*w2*astar_12);

        cnst = (2.0/(w1*wsum))*sqrt(2.0*w1/wsum)*sigratio3;
        q2 = cnst*((2.5 - 1.2*bstar_12)*w2*w2 + 3.0*w1*w1
            + 1.6*w1*w2*astar_12);

        q12 = wmwp*wmwp*15.0*(2.5 - 1.2*bstar_12)
              + 4.0*w1*w2*astar_12*(11.0 - 2.4*bstar_12)/(wsum*wsum)
              +  1.6*wsum*om22_1*om22_2/(om11_12*om11_12*sqw12) 
              * sigratio2 * sigratio3;
            
        cnst = 6.0*cstar_12 - 5.0;
        fkj = 1.0 + 0.1*cnst*cnst *
                    (p1*xk*xk + p2*xj*xj + p12*xk*xj)/
                    (q1*xk*xk + q2*xj*xj + q12*xk*xj);
        fjk = 1.0 + 0.1*cnst*cnst *
                    (p2*xk*xk + p1*xj*xj + p12*xk*xj)/
                    (q2*xk*xk + q1*xj*xj + q12*xk*xj);
    }


    /** 
     * Calculate corrections to the well depth parameter and the
     * diamter for use in computing the binary diffusion coefficient
     * of polar-nonpolar pairs. For more information about this
     * correction, see Dixon-Lewis, Proc. Royal Society (1968).
     */
    void TransportFactory::makePolarCorrections(int i, int j, 
        const TransportParams& tr, doublereal& f_eps, doublereal& f_sigma) {

        if (tr.polar[i] == tr.polar[j]) {
            f_eps = 1.0; f_sigma = 1.0; return;
        }

        // corrections to the effective diameter and well depth 
        // if one is polar and one is non-polar

        int kp = (tr.polar[i] ? i : j);     // the polar one
        int knp = (i == kp ? j : i);        // the nonpolar one

        doublereal d3np, d3p, alpha_star, mu_p_star, xi;
        d3np = pow(tr.sigma[knp],3);
        d3p  = pow(tr.sigma[kp],3);
        alpha_star = tr.alpha[knp]/d3np;
        mu_p_star  = tr.dipole(kp,kp)/sqrt(d3p 
            * tr.eps[kp]);
        xi = 1.0 + 0.25 * alpha_star * mu_p_star * mu_p_star * 
             sqrt(tr.eps[kp]/tr.eps[knp]);
        f_sigma = pow(xi, -1.0/6.0);
        f_eps = xi*xi;
    }


    // Constructor does nothing
    TransportFactory::TransportFactory() : m_integrals(0) {
        m_models["Mix"] = cMixtureAveraged;
        m_models["Multi"] = cMulticomponent;
        m_models["Solid"] = cSolidTransport;
        m_models["None"] = 0;
    }

    TransportFactory::~TransportFactory() {
        delete __factory;
        __factory = 0;
        delete m_integrals;
        m_integrals = 0;
    }

    /**
     *  make one of several transport models, and return a base class
     *  pointer to it.
     */
    Transport* TransportFactory::newTransport(string transportModel,
        thermo_t* phase, int log_level) {

        if (transportModel == "") return new Transport;

        vector_fp state;
        Transport* tr = 0;
        phase->saveState(state);

        switch(m_models[transportModel]) {
        case None:
            tr = new Transport; break;
        case cMulticomponent:
            tr = new MultiTransport;
            initTransport(tr, phase, 0, log_level);
            break;
        case CK_Multicomponent:
            tr = new MultiTransport;
            initTransport(tr, phase, CK_Mode, log_level);
            break;
        case cMixtureAveraged:
            tr = new MixTransport;
            initTransport(tr, phase, 0, log_level);
            break;
        case CK_MixtureAveraged:
            tr = new MixTransport;
            initTransport(tr, phase, CK_Mode, log_level);
            break;
        case cSolidTransport:
            tr = new SolidTransport;
            tr->setThermo(*phase);
            break;
        default:
            throw CanteraError("newTransport","unknown transport model");
        }
        phase->restoreState(state);
        return tr;
    }

//     Transport* TransportFactory::
//     newPowerTransport(const string& transport_database, 
//         Transport::phase_t* mix) {
//         doublereal tref, viscref, lambdaref, viscexp, lambdaexp,
//             diffexp;
//         ifstream f(transport_database.c_str());
//         if (!f) {
//             throw TransportDBError("newPowerTransport: file" + 
//                 transport_database + " not found.");
//         }
//         f >> tref 
//           >> viscref >> viscexp
//           >> lambdaref >> lambdaexp
//           >> diffexp;
//         int nsp = mix->nSpecies();
//         int k;
//         vector_fp dref(nsp);
//         for (k = 0; k < nsp; k++) {
//             f >> dref[k];
//         }
//         PowerTransport* t = new PowerTransport();
//         t->init(mix);
//         t->setExponents(viscexp, lambdaexp, diffexp);
//         t->setRefValues(tref, viscref, lambdaref, nsp,  dref.begin());
//         f.close();
//         return t;
//     }



    /** 
     * Prepare to build a new kinetic-theory-based transport manager
     * for low-density gases. Uses polynomial fits to Monchick & Mason
     * collision integrals.
     */
    void TransportFactory::setupMM(ostream& flog, 
        XML_Node* transport_database, 
        thermo_t* thermo, int mode, int log_level, TransportParams& tr) {
        
        // constant mixture attributes
        //tr.mix = thermo;
        tr.thermo = thermo;
        tr.nsp = tr.thermo->nSpecies();
        int nsp = tr.nsp;

        tr.tmin = thermo->minTemp();
        tr.tmax = thermo->maxTemp();
        tr.mw.resize(nsp);
        tr.log_level = log_level;

        copy(tr.thermo->molecularWeights().begin(), 
            tr.thermo->molecularWeights().end(), tr.mw.begin());

        tr.mode = mode;
        tr.epsilon.resize(nsp, nsp, 0.0);
        tr.delta.resize(nsp, nsp, 0.0);
        tr.reducedMass.resize(nsp, nsp, 0.0);
        tr.dipole.resize(nsp, nsp, 0.0);
        tr.diam.resize(nsp, nsp, 0.0);
        tr.crot.resize(nsp);
        tr.zrot.resize(nsp);
        tr.polar.resize(nsp, false);
        tr.alpha.resize(nsp, 0.0);
        tr.poly.resize(nsp);
        tr.sigma.resize(nsp);
        tr.eps.resize(nsp);


        //readTransportDatabase(flog, transport_database, 
        //     tr.thermo->speciesNames(), tr);

        XML_Node root, log;
        //string infile = findInputFile(transport_database);
        //ifstream fin(infile.c_str());
        //root.build(fin);
        getTransportData(transport_database, log,  
            tr.thermo->speciesNames(), tr);

        int i, j;
        for (i = 0; i < nsp; i++) tr.poly[i].resize(nsp);

        doublereal ts1, ts2, tstar_min = 1.e8, tstar_max = 0.0;
        doublereal f_eps, f_sigma;

        DenseMatrix& diam = tr.diam;
        DenseMatrix& epsilon = tr.epsilon;

        for (i = 0; i < nsp; i++) 
        {
            for (j = i; j < nsp; j++) 
            {
                // the reduced mass
                tr.reducedMass(i,j) = 
                    tr.mw[i] * tr.mw[j] / (Avogadro * (tr.mw[i] + tr.mw[j]));

                // hard-sphere diameter for (i,j) collisions
                diam(i,j) = 0.5*(tr.sigma[i] + tr.sigma[j]);

                // the effective well depth for (i,j) collisions
                epsilon(i,j) = sqrt(tr.eps[i]*tr.eps[j]);

                //  The polynomial fits of collision integrals vs. T* 
                //  will be done for the T* from tstar_min to tstar_max
                ts1 = Boltzmann * tr.tmin/epsilon(i,j);
                ts2 = Boltzmann * tr.tmax/epsilon(i,j);
                if (ts1 < tstar_min) tstar_min = ts1;
                if (ts2 > tstar_max) tstar_max = ts2;

                // the effective dipole moment for (i,j) collisions
                tr.dipole(i,j) = sqrt(tr.dipole(i,i)*tr.dipole(j,j));

                // reduced dipole moment delta* (nondimensional)
                doublereal d = diam(i,j);
                tr.delta(i,j) =  0.5 * tr.dipole(i,j)*tr.dipole(i,j) 
                                / (epsilon(i,j) * d * d * d);

                makePolarCorrections(i, j, tr, f_eps, f_sigma);
                tr.diam(i,j) *= f_sigma;
                epsilon(i,j) *= f_eps;

                // properties are symmetric 
                tr.reducedMass(j,i) = tr.reducedMass(i,j);
                diam(j,i) = diam(i,j);
                epsilon(j,i) = epsilon(i,j);
                tr.dipole(j,i)  = tr.dipole(i,j);
                tr.delta(j,i)   = tr.delta(i,j);
            }
        }


        // Chemkin fits the entire T* range in the Monchick and Mason tables,
        // so modify tstar_min and tstar_max if in Chemkin compatibility mode
        //if (mode == CK_Mode) {
            tstar_min = 0.101; 
            tstar_max = 99.9;
            //}

        // initialize the collision integral calculator for the desired
        // T* range
        tr.xml->XML_open(flog, "collision_integrals");
        m_integrals = new MMCollisionInt;
        m_integrals->init(tr.xml, tstar_min, tstar_max);
        fitCollisionIntegrals(flog, tr);
        tr.xml->XML_close(flog, "collision_integrals");

        // make polynomial fits
        tr.xml->XML_open(flog, "property fits");
        fitProperties(tr,flog);
        tr.xml->XML_close(flog, "property fits");
    }


    void TransportFactory::initTransport(Transport* tran, 
        thermo_t* thermo, int mode, int log_level) { 

        XML_Node* transport_database = thermo->speciesData();
        
        TransportParams tr;
        ofstream flog("transport_log.xml");
        tr.xml = new XML_Writer(flog);
        tr.xml->XML_open(flog, "transport");

        // set up Monchick and Mason collision integrals
        setupMM(flog, transport_database, thermo, mode, log_level, tr);

        // do model-specific initialization
        tran->init(tr);

        tr.xml->XML_close(flog, "transport");
        flog.close();

        // finished with log file
        flog.close();

        return;
    }





    /********************************************************
     *
     *      Collision Integral Fits
     *
     ********************************************************/


    void TransportFactory::fitCollisionIntegrals(ostream& logfile, 
                    TransportParams& tr) {

        doublereal* dptr;
        doublereal dstar;
        int nsp = tr.nsp;
        int mode = tr.mode;
        int i, j;

        // Chemkin fits to sixth order polynomials
        int degree = (mode == CK_Mode ? 6 : COLL_INT_POLY_DEGREE);

        tr.xml->XML_open(logfile, "tstar_fits");
        tr.xml->XML_comment(logfile, "fits to A*, B*, and C* vs. log(T*).\n"
            "These are done only for the required dstar(j,k) values."); 
        if (tr.log_level < 3) 
            tr.xml->XML_comment(logfile, "*** polynomial coefficients not printed (log_level < 3) ***");

        for (i = 0; i < nsp; i++) 
        {
            for (j = i; j < nsp; j++) 
            {
                // Chemkin fits only delta* = 0
                if (mode != CK_Mode) 
                    dstar = tr.delta(i,j);
                else 
                    dstar = 0.0;

                // if a fit has already been generated for 
                // delta* = tr.delta(i,j), then use it. Otherwise,
                // make a new fit, and add tr.delta(i,j) to the list
                // of delta* values for which fits have been done.
 
                // 'find' returns a pointer to end() if not found
                if (dptr = find(tr.fitlist.begin(), tr.fitlist.end(), 
                        dstar), dptr == tr.fitlist.end()) 
                { 
                    vector_fp ca(degree+1), cb(degree+1), cc(degree+1);
                    vector_fp co22(degree+1);
                    m_integrals->fit(logfile, degree, dstar,  
                        ca.begin(), cb.begin(), cc.begin());
                    m_integrals->fit_omega22(degree, dstar,
                        co22.begin());
                    tr.omega22_poly.push_back(co22);
                    tr.astar_poly.push_back(ca);
                    tr.bstar_poly.push_back(cb);
                    tr.cstar_poly.push_back(cc);
                    tr.poly[i][j] = tr.astar_poly.size() - 1;
                    tr.fitlist.push_back(dstar);
                }

                // delta* found in fitlist, so just point to this
                // polynomial
                else {
                    tr.poly[i][j] = (dptr - tr.fitlist.begin());
                }
                tr.poly[j][i] = tr.poly[i][j];
            } 
        }
        tr.xml->XML_close(logfile, "tstar_fits");
    }




    /*********************************************************
     *
     *                Read Transport Database
     *
     *********************************************************/

 
    /** 
     * Read transport property data from a file for a list of species.
     * Given the name of a file containing transport property
     * parameters and a list of species names, this method returns an
     * instance of TransportParams containing the transport data for
     * these species read from the file.
     */
//     void TransportFactory::readTransportDatabase( 
//         ostream& logfile, 
//         const string& database_file,
//         const vector<string>& names, 
//         TransportParams& tr)
//     {

//         string dbase_file = findInputFile(database_file);
//         ifstream dbase(dbase_file.c_str());        

//         if (!dbase) throw CanteraError("readTransportDatabase",
//             "cannot open file " + database_file);

//         string name, rest;
//         int geom, linenum = 0;
//         map<string, TransportData> datatable;
//         doublereal welldepth, diam, dipole, polar, rot;


        
//         // read all entries in database into 'datatable' and check for 
//         // errors. Note that this procedure validates all entries, not 
//         // only those for the species listed in 'names'.

//         while (!dbase.eof()) {
//             dbase >> name;
//             linenum++;
//             if (name[0] != '!' && !dbase.eof()) {
//                 dbase >> geom >> welldepth >> diam 
//                       >> dipole >> polar >> rot;

//                 // get the rest of the line, in case there are comments
//                 getline(dbase, rest);

//                 TransportData data;
//                 data.speciesName = name;
//                 data.geometry = geom;
//                 data.wellDepth = welldepth;
//                 data.diameter = diam;
//                 data.dipoleMoment = dipole;
//                 data.polarizability = polar;
//                 data.rotRelaxNumber = rot;

//                 datatable[name] = data;

//                 if (welldepth >= 0.0) data.wellDepth = welldepth;
//                 else throw TransportDBError(linenum,
//                     "negative well depth");

//                 if (diam > 0.0) data.diameter = diam;
//                 else throw TransportDBError(linenum,
//                     "negative or zero diameter");

//                 if (dipole >= 0.0) data.dipoleMoment = dipole;
//                 else throw TransportDBError(linenum,
//                     "negative dipole moment");

//                 if (polar >= 0.0) data.polarizability = polar;
//                 else throw TransportDBError(linenum,
//                     "negative polarizability");

//                 if (rot >= 0.0) data.rotRelaxNumber = rot;
//                 else throw TransportDBError(linenum,
//                     "negative rotation relaxation number");

//                 datatable[name] = data;
//             }
//         }
//         dbase.close();

//         // look up the entries for the species listed in 'names'.
//         tr.xml->XML_open(logfile,"database_parameters");
//         tr.xml->XML_item(logfile,"file_name",database_file);
//         for (int i = 0; i < tr.nsp; i++) {

//             TransportData& trdat = datatable[names[i]];
            
//             // 'datatable' returns a default TransportData object if
//             // the species name is not one in the transport database.
//             // This can be detected by examining 'geometry'.
//             if (trdat.geometry < 0) {
//                 throw TransportDBError(0,"no transport data found for species " 
//                     + names[i]);
//             }

//             // parameters are converted to SI units before storing            

//             // rotational heat capacity / R
//             switch (trdat.geometry) {
//             case 0:
//                 tr.crot[i] = 0.0;     // monatomic
//                 break;
//             case 1:
//                 tr.crot[i] = 1.0;     // linear
//                 break;
//             default:
//                 tr.crot[i] = 1.5;     // nonlinear
//             }


//             tr.dipole(i,i) = 1.e-25 * SqrtTen * trdat.dipoleMoment;

//             if (trdat.dipoleMoment > 0.0) 
//                 tr.polar[i] = true;
//             else
//                 tr.polar[i] = false;

//             // A^3 -> m^3
//             tr.alpha[i] = 1.e-30 * trdat.polarizability;

//             tr.sigma[i] = 1.e-10 * trdat.diameter;

//             tr.eps[i] = Boltzmann * trdat.wellDepth;
//             tr.zrot[i]  = fmaxx(1.0, trdat.rotRelaxNumber);

//             // write database parameters to log file
//             tr.xml->XML_open(logfile, names[i]);
//             tr.xml->XML_item(logfile, "geom", trdat.geometry);
//             tr.xml->XML_item(logfile, "epsilon", trdat.wellDepth);
//             tr.xml->XML_item(logfile, "sigma", trdat.diameter);
//             tr.xml->XML_item(logfile, "dipole", trdat.dipoleMoment);
//             tr.xml->XML_item(logfile, "alpha", trdat.polarizability);
//             tr.xml->XML_item(logfile, "zrot", trdat.rotRelaxNumber);
//             tr.xml->XML_close(logfile, names[i]);
//         }
//         tr.xml->XML_close(logfile,"database_parameters");
    //    }



    /*********************************************************
     *
     *                Read Transport Database
     *
     *********************************************************/

 
    /** 
     * Read transport property data from a file for a list of species.
     * Given the name of a file containing transport property
     * parameters and a list of species names, this method returns an
     * instance of TransportParams containing the transport data for
     * these species read from the file.
     */
    void TransportFactory::getTransportData(XML_Node* transport_database,  
        XML_Node& log, const vector<string>& names, TransportParams& tr)
    {
        string name;
        int geom;
        map<string, GasTransportData> datatable;
        doublereal welldepth, diam, dipole, polar, rot;

        //XML_Node* sparray = find_XML("", &root, "", "", "speciesData");
        vector<XML_Node*> xspecies;
        transport_database->getChildren("species",xspecies);
        int nsp = xspecies.size();
        
        // read all entries in database into 'datatable' and check for 
        // errors. Note that this procedure validates all entries, not 
        // only those for the species listed in 'names'.

        string val, type;
        map<string, int> gindx;
        gindx["atom"] = 100;
        gindx["linear"] = 101;
        gindx["nonlinear"] = 102;
        int linenum = 0;
        int i;
        for (i = 0; i < nsp; i++) {
            XML_Node& sp = *xspecies[i];
            name = sp["name"];
            XML_Node& tr = sp.child("transport");
            getString(tr, "geometry", val, type);
            geom = gindx[val] - 100;
            map<string, doublereal> fv;

            welldepth = getFloat(tr, "LJ_welldepth");
            diam = getFloat(tr, "LJ_diameter");
            dipole = getFloat(tr, "dipoleMoment");
            polar = getFloat(tr, "polarizability");
            rot = getFloat(tr, "rotRelax");

            //getFloats(tr, fv, false);
            //welldepth = fv["LJ_welldepth"];
            //diam = fv["LJ_diameter"];
            //dipole = fv["dipoleMoment"];
            //polar = fv["polarizability"];
            //rot = fv["rotRelax"];

            GasTransportData data;
            data.speciesName = name;
            data.geometry = geom;
            if (welldepth >= 0.0) data.wellDepth = welldepth;
            else throw TransportDBError(linenum,
                "negative well depth");

            if (diam > 0.0) data.diameter = diam;
            else throw TransportDBError(linenum,
                "negative or zero diameter");

            if (dipole >= 0.0) data.dipoleMoment = dipole;
            else throw TransportDBError(linenum,
                "negative dipole moment");
            
            if (polar >= 0.0) data.polarizability = polar;
            else throw TransportDBError(linenum,
                "negative polarizability");

            if (rot >= 0.0) data.rotRelaxNumber = rot;
            else throw TransportDBError(linenum,
                "negative rotation relaxation number");

            datatable[name] = data;
        }

        // look up the entries for the species listed in 'names'.
        //tr.xml->XML_open(logfile,"database_parameters");
        //tr.xml->XML_item(logfile,"file_name",database_file);

        for (i = 0; i < tr.nsp; i++) {

            GasTransportData& trdat = datatable[names[i]];
            
            // 'datatable' returns a default TransportData object if
            // the species name is not one in the transport database.
            // This can be detected by examining 'geometry'.
            if (trdat.geometry < 0) {
                throw TransportDBError(0,"no transport data found for species " 
                    + names[i]);
            }

            // parameters are converted to SI units before storing            

            // rotational heat capacity / R
            switch (trdat.geometry) {
            case 0:
                tr.crot[i] = 0.0;     // monatomic
                break;
            case 1:
                tr.crot[i] = 1.0;     // linear
                break;
            default:
                tr.crot[i] = 1.5;     // nonlinear
            }


            tr.dipole(i,i) = 1.e-25 * SqrtTen * trdat.dipoleMoment;

            if (trdat.dipoleMoment > 0.0) 
                tr.polar[i] = true;
            else
                tr.polar[i] = false;

            // A^3 -> m^3
            tr.alpha[i] = 1.e-30 * trdat.polarizability;

            tr.sigma[i] = 1.e-10 * trdat.diameter;

            tr.eps[i] = Boltzmann * trdat.wellDepth;
            tr.zrot[i]  = fmaxx(1.0, trdat.rotRelaxNumber);

            // write database parameters to log file
            //tr.xml->XML_open(logfile, names[i]);
            //tr.xml->XML_item(logfile, "geom", trdat.geometry);
            //tr.xml->XML_item(logfile, "epsilon", trdat.wellDepth);
            //tr.xml->XML_item(logfile, "sigma", trdat.diameter);
            //tr.xml->XML_item(logfile, "dipole", trdat.dipoleMoment);
            //tr.xml->XML_item(logfile, "alpha", trdat.polarizability);
            //tr.xml->XML_item(logfile, "zrot", trdat.rotRelaxNumber);
            //tr.xml->XML_close(logfile, names[i]);
        }
        //tr.xml->XML_close(logfile,"database_parameters");
    }


    /*********************************************************
     *
     *                Polynomial fitting 
     *
     *********************************************************/



    /*****************   fitProperties  ***************/

    /**
     * Generate polynomial fits for the pure-species viscosities and
     * for the binary diffusion coefficients. If
     * CK_mode, then the fits are of the
     * form \f[
     * \log(\eta(i)) = \sum_{n = 0}^3 a_n(i) (\log T)^n
     * \f]
     * and \f[
     * \log(D(i,j)) = \sum_{n = 0}^3 a_n(i,j) (\log T)^n
     * \f]
     * Otherwise the fits are of the form
     * \f[
     * \eta(i)/sqrt(k_BT) = \sum_{n = 0}^4 a_n(i) (\log T)^n
     * \f]
     * and \f[
     * D(i,j)/sqrt(k_BT)) = \sum_{n = 0}^4 a_n(i,j) (\log T)^n
     * \f]
     */
    void TransportFactory::fitProperties(TransportParams& tr, 
                ostream& logfile) {
        doublereal tstar;
        int k, j, n, ndeg = 0;
        char s[100];

        // number of points to use in generating fit data
        const int np = 50;

        int mode = tr.mode;
        int degree = (mode == CK_Mode ? 3 : 4);

        doublereal t, om22; 
        doublereal dt = (tr.tmax - tr.tmin)/(np-1);
        vector_fp tlog(np), spvisc(np), spcond(np);
        doublereal val, fit;

        vector_fp w(np), w2(np);
        
        // generate array of log(t) values
        for (n = 0; n < np; n++) {
            t = tr.tmin + dt*n;
            tlog[n] = log(t);
        }

        // vector of polynomial coefficients
        vector_fp c(degree + 1), c2(degree + 1);


        // fit the pure-species viscosity and thermal conductivity for
        // each species

        if (tr.log_level < 2)
            tr.xml->XML_comment(logfile, 
                "*** polynomial coefficients not printed (log_level < 3) ***");

        int ipoly;
        doublereal sqrt_T, visc, err, relerr, 
            mxerr = 0.0, mxrelerr = 0.0, mxerr_cond = 0.0, mxrelerr_cond = 0.0;
        tr.xml->XML_open(logfile, "viscosity");
        tr.xml->XML_comment(logfile,"Polynomial fits for viscosity");
        if (mode == CK_Mode) {
            tr.xml->XML_comment(logfile,"log(viscosity) fit to cubic "
                "polynomial in log(T)");
        }
        else {
            sprintf(s, "viscosity/sqrt(T) fit to "
                "polynomial of degree %d in log(T)",degree);
            tr.xml->XML_comment(logfile,s);
        }

        //        const vector_fp& cp_R = tr.mix->cp_R();

        doublereal cp_R, cond, w_RT, f_int, A_factor, B_factor,
            c1, cv_rot, cv_int, f_rot, f_trans, om11;
        doublereal diffcoeff;

        for (k = 0; k < tr.nsp; k++) 
        {
            for (n = 0; n < np; n++) {
                t = tr.tmin + dt*n;

                tr.thermo->setTemperature(t);
                cp_R = ((IdealGasPhase*)tr.thermo)->cp_R()[k];

                tstar = Boltzmann * t/ tr.eps[k];
                sqrt_T = sqrt(t);
                om22 = m_integrals->omega22(tstar, tr.delta(k,k));
                om11 = m_integrals->omega11(tstar, tr.delta(k,k));

                // self-diffusion coefficient, without polar
                // corrections
                diffcoeff = ThreeSixteenths * 
                            sqrt( 2.0 * Pi/tr.reducedMass(k,k) ) *
                            pow((Boltzmann * t), 1.5)/ 
                            (Pi * tr.sigma[k] * tr.sigma[k] * om11);

                // viscosity
                visc = FiveSixteenths 
                       * sqrt(Pi * tr.mw[k] * Boltzmann * t / Avogadro) / 
                       (om22 * Pi * tr.sigma[k]*tr.sigma[k]);

                // thermal conductivity 
                w_RT = tr.mw[k]/(GasConstant * t);
                f_int = w_RT * diffcoeff/visc;
                cv_rot = tr.crot[k];

                A_factor = 2.5 - f_int;
                B_factor = tr.zrot[k] + TwoOverPi 
                           *(FiveThirds * cv_rot + f_int);
                c1 = TwoOverPi * A_factor/B_factor;
                cv_int = cp_R - 2.5 - cv_rot;

                f_rot = f_int * (1.0 + c1);
                f_trans = 2.5 * (1.0 - c1 * cv_rot/1.5);

                cond = (visc/tr.mw[k])*GasConstant*(f_trans * 1.5 
                    + f_rot * cv_rot + f_int * cv_int);

                if (mode == CK_Mode) {
                    spvisc[n] = log(visc);
                    spcond[n] = log(cond);
                    w[n] = -1.0;
                    w2[n] = -1.0;
                }
                else {
                    spvisc[n] = visc/sqrt_T;
                    spcond[n] = cond/sqrt_T;
                    w[n] = 1.0/(spvisc[n]*spvisc[n]);
                    w2[n] = 1.0/(spcond[n]*spcond[n]);
                }
            }
            polyfit(np, tlog.begin(), spvisc.begin(), 
                w.begin(), degree, ndeg, 0.0, c.begin());
            polyfit(np, tlog.begin(), spcond.begin(), 
                w.begin(), degree, ndeg, 0.0, c2.begin());

            // evaluate max fit errors for viscosity
            for (n = 0; n < np; n++) {
                if (mode == CK_Mode) {
                    val = exp(spvisc[n]);
                    fit = exp(poly3(tlog[n], c.begin())); 
                }
                else {
                    sqrt_T = exp(0.5*tlog[n]);
                    val = sqrt_T * spvisc[n];
                    fit = sqrt_T * poly4(tlog[n], c.begin());
                }
                err = fit - val;
                relerr = err/val;
                if (fabs(err) > mxerr) mxerr = fabs(err);
                if (fabs(relerr) > mxrelerr) mxrelerr = fabs(relerr);               
            }

            // evaluate max fit errors for conductivity
            for (n = 0; n < np; n++) {
                if (mode == CK_Mode) {
                    val = exp(spcond[n]);
                    fit = exp(poly3(tlog[n], c2.begin())); 
                }
                else {
                    sqrt_T = exp(0.5*tlog[n]);
                    val = sqrt_T * spcond[n];
                    fit = sqrt_T * poly4(tlog[n], c2.begin());
                }
                err = fit - val;
                relerr = err/val;
                if (fabs(err) > mxerr_cond) mxerr_cond = fabs(err);
                if (fabs(relerr) > mxrelerr_cond) mxrelerr_cond = fabs(relerr);
            }
            tr.visccoeffs.push_back(c);
            tr.condcoeffs.push_back(c2);

            if (tr.log_level >= 2) {
                tr.xml->XML_writeVector(logfile, "    ", tr.thermo->speciesName(k), 
                    c.size(), c.begin());
            }
        }

        sprintf(s, "Maximum viscosity absolute error:  %12.6g", mxerr);
        tr.xml->XML_comment(logfile,s);
        sprintf(s, "Maximum viscosity relative error:  %12.6g", mxrelerr);
        tr.xml->XML_comment(logfile,s);        
        tr.xml->XML_close(logfile, "viscosity");


        tr.xml->XML_open(logfile, "conductivity");
        tr.xml->XML_comment(logfile,"Polynomial fits for conductivity");
        if (mode == CK_Mode) 
            tr.xml->XML_comment(logfile,"log(conductivity) fit to cubic "
                "polynomial in log(T)");
        else {
            sprintf(s, "conductivity/sqrt(T) fit to "
                "polynomial of degree %d in log(T)",degree);
            tr.xml->XML_comment(logfile,s);
        }
        if (tr.log_level >= 2) 
            for (k = 0; k < tr.nsp; k++) {
                tr.xml->XML_writeVector(logfile, "    ", tr.thermo->speciesName(k), 
                    degree+1, tr.condcoeffs[k].begin());
            }            
        sprintf(s, "Maximum conductivity absolute error:  %12.6g", mxerr_cond);
        tr.xml->XML_comment(logfile,s);
        sprintf(s, "Maximum conductivity relative error:  %12.6g", mxrelerr_cond);
        tr.xml->XML_comment(logfile,s);        
        tr.xml->XML_close(logfile, "conductivity");

        // fit the binary diffusion coefficients for each species pair

        tr.xml->XML_open(logfile, "binary_diffusion_coefficients");
        tr.xml->XML_comment(logfile, "binary diffusion coefficients");
        if (mode == CK_Mode) 
            tr.xml->XML_comment(logfile,"log(D) fit to cubic "
                "polynomial in log(T)");
        else {
            sprintf(s, "D/T**(3/2) fit to "
                "polynomial of degree %d in log(T)",degree);
            tr.xml->XML_comment(logfile,s);
        }
        mxerr = 0.0, mxrelerr = 0.0;
        vector_fp diff(np + 1);
        doublereal eps, sigma;
        for (k = 0; k < tr.nsp; k++) 
        {            
            for (j = k; j < tr.nsp; j++) {

                ipoly = tr.poly[k][j];
                for (n = 0; n < np; n++) {

                    t = tr.tmin + dt*n;
                    
                    eps = tr.epsilon(j,k);
                    tstar = Boltzmann * t/eps;
                    sigma = tr.diam(j,k);
                    om11 = m_integrals->omega11(tstar, tr.delta(j,k));

                    diffcoeff = ThreeSixteenths * 
                                sqrt( 2.0 * Pi/tr.reducedMass(k,j) ) *
                                pow((Boltzmann * t), 1.5)/ 
                                (Pi * sigma * sigma * om11);

                    // 2nd order correction
                    doublereal fkj, fjk;
                    getBinDiffCorrection(t, tr, k, j, 1.0, 1.0, fkj, fjk);
                    //diffcoeff *= fkj;

                    if (mode == CK_Mode) {

                        diff[n] = log(diffcoeff);
                        w[n] = -1.0;
                    }
                    else {
                        diff[n] = diffcoeff/pow(t, 1.5);
                        w[n] = 1.0/(diff[n]*diff[n]);
                    }
                }
                polyfit(np, tlog.begin(), diff.begin(), 
                    w.begin(), degree, ndeg, 0.0, c.begin());

                doublereal pre;
                for (n = 0; n < np; n++) {
                    if (mode == CK_Mode) {
                        val = exp(diff[n]);
                        fit = exp(poly3(tlog[n], c.begin())); 
                    }
                    else {
                        t = exp(tlog[n]);
                        pre = pow(t, 1.5);
                        val = pre * diff[n];
                        fit = pre * poly4(tlog[n], c.begin());
                    }
                    err = fit - val;
                    relerr = err/val;
                    if (fabs(err) > mxerr) mxerr = fabs(err);
                    if (fabs(relerr) > mxrelerr) mxrelerr = fabs(relerr);               
                }
                tr.diffcoeffs.push_back(c);
                if (tr.log_level >= 2)
                    tr.xml->XML_writeVector(logfile, "    ", tr.thermo->speciesName(k)
                        + "__"+tr.thermo->speciesName(j), c.size(), c.begin());
            }
        }
        sprintf(s,"Maximum binary diffusion coefficient absolute error:"
            "  %12.6g", mxerr);
        tr.xml->XML_comment(logfile,s);
        sprintf(s, "Maximum binary diffusion coefficient relative error:"
            "%12.6g", mxrelerr);
        tr.xml->XML_comment(logfile,s);
        tr.xml->XML_close(logfile, "binary_diffusion_coefficients");
    }
}
