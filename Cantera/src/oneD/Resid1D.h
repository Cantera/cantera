/**
 *  @file Resid1D.h
 *
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#ifndef CT_RESID1D_H
#define CT_RESID1D_H


//#include "stringUtils.h"
#include "../ctexceptions.h"
#include "../xml.h"

namespace Cantera {

    // domain types
    const int cFlowType         = 101;
    const int cSurfType         = 102;
    const int cConnectorType    = 103;
    const int cInletType        = 104;
    const int cSymmType         = 105;
    const int cOutletType       = 106;

    class MultiJac;
    class OneDim;


    /**
     * Base class for single-domain, one-dimensional residual function
     * evaluators.
     */
    class Resid1D {
    public:

        /**
         * Constructor.
         * @param nv Number of variables at each grid point.
         * @param points Number of grid points.
         */
        Resid1D(int nv=1, int points=1, 
            doublereal time = 0.0) : 
            m_time(time),
            m_container(0), 
            m_index(-1),
            m_type(0),
            m_iloc(0),
            m_jstart(0),
            m_left(0), 
            m_right(0) {
            resize(nv, points);
        }

        /// Destructor.
        virtual ~Resid1D(){}

        /// Domain type flag.
        const int domainType() { return m_type; }

        const OneDim& container() const{ return *m_container; }

        /**
         * Specify the container object for this domain, and the
         * position of this domain in the list.
         */
        void setContainer(OneDim* c, int index){
            m_container = c;
            m_index = index;
        }

        /** Initialize. Base class method does nothing, but may be
         * overloaded.
         */
        virtual void init(){}

        /**
         * Resize the domain to have nv components and np grid points.
         */
        virtual void resize(int nv, int np) {
            m_nv = nv;
            m_max.resize(m_nv, 0.0);
            m_min.resize(m_nv, 0.0);
            m_rtol.resize(m_nv, 0.0);
            m_atol.resize(m_nv, 0.0);
            m_points = np;
            m_slast.resize(m_nv * m_points, 0.0);
            locate();
        }

        /// Number of components at each grid point.
        int nComponents() const { return m_nv; }

        /// Number of grid points in this domain.
        int nPoints() const { return m_points; }

        /// Name of the nth component. May be overloaded.
        virtual string componentName(int n) const { 
            return "component " + int2str(n); }

        /**
         * Set the lower and upper bounds for each solution component.
         */
        void setBounds(int nl, const doublereal* lower, 
            int nu, const doublereal* upper) {
            if (nl != m_nv || nu != m_nv)
                throw CanteraError("Resid1D::setBounds",
                    "wrong array size for solution bounds");
            copy(upper, upper + m_nv, m_max.begin());
            copy(lower, lower + m_nv, m_min.begin());
        }

        void setTolerances(int nr, const doublereal* rtol, 
            int na, const doublereal* atol) {
            if (nr != m_nv || na != m_nv)
                throw CanteraError("Resid1D::setTolerances",
                    "wrong array size for solution error tolerances. Size should be "+int2str(m_nv));
            copy(rtol, rtol + m_nv, m_rtol.begin());
            copy(atol, atol + m_nv, m_atol.begin());
        }

        doublereal rtol(int n) { return m_rtol[n]; }
        doublereal atol(int n) { return m_atol[n]; }

        doublereal upperBound(int n) const { return m_max[n]; }
        doublereal lowerBound(int n) const { return m_min[n]; }

        void initTimeInteg(doublereal dt, const doublereal* x0) {
            copy(x0 + loc(), x0 + loc() + size(), m_slast.begin());
            m_rdt = 1.0/dt;
        } 
        
        void setSteadyMode() { m_rdt = 0.0; }

        bool steady() { return (m_rdt == 0.0); }
        bool transient() { return (m_rdt != 0.0); }

        void needJacUpdate();

        void evalss(doublereal* x, doublereal* r, integer* mask) {
            eval(-1,x,r,mask,0.0);
        }

        /**
         * Evaluate the residual function at point j. If j < 0,
         * evaluate the residual function at all points.
         */         
        virtual void eval(int j, doublereal* x, doublereal* r, 
            integer* mask, doublereal rdt=0.0) {
            throw CanteraError("Resid1D::eval",
                "residual function not defined.");
        }

        virtual void update(doublereal* x) {}

        doublereal time() { return m_time;}
        void incrementTime(doublereal dt) { m_time += dt; }
        size_t index(int n, int j) const { return m_nv*j + n; }

        virtual void setJac(MultiJac* jac){}
        virtual void save(XML_Node& o, doublereal* sol) {
            throw CanteraError("Resid1D::save","base class method called");
        }

        int size() { return m_nv*m_points; }

        void locate() {
            if (m_left) {
                m_jstart = m_left->lastPoint() + 1;
                m_iloc = m_left->loc() + m_left->size();
            }
            else {
                m_jstart = 0;
                m_iloc = 0;
            }
            if (m_right) m_right->locate();
        }

        virtual int loc(int j = 0) { return m_iloc; }

        int firstPoint() { return m_jstart; }
        int lastPoint() { return m_jstart + m_points - 1; }

        void append(Resid1D* right) {
            linkRight(right);
            right->linkLeft(this);
        }

        void linkLeft(Resid1D* left) { 
            m_left = left;
            locate();
        }
        void linkRight(Resid1D* right) { m_right = right; }

        Resid1D* left() { return m_left; }
        Resid1D* right() { return m_right; }

        double prevSoln(int n, int j) const{
            return m_slast[m_nv*j + n];
        }
        
        void setID(const string& s) {m_id = s;}
        void setDesc(const string& s) {m_desc = s;}

        virtual void getTransientMask(integer* mask){}

    protected:

        doublereal m_rdt;
        int m_nv;
        int m_points;
        vector_fp m_slast;
        doublereal m_time;
        vector_fp m_max;
        vector_fp m_min;
        vector_fp m_rtol;
        vector_fp m_atol;
        OneDim* m_container;
        int m_index;
        int m_type;
        int m_iloc;
        int m_jstart;
        Resid1D *m_left, *m_right;
        string m_id, m_desc;

    private:

    };
}

#endif


