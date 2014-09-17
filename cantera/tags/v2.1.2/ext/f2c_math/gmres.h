//  -*- C++ -*-

#ifndef GMRES_BLAS_H
#define GMRES_BLAS_H

// ============================================================================
//
//  GMRES nach Saad, Schultz
//     GMRES: a generalized minimal residual algorithm for solving nonsymmetric
//     linear systems
//     SIAM J Sci Stat Comput 7, 856-869 (1986)
//
//                                                 ----------------------------
//                                                 Christian Badura, Mai 1998
//
// ============================================================================


template< class Matrix >
inline int
gmres( int m, int N, const Matrix &A, const doublereal *b, doublereal *x, doublereal eps );


template< class Matrix >
inline int
gmres( int m, int N, const Matrix &A, const doublereal *b, doublereal *x, doublereal eps,
       bool detailed );


// ============================================================================

#include "cblas.h"
#include "../../Cantera/src/ctlapack.h"

template< class Matrix >
inline int
gmres( int m, int n, const Matrix &A, const doublereal *b, doublereal *x, doublereal eps,
    bool detailed ) {
    if ( n<=0 )
        return -1;
    typedef doublereal *doublerealP;
    doublereal *V  = new doublereal[n*(m+1)];
    doublereal *U  = new doublereal[m*(m+1)/2];
    doublereal *r  = new doublereal[n];
    doublereal *y  = new doublereal[m+1];
    doublereal *c  = new doublereal[m];
    doublereal *s  = new doublereal[m];
    doublereal **v = new doublerealP[m+1];
    for ( int i=0; i<=m; ++i ) v[i]=V+i*n;
    int its=-1;
    {
        doublereal beta, h, rd, dd, nrm2b;
        int j, io, uij, u0j;
        nrm2b=dnrm2(n,b,1);
        cout << " norm = " << nrm2b << endl;
        io=0;
        do  { // "aussere Iteration
            ++io;
            //mult(A,x,r);
            A.mult(x,r);
            daxpy(n,-1.,b,1,r,1);
            beta=dnrm2(n,r,1);
            dcopy(n,r,1,v[0],1);
            dscal(n,1./beta,v[0],1);

            y[0]=beta;
            j=0;
            uij=0;
            do { // innere Iteration j=0,...,m-1
                u0j=uij;
                //mult(A,v[j],v[j+1]);
                A.mult(v[j],v[j+1]);

                ct_dgemv(ctlapack::ColMajor, ctlapack::Transpose, n, j+1, 1.0, V, n, v[j+1], 1, 0.0, U+u0j, 1);
                ct_dgemv(ctlapack::ColMajor, ctlapack::NoTranspose, n, j+1, -1.0, V, n, U+u0j, 1, 1.0, v[j+1], 1);

                //dgemv(Transpose,n,j+1,1.,V,n,v[j+1],1,0.,U+u0j,1);
                //dgemv(NoTranspose,n,j+1,-1.,V,n,U+u0j,1,1.,v[j+1],1);

                h=dnrm2(n,v[j+1],1);

                dscal(n,1./h,v[j+1],1);

                for ( int i=0; i<j; ++i ) { // rotiere neue Spalte
                    doublereal tmp = c[i]*U[uij]-s[i]*U[uij+1];
                    U[uij+1]   = s[i]*U[uij]+c[i]*U[uij+1];
                    U[uij]     = tmp;
                    ++uij;
                }
                { // berechne neue Rotation
                    rd     = U[uij];
                    dd     = sqrt(rd*rd+h*h);
                    c[j]   = rd/dd;
                    s[j]   = -h/dd;
                    U[uij] = dd;
                    ++uij;
                }
                { // rotiere rechte Seite y (vorher: y[j+1]=0)
                    y[j+1] = s[j]*y[j];
                    y[j]   = c[j]*y[j];
                }
                  ++j;
                  if ( detailed ) {
                      cout<<"gmres("<<m<<")\t"<<io<<"\t"<<j<<"\t"<<y[j]<< endl;
                  }
            } while ( j<m && fabs(y[j])>=eps*nrm2b );
            { // minimiere bzgl Y
                dtpsv(UpperTriangle,NoTranspose,NotUnitTriangular,j,U,y,1);
                // korrigiere X
                dgemv(NoTranspose,n,j,-1.,V,n,y,1,1.,x,1);
            }
        } while ( fabs(y[j])>=eps*nrm2b );

        // R"uckgabe: Zahl der inneren Iterationen
        its = m*(io-1)+j;
    }

    delete[] V;
    delete[] U;
    delete[] r;
    delete[] y;
    delete[] c;
    delete[] s;
    delete[] v;
    return its;
}


// ============================================================================


template< class Matrix >
inline int
gmres( int m, int n, const Matrix &A, const doublereal *b, doublereal *x, doublereal eps ){
    return gmres(m,n,A,b,x,eps,false);
}

// ============================================================================


#endif // GMRES_BLAS_H
