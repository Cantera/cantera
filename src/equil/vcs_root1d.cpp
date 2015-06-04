/**
 * @file vcs_root1d.cpp
 *  Code for a one dimensional root finder program.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_defs.h"
#include "cantera/numerics/ctlapack.h"

#include <cstdio>

namespace Cantera
{

#define TOL_CONV 1.0E-5

static void print_funcEval(FILE* fp, double xval, double fval, int its)
{
    fprintf(fp,"\n");
    fprintf(fp,"...............................................................\n");
    fprintf(fp,".................. vcs_root1d Function Evaluation .............\n");
    fprintf(fp,"..................  iteration = %5d ........................\n", its);
    fprintf(fp,"..................  value = %12.5g ......................\n", xval);
    fprintf(fp,"..................  funct = %12.5g ......................\n", fval);
    fprintf(fp,"...............................................................\n");
    fprintf(fp,"\n");
}

int vcsUtil_root1d(double xmin, double xmax, size_t itmax,
                   VCS_FUNC_PTR func, void* fptrPassthrough,
                   double FuncTargVal, int varID,
                   double* xbest, int printLvl)
{
    warn_deprecated("vcsUtil_root1d", "To be removed after Cantera 2.2.");
    static int callNum = 0;
    const char* stre = "vcs_root1d ERROR: ";
    const char* strw = "vcs_root1d WARNING: ";
    bool converged = false;
    int err = 0;

#ifdef DEBUG_MODE
    char fileName[80];
#else
    char* fileName;
#endif
    FILE* fp = 0;

    double x1, x2, xnew, f1, f2, fnew, slope;
    size_t its = 0;
    int posStraddle = 0;
    int retn = VCS_SUCCESS;
    bool foundPosF = false;
    bool foundNegF = false;
    bool foundStraddle = false;
    double xPosF = 0.0;
    double xNegF = 0.0;
    double fnorm;   /* A valid norm for the making the function value
            * dimensionless */
    double c[9], f[3], xn1, xn2, x0 = 0.0, f0 = 0.0, root, theta, xquad;

    callNum++;
    if (DEBUG_MODE_ENABLED && printLvl >= 3) {
        sprintf(fileName, "rootfd_%d.log", callNum);
        fp = fopen(fileName, "w");
        fprintf(fp, " Iter   TP_its  xval   Func_val  |  Reasoning\n");
        fprintf(fp, "-----------------------------------------------------"
                "-------------------------------\n");
    } else if (printLvl >= 3) {
        plogf("WARNING: vcsUtil_root1d: printlvl >= 3, but debug mode not turned on\n");
    }
    if (xmax <= xmin) {
        plogf("%sxmin and xmax are bad: %g %g\n", stre, xmin, xmax);
        return VCS_PUB_BAD;
    }
    x1 = *xbest;
    if (x1 < xmin || x1 > xmax) {
        x1 = (xmin + xmax) / 2.0;
    }
    f1 = func(x1, FuncTargVal, varID, fptrPassthrough, &err);
    if (DEBUG_MODE_ENABLED && printLvl >= 3) {
        print_funcEval(fp, x1, f1, static_cast<int>(its));
        fprintf(fp, "%-5d  %-5d  %-15.5E %-15.5E\n", -2, 0, x1, f1);
    }
    if (f1 == 0.0) {
        *xbest = x1;
        return VCS_SUCCESS;
    } else if (f1 > 0.0) {
        foundPosF = true;
        xPosF = x1;
    } else {
        foundNegF = true;
        xNegF = x1;
    }
    x2 = x1 * 1.1;
    if (x2 > xmax) {
        x2 = x1 - (xmax - xmin) / 100.;
    }
    f2 = func(x2, FuncTargVal, varID, fptrPassthrough, &err);
    if (DEBUG_MODE_ENABLED && printLvl >= 3) {
        print_funcEval(fp, x2, f2, static_cast<int>(its));
        fprintf(fp, "%-5d  %-5d  %-15.5E %-15.5E", -1, 0, x2, f2);
    }

    if (FuncTargVal != 0.0) {
        fnorm = fabs(FuncTargVal) + 1.0E-13;
    } else {
        fnorm = 0.5*(fabs(f1) + fabs(f2)) + fabs(FuncTargVal);
    }

    if (f2 == 0.0) {
        return retn;
    } else if (f2 > 0.0) {
        if (!foundPosF) {
            foundPosF = true;
            xPosF = x2;
        }
    } else {
        if (!foundNegF) {
            foundNegF = true;
            xNegF = x2;
        }
    }
    foundStraddle = foundPosF && foundNegF;
    if (foundStraddle) {
        if (xPosF > xNegF) {
            posStraddle = true;
        } else {
            posStraddle = false;
        }
    }
    int ipiv[3];
    int info;

    do {
        /*
        *    Find an estimate of the next point to try based on
        *    a linear approximation.
        */
        slope = (f2 - f1) / (x2 - x1);
        if (slope == 0.0) {
            plogf("%s functions evals produced the same result, %g, at %g and %g\n",
                  strw, f2, x1, x2);
            xnew = 2*x2 - x1 + 1.0E-3;
        } else {
            xnew = x2 - f2 / slope;
        }
        if (DEBUG_MODE_ENABLED && printLvl >= 3) {
            fprintf(fp, " | xlin = %-9.4g", xnew);
        }

        /*
        *  Do a quadratic fit -> Note this algorithm seems
        *  to work OK. The quadratic approximation doesn't kick in until
        *  the end of the run, when it becomes reliable.
        */
        if (its > 0) {
            c[0] = 1.;
            c[1] = 1.;
            c[2] = 1.;
            c[3] = x0;
            c[4] = x1;
            c[5] = x2;
            c[6] = pow(x0, 2);
            c[7] = pow(x1, 2);
            c[8] = pow(x2, 2);
            f[0] = f0;
            f[1] = f1;
            f[2] = f2;

            ct_dgetrf(3, 3, c, 3, ipiv, info);
            if (info) {
                goto QUAD_BAIL;
            }
            ct_dgetrs(ctlapack::NoTranspose, 3, 1, c, 3, ipiv, f, 3, info);
            root = f[1]* f[1] - 4.0 * f[0] * f[2];
            if (root >= 0.0) {
                xn1 = (- f[1] + sqrt(root)) / (2.0 * f[2]);
                xn2 = (- f[1] - sqrt(root)) / (2.0 * f[2]);
                if (fabs(xn2 - x2) < fabs(xn1 - x2) && xn2 > 0.0) {
                    xquad = xn2;
                } else {
                    xquad = xn1;
                }
                theta = fabs(xquad - xnew) / fabs(xnew - x2);
                theta = std::min(1.0, theta);
                xnew = theta * xnew + (1.0 - theta) * xquad;
                if (DEBUG_MODE_ENABLED && printLvl >= 3) {
                    if (theta != 1.0) {
                        fprintf(fp, " | xquad = %-9.4g", xnew);
                    }
                }
            } else {
                /*
                *   Pick out situations where the convergence may be
                *   accelerated.
                */
                if ((sign(xnew - x2) == sign(x2 - x1)) &&
                        (sign(x2   - x1) == sign(x1 - x0))) {
                    xnew += xnew - x2;
                    if (DEBUG_MODE_ENABLED && printLvl >= 3) {
                        fprintf(fp, " | xquada = %-9.4g", xnew);
                    }
                }
            }
        }
QUAD_BAIL:
        ;


        /*
        *
        *  Put heuristic bounds on the step jump
        */
        if ((xnew > x1 && xnew < x2) || (xnew < x1 && xnew > x2)) {
            /*
            *
            *   If we are doing a jump in between two points, make sure
            *   the new trial is between 10% and 90% of the distance
            *   between the old points.
            */
            slope = fabs(x2 - x1) / 10.;
            if (fabs(xnew - x1) < slope) {
                xnew = x1 + sign(xnew-x1) * slope;
                if (DEBUG_MODE_ENABLED && printLvl >= 3) {
                    fprintf(fp, " | x10%% = %-9.4g", xnew);
                }
            }
            if (fabs(xnew - x2) < slope) {
                xnew = x2 + sign(xnew-x2) * slope;
                if (DEBUG_MODE_ENABLED && printLvl >= 3) {
                    fprintf(fp, " | x10%% = %-9.4g", xnew);
                }
            }
        } else {
            /*
            *   If we are venturing into new ground, only allow the step jump
            *   to increase by 100% at each iteration
            */
            slope = 2.0 * fabs(x2 - x1);
            if (fabs(slope) < fabs(xnew - x2)) {
                xnew = x2 + sign(xnew-x2) * slope;
                if (DEBUG_MODE_ENABLED && printLvl >= 3) {
                    fprintf(fp, " | xlimitsize = %-9.4g", xnew);
                }
            }
        }


        if (xnew > xmax) {
            xnew = x2 + (xmax - x2) / 2.0;
            if (DEBUG_MODE_ENABLED && printLvl >= 3) {
                fprintf(fp, " | xlimitmax = %-9.4g", xnew);
            }
        }
        if (xnew < xmin) {
            xnew = x2 + (x2 - xmin) / 2.0;
            if (DEBUG_MODE_ENABLED && printLvl >= 3) {
                fprintf(fp, " | xlimitmin = %-9.4g", xnew);
            }
        }
        if (foundStraddle) {
#ifdef DEBUG_MODE
            slope = xnew;
#endif
            if (posStraddle) {
                if (f2 > 0.0) {
                    if (xnew > x2) {
                        xnew = (xNegF + x2)/2;
                    }
                    if (xnew < xNegF) {
                        xnew = (xNegF + x2)/2;
                    }
                } else {
                    if (xnew < x2) {
                        xnew = (xPosF + x2)/2;
                    }
                    if (xnew > xPosF) {
                        xnew = (xPosF + x2)/2;
                    }
                }
            } else {
                if (f2 > 0.0) {
                    if (xnew < x2) {
                        xnew = (xNegF + x2)/2;
                    }
                    if (xnew > xNegF) {
                        xnew = (xNegF + x2)/2;
                    }
                } else {
                    if (xnew > x2) {
                        xnew = (xPosF + x2)/2;
                    }
                    if (xnew < xPosF) {
                        xnew = (xPosF + x2)/2;
                    }
                }
            }
            if (DEBUG_MODE_ENABLED && printLvl >= 3) {
                if (slope != xnew) {
                    fprintf(fp, " | xstraddle = %-9.4g", xnew);
                }
            }
        }

        fnew = func(xnew, FuncTargVal, varID, fptrPassthrough, &err);
        if (DEBUG_MODE_ENABLED && printLvl >= 3) {
            fprintf(fp,"\n");
            print_funcEval(fp, xnew, fnew, static_cast<int>(its));
            fprintf(fp, "%-5d  %-5d  %-15.5E %-15.5E", (int) its, 0, xnew, fnew);
        }

        if (foundStraddle) {
            if (posStraddle) {
                if (fnew > 0.0) {
                    if (xnew < xPosF) {
                        xPosF = xnew;
                    }
                } else {
                    if (xnew > xNegF) {
                        xNegF = xnew;
                    }
                }
            } else {
                if (fnew > 0.0) {
                    if (xnew > xPosF) {
                        xPosF = xnew;
                    }
                } else {
                    if (xnew < xNegF) {
                        xNegF = xnew;
                    }
                }
            }
        }

        if (! foundStraddle) {
            if (fnew > 0.0) {
                if (!foundPosF) {
                    foundPosF = true;
                    xPosF = xnew;
                    foundStraddle = true;
                    posStraddle = (xPosF > xNegF);
                }
            } else {
                if (!foundNegF) {
                    foundNegF = true;
                    xNegF = xnew;
                    foundStraddle = true;
                    posStraddle = (xPosF > xNegF);
                }
            }
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
        x2 = xnew;
        f2 = fnew;
        if (fabs(fnew / fnorm) < 1.0E-5) {
            converged = true;
        }
        its++;
    } while (! converged && its < itmax);
    if (converged) {
        if (printLvl >= 1) {
            plogf("vcs_root1d success: convergence achieved\n");
        }
        if (DEBUG_MODE_ENABLED && printLvl >= 3) {
            fprintf(fp, " | vcs_root1d success in %d its, fnorm = %g\n", (int) its, fnorm);
        }
    } else {
        retn = VCS_FAILED_CONVERGENCE;
        if (printLvl >= 1) {
            plogf("vcs_root1d ERROR: maximum iterations exceeded without convergence\n");
        }
        if (DEBUG_MODE_ENABLED && printLvl >= 3) {
            fprintf(fp, "\nvcs_root1d failure in %lu its\n", its);
        }
    }
    *xbest = x2;
    if (DEBUG_MODE_ENABLED && printLvl >= 3) {
        fclose(fp);
    }
    return retn;
}

}
