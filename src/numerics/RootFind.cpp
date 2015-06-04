/**
 * @file: RootFind.cpp  root finder for 1D problems
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "cantera/numerics/RootFind.h"

// turn on debugging for now
#ifndef DEBUG_MODE
#define DEBUG_MODE
#undef DEBUG_MODE_ENABLED
#define DEBUG_MODE_ENABLED 1
#endif

#include "cantera/base/utilities.h"
#include "cantera/base/stringUtils.h"

using namespace std;
namespace Cantera
{

//!  Print out a form for the current function evaluation
/*!
 *  @param fp     Pointer to the FILE object
 *  @param xval   Current value of x
 *  @param fval   Current value of f
 *  @param its    Current iteration value
 */
static void print_funcEval(FILE* fp, doublereal xval, doublereal fval, int its)
{
    fprintf(fp,"\n");
    fprintf(fp,"...............................................................\n");
    fprintf(fp,".................. RootFind Function Evaluation ...............\n");
    fprintf(fp,"..................  iteration = %5d ........................\n", its);
    fprintf(fp,"..................  value = %12.5g ......................\n", xval);
    fprintf(fp,"..................  funct = %12.5g ......................\n", fval);
    fprintf(fp,"...............................................................\n");
    fprintf(fp,"\n");
}

RootFind::RootFind(ResidEval* resid) :
    m_residFunc(resid),
    m_funcTargetValue(0.0),
    m_atolf(1.0E-11),
    m_atolx(1.0E-11),
    m_rtolf(1.0E-5),
    m_rtolx(1.0E-5),
    m_maxstep(1000),
    printLvl(0),
    writeLogAllowed_(false),
    DeltaXnorm_(0.01),
    specifiedDeltaXnorm_(0),
    DeltaXMax_(1.0E6),
    specifiedDeltaXMax_(0),
    FuncIsGenerallyIncreasing_(false),
    FuncIsGenerallyDecreasing_(false),
    deltaXConverged_(0.0),
    x_maxTried_(-1.0E300),
    fx_maxTried_(0.0),
    x_minTried_(1.0E300),
    fx_minTried_(0.0)
{
}

RootFind::RootFind(const RootFind& r) :
    m_residFunc(r.m_residFunc),
    m_funcTargetValue(0.0),
    m_atolf(1.0E-11),
    m_atolx(1.0E-11),
    m_rtolf(1.0E-5),
    m_rtolx(1.0E-5),
    m_maxstep(1000),
    printLvl(0),
    writeLogAllowed_(false),
    DeltaXnorm_(0.01),
    specifiedDeltaXnorm_(0),
    DeltaXMax_(1.0E6),
    specifiedDeltaXMax_(0),
    FuncIsGenerallyIncreasing_(false),
    FuncIsGenerallyDecreasing_(false),
    deltaXConverged_(0.0),
    x_maxTried_(-1.0E300),
    fx_maxTried_(0.0),
    x_minTried_(1.0E300),
    fx_minTried_(0.0)
{
    *this = r;
}

RootFind& RootFind::operator=(const RootFind& right)
{
    if (this == &right) {
        return *this;
    }
    m_residFunc = right.m_residFunc;
    m_funcTargetValue = right.m_funcTargetValue;
    m_atolf = right.m_atolf;
    m_atolx = right.m_atolx;
    m_rtolf = right.m_rtolf;
    m_rtolx = right.m_rtolx;
    m_maxstep = right.m_maxstep;
    printLvl  = right.printLvl;
    writeLogAllowed_ = right.writeLogAllowed_;
    DeltaXnorm_  = right.DeltaXnorm_;
    specifiedDeltaXnorm_ = right.specifiedDeltaXnorm_;
    DeltaXMax_  = right.DeltaXMax_;
    specifiedDeltaXMax_ = right.specifiedDeltaXMax_;
    FuncIsGenerallyIncreasing_   = right.FuncIsGenerallyIncreasing_;
    FuncIsGenerallyDecreasing_  = right.FuncIsGenerallyDecreasing_;
    deltaXConverged_  = right.deltaXConverged_;
    x_maxTried_   = right.x_maxTried_;
    fx_maxTried_  = right.fx_maxTried_;
    x_minTried_   = right.x_minTried_;
    fx_minTried_  = right.fx_minTried_;

    return *this;
}

doublereal RootFind::delXNonzero(doublereal x1) const
{
    doublereal deltaX = 1.0E-14 * fabs(x1);
    doublereal delmin = DeltaXnorm_ * 1.0E-14;
    if (delmin > deltaX) {
        return delmin;
    }
    return deltaX;
}

doublereal RootFind::delXMeaningful(doublereal x1) const
{
    doublereal del = delXNonzero(x1);
    if (deltaXConverged_ > del) {
        return deltaXConverged_;
    }
    return del;
}

double RootFind::deltaXControlled(doublereal x2, doublereal x1) const
{
    doublereal sgnn = 1.0;
    if (x1 > x2) {
        sgnn = -1.0;
    }
    doublereal deltaX = x2 - x1;
    doublereal x = fabs(x2) + fabs(x1);
    doublereal deltaXm = delXNonzero(x);
    if (fabs(deltaX) < deltaXm) {
        deltaX = sgnn * deltaXm;
    }
    return deltaX;
}

bool RootFind::theSame(doublereal x2, doublereal x1, doublereal factor) const
{
    doublereal x = fabs(x2) + fabs(x1);
    doublereal deltaX = delXMeaningful(x);
    doublereal deltaXSmall = factor * deltaX;
    deltaXSmall = std::max(deltaXSmall , x * 1.0E-15);
    if (fabs(x2 - x1) < deltaXSmall) {
        return true;
    }
    return false;
}

int RootFind::solve(doublereal xmin, doublereal xmax, int itmax, doublereal& funcTargetValue, doublereal* xbest)
{
    /*
     *   We store the function target and then actually calculate a modified functional
     *
     *       func = eval(x1) -  m_funcTargetValue = 0
     */
    m_funcTargetValue = funcTargetValue;

    static int callNum = 0;
    const char* stre = "RootFind ERROR: ";
    const char* strw = "RootFind WARNING: ";
    int converged = 0;
    int bottomBump = 0;
    int topBump = 0;
    FILE* fp = 0;
    int doFinalFuncCall = 0;
    doublereal x1, x2, xnew, f1, f2, fnew, slope;
    doublereal deltaX2 = 0.0, deltaXnew = 0.0;

    int posStraddle = 0;
    int retn = ROOTFIND_FAILEDCONVERGENCE;
    int foundPosF = 0;
    int foundNegF = 0;
    int foundStraddle = 0;
    doublereal xPosF = 0.0;
    doublereal fPosF = 1.0E300;
    doublereal xNegF = 0.0;
    doublereal fNegF = -1.0E300;
    doublereal fnorm;   /* A valid norm for the making the function value  dimensionless */
    doublereal xDelMin;
    doublereal sgn;
    doublereal fnoise = 0.0;
    rfHistory_.clear();
    rfTable rfT;
    rfT.clear();
    rfT.reasoning = "First Point: ";

    callNum++;
    if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
        char fileName[80];
        sprintf(fileName, "RootFind_%d.log", callNum);
        fp = fopen(fileName, "w");
        fprintf(fp, " Iter   TP_its  xval   Func_val  |  Reasoning\n");
        fprintf(fp, "-----------------------------------------------------"
                "-------------------------------\n");
    } else if (printLvl >= 3) {
        writelog("WARNING: RootFind: printlvl >= 3, but debug mode not turned on\n");
    }
    if (xmax <= xmin) {
        writelogf("%sxmin and xmax are bad: %g %g\n", stre, xmin, xmax);
        funcTargetValue = func(*xbest);
        return ROOTFIND_BADINPUT;
    }

    /*
     *  If the maximum step size has not been specified, set it here to 1/5 of the
     *  domain range of x.
     */
    if (!specifiedDeltaXMax_) {
        DeltaXMax_ = 0.2 *(xmax - xmin);
    }

    if (!specifiedDeltaXnorm_) {
        DeltaXnorm_ = 0.2 * DeltaXMax_;
    } else {
        if (DeltaXnorm_ > DeltaXMax_) {
            if (specifiedDeltaXnorm_) {
                DeltaXMax_ = DeltaXnorm_;
            } else {
                DeltaXnorm_ = 0.5 * DeltaXMax_;
            }
        }
    }

    /*
     *  Calculate an initial value of deltaXConverged_
     */
    deltaXConverged_ = m_rtolx * (*xbest) + m_atolx;
    if (DeltaXnorm_ < deltaXConverged_) {
        writelogf("%s DeltaXnorm_, %g, is too small compared to tols, increasing to %g\n",
                  stre, DeltaXnorm_,  deltaXConverged_);
        DeltaXnorm_ =  deltaXConverged_;
    }

    /*
     *  Find the first function value f1 = func(x1), by using the value entered into xbest.
     *  Process it
     */
    x1 = *xbest;
    if (x1 < xmin || x1 > xmax) {
        x1 = (xmin + xmax) / 2.0;
        rfT.reasoning += " x1 set middle between xmin and xmax because entrance is outside bounds.";
    } else {
        rfT.reasoning += " x1 set to entrance x.";
    }

    x_maxTried_ = x1;
    x_minTried_ = x1;
    int its = 1;
    f1 = func(x1);

    if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
        print_funcEval(fp, x1, f1, its);
        fprintf(fp, "%-5d  %-5d  %-15.5E %-15.5E\n", -2, 0, x1, f1);
    }

    if (f1 == 0.0) {
        *xbest = x1;
        return 0;
    } else if (f1 > fnoise) {
        foundPosF = 1;
        xPosF = x1;
        fPosF = f1;
    } else if (f1 < -fnoise) {
        foundNegF = 1;
        xNegF = x1;
        fNegF = f1;
    }
    rfT.its = its;
    rfT.TP_its = 0;
    rfT.xval = x1;
    rfT.fval = f1;
    rfT.foundPos = foundPosF;
    rfT.foundNeg = foundNegF;
    rfT.deltaXConverged = m_rtolx * (fabs(x1) + 0.001);
    rfT.deltaFConverged = fabs(f1) * m_rtolf;
    rfT.delX = xmax - xmin;
    rfHistory_.push_back(rfT);
    rfT.clear();

    /*
     * Now, this is actually a tricky part of the algorithm - Find the x value for
     * the second point. It's tricky because we don't have a valid idea of the scale of x yet
     *
     */
    rfT.reasoning = "Second Point: ";
    if (x1 == 0.0) {
        x2 = x1 +  0.01 * DeltaXnorm_;
        rfT.reasoning += "Set by DeltaXnorm_";
    } else {
        x2 = x1 * 1.0001;
        rfT.reasoning += "Set slightly higher.";
    }
    if (x2 > xmax) {
        x2 = x1 - 0.01 * DeltaXnorm_;
        rfT.reasoning += " - But adjusted to be within bounds";
    }

    /*
     *  Find the second function value f2 = func(x2),  Process it
     */
    deltaX2 = x2 - x1;
    its++;
    f2 = func(x2);
    if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
        print_funcEval(fp, x2, f2, its);
        fprintf(fp, "%-5d  %-5d  %-15.5E %-15.5E", -1, 0, x2, f2);
    }

    /*
     * Calculate the norm of the function, this is the nominal value of f. We try
     * to reduce the nominal value of f by rtolf, this is the main convergence requirement.
     */
    if (m_funcTargetValue != 0.0) {
        fnorm = m_atolf + fabs(m_funcTargetValue);
    } else {
        fnorm = 0.5*(fabs(f1) + fabs(f2)) + fabs(m_funcTargetValue) + m_atolf;
    }
    fnoise = 1.0E-100;


    if (f2 > fnoise) {
        if (!foundPosF) {
            foundPosF = 1;
            xPosF = x2;
            fPosF = f2;
        }
    } else if (f2 < - fnoise) {
        if (!foundNegF) {
            foundNegF = 1;
            xNegF = x2;
            fNegF = f2;
        }
    } else if (f2 == 0.0) {
        *xbest = x2;
        return ROOTFIND_SUCCESS;
    }
    rfT.its = its;
    rfT.TP_its = 0;
    rfT.xval = x2;
    rfT.fval = f2;
    rfT.foundPos = foundPosF;
    rfT.foundNeg = foundNegF;


    /*
     *  See if we have already achieved a straddle
     */
    foundStraddle = foundPosF && foundNegF;
    if (foundStraddle) {
        if (xPosF > xNegF) {
            posStraddle = 1;
        } else {
            posStraddle = 0;
        }
    }

    bool useNextStrat = false;
    bool slopePointingToHigher = true;
    // ---------------------------------------------------------------------------------------------
    //                MAIN LOOP
    // ---------------------------------------------------------------------------------------------
    do {
        /*
         *    Find an estimate of the next point, xnew, to try based on
         *    a linear approximation from the last two points.
         */
        if (DEBUG_MODE_ENABLED && fabs(x2 - x1) < 1.0E-14) {
            printf(" RootFind: we are here x2 = %g x1 = %g\n", x2, x1);
        }
        doublereal delXtmp = deltaXControlled(x2, x1);
        slope = (f2 - f1) / delXtmp;
        rfT.slope = slope;
        rfHistory_.push_back(rfT);
        rfT.clear();
        rfT.reasoning = "";
        if (fabs(slope) <= 1.0E-100) {
            if (printLvl >= 2) {
                writelogf("%s functions evals produced the same result, %g, at %g and %g\n",
                          strw, f2, x1, x2);
            }
            xnew = x2 +  DeltaXnorm_;
            slopePointingToHigher = true;
            useNextStrat = true;
            rfT.reasoning += "Slope is close to zero. ";
        } else {
            useNextStrat = false;
            xnew = x2 - f2 / slope;
            if (xnew > x2) {
                slopePointingToHigher = true;
            } else {
                slopePointingToHigher = false;
            }
            rfT.reasoning += "Slope is good. ";
        }
        if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
            fprintf(fp, " | xlin = %-11.5E", xnew);
        }
        deltaXnew = xnew - x2;
        /*
         * If the suggested step size is too big, throw out step
         */
        if (!foundStraddle) {
            if (fabs(xnew - x2) > DeltaXMax_) {
                useNextStrat = true;
                rfT.reasoning += "Too large change in xnew from slope. ";
            }
            if (fabs(deltaXnew) < fabs(deltaX2)) {
                deltaXnew = 1.2 * deltaXnew;
                xnew = x2 + deltaXnew;
            }
        }
        /*
         * If the slope can't be trusted using a different strategy for picking the next point
         */
        if (useNextStrat) {
            rfT.reasoning += "Using DeltaXnorm, " + fp2str(DeltaXnorm_) + " and FuncIsGenerallyIncreasing hints. ";
            if (f2 < 0.0) {
                if (FuncIsGenerallyIncreasing_) {
                    if (slopePointingToHigher) {
                        xnew = std::min(x2 + 3.0*DeltaXnorm_, xnew);
                    } else {
                        xnew = x2 + DeltaXnorm_;
                    }
                } else if (FuncIsGenerallyDecreasing_) {
                    if (!slopePointingToHigher) {
                        xnew = std::max(x2 - 3.0*DeltaXnorm_, xnew);
                    } else {
                        xnew = x2 - DeltaXnorm_;
                    }
                } else {
                    if (slopePointingToHigher) {
                        xnew = x2 + DeltaXnorm_;
                    } else {
                        xnew = x2 - DeltaXnorm_;
                    }
                }
            } else {
                if (FuncIsGenerallyDecreasing_) {
                    if (!slopePointingToHigher) {
                        xnew = std::max(x2 + 3.0*DeltaXnorm_, xnew);
                    } else {
                        xnew = x2 + DeltaXnorm_;
                    }
                } else if (FuncIsGenerallyIncreasing_) {
                    if (! slopePointingToHigher) {
                        xnew = std::min(x2 - 3.0*DeltaXnorm_, xnew);
                    } else {
                        xnew = x2 - DeltaXnorm_;
                    }
                } else {
                    if (slopePointingToHigher) {
                        xnew = x2 + DeltaXnorm_;
                    } else {
                        xnew = x2 - DeltaXnorm_;
                    }
                }
            }
        }

        /*
         *  Here, if we have a straddle, we purposefully overshoot the smaller side by 5%. Yes it does lead to
         *  more iterations. However, we're interested in bounding x, and not just doing Newton's method.
         */
        if (foundStraddle) {
            double delta = fabs(x2 - x1);
            if (fabs(xnew - x1) < .01 * delta) {
                xnew = x1 + 0.01 * (x2 - x1);
            } else if (fabs(xnew - x2) < .01 * delta) {
                xnew = x1 + 0.01 * (x2 - x1);
            } else if ((xnew > x1 && xnew < x2) || (xnew < x1 && xnew > x2)) {
                if (fabs(xnew - x1) < fabs(x2 - xnew)) {
                    xnew = x1 + 20./19. * (xnew - x1);
                } else {
                    xnew = x2 + 20./19. * (xnew - x2);
                }
            }
        }
        /*
         *  OK, we have an estimate xnew.
         *
         *
         *  Put heuristic bounds on the step jump
         */
        if ((xnew > x1 && xnew < x2) || (xnew < x1 && xnew > x2)) {
            /*
             *   If we are doing a jump in between the two previous points, make sure
             *   the new trial is no closer that 10% of the distances between x2-x1 to
             *   any of the original points. This is an important part of finding a good bound.
             */
            xDelMin = fabs(x2 - x1) / 10.;
            if (fabs(xnew - x1) < xDelMin) {
                xnew = x1 + sign(xnew-x1) * xDelMin;
                if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
                    fprintf(fp, " | x10%% = %-11.5E", xnew);
                }
            }
            if (fabs(xnew - x2) < 0.1 * xDelMin) {
                xnew = x2 + sign(xnew-x2) * 0.1 *  xDelMin;
                if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
                    fprintf(fp, " | x10%% = %-11.5E", xnew);
                }
            }
        } else {
            /*
             *   If we are venturing into new ground, only allow the step jump
             *   to increase by 50% at each iteration, unless the step jump is less than
             *   the user has said that it is ok to take
             */
            doublereal xDelMax = 1.5 * fabs(x2 - x1);
            if (specifiedDeltaXnorm_) {
                if (0.5 * DeltaXnorm_ >  xDelMax) {
                    xDelMax = 0.5 *DeltaXnorm_ ;
                }
            }
            if (fabs(xDelMax) < fabs(xnew - x2)) {
                xnew = x2 + sign(xnew-x2) * xDelMax;
                if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
                    fprintf(fp, " | xlimitsize = %-11.5E", xnew);
                }
            }
            /*
             *   If we are doing a jump outside the two previous points, make sure
             *   the new trial is no closer that 10% of the distances between x2-x1 to
             *   any of the original points. This is an important part of finding a good bound.
             */
            xDelMin = 0.1 * fabs(x2 - x1);
            if (fabs(xnew - x2) < xDelMin) {
                xnew = x2 + sign(xnew - x2) * xDelMin;
                if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
                    fprintf(fp, " | x10%% = %-11.5E", xnew);
                }
            }
            if (fabs(xnew - x1) < xDelMin) {
                xnew = x1 + sign(xnew - x1) * xDelMin;
                if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
                    fprintf(fp, " | x10%% = %-11.5E", xnew);
                }
            }
        }
        /*
         *  HKM -> Not sure this section is needed
         */
        if (foundStraddle) {
            double xorig = xnew;
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
                    if (xnew < x2)  {
                        xnew = (xNegF + x2)/2;
                    }
                    if (xnew > xNegF) {
                        xnew = (xNegF + x2)/2;
                    }
                } else {
                    if (xnew > x2)  {
                        xnew = (xPosF + x2)/2;
                    }
                    if (xnew < xPosF) {
                        xnew = (xPosF + x2)/2;
                    }
                }
            }
            if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
                if (xorig != xnew) {
                    fprintf(fp, " | xstraddle = %-11.5E", xnew);
                }
            }
        }

        /*
         *  Enforce a minimum stepsize if we haven't found a straddle.
         */
        deltaXnew = xnew - x2;
        if (fabs(deltaXnew) < 1.2 * delXMeaningful(xnew)) {
            if (!foundStraddle) {
                sgn = 1.0;
                if (x2 > xnew) {
                    sgn = -1.0;
                }
                deltaXnew = 1.2 * delXMeaningful(xnew) * sgn;
                rfT.reasoning += "Enforcing minimum stepsize from " + fp2str(xnew - x2) +
                                 " to " + fp2str(deltaXnew);
                xnew = x2 + deltaXnew;
            }
        }

        /*
         *  Guard against going above xmax or below xmin
         */
        if (xnew > xmax) {
            topBump++;
            if (topBump < 3) {
                xnew = x2 + (xmax - x2) / 2.0;
                rfT.reasoning += ("xval reduced to " + fp2str(xnew) + " because predicted xnew was above max value of " + fp2str(xmax));
            } else {
                if (x2 == xmax || x1 == xmax) {
                    // we are here when we are bumping against the top limit.
                    // No further action is possible
                    retn = ROOTFIND_SOLNHIGHERTHANXMAX;
                    *xbest = xnew;
                    rfT.slope = slope;
                    rfT.reasoning += "Giving up because we're at xmax and xnew point higher: " + fp2str(xnew);
                    goto done;
                } else {
                    rfT.reasoning += "xval reduced from " + fp2str(xnew) + " to the max value, " + fp2str(xmax);
                    xnew = xmax;
                }
            }
            if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
                fprintf(fp, " | xlimitmax = %-11.5E", xnew);
            }
        }
        if (xnew < xmin) {
            bottomBump++;
            if (bottomBump < 3) {
                rfT.reasoning += ("xnew increased from " + fp2str(xnew) +"  to " + fp2str(x2 - (x2 - xmin) / 2.0) +
                                  " because above min value of " + fp2str(xmin));
                xnew = x2 - (x2 - xmin) / 2.0;
            } else {
                if (x2 == xmin || x1 == xmin) {
                    // we are here when we are bumping against the bottom limit.
                    // No further action is possible
                    retn = ROOTFIND_SOLNLOWERTHANXMIN;
                    *xbest = xnew;
                    rfT.slope = slope;
                    rfT.reasoning = "Giving up because we're already at xmin and xnew points lower: " + fp2str(xnew);
                    goto done;
                } else {
                    rfT.reasoning += "xval increased from " + fp2str(xnew) + " to the min value, " + fp2str(xmin);
                    xnew = xmin;
                }
            }
            if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
                fprintf(fp, " | xlimitmin = %-11.5E", xnew);
            }
        }

        its++;
        fnew = func(xnew);

        if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
            fprintf(fp,"\n");
            print_funcEval(fp, xnew, fnew, its);
            fprintf(fp, "%-5d  %-5d  %-15.5E %-15.5E", its, 0, xnew, fnew);
        }
        rfT.xval = xnew;
        rfT.fval = fnew;
        rfT.its = its;
        if (foundStraddle) {
            if (posStraddle) {
                if (fnew > 0.0) {
                    if (xnew < xPosF) {
                        xPosF = xnew;
                        fPosF = fnew;
                    }
                } else {
                    if (xnew > xNegF) {
                        xNegF = xnew;
                        fNegF = fnew;
                    }
                }
            } else {
                if (fnew > 0.0) {
                    if (xnew > xPosF) {
                        xPosF = xnew;
                        fPosF = fnew;
                    }
                } else {
                    if (xnew < xNegF) {
                        xNegF = xnew;
                        fNegF = fnew;
                    }
                }
            }
        }

        if (! foundStraddle) {
            if (fnew > fnoise) {
                if (!foundPosF) {
                    foundPosF = 1;
                    rfT.foundPos = 1;
                    xPosF = xnew;
                    fPosF = fnew;
                    foundStraddle = 1;
                    if (xPosF > xNegF) {
                        posStraddle = 1;
                    } else {
                        posStraddle = 0;
                    }
                }
            } else if (fnew < - fnoise) {
                if (!foundNegF) {
                    foundNegF = 1;
                    rfT.foundNeg = 1;
                    xNegF = xnew;
                    fNegF = fnew;
                    foundStraddle = 1;
                    if (xPosF > xNegF) {
                        posStraddle = 1;
                    } else {
                        posStraddle = 0;
                    }
                }
            }
        }

        x1 = x2;
        f1 = f2;

        x2 = xnew;
        f2 = fnew;

        /*
         * As we go on to new data points, we make sure that
         * we have the best straddle of the solution with the choice of F1 and F2 when
         * we do have a straddle to work with.
         */
        if (foundStraddle) {
            bool foundBetterPos = false;
            bool foundBetterNeg = false;
            if (posStraddle) {
                if (f2 > 0.0) {
                    if (xPosF < x2) {
                        foundBetterPos = false;
                        x2 = xPosF;
                        f2 = fPosF;
                    }
                    if (f1 > 0.0) {
                        if (foundBetterPos) {
                            x1 = xNegF;
                            f1 = fNegF;
                        } else {
                            if (x1 >= x2) {
                                x1 = xNegF;
                                f1 = fNegF;
                            }
                        }
                    }
                } else {
                    if (xNegF > x2) {
                        foundBetterNeg = false;
                        x2 = xNegF;
                        f2 = fNegF;
                    }
                    if (f1 < 0.0) {
                        if (foundBetterNeg) {
                            x1 = xPosF;
                            f1 = fPosF;
                        } else {
                            if (x1 <= x2) {
                                x1 = xPosF;
                                f1 = fPosF;
                            }
                        }
                    }
                }
            } else {
                if (f2 < 0.0) {
                    if (xNegF < x2) {
                        foundBetterNeg = false;
                        x2 = xNegF;
                        f2 = fNegF;
                    }
                    if (f1 < 0.0) {
                        if (foundBetterNeg) {
                            x1 = xPosF;
                            f1 = fPosF;
                        } else {
                            if (x1 >= x2) {
                                x1 = xPosF;
                                f1 = fPosF;
                            }
                        }
                    }
                } else {
                    if (xPosF > x2) {
                        foundBetterPos = true;
                        x2 = xPosF;
                        f2 = fPosF;
                    }
                    if (f1 > 0.0) {
                        if (foundBetterNeg) {
                            x1 = xNegF;
                            f1 = fNegF;
                        } else {
                            if (x1 <= x2) {
                                x1 = xNegF;
                                f1 = fNegF;
                            }
                        }
                    }
                }
            }
            AssertThrow((f1* f2 <= 0.0), "F1 and F2 aren't bounding");
        }

        deltaX2 = deltaXnew;
        deltaXnew = x2 - x1;
        deltaXConverged_ = 0.5 * deltaXConverged_ + 0.5 * (m_rtolx * 0.5 * (fabs(x2) + fabs(x1)) + m_atolx);
        rfT.deltaXConverged =  deltaXConverged_;
        rfT.deltaFConverged =  fnorm * m_rtolf;
        if (foundStraddle) {
            rfT.delX = std::max(fabs(deltaX2), fabs(deltaXnew));
        } else {
            rfT.delX = std::max(fabs(deltaX2), fabs(deltaXnew));
            if (x2 < x1) {
                rfT.delX = std::max(rfT.delX, x2 - xmin);
            } else {
                rfT.delX = std::max(rfT.delX, xmax - x2);
            }
        }
        /*
         *     Section To Determine CONVERGENCE criteria
         */
        doFinalFuncCall = 0;
        if ((fabs(fnew / fnorm) < m_rtolf) && foundStraddle) {
            if (fabs(deltaX2) < deltaXConverged_ && fabs(deltaXnew) < deltaXConverged_) {
                converged = 1;
                rfT.reasoning += "NormalConvergence";
                retn = ROOTFIND_SUCCESS;
            }

            else if (fabs(slope) > 1.0E-100) {
                double xdels = fabs(fnew / slope);
                if (xdels < deltaXConverged_ * 0.3) {
                    converged = 1;
                    rfT.reasoning += "NormalConvergence-SlopelimitsDelX";
                    doFinalFuncCall = 1;
                    retn = ROOTFIND_SUCCESS;
                }
            }


            /*
             *  Check for excess convergence in the x coordinate
             */
            if (!converged) {
                if (foundStraddle) {
                    doublereal denom = fabs(x1 - x2);
                    if (denom < 1.0E-200) {
                        retn = ROOTFIND_FAILEDCONVERGENCE;
                        converged = true;
                        rfT.reasoning += "ConvergenceFZero but X1X2Identical";
                    }
                    if (theSame(x2, x1, 1.0E-2)) {
                        converged = true;
                        rfT.reasoning += " ConvergenceF and XSame";
                        retn = ROOTFIND_SUCCESS;
                    }
                }
            }
        } else {
            /*
             *  We are here when F is not converged, but we may want to end anyway
             */
            if (!converged) {
                if (foundStraddle) {
                    doublereal denom = fabs(x1 - x2);
                    if (denom < 1.0E-200) {
                        retn = ROOTFIND_FAILEDCONVERGENCE;
                        converged = true;
                        rfT.reasoning += "FNotConverged but X1X2Identical";
                    }
                    /*
                     *  The premise here is that if x1 and x2 get close to one another,
                     *  then the accuracy of the calculation gets destroyed.
                     */
                    if (theSame(x2, x1, 1.0E-5)) {
                        converged = true;
                        retn = ROOTFIND_SUCCESS_XCONVERGENCEONLY;
                        rfT.reasoning += "FNotConverged but XSame";
                    }
                }
            }
        }
    } while (! converged && its < itmax);

done:
    if (converged) {
        rfT.slope = slope;
        rfHistory_.push_back(rfT);
        rfT.clear();
        rfT.its = its;
        AssertThrow((f1* f2 <= 0.0), "F1 and F2 aren't bounding");

        double x_fpos = x2;
        double x_fneg = x1;
        if (f2 < 0.0) {
            x_fpos = x1;
            x_fneg = x2;
        }
        rfT.delX = fabs(x_fpos - x_fneg);
        if (doFinalFuncCall || (fabs(f1) < 2.0 * fabs(f2))) {
            double delXtmp = deltaXControlled(x2, x1);
            slope = (f2 - f1) / delXtmp;
            xnew = x2 - f2 / slope;
            its++;
            fnew = func(xnew);
            if (fnew > 0.0) {
                if (fabs(xnew - x_fneg) < fabs(x_fpos - x_fneg)) {
                    x_fpos = xnew;
                    rfT.delX = fabs(xnew - x_fneg);
                }
            } else {
                if (fabs(xnew - x_fpos) < fabs(x_fpos - x_fneg)) {
                    x_fneg = xnew;
                    rfT.delX = fabs(xnew - x_fpos);
                }
            }
            rfT.its = its;
            if (fabs(fnew) < fabs(f2) && (fabs(fnew) < fabs(f1))) {
                *xbest = xnew;
                if (doFinalFuncCall) {
                    rfT.reasoning += "CONVERGENCE: Another Evaluation Requested";
                    rfT.delX = fabs(xnew - x2);
                } else {
                    rfT.reasoning += "CONVERGENCE: Another Evaluation done because f1 < f2";
                    rfT.delX = fabs(xnew - x1);
                }
                rfT.fval = fnew;
                rfT.xval = xnew;
                x2 = xnew;
                f2 = fnew;
            } else if (fabs(f1) < fabs(f2)) {
                rfT.its = its;
                rfT.xval = xnew;
                rfT.fval = fnew;

                rfT.slope = slope;
                rfT.reasoning += "CONVERGENCE: Another Evaluation not as good as Second Point ";
                rfHistory_.push_back(rfT);
                rfT.clear();
                rfT.its = its;
                std::swap(f1, f2);
                std::swap(x1, x2);
                *xbest = x2;
                if (fabs(fnew) < fabs(f1)) {
                    if (f1 * fnew > 0.0) {
                        std::swap(f1, fnew);
                        std::swap(x1, xnew);
                    }
                }

                rfT.its = its;
                rfT.xval = *xbest;
                rfT.fval = f2;
                rfT.delX = fabs(x_fpos - x_fneg);
                rfT.reasoning += "CONVERGENCE: NormalEnding -> Second point used";
            } else {
                rfT.its = its;
                rfT.xval = xnew;
                rfT.fval = fnew;

                rfT.slope = slope;
                rfT.reasoning += "CONVERGENCE: Another Evaluation not as good as First Point ";
                rfHistory_.push_back(rfT);
                rfT.clear();
                rfT.its = its;
                *xbest = x2;
                rfT.xval = *xbest;
                rfT.fval = f2;
                rfT.delX = fabs(x_fpos - x_fneg);
                rfT.reasoning += "CONVERGENCE: NormalEnding -> Last point used";
            }
        } else {

            *xbest = x2;

            rfT.xval = *xbest;
            rfT.fval = f2;
            rfT.delX = fabs(x2 - x1);
            rfT.reasoning += "CONVERGENCE: NormalEnding -> Last point used";
        }
        funcTargetValue = f2 + m_funcTargetValue;
        rfT.slope = slope;

        if (printLvl >= 1) {
            writelogf("RootFind success: convergence achieved\n");
        }
        if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
            fprintf(fp, " | RootFind success in %d its, fnorm = %g\n", its, fnorm);
        }
        rfHistory_.push_back(rfT);
    } else {
        rfT.reasoning = "FAILED CONVERGENCE ";
        rfT.slope = slope;
        rfT.its = its;
        if (retn == ROOTFIND_SOLNHIGHERTHANXMAX) {
            if (printLvl >= 1) {
                writelogf("RootFind ERROR: Soln probably lies higher than xmax, %g: best guess = %g\n", xmax, *xbest);
            }
            rfT.reasoning += "Soln probably lies higher than xmax, " + fp2str(xmax) + ": best guess = " + fp2str(*xbest);
        } else   if (retn == ROOTFIND_SOLNLOWERTHANXMIN) {
            if (printLvl >= 1) {
                writelogf("RootFind ERROR: Soln probably lies lower than xmin, %g: best guess = %g\n", xmin, *xbest);
            }
            rfT.reasoning += "Soln probably lies lower than xmin, " + fp2str(xmin) + ": best guess = " + fp2str(*xbest);
        } else {
            retn = ROOTFIND_FAILEDCONVERGENCE;
            if (printLvl >= 1) {
                writelogf("RootFind ERROR: maximum iterations exceeded without convergence, cause unknown\n");
            }
            rfT.reasoning += "Maximum iterations exceeded without convergence, cause unknown";
        }
        if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
            fprintf(fp, "\nRootFind failure in %d its\n", its);
        }

        *xbest = x2;
        funcTargetValue = f2 + m_funcTargetValue;
        rfT.xval = *xbest;
        rfT.fval = f2;
        rfHistory_.push_back(rfT);
    }
    if (DEBUG_MODE_ENABLED && printLvl >= 3 && writeLogAllowed_) {
        fclose(fp);
    }

    if (printLvl >= 2) {
        printTable();
    }

    return retn;
}

doublereal RootFind::func(doublereal x)
{
    doublereal r;
    if (DEBUG_MODE_ENABLED) {
        checkFinite(x);
    }
    m_residFunc->evalSS(0.0, &x, &r);
    if (DEBUG_MODE_ENABLED) {
        checkFinite(r);
    }
    doublereal ff = r  - m_funcTargetValue;
    if (x >= x_maxTried_) {
        x_maxTried_ = x;
        fx_maxTried_ = ff;
    }
    if (x <= x_minTried_) {
        x_minTried_ = x;
        fx_minTried_ = ff;
    }
    return ff;
}

void RootFind::setTol(doublereal rtolf, doublereal atolf, doublereal rtolx, doublereal atolx)
{
    m_atolf = atolf;
    m_rtolf = rtolf;
    if (rtolx <= 0.0) {
        m_rtolx = atolf;
    } else {
        m_rtolx = rtolx;
    }
    if (atolx <= 0.0) {
        m_atolx = atolf;
    } else {
        m_atolx = atolx;
    }
}

void RootFind::setPrintLvl(int printlvl)
{
    printLvl = printlvl;
}

void RootFind::setFuncIsGenerallyIncreasing(bool value)
{
    if (value) {
        FuncIsGenerallyDecreasing_ = false;
    }
    FuncIsGenerallyIncreasing_ = value;
}

void RootFind::setFuncIsGenerallyDecreasing(bool value)
{
    if (value) {
        FuncIsGenerallyIncreasing_ = false;
    }
    FuncIsGenerallyDecreasing_ = value;
}

void RootFind::setDeltaX(doublereal deltaXNorm)
{
    DeltaXnorm_ = deltaXNorm;
    specifiedDeltaXnorm_ = 1;
}

void RootFind::setDeltaXMax(doublereal deltaX)
{
    DeltaXMax_ = deltaX;
    specifiedDeltaXMax_ = 1;
}

void RootFind::printTable()
{
    printf("\t----------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("\t  RootFinder Summary table: \n");
    printf("\t         FTarget = %g\n", m_funcTargetValue);
    printf("\t Iter |       xval             delX        deltaXConv    |    slope    | foundP foundN|   F - F_targ  deltaFConv  |   Reasoning\n");
    printf("\t----------------------------------------------------------------------------------------------------------------------------------------\n");
    for (int i = 0; i < (int) rfHistory_.size(); i++) {
        struct rfTable rfT = rfHistory_[i];
        printf("\t  %3d |%- 17.11E %- 13.7E  %- 13.7E |%- 13.5E|   %3d   %3d  | %- 12.5E %- 12.5E | %s \n",
               rfT.its, rfT.xval, rfT.delX, rfT.deltaXConverged, rfT.slope, rfT.foundPos, rfT.foundNeg, rfT.fval,
               rfT.deltaFConverged, (rfT.reasoning).c_str());
    }
    printf("\t----------------------------------------------------------------------------------------------------------------------------------------\n");
}

}
