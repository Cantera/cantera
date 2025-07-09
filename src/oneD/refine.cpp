//! @file refine.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/refine.h"
#include "cantera/oneD/Flow1D.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{
Refiner::Refiner(Domain1D& domain) :
    m_domain(&domain)
{
    m_nv = m_domain->nComponents();
    m_active.resize(m_nv, true);
}

void Refiner::setCriteria(double ratio, double slope, double curve, double prune)
{
    if (ratio < 2.0) {
        throw CanteraError("Refiner::setCriteria",
            "'ratio' must be greater than 2.0 ({} was specified).", ratio);
    } else if (slope < 0.0 || slope > 1.0) {
        throw CanteraError("Refiner::setCriteria",
            "'slope' must be between 0.0 and 1.0 ({} was specified).", slope);
    } else if (curve < 0.0 || curve > 1.0) {
        throw CanteraError("Refiner::setCriteria",
            "'curve' must be between 0.0 and 1.0 ({} was specified).", curve);
    } else if (prune > curve || prune > slope) {
        throw CanteraError("Refiner::setCriteria",
            "'prune' must be less than 'curve' and 'slope' ({} was specified).",
            prune);
    }
    m_ratio = ratio;
    m_slope = slope;
    m_curve = curve;
    m_prune = prune;
}

int Refiner::analyze(size_t n, const double* z, const double* x)
{
    if (n >= m_npmax) {
        throw CanteraError("Refiner::analyze", "max number of grid points reached ({}).", m_npmax);
    }

    if (m_domain->nPoints() <= 1) {
        return 0;
    }

    // check consistency
    if (n != m_domain->nPoints()) {
        throw CanteraError("Refiner::analyze", "number of grid points provided does not match domain size.");
    }

    // Reset the state of the refiner
    m_insertPts.clear();
    m_componentNames.clear();
    m_keep.clear();

    // Keep the first and last grid points
    m_keep[0] = KEEP;
    m_keep[n-1] = KEEP;

    m_nv = m_domain->nComponents();

    // find locations where cell size ratio is too large.
    vector<double> val(n);
    vector<double> slope(n-1);

    vector<double> dz(n-1); // Store the right-looking grid spacings
    for (size_t j = 0; j < n-1; j++) {
        dz[j] = z[j+1] - z[j];
    }

    for (size_t i = 0; i < m_nv; i++) {
        if (m_active[i]) {
            string name = m_domain->componentName(i);
            // get component i at all points
            for (size_t j = 0; j < n; j++) {
                val[j] = value(x, i, j);
            }

            // slope of component i (using forward difference)
            for (size_t j = 0; j < n-1; j++) {
                slope[j] = (val[j+1] - val[j]) / dz[j];
            }

            // find the range of values and slopes of component i over the domain
            double valMin = *min_element(val.begin(), val.end());
            double valMax = *max_element(val.begin(), val.end());
            double slopeMin = *min_element(slope.begin(), slope.end());
            double slopeMax = *max_element(slope.begin(), slope.end());

            // max absolute values of val and slope
            double valMagnitude = std::max(fabs(valMax), fabs(valMin));
            double slopeMagnitude = std::max(fabs(slopeMax), fabs(slopeMin));

            // refine based on component i only if the range of val is greater than a
            // fraction 'min_range' of max |val|. This eliminates components that
            // consist of small fluctuations around a constant value.
            if (valMax - valMin > m_min_range*valMagnitude) {
                // maximum allowable difference in value between adjacent points. Based
                // on the global min and max values of the component over the domain.
                double max_change = m_slope*(valMax - valMin);
                for (size_t j = 0; j < n-1; j++) {
                    double ratio = fabs(val[j+1] - val[j]) / (max_change + m_thresh);
                    if (ratio > 1.0 && dz[j] >= 2 * m_gridmin) {
                        m_insertPts.insert(j);
                        m_componentNames.insert(name);
                    }
                    if (ratio >= m_prune) {
                        m_keep[j] = KEEP;
                        m_keep[j+1] = KEEP;
                    } else if (m_keep[j] == UNSET) {
                        m_keep[j] = REMOVE;
                    }
                }
            }

            // refine based on the slope of component i only if the range of s is
            // greater than a fraction 'min_range' of max|s|. This eliminates
            // components that consist of small fluctuations on a constant slope
            // background.
            if (slopeMax - slopeMin > m_min_range*slopeMagnitude) {
                // maximum allowable difference in slope between adjacent points.
                double max_change = m_curve*(slopeMax - slopeMin);
                for (size_t j = 0; j < n-2; j++) {
                    // Using the solution component absolute tolerance (m_thresh),
                    // an absolute tolerance for the change in slope can be estimated
                    // for an interval dz as m_thresh/dz.
                    double ratio = fabs(slope[j+1] - slope[j]) / (max_change + m_thresh/dz[j]);
                    if (ratio > 1.0 && dz[j] >= 2*m_gridmin && dz[j+1] >= 2*m_gridmin) {
                        m_componentNames.insert(name);
                        m_insertPts.insert(j);
                        m_insertPts.insert(j+1);
                    }
                    if (ratio >= m_prune) {
                        m_keep[j+1] = KEEP;
                    } else if (m_keep[j+1] == UNSET) {
                        m_keep[j+1] = REMOVE;
                    }
                }
            }
        }
    }

    // Refine based on properties of the grid itself
    for (size_t j = 1; j < n-1; j++) {
        // Add a new point if the ratio with left interval is too large.
        // Extra points around the interval set under consideration are kept.
        if (dz[j] > m_ratio*dz[j-1]) {
            m_insertPts.insert(j);
            m_componentNames.insert(fmt::format("point {}", j));
            m_keep[j-1] = KEEP;
            m_keep[j] = KEEP;
            m_keep[j+1] = KEEP;
            m_keep[j+2] = KEEP;
        }

        // Add a point if the ratio with right interval is too large
        if (dz[j-1] > m_ratio*dz[j]) {
            m_insertPts.insert(j-1);
            m_componentNames.insert(fmt::format("point {}", j-1));
            m_keep[j-2] = KEEP;
            m_keep[j-1] = KEEP;
            m_keep[j] = KEEP;
            m_keep[j+1] = KEEP;
        }

        // Keep the point if removing would make the ratio with the left interval too
        // large.
        if (j > 1 && z[j+1]-z[j-1] > m_ratio * dz[j-2]) {
            m_keep[j] = KEEP;
        }

        // Keep the point if removing would make the ratio with the right interval too
        // large.
        if (j < n-2 && z[j+1]-z[j-1] > m_ratio * dz[j+1]) {
            m_keep[j] = KEEP;
        }

        Flow1D* fflame = dynamic_cast<Flow1D*>(m_domain);
        // Keep the point where the temperature is fixed
        if (fflame && fflame->isFree() && z[j] == fflame->m_zfixed) {
            m_keep[j] = KEEP;
        }

        // Keep the point if it is a control point used for two-point flame control
        if (fflame && fflame->twoPointControlEnabled() &&
            (z[j] == fflame->leftControlPointCoordinate() ||
             z[j] == fflame->rightControlPointCoordinate()))
        {
            m_keep[j] = KEEP;
        }
    }

    // Don't allow pruning to remove multiple adjacent grid points
    // in a single pass.
    for (size_t j = 2; j < n-1; j++) {
        if (m_keep[j] == REMOVE && m_keep[j-1] == REMOVE) {
            m_keep[j] = KEEP;
        }
    }

    return int(m_insertPts.size());
}

double Refiner::value(const double* x, size_t n, size_t j)
{
    return x[m_domain->index(n,j)];
}

void Refiner::show()
{
    if (!m_insertPts.empty()) {
        writeline('#', 78);
        writelog(string("Refining grid in ") +
                 m_domain->id()+".\n"
                 +"    New points inserted after grid points ");
        for (const auto& loc : m_insertPts) {
            writelog("{} ", loc);
        }
        writelog("\n");
        writelog("    to resolve ");
        for (const auto& c : m_componentNames) {
            writelog(c + " ");
        }
        writelog("\n");
        writeline('#', 78);
    } else if (m_domain->nPoints() > 1) {
        writelog("no new points needed in "+m_domain->id()+"\n");
    }
}

int Refiner::getNewGrid(int n, const double* z, int nn, double* zn)
{
    warn_deprecated(
        "Refiner::getNewGrid",
        "Deprecated in Cantera 3.1; unused function that will be removed.");

    int nnew = static_cast<int>(m_insertPts.size());
    if (nnew + n > nn) {
        throw CanteraError("Refine::getNewGrid", "array size too small.");
    }

    if (m_insertPts.empty()) {
        copy(z, z + n, zn);
        return 0;
    }

    int jn = 0;
    for (int j = 0; j < n - 1; j++) {
        zn[jn] = z[j];
        jn++;
        if (m_insertPts.count(j)) {
            zn[jn] = 0.5*(z[j] + z[j+1]);
            jn++;
        }
    }
    zn[jn] = z[n-1];
    return 0;
}


double Refiner::remeshFromSolution(int np, const doublereal* z, const doublereal* x, const double dist_min, const double domain_size) //from MUTAGEN
// int np;   : number of points in the old mesh
// double* z : array of the old mesh points
// double* x : array of variables on the old mesh
{

  double distance;

  if (m_domain->nPoints() == 1)
  {
      m_z_new.clear();
      m_z_new.push_back(z[0]);
      distance=0.0;
      return distance;
  }

  if (m_domain->nPoints() <= 0)
  {
      distance=0.0;
      return distance;
  }

  // a few defintions
  double z_start     = z[0];
  double z_end       = z[np-1];

  // compute w_old
  // [w_old] = [m^-1]
  vector<double> w_old(np-1);
  for (int i=0;i<np-1;i++)
  {
    w_old[i] = 1.0/(z[i+1]-z[i]);
  }

  // find flame center (minimum dz)
  double dz_min = 1e9;
  double z_center = 0.0;
  for (int i=0;i<np-1;i++)
  {
    if (dz_min>z[i+1]-z[i])
    {
      dz_min = z[i+1]-z[i];
      z_center = z[i];
    }
  }

  // init w1 with a very small (but non-zero) value
  // [w1] = [m^-1]
  vector<double> w1(np-1);
  for (int i=0;i<np-1;i++)
  {
    w1[i] = 1.0e-6/(z_end - z_start);
  }

  // get gradient and curvature criterions in w1
  vector<double> vn(np);
  vector<double> gn(np);
  vector<double> gc(np-1);
  vector<double> cc(np-1);
  // loop over active components
  for (int icomp=0;icomp<m_nv;icomp++)
  {
    if (m_active[icomp])
    {

      // node values of component icomp
      for (int i=0;i<np;i++)
      {
        vn[i] = value(x,icomp,i);
      }

      // cell centered gradient of component icomp
      for (int i=0;i<np-1;i++)
      {
        gc[i] = (vn[i+1]-vn[i])/(z[i+1]-z[i]);
      }

      // node centered gradient of component icomp
      gn[0] = gc[0];
      for (int i=1;i<np-1;i++)
      {
        double zc_1  = 0.5*(z[i  ]+z[i-1]); // center of leftmost cell
        double zc_2  = 0.5*(z[i+1]+z[i  ]); // center of rightmost cell
        double alpha = (z[i]-zc_1)/(zc_2-zc_1); // interpolation coeff
        gn[i] = (1.0-alpha)*gc[i-1]+alpha*gc[i];
      }
      gn[np-1] = gc[np-2];

      // cell centered curvature of component icomp
      for (int i=0;i<np-1;i++)
      {
        cc[i] = (gn[i+1]-gn[i])/(z[i+1]-z[i]);
      }

      // max/min of value/gradient of component icomp
      double vmin = *min_element(vn.begin(), vn.end());
      double vmax = *max_element(vn.begin(), vn.end());
      double gmin = *min_element(gn.begin(), gn.end());
      double gmax = *max_element(gn.begin(), gn.end());

      // max absolute values of v and g
      double vv = std::max(fabs(vmax), fabs(vmin));
      double gg = std::max(fabs(gmax), fabs(gmin));

      if ((vmax - vmin) > m_min_range*vv)
      {
        // compute w1 for gradient
        for (int i=0;i<np-1;i++)
        {
          double w1_tmp = gc[i]/((vmax-vmin+1.0e-6)*m_slope);
          w1[i] = std::max(w1[i],w1_tmp);
        }
      }

      if ((gmax - gmin) > m_min_range*gg)
      {
        // compute w1 for curvature
        for (int i=0;i<np-1;i++)
        {
          double w1_tmp = cc[i]/((gmax-gmin+1.0e-6)*m_curve);
          w1[i] = std::max(w1[i],w1_tmp);
        }
      }

    }
  }

  // correct w1 with a broadening factor which takes into accounts the neighbours
  // the target cell ratio is equal to bndr
  // [w2] = [m^-1]
  vector<double> w2(np-1);
  // first, copy w1 in w2
  copy(w1.begin(),w1.end(),w2.begin());
  // parameters
  double bndr = m_ratio;
  double bndq = 2.0*log10(bndr);
  for (int i=0;i<np-1;i++)
  {
    double zc_i  = 0.5*(z[i+1]+z[i  ]); // center of cell i
    for (int j=0;j<np-1;j++)
    {
      double zc_j  = 0.5*(z[j+1]+z[j  ]); // center of cell j
      double bf = 1.0/(1.0/w1[j] + bndq*fabs(zc_i-zc_j));
      w2[i] = std::max(w2[i],bf);
    }
  }


  // regularization step : w3 is the convolution of w2 with a gaussian
  // w3(z) = (\int_{-alpha.sigma}^{+alpha.sigma} w2(z+c).exp(-c^2/sigma^2)dc) / (\int_{-alpha.sigma}^{+alpha.sigma}exp(-c^2/sigma^2)dc)
  // the idea is to chose alpha big enough so that the integration encompass most of the integral
  // on the other hand, alpha should be small enough so that precision is preserved
  // [w3] = [m^-1]
  vector<double> w3(np-1);
  double rwidth  = 5.0;                                   // parameter for the gaussian width
  double alpha   = 2.5;                                   // the integration interval is [-alpha.sigma;+alpha.sigma]
  int    nint    = 100;                                   // half the number of integration intervals
  double w2_max  = *max_element(w2.begin(), w2.end());    // max value of w2
  double sigma   = rwidth/w2_max;                         // width of the gaussian [m]
  double gintg   = sigma*sqrt(Pi);                        // integral of the gaussian [m]
  double dz_int  = alpha*sigma/static_cast<double>(nint); // integration step size [m]
  for (int i=0;i<np-1;i++)
  {
    // set the integral to zero before the integration loop
    w3[i] = 0.0;
    // loop over integration points
    for (int j=-nint; j<nint; j++)
    {
      // c1 and c2 are the bounds of the current integration interval
      double c  = (static_cast<double>(j)+0.5)*dz_int;
      double zc = 0.5*(z[i+1]+z[i  ]);
      double zint = zc+c;
      double f;
      // value at zc+c
      if (zint<=z_start)
      {
        f = 0.0;
      }
      else if (zint>=z_end)
      {
        f = 0.0;
      }
      else
      {
        int i1=indxtp(np,zint,z);
        f = w2[i1];
      }
      // convolution product
      f = f * exp(-c*c/(sigma*sigma));
      // integration and normalization
      w3[i] = w3[i] + f*dz_int/gintg;
    }
  }

  // compute distance between the two meshes
  distance=0.0;
  for (int i=0;i<np-1;i++)
  {
    double this_dist = abs(w3[i]-w_old[i])/(0.5*(w3[i]+w_old[i]));
    distance += this_dist;
  }
  distance = distance/np;

  if (distance <= dist_min)
  {
    distance=0.0;
    m_z_new.clear();
    for (int i=0;i<np;i++)
    {
      m_z_new.push_back(z[i]);
    }
    return distance;
  }


  // w4 = integral of w3
  // [w4] = [-]
  vector<double> w4(np);
  w4[0] = 0.0;
  for (int i=0;i<np-1;i++)
  {
    double dz  = z[i+1]-z[i];
    w4[i+1] = w4[i] + dz*w3[i];
  }


  // new size
  double w4_last = w4[np-1];
  int np_new = static_cast<int>(1.0+w4_last + 0.5);

  // final remeshing
  m_z_new.clear();
  for (int i_new=0;i_new<np_new;i_new++)
  {
    double w      = static_cast<double>(i_new);
    int    i      = indxtp(np,w,&w4[0]);
    double slope  = (z[i+1] -z[i])/(w4[i+1]-w4[i]);
    double dz     = slope*(w-w4[i]);
    double z_add = z[i] + dz;
    if (z_add<=z_center+domain_size)
    {
      m_z_new.push_back(z_add);
    }
  }

  // New size
  np_new = m_z_new.size();

  // Add points if necessary
  bool done=false;
  if (m_z_new[np_new-1]>=z_center+domain_size)
  {
    done = true;
  }
  while(!done)
  {
    double ratio_last = (m_z_new[np_new-1] - m_z_new[np_new-2])/(m_z_new[np_new-2] - m_z_new[np_new-3]);
    double z_add = m_z_new[np_new-1] + ratio_last*(m_z_new[np_new-1] - m_z_new[np_new-2]);
    m_z_new.push_back(z_add);
    np_new++;
    if (m_z_new[np_new-1]>=z_center+domain_size)
    {
      done = true;
    }
  }

/* DO NOTHING...
  m_z_new.clear();
  for (int i=0;i<np;i++)
  {
    m_z_new.push_back(z[i]);
  }
*/

  return distance;
}

void Refiner::getRatio(int np, const doublereal* z, doublereal* x) //from MUTAGEN
// int np;   : number of points in the old mesh
// double* z : array of the old mesh points
// double* x : array of variables on the old mesh
{

  if (m_domain->nPoints() == 1)
  {
      return;
  }

  if (m_domain->nPoints() <= 0)
  {
      return;
  }

  // get gradient and curvature criterions in w1
  vector<double> vn(np);
  vector<double> gn(np);
  vector<double> gc(np-1);
  vector<double> cc(np-1);
  vector<double> grad_max(np-1);
  vector<double> curve_max(np-1);
  // loop over active components
  for (int icomp=0;icomp<m_nv;icomp++)
  {
    if (m_active[icomp])
    {

      // node values of component icomp
      for (int i=0;i<np;i++)
      {
        vn[i] = value(x,icomp,i);
      }

      // cell centered gradient of component icomp
      for (int i=0;i<np-1;i++)
      {
        gc[i] = (vn[i+1]-vn[i])/(z[i+1]-z[i]);
      }

      // node centered gradient of component icomp
      gn[0] = gc[0];
      for (int i=1;i<np-1;i++)
      {
        double zc_1  = 0.5*(z[i  ]+z[i-1]); // center of leftmost cell
        double zc_2  = 0.5*(z[i+1]+z[i  ]); // center of rightmost cell
        double alpha = (z[i]-zc_1)/(zc_2-zc_1); // interpolation coeff
        gn[i] = (1.0-alpha)*gc[i-1]+alpha*gc[i];
      }
      gn[np-1] = gc[np-2];

      // cell centered curvature of component icomp
      for (int i=0;i<np-1;i++)
      {
        cc[i] = (gn[i+1]-gn[i])/(z[i+1]-z[i]);
      }

      // max/min of value/gradient of component icomp
      double vmin = *min_element(vn.begin(), vn.end());
      double vmax = *max_element(vn.begin(), vn.end());
      double gmin = *min_element(gn.begin(), gn.end());
      double gmax = *max_element(gn.begin(), gn.end());

      // max absolute values of v and g
      double vv = std::max(fabs(vmax), fabs(vmin));
      double gg = std::max(fabs(gmax), fabs(gmin));

      if ((vmax - vmin) > m_min_range*vv)
      {
        // compute w1 for gradient
        for (int i=0;i<np-1;i++)
        {
          double w1_tmp = (vn[i+1]-vn[i])/(vmax-vmin+1.0e-6);
          grad_max[i] = std::max(grad_max[i],w1_tmp);
        }
      }

      if ((gmax - gmin) > m_min_range*gg)
      {
        // compute w1 for curvature
        for (int i=0;i<np-1;i++)
        {
          double w1_tmp = (gn[i+1]-gn[i])/(gmax-gmin+1.0e-6);
          curve_max[i] = std::max(curve_max[i],w1_tmp);
        }
      }

    }
  }

  m_grad_max = *max_element(grad_max.begin(), grad_max.end());
  m_curve_max = *max_element(curve_max.begin(), curve_max.end());

}

int Refiner::indxtp(int np,double val,const double* array)
{
  if (val<=array[0])
  {
    return(0);
  }
  for (int i=0; i<np-1;i++)
  {
    if ((array[i]<=val)&&(val<=array[i+1]))
    {
       return(i);
    }
  }
  if (array[np-1]<=val)
  {
    return(np-2);
  }
  // this statement is never executed
  return(0);
}

template<class M>
bool has_key(const M& m, int j) {
    if (m.find(j) != m.end()) return true;
    return false;
}

static void r_drawline() {
    string s(78,'#');
    s += '\n';
    writelog(s.c_str());
}

/**
 * Return the square root of machine precision.
 */
static doublereal eps() {
    doublereal e = 1.0;
    while (1.0 + e != 1.0) e *= 0.5;
    return sqrt(e);
}

}
