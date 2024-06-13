/**
 * @file Sim1D.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/MultiJac.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/oneD/MultiNewton.h"
#include "cantera/oneD/refine.h"
#include "cantera/numerics/funcs.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/SolutionArray.h"
#include "cantera/numerics/Func1.h"
#include <limits>
#include <fstream>

using namespace std;

namespace Cantera
{

Sim1D::Sim1D(vector<shared_ptr<Domain1D>>& domains) :
    OneDim(domains),
    m_steady_callback(0)
{
    // resize the internal solution vector and the work array, and perform
    // domain-specific initialization of the solution vector.
    resize();
    for (size_t n = 0; n < nDomains(); n++) {
        domain(n)._getInitialSoln(m_state->data() + start(n));
    }

    // set some defaults
    m_tstep = 1.0e-5;
    m_steps = { 10 };
}

void Sim1D::setInitialGuess(const string& component, vector<double>& locs,
                            vector<double>& vals)
{
    for (size_t dom=0; dom<nDomains(); dom++) {
        Domain1D& d = domain(dom);
        size_t ncomp = d.nComponents();
        for (size_t comp=0; comp<ncomp; comp++) {
            if (d.componentName(comp)==component) {
                setProfile(dom,comp,locs,vals);
            }
        }
    }
}

void Sim1D::setValue(size_t dom, size_t comp, size_t localPoint, double value)
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_state->size(), "Sim1D::setValue",
                   "Index out of bounds: {} > {}", iloc, m_state->size());
    (*m_state)[iloc] = value;
}

double Sim1D::value(size_t dom, size_t comp, size_t localPoint) const
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_state->size(), "Sim1D::value",
                   "Index out of bounds: {} > {}", iloc, m_state->size());
    return (*m_state)[iloc];
}

double Sim1D::workValue(size_t dom, size_t comp, size_t localPoint) const
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_state->size(), "Sim1D::workValue",
                   "Index out of bounds: {} > {}", iloc, m_state->size());
    return m_xnew[iloc];
}

void Sim1D::setProfile(size_t dom, size_t comp,
                       const vector<double>& pos, const vector<double>& values)
{
    if (pos.front() != 0.0 || pos.back() != 1.0) {
        throw CanteraError("Sim1D::setProfile",
            "`pos` vector must span the range [0, 1]. Got a vector spanning "
            "[{}, {}] instead.", pos.front(), pos.back());
    }
    Domain1D& d = domain(dom);
    double z0 = d.zmin();
    double z1 = d.zmax();
    for (size_t n = 0; n < d.nPoints(); n++) {
        double zpt = d.z(n);
        double frac = (zpt - z0)/(z1 - z0);
        double v = linearInterp(frac, pos, values);
        setValue(dom, comp, n, v);
    }
}

void Sim1D::save(const string& fname, const string& name, const string& desc,
                 bool overwrite, int compression, const string& basis)
{
    size_t dot = fname.find_last_of(".");
    string extension = (dot != npos) ? toLowerCopy(fname.substr(dot+1)) : "";
    if (extension == "csv") {
        for (auto dom : m_dom) {
            auto arr = dom->asArray(m_state->data() + dom->loc());
            if (dom->size() > 1) {
                arr->writeEntry(fname, overwrite, basis);
                break;
            }
        }
        return;
    }
    if (basis != "") {
        warn_user("Sim1D::save",
            "Species basis '{}' not implemented for HDF5 or YAML output.", basis);
    }
    if (extension == "h5" || extension == "hdf"  || extension == "hdf5") {
        SolutionArray::writeHeader(fname, name, desc, overwrite);
        for (auto dom : m_dom) {
            auto arr = dom->asArray(m_state->data() + dom->loc());
            arr->writeEntry(fname, name, dom->id(), overwrite, compression);
        }
        return;
    }
    if (extension == "yaml" || extension == "yml") {
        // Check for an existing file and load it if present
        AnyMap data;
        if (std::ifstream(fname).good()) {
            data = AnyMap::fromYamlFile(fname);
        }
        SolutionArray::writeHeader(data, name, desc, overwrite);

        for (auto dom : m_dom) {
            auto arr = dom->asArray(m_state->data() + dom->loc());
            arr->writeEntry(data, name, dom->id(), overwrite);
        }

        // Write the output file and remove the now-outdated cached file
        std::ofstream out(fname);
        out << data.toYamlString();
        AnyMap::clearCachedFile(fname);
        return;
    }
    throw CanteraError("Sim1D::save", "Unsupported file format '{}'.", extension);
}

void Sim1D::saveResidual(const string& fname, const string& name,
                         const string& desc, bool overwrite, int compression)
{
    vector<double> res(m_state->size(), -999);
    OneDim::eval(npos, m_state->data(), &res[0], 0.0);
    // Temporarily put the residual into m_state, since this is the vector that the
    // save() function reads.
    vector<double> backup(*m_state);
    *m_state = res;
    save(fname, name, desc, overwrite, compression);
    *m_state = backup;
}

namespace { // restrict scope of helper function to local translation unit

//! convert data format used by Python h5py export (Cantera < 3.0)
AnyMap legacyH5(shared_ptr<SolutionArray> arr, const AnyMap& header={})
{
    auto meta = arr->meta();
    AnyMap out;

    map<string, string> meta_pairs = {
        {"type", "Domain1D_type"},
        {"name", "name"},
        {"emissivity-left", "emissivity_left"},
        {"emissivity-right", "emissivity_right"},
    };
    for (const auto& [newName, oldName] : meta_pairs) {
        if (meta.hasKey(oldName)) {
            out[newName] = meta[oldName];
        }
    }

    map<string, string> tol_pairs = {
        {"transient-abstol", "transient_abstol"},
        {"steady-abstol", "steady_abstol"},
        {"transient-reltol", "transient_reltol"},
        {"steady-reltol", "steady_reltol"},
    };
    for (const auto& [newName, oldName] : tol_pairs) {
        if (meta.hasKey(oldName)) {
            out["tolerances"][newName] = meta[oldName];
        }
    }

    if (meta.hasKey("phase")) {
        out["phase"]["name"] = meta["phase"]["name"];
        out["phase"]["source"] = meta["phase"]["source"];
    }

    if (arr->size() <= 1) {
        return out;
    }

    map<string, string> header_pairs = {
        {"transport-model", "transport_model"},
        {"radiation-enabled", "radiation_enabled"},
        {"energy-enabled", "energy_enabled"},
        {"Soret-enabled", "soret_enabled"},
        {"species-enabled", "species_enabled"},
    };
    for (const auto& [newName, oldName] : header_pairs) {
        if (header.hasKey(oldName)) {
            out[newName] = header[oldName];
        }
    }

    map<string, string> refiner_pairs = {
        {"ratio", "ratio"},
        {"slope", "slope"},
        {"curve", "curve"},
        {"prune", "prune"},
        // {"grid-min", "???"}, // missing
        {"max-points", "max_grid_points"},
    };
    for (const auto& [newName, oldName] : refiner_pairs) {
        if (header.hasKey(oldName)) {
            out["refine-criteria"][newName] = header[oldName];
        }
    }

    if (header.hasKey("fixed_temperature")) {
        double temp = header.getDouble("fixed_temperature", -1.);
        auto profile = arr->getComponent("T").as<vector<double>>();
        int ix = 0;
        while (profile[ix] <= temp && ix < arr->size()) {
            ix++;
        }
        if (ix != 0) {
            auto grid = arr->getComponent("grid").as<vector<double>>();
            out["fixed-point"]["location"] = grid[ix - 1];
            out["fixed-point"]["temperature"] = temp;
        }
    }

    return out;
}

} // end unnamed namespace

AnyMap Sim1D::restore(const string& fname, const string& name)
{
    size_t dot = fname.find_last_of(".");
    string extension = (dot != npos) ? toLowerCopy(fname.substr(dot+1)) : "";
    if (extension == "xml") {
        throw CanteraError("Sim1D::restore",
                           "Restoring from XML is no longer supported.");
    }
    AnyMap header;
    if (extension == "h5" || extension == "hdf"  || extension == "hdf5") {
        map<string, shared_ptr<SolutionArray>> arrs;
        header = SolutionArray::readHeader(fname, name);

        for (auto dom : m_dom) {
            auto arr = SolutionArray::create(dom->solution());
            try {
                arr->readEntry(fname, name, dom->id());
            } catch (CanteraError& err) {
                throw CanteraError("Sim1D::restore",
                    "Encountered exception when reading entry '{}' from '{}':\n{}",
                    name, fname, err.getMessage());
            }
            dom->resize(dom->nComponents(), arr->size());
            if (!header.hasKey("generator")) {
                arr->meta() = legacyH5(arr, header);
            }
            arrs[dom->id()] = arr;
        }
        resize();
        m_xlast_ts.clear();
        for (auto dom : m_dom) {
            try {
                dom->fromArray(*arrs[dom->id()], m_state->data() + dom->loc());
            } catch (CanteraError& err) {
                throw CanteraError("Sim1D::restore",
                    "Encountered exception when restoring domain '{}' from HDF:\n{}",
                    dom->id(), err.getMessage());
            }
        }
        finalize();
    } else if (extension == "yaml" || extension == "yml") {
        AnyMap root = AnyMap::fromYamlFile(fname);
        map<string, shared_ptr<SolutionArray>> arrs;
        header = SolutionArray::readHeader(root, name);

        for (auto dom : m_dom) {
            auto arr = SolutionArray::create(dom->solution());
            try {
                arr->readEntry(root, name, dom->id());
            } catch (CanteraError& err) {
                throw CanteraError("Sim1D::restore",
                    "Encountered exception when reading entry '{}' from '{}':\n{}",
                    name, fname, err.getMessage());
            }
            dom->resize(dom->nComponents(), arr->size());
            arrs[dom->id()] = arr;
        }
        resize();
        m_xlast_ts.clear();
        for (auto dom : m_dom) {
            try {
                dom->fromArray(*arrs[dom->id()], m_state->data() + dom->loc());
            } catch (CanteraError& err) {
                throw CanteraError("Sim1D::restore",
                    "Encountered exception when restoring domain '{}' from YAML:\n{}",
                    dom->id(), err.getMessage());
            }
        }
        finalize();
    } else {
        throw CanteraError("Sim1D::restore",
                           "Unknown file extension '{}'; supported extensions include "
                           "'h5'/'hdf'/'hdf5' and 'yml'/'yaml'.", extension);
    }
    return header;
}

void Sim1D::setFlatProfile(size_t dom, size_t comp, double v)
{
    size_t np = domain(dom).nPoints();
    for (size_t n = 0; n < np; n++) {
        setValue(dom, comp, n, v);
    }
}

void Sim1D::show(ostream& s)
{
    for (size_t n = 0; n < nDomains(); n++) {
        if (domain(n).type() != "empty") {
            domain(n).show(s, m_state->data() + start(n));
        }
    }
}

void Sim1D::show()
{
    for (size_t n = 0; n < nDomains(); n++) {
        if (domain(n).type() != "empty") {
            writelog("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "+domain(n).id()
                     +" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
            domain(n).show(m_state->data() + start(n));
        }
    }
}

void Sim1D::restoreTimeSteppingSolution()
{
    if (m_xlast_ts.empty()) {
        throw CanteraError("Sim1D::restoreTimeSteppingSolution",
                           "No successful time steps taken on this grid.");
    }
    *m_state = m_xlast_ts;
}

void Sim1D::restoreSteadySolution()
{
    if (m_xlast_ss.empty()) {
        throw CanteraError("Sim1D::restoreSteadySolution",
                           "No successful steady state solution");
    }
    *m_state = m_xlast_ss;
    for (size_t n = 0; n < nDomains(); n++) {
        vector<double>& z = m_grid_last_ss[n];
        domain(n).setupGrid(z.size(), z.data());
    }
}

void Sim1D::getInitialSoln()
{
    for (size_t n = 0; n < nDomains(); n++) {
        domain(n)._getInitialSoln(m_state->data() + start(n));
    }
}

void Sim1D::finalize()
{
    for (size_t n = 0; n < nDomains(); n++) {
        domain(n)._finalize(m_state->data() + start(n));
    }
}

void Sim1D::setTimeStep(double stepsize, size_t n, const int* tsteps)
{
    m_tstep = stepsize;
    m_steps.resize(n);
    for (size_t i = 0; i < n; i++) {
        m_steps[i] = tsteps[i];
    }
}

int Sim1D::newtonSolve(int loglevel)
{
    int m = OneDim::solve(m_state->data(), m_xnew.data(), loglevel);
    if (m >= 0) {
        *m_state = m_xnew;
        return 0;
    } else if (m > -10) {
        return -1;
    } else {
        throw CanteraError("Sim1D::newtonSolve",
                           "ERROR: OneDim::solve returned m = {}", m);
    }
}

void Sim1D::solve(int loglevel, bool refine_grid)
{
    int new_points = 1;
    double dt = m_tstep;
    m_nsteps = 0;
    finalize();

    while (new_points > 0) {
        size_t istep = 0;
        int nsteps = m_steps[istep];

        bool ok = false;
        if (loglevel > 0) {
            writeline('.', 78, true, true);
        }
        while (!ok) {
            // Attempt to solve the steady problem
            setSteadyMode();
            newton().setOptions(m_ss_jac_age);
            debuglog("Attempt Newton solution of steady-state problem...", loglevel);
            int status = newtonSolve(loglevel-1);

            if (status == 0) {
                if (loglevel > 0) {
                    writelog("    success.\n\n");
                    writelog("Problem solved on [");
                    for (size_t mm = 1; mm < nDomains(); mm+=2) {
                        writelog("{}", domain(mm).nPoints());
                        if (mm + 2 < nDomains()) {
                            writelog(", ");
                        }
                    }
                    writelog("] point grid(s).\n");
                }
                if (m_steady_callback) {
                    m_steady_callback->eval(0);
                }

                if (loglevel > 6) {
                    save("debug_sim1d.yaml", "solution",
                         "After successful Newton solve", true);
                }
                if (loglevel > 7) {
                    saveResidual("debug_sim1d.yaml", "residual",
                                 "After successful Newton solve", true);
                }
                ok = true;
            } else {
                debuglog("    failure. \n", loglevel);
                if (loglevel > 6) {
                    save("debug_sim1d.yaml", "solution",
                         "After unsuccessful Newton solve", true);
                }
                if (loglevel > 7) {
                    saveResidual("debug_sim1d.yaml", "residual",
                                 "After unsuccessful Newton solve", true);
                }
                if (loglevel > 0) {
                    writelog("Take {} timesteps   ", nsteps);
                }
                dt = timeStep(nsteps, dt, m_state->data(), m_xnew.data(), loglevel-1);
                m_xlast_ts = *m_state;
                if (loglevel > 6) {
                    save("debug_sim1d.yaml", "solution", "After timestepping", true);
                }
                if (loglevel > 7) {
                    saveResidual("debug_sim1d.yaml", "residual",
                                 "After timestepping", true);
                }

                if (loglevel == 1) {
                    writelog(" {:10.4g} {:10.4g}\n", dt,
                             log10(ssnorm(m_state->data(), m_xnew.data())));
                }
                istep++;
                if (istep >= m_steps.size()) {
                    nsteps = m_steps.back();
                } else {
                    nsteps = m_steps[istep];
                }
                dt = std::min(dt, m_tmax);
            }
        }
        if (loglevel > 0) {
            writeline('.', 78, true, true);
        }
        if (loglevel > 2) {
            show();
        }

        if (refine_grid) {
            new_points = refine(loglevel);
            if (new_points) {
                // If the grid has changed, preemptively reduce the timestep
                // to avoid multiple successive failed time steps.
                dt = m_tstep;
            }
            if (new_points && loglevel > 6) {
                save("debug_sim1d.yaml", "solution", "After regridding", true);
            }
            if (new_points && loglevel > 7) {
                saveResidual("debug_sim1d.yaml", "residual",
                             "After regridding", true);
            }
        } else {
            debuglog("grid refinement disabled.\n", loglevel);
            new_points = 0;
        }
    }
}

int Sim1D::refine(int loglevel)
{
    int ianalyze, np = 0;
    vector<double> znew, xnew;
    vector<size_t> dsize;

    m_xlast_ss = *m_state;
    m_grid_last_ss.clear();

    for (size_t n = 0; n < nDomains(); n++) {
        Domain1D& d = domain(n);
        Refiner& r = d.refiner();

        // Save the old grid corresponding to the converged solution
        m_grid_last_ss.push_back(d.grid());

        // determine where new points are needed
        ianalyze = r.analyze(d.grid().size(), d.grid().data(),
                             m_state->data() + start(n));
        if (ianalyze < 0) {
            return ianalyze;
        }

        if (loglevel > 0) {
            r.show();
        }

        np += r.nNewPoints();
        size_t comp = d.nComponents();

        // loop over points in the current grid
        size_t npnow = d.nPoints();
        size_t nstart = znew.size();
        for (size_t m = 0; m < npnow; m++) {
            if (r.keepPoint(m)) {
                // add the current grid point to the new grid
                znew.push_back(d.grid(m));

                // do the same for the solution at this point
                for (size_t i = 0; i < comp; i++) {
                    xnew.push_back(value(n, i, m));
                }

                // now check whether a new point is needed in the interval to
                // the right of point m, and if so, add entries to znew and xnew
                // for this new point
                if (r.newPointNeeded(m) && m + 1 < npnow) {
                    // add new point at midpoint
                    double zmid = 0.5*(d.grid(m) + d.grid(m+1));
                    znew.push_back(zmid);
                    np++;

                    // for each component, linearly interpolate
                    // the solution to this point
                    for (size_t i = 0; i < comp; i++) {
                        double xmid = 0.5*(value(n, i, m) + value(n, i, m+1));
                        xnew.push_back(xmid);
                    }
                }
            } else {
                if (loglevel > 0) {
                    writelog("refine: discarding point at {}\n", d.grid(m));
                }
            }
        }
        dsize.push_back(znew.size() - nstart);
    }

    // At this point, the new grid znew and the new solution vector xnew have
    // been constructed, but the domains themselves have not yet been modified.
    // Now update each domain with the new grid.

    size_t gridstart = 0, gridsize;
    for (size_t n = 0; n < nDomains(); n++) {
        Domain1D& d = domain(n);
        gridsize = dsize[n];
        if (gridsize != 0) {
            d.setupGrid(gridsize, &znew[gridstart]);
        }
        gridstart += gridsize;
    }

    // Replace the current solution vector with the new one
    *m_state = xnew;
    resize();
    finalize();
    return np;
}

int Sim1D::setFixedTemperature(double t)
{
    int np = 0;
    vector<double> znew, xnew;
    double zfixed = 0.0;
    double z1 = 0.0, z2 = 0.0;
    vector<size_t> dsize;

    for (size_t n = 0; n < nDomains(); n++) {
        Domain1D& d = domain(n);
        size_t comp = d.nComponents();
        size_t mfixed = npos;

        // loop over current grid to determine where new point is needed
        StFlow* d_free = dynamic_cast<StFlow*>(&domain(n));
        size_t npnow = d.nPoints();
        size_t nstart = znew.size();
        if (d_free && d_free->isFree()) {
            for (size_t m = 0; m < npnow - 1; m++) {
                bool fixedpt = false;
                double t1 = value(n, c_offset_T, m);
                double t2 = value(n, c_offset_T, m + 1);
                // threshold to avoid adding new point too close to existing point
                double thresh = min(1., 1.e-1 * (t2 - t1));
                z1 = d.grid(m);
                z2 = d.grid(m + 1);
                if (fabs(t - t1) <= thresh) {
                    zfixed = z1;
                    fixedpt = true;
                } else if (fabs(t2 - t) <= thresh) {
                    zfixed = z2;
                    fixedpt = true;
                } else if ((t1 < t) && (t < t2)) {
                    mfixed = m;
                    zfixed = (z1 - z2) / (t1 - t2) * (t - t2) + z2;
                    fixedpt = true;
                }

                if (fixedpt) {
                    d_free->m_zfixed = zfixed;
                    d_free->m_tfixed = t;
                    break;
                }
            }
        }

        // copy solution domain and push back values
        for (size_t m = 0; m < npnow; m++) {
            // add the current grid point to the new grid
            znew.push_back(d.grid(m));

            // do the same for the solution at this point
            for (size_t i = 0; i < comp; i++) {
                xnew.push_back(value(n, i, m));
            }
            if (m == mfixed) {
                // add new point at zfixed (mfixed is not npos)
                znew.push_back(zfixed);
                np++;
                double interp_factor = (zfixed - z2) / (z1 - z2);
                // for each component, linearly interpolate
                // the solution to this point
                for (size_t i = 0; i < comp; i++) {
                    double xmid = interp_factor*(value(n, i, m) - value(n, i, m+1)) + value(n,i,m+1);
                    xnew.push_back(xmid);
                }
            }
        }
        dsize.push_back(znew.size() - nstart);
    }

    // At this point, the new grid znew and the new solution vector xnew have
    // been constructed, but the domains themselves have not yet been modified.
    // Now update each domain with the new grid.
    size_t gridstart = 0;
    for (size_t n = 0; n < nDomains(); n++) {
        Domain1D& d = domain(n);
        size_t gridsize = dsize[n];
        d.setupGrid(gridsize, &znew[gridstart]);
        gridstart += gridsize;
    }

    // Replace the current solution vector with the new one
    *m_state = xnew;

    resize();
    finalize();
    return np;
}

double Sim1D::fixedTemperature()
{
    double t_fixed = std::numeric_limits<double>::quiet_NaN();
    for (size_t n = 0; n < nDomains(); n++) {
        StFlow* d = dynamic_cast<StFlow*>(&domain(n));
        if (d && d->isFree() && d->m_tfixed > 0) {
            t_fixed = d->m_tfixed;
            break;
        }
    }
    return t_fixed;
}

double Sim1D::fixedTemperatureLocation()
{
    double z_fixed = std::numeric_limits<double>::quiet_NaN();
    for (size_t n = 0; n < nDomains(); n++) {
        StFlow* d = dynamic_cast<StFlow*>(&domain(n));
        if (d && d->isFree() && d->m_tfixed > 0) {
            z_fixed = d->m_zfixed;
            break;
        }
    }
    return z_fixed;
}

void Sim1D::setLeftControlPoint(double temperature)
{
    bool two_point_domain_found = false;
    for (size_t n = 0; n < nDomains(); n++) {
        Domain1D& d = domain(n);

        // Skip if the domain type doesn't match
        if (d.domainType() != "axisymmetric-flow") {
            continue;
        }

        StFlow& d_axis = dynamic_cast<StFlow&>(domain(n));
        size_t np = d_axis.nPoints();

        // Skip if two-point control is not enabled
        if (!d_axis.twoPointControlEnabled()) {
            continue;
        }
        two_point_domain_found = true; // At least one domain with two-point control enabled was found

        for (size_t m = 0; m < np-1; m++) {
            if ((value(n,c_offset_T,m) - temperature) * (value(n,c_offset_T,m+1) - temperature) < 0.0) {
                // Pick the coordinate of the point with the temperature closest to the desired temperature
                int closest_index = 0;
                if (std::abs(value(n,c_offset_T,m) - temperature) < std::abs(value(n,c_offset_T,m+1) - temperature)) {
                    closest_index = m;
                } else {
                    closest_index = m+1;
                }
                d_axis.setLeftControlPointCoordinate(d_axis.grid(closest_index));
                d_axis.setLeftControlPointTemperature(value(n,c_offset_T,closest_index));
                return;
            }
        }
    }

    if (!two_point_domain_found){
        throw CanteraError("Sim1D::setLeftControlPoint",
            "No domain with two-point control enabled was found.");
    } else {
        throw CanteraError("Sim1D::setLeftControlPoint",
            "No control point with temperature {} was able to be found in the flame's temperature range.", temperature);
    }

}

void Sim1D::setRightControlPoint(double temperature)
{
    bool two_point_domain_found = false;
    for (size_t n = 0; n < nDomains(); n++) {
        Domain1D& d = domain(n);

        // Skip if the domain type doesn't match
        if (d.domainType() != "axisymmetric-flow") {
            continue;
        }

        StFlow& d_axis = dynamic_cast<StFlow&>(domain(n));
        size_t np = d_axis.nPoints();

        // Skip if two-point control is not enabled
        if (!d_axis.twoPointControlEnabled()) {
            continue;
        }
        two_point_domain_found = true; // At least one domain with two-point control enabled was found


        for (size_t m = np-1; m > 0; m--) {
            if ((value(n,c_offset_T,m) - temperature) * (value(n,c_offset_T,m-1) - temperature) < 0.0) {
                // Pick the coordinate of the point with the temperature closest to the desired temperature
                int closest_index = 0;
                if (std::abs(value(n,c_offset_T,m) - temperature) < std::abs(value(n,c_offset_T,m-1) - temperature)) {
                    closest_index = m;
                } else {
                    closest_index = m-1;
                }
                d_axis.setRightControlPointCoordinate(d_axis.grid(closest_index));
                d_axis.setRightControlPointTemperature(value(n,c_offset_T,closest_index));
                return;
            }
        }
    }

    if (!two_point_domain_found){
        throw CanteraError("Sim1D::setRightControlPoint",
            "No domain with two-point control enabled was found.");
    } else {
        throw CanteraError("Sim1D::setRightControlPoint",
            "No control point with temperature {} was able to be found in the flame's temperature range.", temperature);
    }

}

void Sim1D::setRefineCriteria(int dom, double ratio,
                              double slope, double curve, double prune)
{
    if (dom >= 0) {
        Refiner& r = domain(dom).refiner();
        r.setCriteria(ratio, slope, curve, prune);
    } else {
        for (size_t n = 0; n < nDomains(); n++) {
            Refiner& r = domain(n).refiner();
            r.setCriteria(ratio, slope, curve, prune);
        }
    }
}

vector<double> Sim1D::getRefineCriteria(int dom)
{
   if (dom >= 0) {
       Refiner& r = domain(dom).refiner();
       return r.getCriteria();
   } else {
       throw CanteraError("Sim1D::getRefineCriteria",
           "Must specify domain to get criteria from");
   }
}

void Sim1D::setGridMin(int dom, double gridmin)
{
    if (dom >= 0) {
        Refiner& r = domain(dom).refiner();
        r.setGridMin(gridmin);
    } else {
        for (size_t n = 0; n < nDomains(); n++) {
            Refiner& r = domain(n).refiner();
            r.setGridMin(gridmin);
        }
    }
}

void Sim1D::setMaxGridPoints(int dom, int npoints)
{
    if (dom >= 0) {
        Refiner& r = domain(dom).refiner();
        r.setMaxPoints(npoints);
    } else {
        for (size_t n = 0; n < nDomains(); n++) {
            Refiner& r = domain(n).refiner();
            r.setMaxPoints(npoints);
        }
    }
}

size_t Sim1D::maxGridPoints(size_t dom)
{
    Refiner& r = domain(dom).refiner();
    return r.maxPoints();
}

double Sim1D::jacobian(int i, int j)
{
    return OneDim::jacobian().value(i,j);
}

void Sim1D::evalSSJacobian()
{
    OneDim::evalSSJacobian(m_state->data(), m_xnew.data());
}

void Sim1D::solveAdjoint(const double* b, double* lambda)
{
    for (auto& D : m_dom) {
        D->forceFullUpdate(true);
    }
    evalSSJacobian();
    for (auto& D : m_dom) {
        D->forceFullUpdate(false);
    }

    // Form J^T
    size_t bw = bandwidth();
    BandMatrix Jt(size(), bw, bw);
    for (size_t i = 0; i < size(); i++) {
        size_t j1 = (i > bw) ? i - bw : 0;
        size_t j2 = (i + bw >= size()) ? size() - 1: i + bw;
        for (size_t j = j1; j <= j2; j++) {
            Jt(j,i) = m_jac->value(i,j);
        }
    }

    Jt.solve(b, lambda);
}

void Sim1D::resize()
{
    OneDim::resize();
    m_xnew.resize(size(), 0.0);
}

}
