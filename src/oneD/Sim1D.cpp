/**
 * @file Sim1D.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/MultiJac.h"
#include "cantera/oneD/Flow1D.h"
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
        domain(n)._getInitialSoln(
            span<double>(m_state->data() + start(n), domain(n).size()));
    }
}

void Sim1D::_setValue(size_t dom, size_t comp, size_t localPoint, double value)
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_state->size(), "Sim1D::setValue",
                   "Index out of bounds: {} > {}", iloc, m_state->size());
    (*m_state)[iloc] = value;
}

double Sim1D::_value(size_t dom, size_t comp, size_t localPoint) const
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_state->size(), "Sim1D::value",
                   "Index out of bounds: {} > {}", iloc, m_state->size());
    return (*m_state)[iloc];
}

double Sim1D::_workValue(size_t dom, size_t comp, size_t localPoint) const
{
    size_t iloc = domain(dom).loc() + domain(dom).index(comp, localPoint);
    AssertThrowMsg(iloc < m_state->size(), "Sim1D::workValue",
                   "Index out of bounds: {} > {}", iloc, m_state->size());
    return m_xnew[iloc];
}

void Sim1D::save(const string& fname, const string& name, const string& desc,
                 bool overwrite, int compression, const string& basis)
{
    size_t dot = fname.find_last_of(".");
    string extension = (dot != npos) ? toLowerCopy(fname.substr(dot+1)) : "";
    if (extension == "csv") {
        for (auto dom : m_dom) {
            auto arr = dom->toArray();
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
            auto arr = dom->toArray();
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
            auto arr = dom->toArray();
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
    OneDim::eval(npos, *m_state, res, 0.0);
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
            auto arr = SolutionArray::create(dom->phase());
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
                dom->fromArray(arrs[dom->id()]);
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
            auto arr = SolutionArray::create(dom->phase());
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
                dom->fromArray(arrs[dom->id()]);
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

void Sim1D::_restore(const string& fname, const string& name)
{
    restore(fname, name);
}

void Sim1D::show()
{
    for (size_t n = 0; n < nDomains(); n++) {
        if (domain(n).domainType() != "empty") {
            writelog("\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "+domain(n).id()
                     +" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
            domain(n).show(
                span<const double>(m_state->data() + start(n), domain(n).size()));
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
        domain(n).setupGrid(z);
    }
}

void Sim1D::getInitialSoln()
{
    for (size_t n = 0; n < nDomains(); n++) {
        domain(n)._getInitialSoln(
            span<double>(m_state->data() + start(n), domain(n).size()));
    }
}

void Sim1D::finalize()
{
    for (size_t n = 0; n < nDomains(); n++) {
        domain(n)._finalize(
            span<const double>(m_state->data() + start(n), domain(n).size()));
    }
}

void Sim1D::solve(int loglevel, bool refine_grid)
{
    int new_points = 1;
    m_attempt_counter = 0;
    finalize();
    if (loglevel > 6) {
        clearDebugFile();
    }

    while (new_points > 0) {
        SteadyStateSystem::solve(loglevel);
        if (loglevel > 0) {
            writelog("\nNewton steady-state solve succeeded.\n\n");
            writelog("Problem solved on [");
            for (size_t mm = 1; mm < nDomains(); mm+=2) {
                writelog("{}", domain(mm).nPoints());
                if (mm + 2 < nDomains()) {
                    writelog(", ");
                }
            }
            writelog("] point grid(s).\n");
            if (loglevel > 3) {
                show();
            }
        }
        if (m_steady_callback) {
            m_steady_callback->eval(0);
        }
        if (refine_grid) {
            new_points = refine(loglevel);
            writeDebugInfo("Regridding", "After regridding", loglevel,
                           m_attempt_counter);
        } else {
            debuglog("grid refinement disabled.\n", loglevel);
            new_points = 0;
        }
    }
    if (new_points < 0) {
        // If the solver finished after removing grid points, do one final evaluation
        // of the governing equations to update internal arrays in each domain that may
        // be used for data saved in output files.
        for (auto dom : m_dom) {
            span<const double> x(m_state->data() + dom->loc(), dom->size());
            span<double> r(m_xnew.data() + dom->loc(), dom->size());
            span<int> mask(m_mask.data() + dom->loc(), dom->size());
            dom->eval(npos, x, r, mask);
        }
    }
}

int Sim1D::refine(int loglevel)
{
    int added = 0;
    int discarded = 0;
    vector<double> znew, xnew;
    vector<size_t> dsize;

    m_xlast_ss = *m_state;
    m_grid_last_ss.clear();

    for (size_t n = 0; n < nDomains(); n++) {
        Domain1D& d = domain(n);
        Refiner& r = d.refiner();

        // Save the old grid corresponding to the converged solution
        auto grid = d.grid();
        m_grid_last_ss.emplace_back(grid.begin(), grid.end());

        // determine where new points are needed
        r.analyze(grid.size(), grid,
                  span<const double>(m_state->data() + start(n), d.size()));

        if (loglevel > 0) {
            r.show();
        }

        added += r.nNewPoints();
        size_t comp = d.nComponents();

        // loop over points in the current grid
        size_t npnow = d.nPoints();
        size_t nstart = znew.size();
        for (size_t m = 0; m < npnow; m++) {
            if (r.keepPoint(m)) {
                // add the current grid point to the new grid
                znew.push_back(d.z(m));

                // do the same for the solution at this point
                for (size_t i = 0; i < comp; i++) {
                    xnew.push_back(_value(n, i, m));
                }

                // now check whether a new point is needed in the interval to
                // the right of point m, and if so, add entries to znew and xnew
                // for this new point
                if (r.newPointNeeded(m) && m + 1 < npnow) {
                    // add new point at midpoint
                    double zmid = 0.5*(d.z(m) + d.z(m+1));
                    znew.push_back(zmid);
                    added++;

                    // for each component, linearly interpolate the solution to this
                    // point
                    for (size_t i = 0; i < comp; i++) {
                        double xmid = 0.5*(_value(n, i, m) + _value(n, i, m+1));
                        xnew.push_back(xmid);
                    }
                }
            } else {
                discarded++;
                if (loglevel > 0) {
                    writelog("refine: discarding point at {}\n", d.z(m));
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
            d.setupGrid(span<const double>(znew.data() + gridstart, gridsize));
        }
        gridstart += gridsize;
    }

    // Replace the current solution vector with the new one
    *m_state = xnew;
    resize();
    finalize();
    return added || -discarded;
}

void Sim1D::clearDebugFile()
{
    std::filesystem::remove("debug_sim1d.yaml");
}

void Sim1D::writeDebugInfo(const string& header_suffix, const string& message,
                           int loglevel, int attempt_counter)
{
    string file_header;
    if (loglevel > 6) {
        file_header = fmt::format("solution_{}_{}", attempt_counter, header_suffix);
        save("debug_sim1d.yaml", file_header, message, true);
    }
    if (loglevel > 7) {
        file_header = fmt::format("residual_{}_{}", attempt_counter, header_suffix);
        saveResidual("debug_sim1d.yaml", file_header, message, true);
    }
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
        Flow1D* d_free = dynamic_cast<Flow1D*>(&domain(n));
        size_t npnow = d.nPoints();
        size_t nstart = znew.size();
        if (d_free && d_free->isFree()) {
            for (size_t m = 0; m < npnow - 1; m++) {
                bool fixedpt = false;
                double t1 = _value(n, c_offset_T, m);
                double t2 = _value(n, c_offset_T, m + 1);
                // threshold to avoid adding new point too close to existing point
                double thresh = min(1., 1.e-1 * (t2 - t1));
                z1 = d.z(m);
                z2 = d.z(m + 1);
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
            znew.push_back(d.z(m));

            // do the same for the solution at this point
            for (size_t i = 0; i < comp; i++) {
                xnew.push_back(_value(n, i, m));
            }
            if (m == mfixed) {
                // add new point at zfixed (mfixed is not npos)
                znew.push_back(zfixed);
                np++;
                double interp_factor = (zfixed - z2) / (z1 - z2);
                // for each component, linearly interpolate
                // the solution to this point
                for (size_t i = 0; i < comp; i++) {
                    double xmid = interp_factor*(
                        _value(n, i, m) - _value(n, i, m+1)) + _value(n,i,m+1);
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
        d.setupGrid(span<const double>(znew.data() + gridstart, gridsize));
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
        Flow1D* d = dynamic_cast<Flow1D*>(&domain(n));
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
        Flow1D* d = dynamic_cast<Flow1D*>(&domain(n));
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

        Flow1D& d_axis = dynamic_cast<Flow1D&>(domain(n));
        size_t np = d_axis.nPoints();

        // Skip if two-point control is not enabled
        if (!d_axis.twoPointControlEnabled()) {
            continue;
        }
        two_point_domain_found = true;

        double current_val, next_val;
        for (size_t m = 0; m < np-1; m++) {
            current_val = _value(n,c_offset_T,m);
            next_val = _value(n,c_offset_T,m+1);
            if ((current_val - temperature) * (next_val - temperature) < 0.0) {
                // Pick the coordinate of the point with the temperature closest
                // to the desired temperature
                size_t index = 0;
                if (std::abs(current_val - temperature) <
                    std::abs(next_val - temperature)) {
                    index = m;
                } else {
                    index = m+1;
                }
                d_axis.setLeftControlPointCoordinate(d_axis.z(index));
                d_axis.setLeftControlPointTemperature(_value(n,c_offset_T,index));
                return;
            }
        }
    }

    if (!two_point_domain_found) {
        throw CanteraError("Sim1D::setLeftControlPoint",
            "No domain with two-point control enabled was found.");
    } else {
        throw CanteraError("Sim1D::setLeftControlPoint",
            "No control point with temperature {} was able to be found in the"
            "flame's temperature range.", temperature);
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

        Flow1D& d_axis = dynamic_cast<Flow1D&>(domain(n));
        size_t np = d_axis.nPoints();

        // Skip if two-point control is not enabled
        if (!d_axis.twoPointControlEnabled()) {
            continue;
        }
        two_point_domain_found = true;

        double current_val, next_val;
        for (size_t m = np-1; m > 0; m--) {
            current_val = _value(n,c_offset_T,m);
            next_val = _value(n,c_offset_T,m-1);
            if ((current_val - temperature) * (next_val - temperature) < 0.0) {
                // Pick the coordinate of the point with the temperature closest
                // to the desired temperature
                size_t index = 0;
                if (std::abs(current_val - temperature) <
                    std::abs(next_val - temperature)) {
                    index = m;
                } else {
                    index = m-1;
                }
                d_axis.setRightControlPointCoordinate(d_axis.z(index));
                d_axis.setRightControlPointTemperature(_value(n,c_offset_T,index));
                return;
            }
        }
    }

    if (!two_point_domain_found) {
        throw CanteraError("Sim1D::setRightControlPoint",
            "No domain with two-point control enabled was found.");
    } else {
        throw CanteraError("Sim1D::setRightControlPoint",
            "No control point with temperature {} was able to be found in the"
            "flame's temperature range.", temperature);
    }

}

void Sim1D::setRefineCriteria(int dom, double ratio,
                              double slope, double curve, double prune)
{
    if (dom >= 0) {
        domain(dom).setRefineCriteria(ratio, slope, curve, prune);
    } else {
        for (size_t n = 0; n < nDomains(); n++) {
            domain(n).setRefineCriteria(ratio, slope, curve, prune);
        }
    }
}

vector<double> Sim1D::getRefineCriteria(int dom)
{
   if (dom >= 0) {
       return domain(dom).getRefineCriteria();
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

void Sim1D::evalSSJacobian()
{
    OneDim::evalSSJacobian(*m_state);
}

void Sim1D::solveAdjoint(span<const double> b, span<double> lambda)
{
    for (auto& D : m_dom) {
        D->forceFullUpdate(true);
    }
    evalSSJacobian();
    for (auto& D : m_dom) {
        D->forceFullUpdate(false);
    }

    auto multijac = dynamic_pointer_cast<MultiJac>(m_jac);
    if (!multijac) {
        throw CanteraError("Sim1D::solveAdjoint", "Banded (MultiJac) required");
    }
    // Form J^T
    size_t bw = bandwidth();
    BandMatrix Jt(size(), bw, bw);
    for (size_t i = 0; i < size(); i++) {
        size_t j1 = (i > bw) ? i - bw : 0;
        size_t j2 = (i + bw >= size()) ? size() - 1: i + bw;
        for (size_t j = j1; j <= j2; j++) {
            Jt(j,i) = multijac->value(i,j);
        }
    }

    Jt.solve(b, lambda);
}

void Sim1D::resize()
{
    OneDim::resize();
    m_xnew.resize(size(), 0.0);
}

shared_ptr<Sim1D> newSim1D(vector<shared_ptr<Domain1D>>& domains)
{
    return make_shared<Sim1D>(domains);
}

}
