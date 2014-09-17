
#include <cstdlib>

static PyObject*
py_domain_clear(PyObject* self, PyObject* args)
{
    int _val;
    _val = domain_clear();
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_del(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:domain_del", &i)) {
        return NULL;
    }

    _val = domain_del(i);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_type(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:domain_type", &i)) {
        return NULL;
    }

    _val = domain_type(i);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_index(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:domain_index", &i)) {
        return NULL;
    }

    _val = int(domain_index(i));
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_nComponents(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:domain_nComponents", &i)) {
        return NULL;
    }

    _val = int(domain_nComponents(i));
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_nPoints(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:domain_nPoints", &i)) {
        return NULL;
    }

    _val = int(domain_nPoints(i));
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_componentName(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int n;
    if (!PyArg_ParseTuple(args, "ii:domain_componentName", &i, &n)) {
        return NULL;
    }

    int nameout_sz = 80;
    char* nameout = new char[nameout_sz];

    _val = domain_componentName(i,n,nameout_sz,nameout);
    PyObject* _ret = Py_BuildValue("s",nameout);
    delete[] nameout;
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return _ret;

}


static PyObject*
py_domain_componentIndex(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    char* name;
    if (!PyArg_ParseTuple(args, "is:domain_componentIndex", &i, &name)) {
        return NULL;
    }

    _val = int(domain_componentIndex(i,name));
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_setBounds(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int n;
    double lower;
    double upper;
    if (!PyArg_ParseTuple(args, "iidd:domain_setBounds", &i, &n, &lower, &upper)) {
        return NULL;
    }

    _val = domain_setBounds(i,n,lower,upper);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_lowerBound(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int n;
    if (!PyArg_ParseTuple(args, "ii:domain_lowerBound", &i, &n)) {
        return NULL;
    }

    _val = domain_lowerBound(i,n);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_domain_upperBound(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int n;
    if (!PyArg_ParseTuple(args, "ii:domain_upperBound", &i, &n)) {
        return NULL;
    }

    _val = domain_upperBound(i,n);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_domain_setTolerances(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int n;
    double rtol, atol;
    int itime;
    if (!PyArg_ParseTuple(args, "iiddi:domain_setTolerances", &i, &n, &rtol, &atol, &itime)) {
        return NULL;
    }

    _val = domain_setTolerances(i,n, rtol, atol,itime);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_rtol(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int n;
    if (!PyArg_ParseTuple(args, "ii:domain_rtol", &i, &n)) {
        return NULL;
    }

    _val = domain_rtol(i,n);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_domain_atol(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int n;
    if (!PyArg_ParseTuple(args, "ii:domain_atol", &i, &n)) {
        return NULL;
    }

    _val = domain_atol(i,n);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_domain_setupGrid(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    PyObject* grid;
    if (!PyArg_ParseTuple(args, "iO:domain_setupGrid", &i, &grid)) {
        return NULL;
    }


    PyArrayObject* grid_array = (PyArrayObject*)
                                PyArray_ContiguousFromObject(grid, PyArray_DOUBLE, 1, 1);
    double* grid_data = (double*)(grid_array->data);
    size_t grid_len = grid_array->dimensions[0];

    _val = domain_setupGrid(i,grid_len,grid_data);
    Py_DECREF(grid_array);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_setID(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    char* id;
    if (!PyArg_ParseTuple(args, "is:domain_setID", &i, &id)) {
        return NULL;
    }

    _val = domain_setID(i,id);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_setDesc(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    char* desc;
    if (!PyArg_ParseTuple(args, "is:domain_setDesc", &i, &desc)) {
        return NULL;
    }

    _val = domain_setDesc(i,desc);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_inlet_setSpreadRate(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double v;
    if (!PyArg_ParseTuple(args, "id:inlet_setSpreadRate", &i, &v)) {
        return NULL;
    }

    _val = inlet_setSpreadRate(i,v);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_domain_grid(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int n;
    if (!PyArg_ParseTuple(args, "ii:domain_grid", &i, &n)) {
        return NULL;
    }

    _val = domain_grid(i,n);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_bdry_setMdot(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double mdot;
    if (!PyArg_ParseTuple(args, "id:bdry_setMdot", &i, &mdot)) {
        return NULL;
    }

    _val = bdry_setMdot(i,mdot);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_bdry_setTemperature(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double t;
    if (!PyArg_ParseTuple(args, "id:bdry_setTemperature", &i, &t)) {
        return NULL;
    }

    _val = bdry_setTemperature(i,t);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_bdry_setMoleFractions(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    char* x;
    if (!PyArg_ParseTuple(args, "is:bdry_setMoleFractions", &i, &x)) {
        return NULL;
    }

    _val = bdry_setMoleFractions(i,x);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_bdry_temperature(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:bdry_temperature", &i)) {
        return NULL;
    }

    _val = bdry_temperature(i);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_bdry_massFraction(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int k;
    if (!PyArg_ParseTuple(args, "ii:bdry_massFraction", &i, &k)) {
        return NULL;
    }

    _val = bdry_massFraction(i,k);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_bdry_mdot(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:bdry_mdot", &i)) {
        return NULL;
    }

    _val = bdry_mdot(i);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_reactingsurf_setkineticsmgr(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int j;
    if (!PyArg_ParseTuple(args, "ii:reactingsurf_setkineticsmgr", &i, &j)) {
        return NULL;
    }

    _val = reactingsurf_setkineticsmgr(i,j);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_reactingsurf_enableCoverageEqs(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int onoff;
    if (!PyArg_ParseTuple(args, "ii:reactingsurf_enableCoverageEqs", &i, &onoff)) {
        return NULL;
    }

    _val = reactingsurf_enableCoverageEqs(i,onoff);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_inlet_new(PyObject* self, PyObject* args)
{
    int _val;
    _val = inlet_new();
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_outlet_new(PyObject* self, PyObject* args)
{
    int _val;
    _val = outlet_new();
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_outletres_new(PyObject* self, PyObject* args)
{
    int _val;
    _val = outletres_new();
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_symm_new(PyObject* self, PyObject* args)
{
    int _val;
    _val = symm_new();
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_surf_new(PyObject* self, PyObject* args)
{
    int _val;
    _val = surf_new();
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_reactingsurf_new(PyObject* self, PyObject* args)
{
    int _val;
    _val = reactingsurf_new();
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_stflow_new(PyObject* self, PyObject* args)
{
    int _val;
    int iph;
    int ikin;
    int itr;
    int itype;
    if (!PyArg_ParseTuple(args, "iiii:stflow_new", &iph, &ikin, &itr, &itype)) {
        return NULL;
    }

    _val = stflow_new(iph,ikin,itr,itype);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_stflow_setTransport(PyObject* self, PyObject* args)
{
    int _val;
    int i, itr, isoret;
    if (!PyArg_ParseTuple(args, "iii:stflow_setTransport", &i, &itr, &isoret)) {
        return NULL;
    }

    _val = stflow_setTransport(i,itr, isoret);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_stflow_enableSoret(PyObject* self, PyObject* args)
{
    int _val;
    int i, isoret;
    if (!PyArg_ParseTuple(args, "ii:stflow_enableSoret", &i, &isoret)) {
        return NULL;
    }

    _val = stflow_enableSoret(i,isoret);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_stflow_setPressure(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double p;
    if (!PyArg_ParseTuple(args, "id:stflow_setPressure", &i, &p)) {
        return NULL;
    }

    _val = stflow_setPressure(i,p);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_stflow_pressure(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int n;
    if (!PyArg_ParseTuple(args, "i:stflow_pressure", &i)) {
        return NULL;
    }

    _val = stflow_pressure(i);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}

static PyObject*
py_stflow_setFixedTempProfile(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    PyObject* pos;
    PyObject* temp;
    if (!PyArg_ParseTuple(args, "iOO:stflow_setFixedTempProfile", &i, &pos, &temp)) {
        return NULL;
    }


    PyArrayObject* pos_array = (PyArrayObject*)
                               PyArray_ContiguousFromObject(pos, PyArray_DOUBLE, 1, 1);
    double* pos_data = (double*)(pos_array->data);
    size_t pos_len = pos_array->dimensions[0];


    PyArrayObject* temp_array = (PyArrayObject*)
                                PyArray_ContiguousFromObject(temp, PyArray_DOUBLE, 1, 1);


    double* temp_data = (double*)(temp_array->data);
    size_t temp_len = temp_array->dimensions[0];

    _val = stflow_setFixedTempProfile(i,pos_len,pos_data,temp_len,temp_data);
    Py_DECREF(pos_array);
    Py_DECREF(temp_array);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_stflow_solveSpeciesEqs(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int flag;
    if (!PyArg_ParseTuple(args, "ii:stflow_solveSpeciesEqs", &i, &flag)) {
        return NULL;
    }

    _val = stflow_solveSpeciesEqs(i,flag);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_stflow_solveEnergyEqn(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int flag;
    if (!PyArg_ParseTuple(args, "ii:stflow_solveEnergyEqn", &i, &flag)) {
        return NULL;
    }

    _val = stflow_solveEnergyEqn(i,flag);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_clear(PyObject* self, PyObject* args)
{
    int _val;
    _val = sim1D_clear();
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_new(PyObject* self, PyObject* args)
{
    int _val;
    PyObject* domains;
    if (!PyArg_ParseTuple(args, "O:sim1D_new", &domains)) {
        return NULL;
    }

    // Assuming this is screwed up too
    //PyArrayObject* domains_array = (PyArrayObject*)
    // PyArray_ContiguousFromObject(domains, PyArray_INT, 1, 1);
    // int* domains_data = (int*)(domains_array->data);
    PyArrayObject* domains_array = (PyArrayObject*)
                                   PyArray_ContiguousFromObject(domains, PyArray_DOUBLE, 1, 1);
    void* nTMPv = (void*)(domains_array->data);
    double* dd_data = (double*) nTMPv;
    size_t domains_len = domains_array->dimensions[0];

    int* domains_data = (int*) malloc(sizeof(int) * domains_len);
    for (size_t i = 0; i <  domains_len; i++) {
        domains_data[i] = (int) dd_data[i];
    }

    _val = sim1D_new(domains_len, domains_data);

    free(domains_data);
    Py_DECREF(domains_array);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_del(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:sim1D_del", &i)) {
        return NULL;
    }

    _val = sim1D_del(i);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_setValue(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int dom;
    int comp;
    int localPoint;
    double value;
    if (!PyArg_ParseTuple(args, "iiiid:sim1D_setValue", &i, &dom, &comp, &localPoint, &value)) {
        return NULL;
    }

    _val = sim1D_setValue(i,dom,comp,localPoint,value);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_setProfile(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int dom;
    int comp;
    PyObject* pos;
    PyObject* v;
    if (!PyArg_ParseTuple(args, "iiiOO:sim1D_setProfile", &i, &dom, &comp, &pos, &v)) {
        return NULL;
    }


    PyArrayObject* pos_array = (PyArrayObject*)
                               PyArray_ContiguousFromObject(pos, PyArray_DOUBLE, 1, 1);
    double* pos_data = (double*)(pos_array->data);
    size_t pos_len = pos_array->dimensions[0];


    PyArrayObject* v_array = (PyArrayObject*)
                             PyArray_ContiguousFromObject(v, PyArray_DOUBLE, 1, 1);
    double* v_data = (double*)(v_array->data);
    size_t v_len = v_array->dimensions[0];

    _val = sim1D_setProfile(i,dom,comp,pos_len,pos_data,v_len,v_data);
    Py_DECREF(pos_array);
    Py_DECREF(v_array);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_setFlatProfile(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int dom;
    int comp;
    double v;
    if (!PyArg_ParseTuple(args, "iiid:sim1D_setFlatProfile", &i, &dom, &comp, &v)) {
        return NULL;
    }

    _val = sim1D_setFlatProfile(i,dom,comp,v);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_showSolution(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    char* fname;
    if (!PyArg_ParseTuple(args, "is:sim1D_showSolution", &i, &fname)) {
        return NULL;
    }

    _val = sim1D_showSolution(i,fname);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_setTimeStep(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double stepsize;
    PyObject* nsteps;
    if (!PyArg_ParseTuple(args, "idO:sim1D_setTimeStep", &i, &stepsize, &nsteps)) {

        return NULL;
    }

    // Found that the PyArray_INT had to be replaced by PyArray_DOUBLE
    // in the call below. This needs exploring.
    // PyArrayObject* nsteps_array = (PyArrayObject*)
    //  PyArray_ContiguousFromObject(nsteps, PyArray_INT, 1, 1);
    PyArrayObject* nsteps_array = (PyArrayObject*)
                                  PyArray_ContiguousFromObject(nsteps, PyArray_DOUBLE, 1, 1);

    void* nTMPv = (void*)(nsteps_array->data);
    double* nsteps_data = (double*) nTMPv;
    size_t nsteps_len = nsteps_array->dimensions[0];

    int* nsteps_datai = (int*) malloc(sizeof(int) * nsteps_len);
    for (size_t i = 0; i <  nsteps_len; i++) {
        nsteps_datai[i] = (int) nsteps_data[i];
    }
    _val = sim1D_setTimeStep(i, stepsize, nsteps_len, nsteps_datai);
    free(nsteps_datai);

    Py_DECREF(nsteps_array);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_solve(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int loglevel;
    int refine_grid;
    if (!PyArg_ParseTuple(args, "iii:sim1D_solve", &i, &loglevel, &refine_grid)) {
        return NULL;
    }

    _val = sim1D_solve(i,loglevel,refine_grid);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_refine(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int loglevel;
    if (!PyArg_ParseTuple(args, "ii:sim1D_refine", &i, &loglevel)) {
        return NULL;
    }

    _val = sim1D_refine(i,loglevel);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_setRefineCriteria(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int dom;
    double ratio;
    double slope;
    double curve;
    double prune;
    if (!PyArg_ParseTuple(args, "iidddd:sim1D_setRefineCriteria", &i, &dom, &ratio, &slope, &curve, &prune)) {
        return NULL;
    }

    _val = sim1D_setRefineCriteria(i,dom,ratio,slope,curve,prune);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_sim1D_setGridMin(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int dom;
    double gridmin;
    if (!PyArg_ParseTuple(args, "iid:sim1D_setGridMin", &i, &dom, &gridmin)) {
        return NULL;
    }

    _val = sim1D_setGridMin(i,dom,gridmin);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_sim1D_getInitialSoln(PyObject* self, PyObject* args)
{
    int i, iok;
    if (!PyArg_ParseTuple(args, "i:sim1D_getInitialSoln", &i)) {
        return NULL;
    }

    iok = sim1D_getInitialSoln(i);
    if (iok == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",iok);
}

static PyObject*
py_sim1D_save(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    char* fname;
    char* id;
    char* desc;
    if (!PyArg_ParseTuple(args, "isss:sim1D_save", &i, &fname, &id, &desc)) {
        return NULL;
    }

    _val = sim1D_save(i,fname,id,desc);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_restore(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    char* fname;
    char* id;
    if (!PyArg_ParseTuple(args, "iss:sim1D_restore", &i, &fname, &id)) {
        return NULL;
    }

    _val = sim1D_restore(i,fname,id);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_writeStats(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int printTime;
    if (!PyArg_ParseTuple(args, "ii:sim1D_writeStats", &i, &printTime)) {
        return NULL;
    }
    _val = sim1D_writeStats(i, printTime);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_domainIndex(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    char* name;
    if (!PyArg_ParseTuple(args, "is:sim1D_domainIndex", &i, &name)) {
        return NULL;
    }

    _val = sim1D_domainIndex(i,name);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_value(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int idom;
    int icomp;
    int localPoint;
    if (!PyArg_ParseTuple(args, "iiii:sim1D_value", &i, &idom, &icomp,
                          &localPoint)) {
        return NULL;
    }

    _val = sim1D_value(i,idom,icomp,localPoint);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_sim1D_workValue(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int idom;
    int icomp;
    int localPoint;
    if (!PyArg_ParseTuple(args, "iiii:sim1D_workValue", &i, &idom, &icomp, &localPoint)) {
        return NULL;
    }

    _val = sim1D_workValue(i,idom,icomp,localPoint);
    if (_val == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_sim1D_eval(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double rdt;
    int count;
    if (!PyArg_ParseTuple(args, "idi:sim1D_eval", &i, &rdt, &count)) {
        return NULL;
    }

    _val = sim1D_eval(i,rdt,count);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_setMaxJacAge(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int ss_age;
    int ts_age;
    if (!PyArg_ParseTuple(args, "iii:sim1D_setMaxJacAge", &i, &ss_age, &ts_age)) {
        return NULL;
    }

    _val = sim1D_setMaxJacAge(i,ss_age,ts_age);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_timeStepFactor(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double tfactor;
    if (!PyArg_ParseTuple(args, "id:sim1D_timeStepFactor", &i, &tfactor)) {
        return NULL;
    }

    _val = sim1D_timeStepFactor(i,tfactor);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_setTimeStepLimits(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double tsmin;
    double tsmax;
    if (!PyArg_ParseTuple(args, "idd:sim1D_setTimeStepLimits", &i, &tsmin, &tsmax)) {
        return NULL;
    }

    _val = sim1D_setTimeStepLimits(i,tsmin,tsmax);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_sim1D_setFixedTemperature(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double temp;
    if (!PyArg_ParseTuple(args, "id:sim1D_setFixedTemperature", &i, &temp)) {
        return NULL;
    }
    _val = sim1D_setFixedTemperature(i,temp);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

