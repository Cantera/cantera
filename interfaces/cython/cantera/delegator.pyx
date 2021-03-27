cdef void callback_v(PyFuncInfo& funcInfo):
    try:
        (<object>funcInfo.func())()
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)


cdef void callback_v_d(PyFuncInfo& funcInfo, double arg):
    try:
        (<object>funcInfo.func())(arg)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)


cdef void callback_v_b(PyFuncInfo& funcInfo, cbool arg):
    try:
        (<object>funcInfo.func())(arg)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)


cdef void callback_v_dp(PyFuncInfo& funcInfo, size_array1 sizes, double* arg):
    cdef double[:] view = <double[:sizes[0]]>arg if sizes[0] else None

    try:
        (<object>funcInfo.func())(view)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)


cdef void callback_v_d_dp(PyFuncInfo& funcInfo, size_array1 sizes, double arg1,
                                double* arg2):
    cdef double[:] view = <double[:sizes[0]]>arg2 if sizes[0] else None

    try:
        (<object>funcInfo.func())(arg1, view)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)


cdef int callback_i_dr_d_dp(PyFuncInfo& funcInfo, double& out,
                            size_array1 sizes, double arg1, double* arg2):
    cdef double[:] view = <double[:sizes[0]]>arg2 if sizes[0] else None
    try:
        ret = (<object>funcInfo.func())(arg1, view)
        if ret is None:
            return 0
        else:
            (&out)[0] = ret
            return 1
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)
    return -1


cdef int callback_i_sr_z(PyFuncInfo& funcInfo, string& out, size_t arg):
    try:
        ret = (<object>funcInfo.func())(arg)
        if ret is None:
            return 0
        else:
            (&out)[0] = stringify(ret)
            return 1
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)
    return -1


cdef int callback_i_zr_csr(PyFuncInfo& funcInfo, size_t& out, const string& arg):
    try:
        ret = (<object>funcInfo.func())(pystr(arg))
        if ret is None:
            return 0
        else:
            (&out)[0] = ret
            return 1
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)
    return -1


cdef void assign_delegates(obj, CxxDelegator* delegator):
    # Find all delegate methods, and make sure there aren't multiple
    # conflicting implementations
    cdef string cxx_name
    cdef string cxx_when
    for name in obj.delegatable_methods:
        when = None
        replace = 'replace_{}'.format(name)
        if hasattr(obj, replace):
            when = 'replace'
            method = getattr(obj, replace)

        before = 'before_{}'.format(name)
        if hasattr(obj, before):
            if when is not None:
                raise CanteraError(
                    "Only one delegate supported for '{}'".format(name))
            when = 'before'
            method = getattr(obj, before)

        after = 'after_{}'.format(name)
        if hasattr(obj, after):
            if when is not None:
                raise CanteraError(
                    "Only one delegate supported for '{}'".format(name))
            when = 'after'
            method = getattr(obj, after)

        if when is None:
            continue

        cxx_name = stringify(obj.delegatable_methods[name][0])
        callback = obj.delegatable_methods[name][1].replace(' ', '')
        cxx_when = stringify(when)
        if callback == 'void()':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v), cxx_when)
        if callback == 'void(double)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_d), cxx_when)
        if callback == 'void(double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_dp), cxx_when)
        if callback == 'void(bool)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_b), cxx_when)
        if callback == 'void(double,double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_d_dp), cxx_when)
        if callback == 'double(double,double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_i_dr_d_dp), cxx_when)
        if callback == 'string(size_t)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_i_sr_z), cxx_when)
        if callback == 'size_t(string)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_i_zr_csr), cxx_when)
