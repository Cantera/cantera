""" Python script to generate a Python extension module from a clib
header file.  """

import sys

_class = ''
_newclass = 1

def getargs(line):
    """Get the function name and arguments."""
    i1 = line.find('(')
    i2 = line.find(')')
    if (i1 < 0 or i2 < 0):
        raise 'syntax error: missing open or close quote'
    nm = line[:i1].split()
    nm = nm[-1]
    argline = line[i1+1:i2]
    args = argline.split(',')
    for n in range(len(args)): args[n] = args[n].split()
    v = []
    for a in args:
        if len(a) == 2: v.append(a)
    return nm, v

_itype = {'int':'i', 'double':'d', 'char*':'s', 'double*':'O', 'int*':'O'}


def isoutput(name):
    if len(name) >= 3 and name[-3:] == 'out':
        return 1
    else:
        return 0
    
def writepyfunc(rtype, name, args):
    """Write the Python extension module function."""
    
    print """
static PyObject *
py_"""+name+"""(PyObject *self, PyObject *args)
{
    """+rtype+""" _val;"""

    global _class, _newclass
    toks = name.split('_')    
    cls = toks[0]
    if len(toks) == 2:
        func = toks[1]
    else:
        func = toks[1] + toks[2]
        
    if cls != _class:
        _class = cls
        _newclass = 1
    else:
        _newclass = 0
        
    na = len(args)
    ain = []
    output = []
    if na > 0:
        vtype = []
        for a in args:
            # if the argument is an array, then the previous argument
            # must have been the array size. The Python argument list
            # will not include the size
            if a[0] == 'double*' or a[0] == 'int*':
                if not isoutput(a[1]):
                    vtype[-1] = 'PyObject*'
                    ain[-1] = a
                else:
                    output.append(a)
            elif a[0] == 'char*' and isoutput(a[1]):
                output.append(a)
                ain.pop()
            else:
                vtype.append(a[0])
                ain.append(a)
        for n in range(len(ain)):            
            print '   ',vtype[n],ain[n][1]+';'
        print '    if (!PyArg_ParseTuple(args,',
        s = '"'
        for a in ain:
                s += _itype[a[0]]
        s += ':'+name+'",'
        for a in ain:
                s += ' &'+a[1]+','
        s = s[:-1]+'))'
        print s,
        print """
        return NULL;
        """
    v = []
    for a in output:
        if a[0] == 'char*':
            print '    int '+a[1]+'_sz = 80;'
            print '    char* '+a[1]+' = new char['+a[1]+'_sz];'
            print
            
    for a in args:
            if a[0] == 'double*' or a[0] == 'int*':
                v[-1] = a[1]+'_len'
                v.append(a[1]+'_data')
                array = a[1]+'_array'
                print
                print '    PyArrayObject* '+array+' = (PyArrayObject*)'+a[1]+';'
                print '    '+a[0]+' '+a[1]+'_data = ('+a[0]+')'+array+'->data;'
                print '    int '+a[1]+'_len = '+array+'->dimensions[0];'
                print
            elif a[0] == 'char*' and isoutput(a[1]):
                v[-1] = a[1]+'_sz'
                v.append(a[1])
            else:
                v.append(a[1])
                
    s = '    _val = '+name+'('
    for a in v:
        s += a+','
    if s[-1] == ',': s = s[:-1]
    s += ');'
    print s,
    if (output):
        print '\n    PyObject* _ret = Py_BuildValue("'+_itype[output[0][0]]+'",'+output[0][1]+');'
        print '    delete '+output[0][1]+';'
        print '    if (int(_val) == -1) return reportCanteraError();'
        print """    return _ret;\n
}
"""
    else:
        print """
    if (int(_val) == -1) return reportCanteraError();
    """+'return Py_BuildValue("'+_itype[rtype]+'",_val);'+"""
}
"""
    return ain


def writepyclass(f, name, args):
    global _newclass
    if _newclass == 1:
        f.write("class "+_class.capitalize()+":\n")
        f.write("    def __init__(self):\n")
        f.write("        pass\n");
        _newclass = 0

    toks = name.split('_')    
    cls = toks[0]
    if len(toks) == 2:
        nm = toks[1]
    else:
        nm = toks[1] + toks[2]
        
    f.write('    def '+nm+'(self')
    for a in args[1:]:
        f.write(', '+a[1])
    f.write('):\n')
    f.write('        return _cantera.'+name+'(self._index')
    for a in args[1:]:
        f.write(', '+a[1])    
    f.write(')\n')
    
fname = sys.argv[1]
base, ext = fname.split('.')
mfile = 'py'+base+'_methods.h'
pfile = base+'.py'

_rtypes = ['int', 'double']

f = open(fname,'r')
fm = open(mfile,'w')
fp = open(pfile,'w')

lines = f.readlines()
f.close()

infunc = 0
funcline = ''
for line in lines:
    toks = line.split()
    if len(toks) > 0:
        if not infunc and toks[0] in _rtypes:
            infunc = 1
            funcline = line
        elif infunc:
            funcline += line
        last = toks[-1]
        if last[-1] == ';':
            infunc = 0
            name, args = getargs(funcline)
            toks = funcline.split()
            a = writepyfunc(toks[0], name, args)
            writepyclass(fp, name, a)
            fm.write('   {"'+name+'", py_'+name+', METH_VARARGS},\n')
            funcline = ''

fm.close()
fp.close()

        
