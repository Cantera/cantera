import sys

def getargs(line):
    i1 = line.find('(')
    i2 = line.find(')')
    if (i1 < 0 or i2 < 0):
        raise 'syntax error: missing open or close quote'
    nm = line[:i1].split()
    nm = nm[-1]
    argline = line[i1+1:i2]
    args = argline.split(',')
    for n in range(len(args)): args[n] = args[n].split()
    return nm, args

_itype = {'int':'i', 'double':'d', 'char*':'s', 'double*':'O', 'int*':'O'}

def writepyfunc(name, args):
    print """
static PyObject *
py_"""+name+"""(PyObject *self, PyObject *args)
{
    int _iok;"""
    for a in args:
        if len(a) == 2:
            print '   ',a[0],a[1]+';'
    print '    if (!PyArg_ParseTuple(args,',
    s = '"'
    for a in args:
        if len(a) == 2:
            s += _itype[a[0]]
    s += ':'+name+'",'
    for a in args:
        if len(a) == 2:
            s += ' &'+a[1]+','
    s = s[:-1]+'))'
    print s,
    print """
        return NULL;
        """
    s = '    _iok = '+name+'('
    for a in args:
        if len(a) == 2:
            s += a[1]+','
    s = s[:-1]+')'
    print s,
    print """
    if (_iok == -1) return reportCanteraError();
    return Py_BuildValue("i",_iok);
}
    """

    
fname = sys.argv[1]

_rtypes = ['int', 'double']

f = open(fname,'r')
lines = f.readlines()

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
            writepyfunc(name, args)
            funcline = ''


        
