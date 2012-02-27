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
        if len(a) <= 1:
            pass
        elif len(a) == 2:
            v.append(a)
        elif len(a) == 3 and a[0] == 'const':
            v.append([a[0]+' '+a[1],a[2]])
        else:
            print 'a = ',a, args
            print 'line = ',line
            raise 'malformed argument: '
    return nm, v

_c2fout = {'int*':'integer', 'integer*':'integer',
           'double*':'double precision', 'doublereal*':'double precision',
           'char*':'character*(*)'}

_c2fin  = {'const int*':'integer', 'const integer*':'integer',
           'const double*':'double precision',
           'const doublereal*':'double precision',
           'const char*':'character*(*)'}

_c2fret = {'int':'integer', 'integer':'integer',
           'double':'double precision',
           'status_t':'integer',
           'doublereal':'double precision'}


def writeinterface(fint, rtype, name, args):
    s = '    '+_c2fret[rtype] + ' function '
    if name[-1] == '_':
        name = name[:-1]
    s += name + '('
    argstr=''
    for a in args:
        if a[0] <> 'ftnlen':
            argstr += a[1] + ', '
    if len(argstr) > 0:
        argstr = argstr[:-2]
    s += argstr+')\n'
    fint.write(s)
    for a in args:
        if a[0][0] == 'd':
            arr = '(*)'
        else:
            arr = ''
        if a[0] in _c2fin:
            fint.write('        '+_c2fin[a[0]]+', intent(in) :: '+a[1]+arr+'\n')
        elif a[0] in _c2fout:
            fint.write('        '+_c2fout[a[0]]+', intent(out) :: '+a[1]+arr+'\n')
    fint.write('    end function '+name+'\n\n')


def writef90(fmod, rtype, otype, hndl, name, args):
    subroutine = 0
    if rtype == 'status_t':
        subroutine = 1
    if subroutine:
        s = '    subroutine '
    else:
        s = '    '+_c2fret[rtype] + ' function '
    if name[-1] == '_':
        name = name[:-1]
    wname = 'ct'+name[1:]

    s += wname + '('
    argstr='self, '
    for a in args[1:]:
        if a[0] <> 'ftnlen':
            argstr += a[1] + ', '
    if len(argstr) > 0:
        argstr = argstr[:-2]
    s += argstr+')\n'
    fmod.write(s)
    fmod.write("""      implicit none
      type("""+otype+'), intent(inout) :: self\n')
    for a in args[1:]:
        if a[0][0] == 'd':
            arr = '(*)'
        else:
            arr = ''
        if a[0] in _c2fin:
            fmod.write('      '+_c2fin[a[0]]+', intent(in) :: '+a[1]+arr+'\n')
        elif a[0] in _c2fout:
            fmod.write('      '+_c2fout[a[0]]+', intent(out) :: '+a[1]+arr+'\n')
    if subroutine:
        s = '      self%err = '+name+'(self%'+hndl+', '
    else:
        s = '      '+wname+' = '+name+'(self%'+hndl+', '
    argstr = ''
    for a in args[1:]:
        if a[0] <> 'ftnlen':
            argstr += a[1] + ', '
    argstr = argstr[:-2]
    s += argstr+')'
    if subroutine:
        fmod.write(s+'\n    end subroutine '+wname+'\n\n')
    else:
        fmod.write(s+'\n    end function '+wname+'\n\n')



fname = sys.argv[1]    # fctxml
otype = sys.argv[2]    # XML_Node
hndl = sys.argv[3]     # xml_id

base, ext = fname.split('.')
modname = base

f = open(fname,'r')
fint = open(base+'_interface.f90','w')
#fmod = open(modname+'.f90','w')
fmod = open('cantera_'+modname+'.f90','w')

lines = f.readlines()
f.close()

#_rtypes = ['int', 'double', 'integer']

infunc = 0
funcline = ''
extern = 0

fint.write('module '+base+'\n')

fint.write('interface\n')
for line in lines:
    toks = line.split()
    if len(toks) > 0 and toks[0][:2] <> '//':
        if toks[0] == 'extern':
            extern = 1
        if extern:
            if not infunc:
                if line.find('CANTERA_CAPI') > 0:
                    infunc = 1
                    funcline = line
            elif infunc:
                funcline += line
            last = toks[-1]
            if infunc and last[-1] == '{':
                infunc = 0
                name, args = getargs(funcline)
                toks = funcline.split()
                a = writeinterface(fint, toks[0], name, args)
                funcline = ''
fint.write('end interface\n')
fint.write('end module '+base+'\n')
fint.close()

fmod.write('module '+modname+'\n')
fmod.write('  use '+base+"""

  type """+otype+"""
    integer :: """+hndl+"""
  end type """+otype+"""

contains

""")

for line in lines:
    toks = line.split()
    if len(toks) > 0 and toks[0][:2] <> '//':
        if toks[0] == 'extern':
            extern = 1
        if extern:
            if not infunc:
                if line.find('CANTERA_CAPI') > 0:
                    infunc = 1
                    funcline = line
            elif infunc:
                funcline += line
            last = toks[-1]
            if infunc and last[-1] == '{':
                infunc = 0
                name, args = getargs(funcline)
                toks = funcline.split()
                a = writef90(fmod, toks[0], otype, hndl, name, args)
                funcline = ''

fmod.write('end module '+base+'\n')
fmod.close()
