import sys, os
from tempfile import mktemp

def process(name):
    parts = name.split('.')
    base = parts[0]
    if len(parts) == 2:
        ext = parts[1]
    fname = mktemp('.py')
    fo = open(fname,'w')
    txt = """from Cantera.ctml_writer import *
import sys, os
f = sys.argv[1]
b = sys.argv[2]
try:
   os.remove(b+'.xml')
except:
   pass
execfile(f)
write()
"""
    fo.write(txt)
    fo.close()
    cmd = sys.executable+' '+fname+' '+name+' '+base
    os.system(cmd)
    os.remove(fname)



