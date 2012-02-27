import os

def writepydoc(pkg, out):
    nm = pkg.__name__
    out.write('pydoc -w '+nm+'\n')
    p = pkg.__path__[0]
    files = os.listdir(p)
    for f in files:
        n = f.split('.')
        if len(n) >= 2 and n[-1] == 'py' and n[0][0] <> '_':
            try:
                exec('import '+nm+'.'+n[0])
                out.write('pydoc -w '+nm+'.'+n[0]+'\n')
            except:
                pass


out = open('mkdoc','w')
import Cantera
writepydoc(Cantera, out)
import Cantera.OneD
writepydoc(Cantera.OneD, out)

out.close()
