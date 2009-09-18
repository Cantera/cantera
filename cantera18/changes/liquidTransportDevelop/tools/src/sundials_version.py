#
# find the version of the installed sundials package
#
import string
import sys
args = sys.argv

def splitversion(s):
    toks = s.split('.')
    if len(toks) <> 3:
        return (0,0,0)
    return (string.atoi(toks[0]), string.atoi(toks[1]), string.atoi(toks[2]))

sundials_home = args[1]
vinst = (0,0,0)

try:
    readme = open(sundials_home+'/README','r')
except:
    print "can't open "+sundials_home+"/README"
    sys.exit(-1)

try:
    lines = readme.readlines()
    for line in lines:
        toks = line.split()
        if 'Release' in toks:
            n = toks.index('Release')
            version = toks[n+1]
            if version[-1] == ',':
                version = version[:-1]
                vinst = splitversion(version)
                break
except:
    vinst = (-1,-1,-1)

fout = open('sundials_includes.h','w')
fout.write('#define SUNDIALS_MAJOR_VERSION = '+`vinst[0]`+'\n#define SUNDIALS_MINOR_VERSION = '+`vinst[1]`+'\n#define SUNDIALS_REVISION = '+`vinst[2]`)

fout.write('\n#define SUNDIALS_VERSION_'+`vinst[0]`+'\n')

fout.close()



