f = open('modules','r')
lines = f.readlines()
gen = {}
mods = []
for m in lines:
    file, prefix = m.split()
    fm = open(file,'r')
    mlines = fm.readlines()
    for line in mlines:
        toks = line.split()
        if len(toks) == 2:
            if toks[0] == 'module':
                mods.append(toks[1])

        le = line.find('end')
        lf = line.find('function')
        ls = line.find('subroutine')
        loc = line.find(prefix+'_')
        if le < 0 and (lf > 0 or ls > 0):
            if loc > 0:
                sline = line[loc:]
                n = len(prefix)
                p = sline.find('(')
                nm = sline[n+1:p]
                if nm <> '':
                    if gen.has_key(nm):
                        gen[nm].append(prefix+'_'+nm)
                    else:
                        gen[nm] = [prefix+'_'+nm]

fout = open('canteramod.f90','w')
fout.write('MODULE CANTERA\n\n')
for m in mods:
    fout.write('  USE '+m+'\n')

funcs = gen.keys()
funcs.sort()
for fn in funcs:
    fout.write('\n  INTERFACE '+fn+'\n')
    for cf in gen[fn]:
        fout.write('     MODULE PROCEDURE '+cf+'\n')
    fout.write('  END INTERFACE '+fn+'\n')

fout.write('\nEND MODULE CANTERA\n\n')
