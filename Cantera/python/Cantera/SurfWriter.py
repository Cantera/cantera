"""
Write C functions implementing a surface reaction mechanism.
NOT CURRENTLY FUNCTIONAL
"""

from Cantera import CanteraError
from Cantera import units
from Numeric import array
from Cantera import constants
import math
import types

hdr = """
/*

The identity of the species in each phase is as listed below.
"""

ropheader = """

  /*
  Get the reaction rates of progress given the concentrations.

  conc --- concentrations in kmol/m^3 or kmol/m^2. Array c
           must be dimensioned at least KB1 + KB2 + KS. The first
           KB1 entries are the concentrations of species in bulk phase 1,
           the next KB2 entries are the concentrations of species in bulk
           phase 2, and the final KS entries are the concentrations of the
           surface species.

  ropf --- rates of progress of the surface reactions in kmol/m^2/s.
           Must be dimensioned at least NSR, the number of surface reactions.

  */
  
  void get_rop(double* c, double* kf, double* ropf) {
"""

sdotheader = """

  void get_sdot(double* r, double* sdot) {
"""

rateheader = """

  /* get the forward reaction rates */
  void get_kf(double T, double* kf) {
    double logT = log(T);
    double rt = 1.0/T;
"""

finaltxt = """
 }
    """

class SurfWriter:

    def __init__(self, iface, mechname):
        self.bulk1 = iface.p1
        self.bulk2 = iface.p2
        self.surf = iface
        self.sp = {}        
        self.f = open(mechname+'.c', 'w')
        self.write_top()
        self.f.write('\n#ifdef __cplusplus\nextern "C" {\n#endif')
        self.rop = ropheader
        self.rate = rateheader
        self.nrxns = 0
        self.sdot = {}

    def write(self):
        self.f.write(self.rop)
        self.f.write('  }\n\n')
        self.f.write(sdotheader)
        for k in self.sdot.keys():
            self.f.write('\n    /*  '+self.sp[k]+'  */\n')
            self.f.write('    sdot['+`k`+'] = '+self.sdot[k]+';\n')
        self.f.write('  }\n\n')
        self.f.write(self.rate)
        self.f.write('  }\n')
        self.f.write('#ifdef __cplusplus\n}\n#endif\n')        
        self.f.close()

    def write_top(self):
        self.f.write(hdr)
        n1 = self.bulk1.nSpecies()
        self.f.write('Bulk phase 1 species: '+`n1`+' total.\n')

        sp = self.bulk1.speciesNames()
        i = 0
        k = 0
        for s in sp:
            self.f.write('%3d  %s\n' % (i, s))
            i += 1
            self.sp[k] = s
            k += 1
        self.f.write('\n\n')
        
        if self.bulk2:
            n2 = self.bulk2.nSpecies()
            self.f.write('Bulk phase 2 species: '+`n2`+' total.\n')
            
            sp = self.bulk2.speciesNames()
            i = 0
            for s in sp:
                self.f.write('%3d  %s\n' % (i, s))
                i += 1
                self.sp[k] = s
                k += 1                
            self.f.write('\n\n')
            
        else:
            self.f.write('No second bulk phase.\n\n')

        ns = self.surf.nSpecies()
        self.f.write('Surface species: '+`ns`+' total.\n')
        sp = self.surf.speciesNames()
        i = 0
        for s in sp:
            self.f.write('%3d  %s\n' % (i, s))
            i += 1
            self.sp[k] = s
            k += 1            
        self.f.write('\n\n*/')
        
        self.f.write('\n\n')
        
        
    def write_update_rate(self, rate):
        i = self.nrxns
        self.rate += '    kf['+`i`+'] = '+'%17.10e' % (rate[0],)
        if rate[1] == 0.0 and rate[2] == 0.0:
            self.rate += ';\n'            
            return
        self.rate += ' * exp('
        if rate[1] <> 0.0:
            self.rate += '%f * logT' % (rate[1],)
        if rate[2] <> 0.0:
            self.rate += '- %f * rt' % (rate[2],)
        self.rate += ');\n'

            
    def write_ROP(self, rindex, rstoich, rorder,
                    pindex, pstoich, rate):
        f = ''
        i = self.nrxns
        f += '    ropf['+`i`+'] = kf['+`i`+']*'
        nr = len(rindex)
        for n in range(nr):
            k = rindex[n]
            if rorder[n] == rstoich[n]:
                for j in range(rstoich[n]):
                    f +='c['+`k`+']*'
            else:
                f += 'pow(c['+`k`+'], '+`rorder[n]`+')*'
        f = f[:-1]+';\n'
        self.rop += f

    def write_sdot(self, rindex, rstoich, pindex, pstoich):
        i = self.nrxns
        nr = len(rindex)
        for n in range(nr):
            k = rindex[n]
            if rstoich[n] == 1: st = ''
            else: st = `rstoich[n]`+'*'
            if not self.sdot.has_key(k):
                self.sdot[k] = ' -'+st+'r['+`i`+']'
            else:
                self.sdot[k] += ' - '+st+'r['+`i`+']'
                
        np = len(pindex)
        for n in range(np):
            k = pindex[n]
            if pstoich[n] == 1: st = ''
            else: st = `pstoich[n]`+'*'            
            if not self.sdot.has_key(k):
                self.sdot[k] = st+'r['+`i`+']'
            else:
                self.sdot[k] += ' + '+st+'r['+`i`+']'            
        
