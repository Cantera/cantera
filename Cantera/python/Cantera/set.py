from exceptions import CanteraError

def set(a, **options):

    pval = None
    hval = None
    uval = None
    sval = None
    vval = None
    np = 0
    
    for o in options.keys():
        val = options[o]
        if o == 'Temperature' or o == 'T':
            a.setTemperature(val)
        elif o == 'Density' or o == 'Rho':
            a.setDensity(val)
            vval = 1.0/val
        elif o == 'V':
            a.setDensity(1.0/val)
            vval = val
        elif o == 'MoleFractions' or o == 'X':
            a.setMoleFractions(val)
        elif o == 'MassFractions' or o == 'Y':
            a.setMassFractions(val)            
        elif o == 'Pressure' or o == 'P':
            pval = val
            np = np + 1
        elif o == 'Enthalpy' or o == 'H':
            hval = val
            np = np + 1
        elif o == 'IntEnergy' or o == 'U':
            uval = val
            np = np + 1             
        elif o == 'Entropy' or o == 'S':
            sval = val
            np = np + 1
        else:
            raise CanteraError('unknown property: '+o)

    if np == 1:
        if pval:
            a.setPressure(pval)

    if np >= 2: 
        if pval and hval:
            a.setState_HP(hval,pval)
        elif uval and vval:
            a.setState_UV(uval,vval)
        elif sval and pval:
            a.setState_SP(sval,pval)   
        elif sval and vval:
            a.setState_SV(sval,vval)      
        else:
            raise CanteraError('unimplemented property pair')
   


