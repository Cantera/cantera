from exceptions import CanteraError

def setByName(a, options):
    """Set properties of phase 'a' by specifying keywords. Either the full
    property name or the short form may be given. The capitalization must
    be exactly as shown here.

    Note: all extensive property values are specified for a unit mass
    - i.e., the *specific* (not molar) property value should be
    specified.

    keyword       short form         property
    --------------------------------------------
    Pressure           P             pressure
    Density            Rho           density
    Temperature        T             temperature
    Volume             V             specific volume
    MoleFractions      X             mole fractions
    MassFracttions     Y             mass fractions
    Enthalpy           H             specific enthalpy
    IntEnergy          U             specific internal energy
    Entropy            S             specific entropy
    Vapor              Vap           vapor fraction in a two-phase mixture
    Liquid             Liq           liquid fraction in a two-phase mixture


    """

    tval = None
    pval = None
    hval = None
    uval = None
    sval = None
    vval = None
    qval = None

    np = 0
    nt = 0
    nv = 0
    nx = 0
    ny = 0
    ns = 0
    nh = 0
    nu = 0
    nq = 0

    for o in options.keys():
        val = options[o]
        if o == 'Temperature' or o == 'T':
            nt += 1
            tval = val
        elif o == 'Density' or o == 'Rho':
            nv += 1
            vval = 1.0/val
        elif o == 'Volume' or o == 'V':
            nv += 1
            vval = val
        elif o == 'MoleFractions' or o == 'X':
            nx += 1
            a.setMoleFractions(val)
        elif o == 'MassFractions' or o == 'Y':
            ny += 1
            a.setMassFractions(val)
        elif o == 'Pressure' or o == 'P':
            pval = val
            np += 1
        elif o == 'Enthalpy' or o == 'H':
            hval = val
            nh += 1
        elif o == 'IntEnergy' or o == 'U':
            uval = val
            nu += 1
        elif o == 'Entropy' or o == 'S':
            sval = val
            ns += 1
        elif o == 'Vapor' or o == 'Vap':
            nq += 1
            qval = val
        elif o == 'Liquid' or o == 'Liq':
            nq += 1
            qval = 1.0 - val

        else:
            raise CanteraError('unknown property: '+o)

    if nx + ny > 1:
        raise CanteraError('composition specified multiple times')

    nn = [nt, np, nv, ns, nh, nu, nq]
    for n in nn:
        if n > 1:
            raise CanteraError('property specified multiple times')

    ntot = nt + np + nv + ns + nh + nu + nq

    # set individual properties
    if ntot == 1:
        if nt == 1:
            a.setTemperature(tval)
        elif nv == 1:
            a.setDensity(1.0/vval)
        elif np == 1:
            a.setPressure(pval)
        else:
            props = options.keys()
            raise CanteraError('property '+props[0]+
                               ' can only be set in combination with '
                               +'another property')
    # set property pairs
    elif ntot == 2:
        if np == 1 and nh == 1:
            a.setState_HP(hval, pval)
        elif nu == 1 and nv == 1:
            a.setState_UV(uval, vval)
        elif ns == 1 and np == 1:
            a.setState_SP(sval, pval)
        elif ns == 1 and nv == 1:
            a.setState_SV(sval, vval)
        elif nt == 1 and np == 1:
            a.setState_TP(tval, pval)
        elif nt == 1 and nv == 1:
            a.setState_TR(tval, 1.0/vval)
        elif nt == 1 and nq == 1:
            a.setState_Tsat(tval, qval)
        elif np == 1 and nq == 1:
            a.setState_Psat(pval, qval)
        else:
            raise CanteraError('unimplemented property pair')


def set(a, **options):
    setByName(a, options)
