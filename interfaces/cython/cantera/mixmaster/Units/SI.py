#!/usr/bin/env python
from .unit import unit, dimensionless

#
# The basic SI units
#
meter = unit(1.0, (1, 0, 0, 0, 0, 0, 0))
kilogram = unit(1.0, (0, 1, 0, 0, 0, 0, 0))
second = unit(1.0, (0, 0, 1, 0, 0, 0, 0))
ampere = unit(1.0, (0, 0, 0, 1, 0, 0, 0))
kelvin = unit(1.0, (0, 0, 0, 0, 1, 0, 0))
mole = unit(1.0, (0, 0, 0, 0, 0, 1, 0))
candela = unit(1.0, (0, 0, 0, 0, 0, 0, 1))

#
# The 21 derived SI units with special names
#
radian = dimensionless                 #  plane angle
steradian = dimensionless              #  solid angle

hertz = 1/second                       #  frequency

newton = meter*kilogram/second**2      #  force
pascal = newton/meter**2               #  pressure
joule = newton*meter                   #  work, heat
watt = joule/second                    #  power, radiant flux

coulomb = ampere*second                #  electric charge
volt = watt/ampere                     #  electric potential difference
farad = coulomb/volt                   #  capacitance
ohm = volt/ampere                      #  electric resistance
siemens = ampere/volt                  #  electric conductance
weber = volt*second                    #  magnetic flux
tesla = weber/meter**2                 #  magnetic flux density
henry = weber/ampere                   #  inductance

celsius = kelvin                       #  Celsius temperature

lumen = candela*steradian              #  luminous flux
lux = lumen/meter**2                   #  illuminance

becquerel = 1/second                   #  radioactivity
gray = joule/kilogram                  #  absorbed dose
sievert = joule/kilogram               #  dose equivalent

#
# The prefixes
#

yotta = 1e24
zetta = 1e21
exa = 1e18
peta = 1e15
tera = 1e12
giga = 1e9
mega = 1e6
kilo = 1000
hecto = 100
deka = 10
deci = .1
centi = .01
milli = .001
micro = 1e-6
nano = 1e-9
pico = 1e-12
femto = 1e-15
atto = 1e-18
zepto = 1e-21
yocto = 1e-24


#
# Test
#

if __name__ == "__main__":

    print("The 7 base SI units:")
    print("             meter: %s" % meter)
    print("          kilogram: %s" % kilogram)
    print("            second: %s" % second)
    print("            ampere: %s" % ampere)
    print("            kelvin: %s" % kelvin)
    print("              mole: %s" % mole)
    print("           candela: %s" % candela)
    print("")
    print("The 21 SI derived units with special names:")
    print("            radian: %s" % radian)
    print("         steradian: %s" % steradian)
    print("             hertz: %s" % hertz)

    print("            newton: %s" % newton)
    print("            pascal: %s" % pascal)
    print("             joule: %s" % joule)
    print("              watt: %s" % watt)

    print("           coulomb: %s" % coulomb)
    print("              volt: %s" % volt)
    print("             farad: %s" % farad)
    print("               ohm: %s" % ohm)
    print("           siemens: %s" % siemens)
    print("             weber: %s" % weber)
    print("             tesla: %s" % tesla)
    print("             henry: %s" % henry)

    print("    degree Celsius: %s" % celsius)

    print("             lumen: %s" % lumen)
    print("               lux: %s" % lux)

    print("         becquerel: %s" % becquerel)
    print("              gray: %s" % gray)
    print("           sievert: %s" % sievert)
