####################################################################
#
#    Tutorial 2: Using your own reaction mechanism files
#
####################################################################

# You can build a gas mixture object by importing element, species,
# and reaction definitions from input files in supported
# formats. Currently, two formats are supported. 

# Importing CK-format files
# -------------------------

# By 'CK format', we mean the input file format developed for use
# with the Chemkin-II software package. [See R. J. Kee,
# F. M. Rupley, and J. A. Miller, Sandia National Laboratories
# Report SAND89-8009 (1989).]

# These files contain no equation of state information, since an
# ideal gas mixture is implicitly assumed. (Chemkin-II does not
# handle non-ideal gases.) Therefore, it is appropriate in Cantera
# to build from them objects that represent ideal gas
# mixtures. This is done using function IdealGasMix:

from Cantera import *
gas1 = IdealGasMix('mech.inp')

# This statement creates a mixture that implements GRI-Mech 3.0,
# much like function GRI30 does. File 'gri30.inp' is in the 'data'
# directory. Under Windows, this directory 
#
# Cantera always looks in the local directory first, however. So if
# you have a file of the same name in the local directory,
# it will be used instead.


# The CK file format specification does not require that all
# species data be contained in the file. Missing species
# definitions (usually called 'thermo' data but in fact defining
# all properties of the species, including name, phase, and
# elemental composition) are to be looked up in a second
# 'thermodynamic database' file. To create the object from an
# incomplete input file, give both file names as arguments:

gas2 = IdealGasMix('air.inp','nasathermo.dat')

# The CK file specification does not include transport data for the
# species. These too are taken from an external database. If you
# need transport properties, include the transport file name and
# transport model to implement as follows:

gas3 = IdealGasMix(src = 'gri30.inp',
                   transport_db = 'gri30_tran.dat',
                   transport = 'Multi')

# Allowable values for the transport model are 'Multi' and
# 'Mix'. If the model is omitted, 'Mix' is assumed.


# Importing CTML files
# -------------------

# Cantera can also read input files in an XML-based format called
# 'CTML' (Cantera Markup Language).  These input files are complete,
# and do not require auxiliary database files for either thermodynamic
# or transport properties. To import a CTML file, simply give the file
# name in the call to IdealGasMix:

gxml = IdealGasMix('gri30.xml')

# Cantera determines the file format by examining its contents.  A
# conversion utility is available that converts CK-format files into
# CTML.



