####################################################################
#
#    Tutorial 2: Using your own reaction mechanism files
#
####################################################################

# You can build a gas mixture object by importing element, species,
# and reaction definitions from input files in the format described in
# the document "Defining Phases and Interfaces"). A set of input files
# in this format is contained in the data folder. 

# Many existing reaction mechanism files are in "CK format," by
# which we mean the input file format developed for use with the
# Chemkin-II software package. [See R. J. Kee, F. M. Rupley, and
# J. A. Miller, Sandia National Laboratories Report SAND89-8009
# (1989).]

# Cantera comes with a converter utility program 'ck2cti' (or 'ck2cti.exe')
# that converts CK format into Cantera format. This program should be run
# from the command line first to convert any CK files you plan to use into
# Cantera format.

from Cantera import *
gas1 = IdealGasMix('gri30.cti')

# This statement creates a mixture that implements GRI-Mech 3.0, much
# like function GRI30 does. File 'gri30.cti' is in the 'data'
# directory. Under Windows, this directory is in
# C:\Program Files\Common Files\Cantera.

# A Cantera input file may contain more than one phase specification, or may
# contain specifications of interfaces (surfaces).

# Use importPhase to import a phase:
gas2 = importPhase('diamond.cti', 'gas')        # a gas
diamond = importPhase('diamond.cti','diamond')  # bulk diamond

# Use importInterface to import a surface:
diamonnd_surf = importInterface('diamond.cti','diamond_100',
                                phases = [gas2, diamond])
# Note that the bulk (i.e., 3D) phases that participate in the surface reactions
# must also be passed as arguments to importInterface.

# How does Cantera find input files like diamond.cti?  Cantera always
# looks in the local directory first. If it is not there, Cantera
# looks for it on its search path. It looks for it in the data
# directory specified when Cantera was built (by default this is
# /usr/local/cantera/data on unix systems). If you define environment
# variable CANTERA_DATA_DIR, it will also look there, or else you can
# call function addDirectory to add a directory to the search path.

# Warning: when Cantera reads a .cti input file, wherever it is
# located, it always writes a file of the same name but with extension
# .xml *in the local directory*. If you happen to have some other file
# by that name, it will be overwritten. Once the XML file is created,
# you can use it instead of the .cti file, which will result in
# somewhat faster startup.



