####################################################################
print """

    Tutorial 2: Using your own reaction mechanism files

"""
####################################################################
from Cantera import *
from time import clock

# In the last tutorial, we used function GRI30 to create an object
# that models an ideal gas mixture with the species and reactions of
# GRI-Mech 3.0. Another way to do this is shown here:

gas = importPhase('gri30.cti', 'gri30')

# Function 'importPhase' constructs an object representing a phase of
# matter by reading in attributes of the phase from a file, which in
# this case is 'gri30.cti'. This file contains a complete
# specification of the GRI-Mech 3.0 reaction mechanism, including
# element data (name, atomic weight), species data (name, elemental
# composition, coefficients to compute thermodynamic and transport
# properties), and reaction data (stoichiometry, rate coefficient
# parameters). The file is written in a format understood by Cantera,
# which is described in the document "Defining Phases and Interfaces."


# CTI files distributed with Cantera
#---------------------------------

# Several reaction mechanism files in this format are included in the
# Cantera distribution, including ones that model high-temperature
# air, a hydrogen/oxygen reaction mechanism, and a few surface
# reaction mechanisms. Under Windows, the installation program puts
# these files in 'C:\Program File\Common Files\Cantera.'  On a
# unix/linux/Mac OSX machine, they are usually kept in the 'data'
# subdirectory within the Cantera installation directory.

# If for some reason Cantera has difficulty finding where these files
# are on your system, set environment variable CANTERA_DATA to the
# directory where they are located. Alternatively, you can call function
# addDirectory to add a directory to the Cantera search path:
addDirectory('/usr/local/data')
ggg = importPhase('dummy.cti')

# Cantera input files are plain text files, and can be created with
# any text editor. See the document 'Defining Phases and Interfaces'
# for more information.

from Cantera import *
t0 = clock()
gas1 = importPhase('gri30.cti')
print 'time to create gas1 = ',clock() - t0 

# This statement creates a mixture that implements GRI-Mech 3.0, much
# like function GRI30 does. File 'gri30.cti' is in the 'data'
# directory. Under Windows, this directory is in C:\Program
# Files\Common Files\Cantera and/or C:\CANTERA\DATA. On most other
# platforms, it is usually in /usr/local/cantera/data.


# A Cantera input file may contain more than one phase specification, or may
# contain specifications of interfaces (surfaces).

# Use importPhase to import a phase:
t0 = clock()
gas2 = importPhase('diamond.cti', 'gas')        # a gas
print 'time to create gas2 = ',clock() - t0

t0 = clock()
diamond = importPhase('diamond.cti','diamond')  # bulk diamond
print 'time to create diamond = ',clock() - t0

# Use importInterface to import a surface:
t0 = clock()
diamonnd_surf = importInterface('diamond.cti','diamond_100',
                                phases = [gas2, diamond])
print 'time to create diamond_surf = ',clock() - t0

# Note that the bulk (i.e., 3D) phases that participate in the surface
# reactions must also be passed as arguments to importInterface.

# Multiple phases defined in the same input file can be imported with
# one statement:
t0 = clock()
[gas3, diamond2] = importPhases('diamond.cti', ['gas','diamond'])
print 'time to create both gas3 and diamond2 = ',clock() - t0

# Note that importing from a file is much faster the second time. This
# is because the file is only read and converted to XML once. The XML
# tree is kept in memory once it is read in case it is needed later.

# How does Cantera find input files like diamond.cti?  Cantera always
# looks in the local directory first. If it is not there, Cantera
# looks for it on its search path. It looks for it in the data
# directory specified when Cantera was built (by default this is
# /usr/local/cantera/data on unix systems). If you define environment
# variable CANTERA_DATA, it will also look there, or else you can
# call function addDirectory to add a directory to the search path.

# Warning: when Cantera reads a .cti input file, wherever it is
# located, it always writes a file of the same name but with extension
# .xml *in the local directory*. If you happen to have some other file
# by that name, it will be overwritten. Once the XML file is created,
# you can use it instead of the .cti file, which will result in
# somewhat faster startup.

gas4 = IdealGasMix('gri30.xml')
# Note that the function 'IdealGasMix' simply calls 'importPhase', and
# checks that the phase represents an ideal gas mixture

# Interfaces can be imported from XML files too.
diamonnd_surf2 = importInterface('diamond.xml','diamond_100',
                                 phases = [gas2, diamond])


# Converting CK-format files
# --------------------------

# Many existing reaction mechanism files are in "CK format," by
# which we mean the input file format developed for use with the
# Chemkin-II software package. [See R. J. Kee, F. M. Rupley, and
# J. A. Miller, Sandia National Laboratories Report SAND89-8009
# (1989).]

# Cantera comes with a converter utility program 'ck2cti' (or
# 'ck2cti.exe') that converts CK format into Cantera format. This
# program should be run from the command line first to convert any CK
# files you plan to use into Cantera format. This utility program can
# also be downloaded from the Cantera User's Group web site.



