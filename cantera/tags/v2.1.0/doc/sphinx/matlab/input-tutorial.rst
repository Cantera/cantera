
**********************************
Tutorial: Working with input files
**********************************

.. highlight:: matlab

CTI files
---------

This is the typical way to create a Cantera "phase" object in Matlab::

    gas1 = importPhase('gri30.cti', 'gri30');

Function ``importPhase`` constructs an object representing a phase of matter by
reading in attributes of the phase from a file, which in this case is
``gri30.cti``. This file contains several phase specifications; the one we want
here is ``gri30``, which is specified by the second argument. This file contains
a complete specification of the GRI-Mech 3.0 reaction mechanism, including
element data (name, atomic weight), species data (name, elemental composition,
coefficients to compute thermodynamic and transport properties), and reaction
data (stoichiometry, rate coefficient parameters). The file is written in a
format understood by Cantera, which is described in :ref:`sec-defining-phases`.

CTI files distributed with Cantera
----------------------------------

Several reaction mechanism files in this format are included in the Cantera
distribution, including ones that model high-temperature air, a hydrogen/oxygen
reaction mechanism, and a few surface reaction mechanisms.  These files are kept
in the ``data`` subdirectory within the Cantera installation directory.

If for some reason Cantera has difficulty finding where these files are on your
system, set environment variable ``CANTERA_DATA`` to the directory where they
are located. Alternatively, you can call function ``addDirectory`` to add a
directory to the Cantera search path::

    addDirectory('/usr/local/cantera/my_data_files');

Cantera input files are plain text files, and can be created with any text
editor. See :ref:`sec-defining-phases` for more information.

Importing multiple phases or interfaces
---------------------------------------

A Cantera input file may contain more than one phase specification, or may
contain specifications of interfaces (surfaces). Here we import definitions of
two bulk phases and the interface between them from file ``diamond.cti``::

    gas2 = importPhase('diamond.cti', 'gas');                   % a gas
    diamond = importPhase('diamond.cti', 'diamond');             % bulk diamond
    diamond_surf = importInterface('diamond.cti', 'diamond_100', ...
                                   gas2, diamond);

Note that the bulk (i.e., 3D) phases that participate in the surface reactions
must also be passed as arguments to importInterface.

CTML files
----------

Note that when Cantera reads a ``.cti`` input file, wherever it is located, it
always writes a file of the same name but with extension ``.xml`` *in the local
directory*. If you happen to have some other file by that name, it will be
overwritten. Once the XML file is created, you can use it instead of the
``.cti`` file, which will result in somewhat faster startup::

    gas4 = importPhase('gri30.xml','gri30');

Interfaces can be imported from XML files too::

    diamonnd_surf2 = importInterface('diamond.xml','diamond_100',...
                                     gas2, diamond);

Converting CK-format files
--------------------------

Many existing reaction mechanism files are in "CK format," by which we mean the
input file format developed for use with the Chemkin-II software package.  [See
R. J. Kee, F. M. Rupley, and J. A. Miller, Sandia National Laboratories Report
SAND89-8009 (1989).]

Cantera comes with a converter utility program ``ck2cti`` (or ``ck2cti.exe``)
that converts CK format into Cantera format. This program should be run from the
command line first to convert any CK files you plan to use into Cantera format.

Here's an example of how to use it:

.. code-block:: bash

    ck2cti -i mech.inp -t therm.dat -tr tran.dat -id mymech > mech.cti
