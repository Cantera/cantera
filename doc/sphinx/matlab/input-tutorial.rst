
**********************************
Tutorial: Working with input files
**********************************

.. highlight:: matlab

CTI files
---------

This is the typical way to create a Cantera "phase" object in Matlab::

    gas1 = Solution('gri30.cti', 'gri30');

This statement constructs a ``Solution`` object representing a phase of matter by
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
system, set environment variable ``CANTERA_DATA`` to the directory or
directories (separated using ``;`` on Windows or ``:`` on other operating
systems) where they are located. Alternatively, you can call function
`add_directory` to add a directory to the Cantera search path::

    addDirectory('/usr/local/cantera/my_data_files');

Cantera input files are plain text files, and can be created with any text
editor. See :ref:`sec-defining-phases` for more information.

Importing multiple phases or interfaces
---------------------------------------

A Cantera input file may contain more than one phase specification, or may
contain specifications of interfaces (surfaces). Here we import definitions of
two bulk phases and the interface between them from file ``diamond.cti``::

    gas2 = Solution('diamond.cti', 'gas');                   % a gas
    diamond = Solution('diamond.cti', 'diamond');             % bulk diamond
    diamond_surf = importInterface('diamond.cti', 'diamond_100', ...
                                   gas2, diamond);

Note that the bulk (i.e., 3D) phases that participate in the surface reactions
must also be passed as arguments to importInterface.

Converting CK-format files
--------------------------

See :ref:`sec-ck-format-conversion` in the :ref:`sec-input-files` documentation.
