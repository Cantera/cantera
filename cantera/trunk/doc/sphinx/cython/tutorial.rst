.. py:currentmodule:: cantera

Tutorial
========

Getting Started
---------------

Start by creating a new Python script (filename ending with ".py") in the text editor of your choice.

Put the following statement at the top of your script to import the Cantera Python module::

    >>> import cantera as ct

When using Cantera, the first thing you usually need is an object representing
some phase of matter. Here, we'll create a gas mixture::

    >>> gas1 = ct.Solution('gri30.xml')

To view the state of the mixture, *call* the `gas1` object as if it were a
function::

    >>> gas1()

You should see something like this::

     gri30:

          temperature             300  K
             pressure          101325  Pa
              density       0.0818891  kg/m^3
     mean mol. weight         2.01588  amu

                             1 kg            1 kmol
                          -----------      ------------
             enthalpy         26470.1        5.336e+04     J
      internal energy    -1.21087e+06       -2.441e+06     J
              entropy         64913.9        1.309e+05     J/K
       Gibbs function    -1.94477e+07        -3.92e+07     J
    heat capacity c_p         14311.8        2.885e+04     J/K
    heat capacity c_v         10187.3        2.054e+04     J/K

                              X                 Y          Chem. Pot. / RT
                        -------------     ------------     ------------
                   H2              1                1         -15.7173
                    H              0                0
                    O              0                0
                   O2              0                0
                   OH              0                0
                  H2O              0                0
                  HO2              0                0
                 H2O2              0                0
                    C              0                0
                   CH              0                0
                  CH2              0                0
               CH2(S)              0                0
                  CH3              0                0
                  CH4              0                0
                   CO              0                0
                  CO2              0                0
                  HCO              0                0
                 CH2O              0                0
                CH2OH              0                0
                 CH3O              0                0
                CH3OH              0                0
                  C2H              0                0
                 C2H2              0                0
                 C2H3              0                0
                 C2H4              0                0
                 C2H5              0                0
                 C2H6              0                0
                 HCCO              0                0
                CH2CO              0                0
                HCCOH              0                0
                    N              0                0
                   NH              0                0
                  NH2              0                0
                  NH3              0                0
                  NNH              0                0
                   NO              0                0
                  NO2              0                0
                  N2O              0                0
                  HNO              0                0
                   CN              0                0
                  HCN              0                0
                 H2CN              0                0
                 HCNN              0                0
                 HCNO              0                0
                 HOCN              0                0
                 HNCO              0                0
                  NCO              0                0
                   N2              0                0
                   AR              0                0
                 C3H7              0                0
                 C3H8              0                0
               CH2CHO              0                0
               CH3CHO              0                0

What you have just done is to create an object, `gas1` that implements GRI-
Mech 3.0, the 53-species, 325-reaction natural gas combustion mechanism
developed by Gregory P. Smith, David M. Golden, Michael Frenklach, Nigel W.
Moriarty, Boris Eiteneer, Mikhail Goldenberg, C. Thomas Bowman, Ronald K.
Hanson, Soonho Song, William C. Gardiner, Jr., Vitali V. Lissianski, and
Zhiwei Qin. See http://www.me.berkeley.edu/gri_mech/ for more information.

The `gas1` object has properties you would expect for a gas mixture - it has a
temperature, a pressure, species mole and mass fractions, etc. As we'll soon
see, it has many more properties.

The summary of the state of `gas1` printed above shows that new objects
created from the `gri30.xml` input file start out with a temperature of 300 K,
a pressure of 1 atm, and have a composition that consists of only one species,
in this case hydrogen. There is nothing special about H2 - it just happens to
be the first species listed in the input file defining GRI-Mech 3.0. In
general, whichever species is listed first will initially have a mole fraction
of 1.0, and all of the others will be zero.

Setting the State
~~~~~~~~~~~~~~~~~

The state of the object can easily be changed. For example::

    >>> gas1.TP = 1200, 101325

sets the temperature to 1200 K and the pressure to 101325 Pa (Cantera always
uses SI units). After this statement, calling ``gas1()`` results in::

     gri30:

          temperature            1200  K
             pressure          101325  Pa
              density       0.0204723  kg/m^3
     mean mol. weight         2.01588  amu

                             1 kg            1 kmol
                          -----------      ------------
             enthalpy     1.32956e+07         2.68e+07     J
      internal energy     8.34619e+06        1.682e+07     J
              entropy         85227.6        1.718e+05     J/K
       Gibbs function    -8.89775e+07       -1.794e+08     J
    heat capacity c_p         15377.9          3.1e+04     J/K
    heat capacity c_v         11253.4        2.269e+04     J/K

                              X                 Y          Chem. Pot. / RT
                        -------------     ------------     ------------
                   H2              1                1         -17.9775
                        (other species not shown)

Notice that the temperature has been changed as requested, but the pressure
has changed too. The density and composition have not.

Thermodynamics generally requires that *two* properties in addition to
composition information be specified to fix the intensive state of a substance
(or mixture). The state of the mixture can be set using several combinations
of two properties. The following are all equivalent::

    >>> gas1.TP = 1200, 101325           # temperature, pressure
    >>> gas1.TD = 1200, 0.0204723        # temperature, density
    >>> gas1.HP = 1.32956e7, 101325      # specific enthalpy, pressure
    >>> gas1.UV = 8.34619e6, 1/0.0204723 # specific internal energy, specific volume
    >>> gas1.SP = 85227.6, 101325        # specific entropy, pressure
    >>> gas1.SV = 85227.6, 1/0.0204723   # specific entropy, specific volume

In each case, the values of the extensive properties must be entered *per unit
mass*.

Properties may be read independently or together::

    >>> gas1.T
    1200.0
    >>> gas1.h
    13295567.68
    >>> gas1.UV
    (8346188.494954427, 48.8465747765848)

The composition can be set in terms of either mole fractions (``X``) or mass
fractions (``Y``)::

    >>> gas1.X = 'CH4:1, O2:2, N2:7.52'

When the composition alone is changed, the temperature and density are held
constant. This means that the pressure and other intensive properties will
change. The composition can also be set in conjunction with the intensive
properties of the mixture::

    >>> gas1.TPX = 1200, 101325, 'CH4:1, O2:2, N2:7.52'
    >>> gas1()

results in::

     gri30:

          temperature            1200  K
             pressure          101325  Pa
              density        0.280629  kg/m^3
     mean mol. weight         27.6332  amu

                             1 kg            1 kmol
                          -----------      ------------
             enthalpy          861943        2.382e+07     J
      internal energy          500879        1.384e+07     J
              entropy          8914.3        2.463e+05     J/K
       Gibbs function    -9.83522e+06       -2.718e+08     J
    heat capacity c_p         1397.26        3.861e+04     J/K
    heat capacity c_v         1096.38         3.03e+04     J/K

                              X                 Y          Chem. Pot. / RT
                        -------------     ------------     ------------
                   H2              0                0
                    H              0                0
                    O              0                0
                   O2       0.190114         0.220149         -28.7472
                   OH              0                0
                  H2O              0                0
                  HO2              0                0
                 H2O2              0                0
                    C              0                0
                   CH              0                0
                  CH2              0                0
               CH2(S)              0                0
                  CH3              0                0
                  CH4       0.095057        0.0551863          -35.961
                   CO              0                0
                  CO2              0                0
                  HCO              0                0
                 CH2O              0                0
                CH2OH              0                0
                 CH3O              0                0
                CH3OH              0                0
                  C2H              0                0
                 C2H2              0                0
                 C2H3              0                0
                 C2H4              0                0
                 C2H5              0                0
                 C2H6              0                0
                 HCCO              0                0
                CH2CO              0                0
                HCCOH              0                0
                    N              0                0
                   NH              0                0
                  NH2              0                0
                  NH3              0                0
                  NNH              0                0
                   NO              0                0
                  NO2              0                0
                  N2O              0                0
                  HNO              0                0
                   CN              0                0
                  HCN              0                0
                 H2CN              0                0
                 HCNN              0                0
                 HCNO              0                0
                 HOCN              0                0
                 HNCO              0                0
                  NCO              0                0
                   N2       0.714829         0.724665         -25.6789
                   AR              0                0
                 C3H7              0                0
                 C3H8              0                0
               CH2CHO              0                0
               CH3CHO              0                0

The composition above was specified using a string. The format is a comma-
separated list of ``<species name>:<relative mole numbers>`` pairs. The mole
numbers will be normalized to produce the mole fractions, and therefore they
are "relative" mole numbers. Mass fractions can be set in this way too by
changing ``X`` to ``Y`` in the above statements.

The composition can also be set using an array, which must have the same size
as the number of species. For example, to set all 53 mole fractions to the
same value, do this::

    >>> gas1.X = np.ones(53) # NumPy array of 53 ones

Or, to set all the mass fractions to equal values::

    >>> gas1.Y = np.ones(53)

When setting the state, you can control what properties are held constant by
passing the special value `None` to the property setter. For example, to
change the specific volume to 2.1 m^3/kg while holding entropy constant::

    >>> gas.SV = None, 2.1

Or to set the mass fractions while holding temperature and pressure constant::

    >>> gas.TPX = None, None, 'CH4:1.0, O2:0.5'

Working With Mechanism Files
----------------------------

In previous example, we created an object that models an ideal gas mixture
with the species and reactions of GRI-Mech 3.0, using the ``gri30.xml`` input
file included with Cantera. This is a "pre-processed" XML input file written
in a format that is easy for Cantera to parse. Cantera also supports an input
file format that is easier to write, called *CTI*. Several reaction mechanism
files in this format are included with Cantera, including ones that model
high- temperature air, a hydrogen/oxygen reaction mechanism, and a few surface
reaction mechanisms. These files are usually located in the ``data``
subdirectory of the Cantera installation directory, e.g. ``C:\\Program
Files\\Cantera\\data`` on Windows or ``/usr/local/cantera/data/`` on
Unix/Linux/Mac OS X machines, depending on how you installed Cantera and the
options you specified.

If for some reason Cantera has difficulty finding where these files are on
your system, set environment variable ``CANTERA_DATA`` to the directory where
they are located. Alternatively, you can call function `add_directory` to add
a directory to the Cantera search path::

    >>> ct.add_directory('/usr/local/cantera/my_data_files')

Cantera input files are plain text files, and can be created with any text
editor. See the document :ref:`sec-defining-phases` for more information.

A Cantera input file may contain more than one phase specification, or may
contain specifications of interfaces (surfaces). Here we import definitions of
two bulk phases and the interface between them from file ``diamond.cti``::

    >>> gas2 = ct.Solution('diamond.cti', 'gas')
    >>> diamond = ct.Solution('diamond.cti', 'diamond')
    >>> diamond_surf = ct.Interface('diamond.cti' , 'diamond_100',
                                    [gas2, diamond])

Note that the bulk (i.e., 3D or homogeneous) phases that participate in the
surface reactions must also be passed as arguments to `Interface`.

When Cantera reads a ``.cti`` input file, wherever it is located, it always
writes a file of the same name but with extension ``.xml`` *in the local
directory*. If you happen to have some other file by that name, it will be
overwritten. Once the XML file is created, you can use it instead of the
``.cti`` file, which will result in somewhat faster startup.

Converting CK-format files
~~~~~~~~~~~~~~~~~~~~~~~~~~

Many existing reaction mechanism files are in "CK format," by which we mean
the input file format developed for use with the Chemkin-II software package.
[See R. J. Kee, F. M. Rupley, and J. A. Miller, Sandia National Laboratories
Report SAND89-8009 (1989).]

Cantera comes with a converter utility program ``ck2cti`` (or ``ck2cti.py``)
that converts CK format into Cantera format. This program should be run from
the command line first to convert any CK files you plan to use into Cantera
format. Here's an example of how to use it. The command::

    $python ck2cti.py --input=mech.inp --thermo=therm.dat --transport=tran.dat

will produce the file ``mech.cti`` in the current directory.


Getting Help
------------

In addition to the Sphinx-generated :ref:`sec-cython-documentation`,
documentation of the Python classes and their methods can be accessed from
within the Python interpreter as well.

Suppose you have created a Cantera object and want to know what methods are
available for it, and get help on using the methods::

    >>> g = Cantera.Solution('gri30.xml')

To get help on the Python class that this object is an instance of::

    >>> help(g)

For a simple list of the properties and methods of this object::

    >>> dir(g)

To get help on a specific method, e.g. the ``species_index`` method::

    >>> help(g.species_index)

For properties, getting the documentation is slightly trickier, as the usual
method will give you the help for the *result*, e.g.::

    >>> help(g.T)

will provide help on Python's ``float`` class. To get the help for the
temperature property, ask for the attribute of the class object itself::

    >>> help(g.__class__.T)

If you are using the IPython shell, help can also be obtained using the `?`
syntax::

    In[1]: g.species_index?

Chemical Equilibrium
--------------------

To set a gas mixture to a state of chemical equilibrium, use the equilibrate
method::

    >>> import cantera as ct
    >>> g = ct.Solution('gri30.xml')
    >>> g.TPX = 300.0, ct.one_atm, 'CH4:0.95,O2:2,N2:7.52'
    >>> g.equilibrate('TP')

The above statement sets the state of object ``g`` to the state of chemical
equilibrium holding temperature and pressure fixed. Alternatively, the
specific enthalpy and pressure can be held fixed::

    >>> g.TPX = 300.0, ct.one_atm, 'CH4:0.95,O2:2,N2:7.52'
    >>> g.equilibrate('HP')

Other options are:

    - 'UV'   fixed specific internal energy and specific volume
    - 'SV'   fixed specific entropy and specific volume
    - 'SP'   fixed specific entropy and pressure

How can you tell if ``equilibrate`` has correctly found the chemical equilibrium
state? One way is verify that the net rates of progress of all reversible
reactions are zero. Here is the code to do this:

    >>> g.TPX = 300.0, ct.one_atm, 'CH4:0.95,O2:2,N2:7.52'
    >>> g.equilibrate('HP')

    >>> rf = g.forward_rates_of_progress
    >>> rr = g.reverse_rates_of_progress
    >>> for i in range(g.n_reactions):
    >>>     if g.is_reversible(i) and rf[i] != 0.0:
    >>>         print(' %4i  %10.4g  ' % (i, (rf[i] - rr[i])/rf[i]))

If the magnitudes of the numbers in this list are all very small, then each
reversible reaction is very nearly equilibrated, which only occurs if the gas
is in chemical equilibrium.

You might be wondering how ``equilibrate`` works. (Then again, you might not).
Method ``equilibrate`` invokes Cantera's chemical equilibrium solver, which uses
an element potential method. The element potential method is one of a class of
equivalent *nonstoichiometric* methods that all have the characteristic that
the problem reduces to solving a set of M nonlinear algebraic equations, where
M is the number of elements (not species). The so-called *stoichiometric*
methods, on the other hand, (including Gibbs minimization), require solving K
nonlinear equations, where K is the number of species (usually K >> M). See
Smith and Missen, "Chemical Reaction Equilibrium Analysis" for more
information on the various algorithms and their characteristics.

Cantera uses a damped Newton method to solve these equations, and does a few
other things to generate a good starting guess and to produce a reasonably
robust algorithm. If you want to know more about the details, look at the on-
line documented source code of Cantera C++ class 'ChemEquil.h'.
