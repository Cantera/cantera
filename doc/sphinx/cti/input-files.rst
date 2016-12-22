.. py:currentmodule:: cantera.ctml_writer

.. _sec-input-files:

************************
Working with Input Files
************************

Before we can describe how to define phases, interfaces, and their components
(elements, species, and reactions), we need to go over a few points about the
mechanics of writing and processing input files.

Input File Syntax
=================

An input file consists of *entries* and *directives*, both of which have a
syntax much like functions. An entry defines an object---for example, a
reaction, or a species, or a phase. A directive sets options that affect how the
entry parameters are interpreted, such as the default unit system, or how
certain errors should be handled.

Cantera's input files follow the syntax rules for Python, so if you're familiar
with Python syntax you already understand many of the details and can probably
skip ahead to :ref:`sec-dimensions`.

Entries have fields that can be assigned values. A species entry is shown below
that has fields *name* and *atoms* (plus several others)::

    species(name='C60', atoms='C:60')

Most entries have some fields that are required; these must be assigned values,
or else processing of the file will abort and an error message will be
printed. Other fields may be optional, and take default values if not assigned.

An entry may be either a *top-level entry* or an *embedded entry*. Top-level
entries specify a phase, an interface, an element, a species, or a reaction, and
begin in the first (leftmost) column. Embedded entries specify a model, or a
group of parameters for a top-level entry, and are usually embedded in a field
of another entry.

The fields of an entry are specified in the form ``<field_name> = <value>``, and may
be listed on one line, or extend across several. For example, two entries for
graphite are shown below. The first is compact::

    stoichiometric_solid(name='graphite', species='C(gr)', elements='C', density=(2.2, 'g/cm3'))

and the second is formatted to be easier to read::

    stoichiometric_solid(
        name     = 'graphite',
        elements = 'C',
        species  = 'C(gr)',
        density  = (2.2, 'g/cm3')
    )

Both are completely equivalent.

The species ``C(gr)`` that appears in the definition of the graphite phase is
also defined by a top-level entry. If the heat capacity of graphite is
approximated as constant, then the following definition could be used::

    species(name='C(gr)',
            atoms='C:1',
            thermo=const_cp(t0=298.15,
                            h0=0.0,
                            s0=(5.6, 'J/mol/K'), # NIST
                            cp0=(8.43, 'J/mol/K'))) # Taylor and Groot (1980)

Note that the thermo field is assigned an embedded entry of type
:class:`const_cp`. Entries are stored as they are encountered when the file is
read, and only processed once the end of the file has been reached. Therefore,
the order in which they appear is unimportant.

Comments
--------

The character ``#`` is the comment character. Everything to the right of this
character on a line is ignored::

    # set the default units
    units(length = 'cm',    # use centimeters for length
          quantity = 'mol') # use moles for quantity

Strings
-------

Strings may be enclosed in single quotes or double quotes, but they must
match. To create a string containing single quotes, enclose it in double quotes,
and vice versa. If you want to create a string to extend over multiple lines,
enclose it in triple quotes::

    string1 = 'A string.'
    string2 = "Also a 'string'"
    string3 = """This is
    a
    string too."""

The multi-line form is useful when specifying a phase containing a large number
of species::

    species = """ H2 H O O2 OH H2O HO2 H2O2 C CH
                  CH2 CH2(S) CH3 CH4 CO CO2 HCO CH2O CH2OH CH3O
                  CH3OH C2H C2H2 C2H3 C2H4 C2H5 C2H6 HCCO CH2CO HCCOH
                  N NH NH2 NH3 NNH NO NO2 N2O HNO CN
                  HCN H2CN HCNN HCNO HOCN HNCO NCO N2 AR C3H7
                  C3H8 CH2CHO CH3CHO """

Sequences
---------

A sequence of multiple items is specified by separating the items by commas and
enclosing them in square brackets or parentheses. The individual items can have
any type---strings, integers, floating-point numbers (or even entries or other
lists). Square brackets are often preferred, since parentheses are also used for
other purposes in the input file, but either can be used::

    s0 = (3.5, 'J/mol/K') # these are
    s0 = [3.5, 'J/mol/K'] # equivalent

Variables
---------

Another way to specify the species C(gr) is shown here::

    graphite_thermo = const_cp(t0=298.15,
                               h0=0.0,
                               s0=(5.6, 'J/mol/K'), # NIST
                               cp0=(8.43, 'J/mol/K')) # Taylor and Groot (1980)

    species(name='C(gr)', atoms='C:1', thermo=graphite_thermo)

In this form, the ``const_cp`` entry is stored in a variable, instead of being
directly embedded within the species entry.  The *thermo* field is assigned this
variable.

Variables can also be used for any other parameter type. For example, if you are
defining several phases in the file, and you want to set them all to the same
initial pressure, you could define a pressure variable::

    P_initial = (2.0, 'atm')

and then set the pressure field in each embedded state entry to this variable.

Omitting Field Names
--------------------

Field names may be omitted if the values are entered in the order specified in
the entry declaration. (Entry declarations are the text printed on a colored
background in the following chapters.) It is also possible to omit only some of
the field names, as long as these fields are listed first, in order, before any
named fields.

For example, The first four entries below are equivalent, while the last two are
incorrect and would generate an error when processed::

    element(symbol="Ar", atomic_mass=39.948) # OK
    element(atomic_mass=39.948, symbol='Ar') # OK
    element('Ar', atomic_mass=39.948)        # OK
    element("Ar", 39.948)                    # OK

    element(39.948, "Ar")                    # error
    element(symbol="Ar", 39.948)             # error

Validation
----------

Normally, Cantera will make some checks for errors in the definitions of species
and reactions, such as checking for duplicate reactions. To slightly speed up
processing (if a mechanism has previously been validated), or in case of
spurious validation errors, validation can be disabled using the
:func:`validate` function. For example, to disable validation of reactions, add
the following to the CTI file::

    validate(reactions='no')

.. _sec-dimensions:

Dimensional Values
==================

Many fields have numerical values that represent dimensional quantities---a
pressure, or a density, for example. If these are entered without specifying the
units, the default units (set by the :class:`units` directive described in
:ref:`sec-default-units`) will be used. However, it is also possible to specify
the units for each individual dimensional quantity (unless stated
otherwise). All that is required is to group the value in parentheses or square
brackets with a string specifying the units::

    pressure = 1.0e5 # default is Pascals
    pressure = (1.0, 'bar') # this is equivalent
    density = (4.0, 'g/cm3')
    density = 4000.0 # kg/m3

Compound unit strings may be used, as long as a few rules are followed:

1. Units in the denominator follow ``/``.
2. Units in the numerator follow ``-``, except for the first one.
3. Numerical exponents follow the unit string without a ``^`` character, and must
   be in the range 2--6. Negative values are not allowed.

Examples of compound units::

    A = (1.0e20, 'cm6/mol2/s') # OK
    h = (6.626e-34, 'J-s')     # OK
    density = (3.0, 'g/cm3')   # OK
    A = (1.0e20, 'cm^6/mol/s') # error (^)
    A = (1.0e20, 'cm6/mol2-s') # error ('s' should be in denominator)
    density = (3.0, 'g-cm-3')  # error (negative exponent)

.. _sec-default-units:

Setting the Default Units
-------------------------

The default unit system may be set with the :func:`units` directive. Note
that unit conversions are not done until the entire file has been read. Only one
units directive should be present in a file, and the defaults it specifies apply
to the entire file.  If the file does not contain a units directive, the default
units are meters, kilograms, kilomoles, and seconds.

Shown below are two equivalent ways of specifying the site density for an
interface. In the first version, the site density is specified without a units
string, and so its units are constructed from the default units for quantity and
length, which are set with a units directive::

    units(length = 'cm', quantity = 'molec')
    interface(name = 'Si-100',
              site_density = 1.0e15, # molecules/cm2 (default units)
              # ...
              )

The second version uses a different default unit system, but overrides the
default units by specifying an explicit units string for the site density::

    units(length = 'cm', quantity = 'mol')
    interface(name = 'Si-100',
              site_density = (1.0e15, 'molec/cm2') # override default units
              # ...
              )

The second version is equivalent to the first, but would be very different if
the units of the site density were not specified!

The *length*, *quantity* and *time* units are used to construct the units for
reaction pre-exponential factors. The *energy* units are used for molar
thermodynamic properties, in combination with the units for *quantity*.

Since activation energies are often specified in units other than those used for
thermodynamic properties, a separate field is devoted to the default units for
activation energies::

    units(length = 'cm', quantity = 'mol', act_energy = 'kcal/mol')
    kf = Arrhenius(A = 1.0e14, b = 0.0, E = 54.0) # E is 54 kcal/mol

See :func:`units` for the declaration of the units directive.

Recognized Units
----------------

Cantera recognizes the following units in various contexts:

===========  ==============
field        allowed values
===========  ==============
length       ``'cm', 'm', 'mm'``
quantity     ``'mol', 'kmol', 'molec'``
time         ``'s', 'min', 'hr', 'ms'``
energy       ``'J', 'kJ', 'cal', 'kcal'``
act_energy   ``'kJ/mol', 'J/mol', 'J/kmol', 'kcal/mol', 'cal/mol', 'eV', 'K'``
pressure     ``'Pa', 'atm', 'bar'``
===========  ==============

Processing Input Files
======================

A Two-step Process
------------------

From the point of view of the user, it appears that a Cantera application that
imports a phase definition reads the input file, and uses the information there
to construct the object representing the phase or interface in the
application. While this is the net effect, it is actually a two-step
process. When a constructor like ``Solution`` is called to import a phase definition
from a file, a preprocessor runs automatically to read the input file and create
a string that contains the same information but in an XML-based format called
CTML. After the preprocessor finishes, Cantera imports the phase definition from
this CTML data.

Two File Formats
----------------

Why two file formats? There are several reasons. XML is a widely-used standard
for data files, and it is designed to be relatively easy to parse. This makes it
possible for other applications to use Cantera CTML data files, without
requiring the substantial chemical knowledge that would be required to use .cti
files. For example, "web services" (small applications that run remotely over a
network) are often designed to accept XML input data over the network, perform a
calculation, and send the output in XML back across the network. Supporting an
XML-based data file format facilitates using Cantera in web services or other
network computing applications.

The difference between the high-level description in a .cti input file and the
lower-level description in the CTML file may be illustrated by how reactions are
handled. In the input file, the reaction stoichiometry and its reversibility or
irreversibility are determined from the reaction equation. For example::

    O + HCCO <=> H + 2 CO

specifies a reversible reaction between an oxygen atom and the ketenyl radical
HCCO to produce one hydrogen atom and two carbon monoxide molecules. If ``<=>``
were replaced with ``=>``, then it would specify that the reaction should be
treated as irreversible.

Of course, this convention is not spelled out in the input file---the parser
simply has to know it, and has to also know that a "reactant" appears on the
left side of the equation, a "product" on the right, that the optional number in
front of a species name is its stoichiometric coefficient (but if missing the
value is one), etc. The preprocessor does know all this, but we cannot expect
the same level of knowledge of chemical conventions by a generic XML parser.

Therefore, in the CTML file, reactions are explicitly specified to be reversible
or irreversible, and the reactants and products are explicitly listed with their
stoichiometric coefficients. The XML file is, in a sense, a "dumbed-down"
version of the input file, spelling out explicitly things that are only implied
in the input file syntax, so that "dumb" (i.e., easy to write) parsers can be
used to read the data with minimal risk of misinterpretation.

The reaction definition::

    reaction( "O + HCCO <=> H + 2 CO", [1.00000E+14, 0, 0])

in the input file is translated by the preprocessor to the following CTML text:

.. code-block:: xml

    <reaction id="0028" reversible="yes">
      <equation>O + HCCO [=] H + 2 CO</equation>
      <rateCoeff>
        <Arrhenius>
          <A units="cm3/mol/s"> 1.000000E+14</A>
          <b>0</b>
          <E units="cal/mol">0.000000</E>
        </Arrhenius>
      </rateCoeff>
      <reactants>HCCO:1 O:1</reactants>
      <products>H:1 CO:2</products>
    </reaction>

The CTML version is much more verbose, and would be much more tedious to write
by hand, but is much easier to parse, particularly since it is not necessary to
write a custom parser---virtually any standard XML parser, of which there are
many, can be used to read the CTML data.

So in general files that are easy for knowledgeable users (you) to write are more
difficult for machines to parse, because they make use of high-level
application-specific knowledge and conventions to simplify the
notation. Conversely, files that are designed to be easily parsed are tedious to
write because so much has to be spelled out explicitly. A natural solution is to
use two formats, one designed for writing by humans, the other for reading by
machines, and provide a preprocessor to convert the human-friendly format to the
machine-friendly one.

Preprocessor Internals: the ``ctml_writer`` Module
--------------------------------------------------

If you are interested in seeing the internals of how the preprocessing works,
take a look at file ``ctml_writer.py`` in the Cantera Python package. Or simply
start Python, and type::

    >>> import cantera.ctml_writer
    >>> help(cantera.ctml_writer)

The ``ctml_writer.py`` module can also be run as a script to convert input .cti
files to CTML. For example, if you have an input file ``phasedefs.cti``, then
simply type at the command line::

    python -m cantera.ctml_writer phasedefs.cti

to create CTML file ``phasedefs.xml``. On systems which support running Python
scripts directly, a script to run ``ctml_writer`` directly is also installed. If
the Cantera ``bin`` directory is on your ``PATH``, you can also do the
conversion by running::

    ctml_writer phasedefs.cti

This can be used to generate XML input files for use on systems where the
Cantera Python package is not installed. Of course, most of the time creation of
the CTML file will happen behind the scenes, and you will not need to be
concerned with CTML files at all.

Error Handling
==============

During processing of an input file, errors may be encountered. These could be
syntax errors, or could be ones that are flagged as errors by Cantera due to
some apparent inconsistency in the data---an unphysical value, a species that
contains an undeclared element, a reaction that contains an undeclared species,
missing species or element definitions, multiple definitions of elements,
species, or reactions, and so on.

Syntax Errors
-------------

Syntax errors are caught by the Python preprocessor, not by Cantera, and must be
corrected before proceeding further.  Python prints a "traceback" that allows
you to find the line that contains the error. For example, consider the
following input file, which is intended to create a gas with the species and
reactions of GRI-Mech 3.0, but has a misspelled the field name ``reactions``::

    ideal_gas(name = 'gas',
              elements = 'H O',
              species = 'gri30: all',
              reactionss = 'gri30: all')

When this definition is imported into an application, an error message like the
following would be printed to the screen, and execution of the program or script
would terminate. ::

    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/some/path/Cantera/importFromFile.py", line 18, in importPhase
        return importPhases(file, [name], loglevel, debug)[0]
      File "/some/path/Cantera/importFromFile.py", line 25, in importPhases
        s.append(solution.Solution(src=file,id=nm,loglevel=loglevel,debug=debug))
      File "/some/path/solution.py", line 39, in __init__
        preprocess = 1, debug = debug)
      File "/some/path/Cantera/XML.py", line 35, in __init__
        self._xml_id = _cantera.xml_get_XML_File(src, debug)
    cantera.error:

    ************************************************
                    Cantera Error!
    ************************************************

    Procedure: ct2ctml
    Error:   Error converting input file "./gas.cti" to CTML.
    Python command was: '/usr/bin/python'
    The exit code was: 4
    -------------- start of converter log --------------
    TypeError on line 4 of './gas.cti':
    __init__() got an unexpected keyword argument 'reactionss'

    | Line |
    |    1 | ideal_gas(name = 'gas',
    |    2 |           elements = 'H O',
    |    3 |           species = 'gri30: all',
    >    4 >           reactionss = 'gri30: all')
    |    5 |
    --------------- end of converter log ---------------

The top part of the error message shows the chain of functions that were called
before the error was encountered. For the most part, these are internal Cantera
functions not of direct concern here. The relevant part of this error message is
the part starting with the "Cantera Error" heading, and specifically the
contents of the *converter log* section. This message says that that on line 4
of ``gas.cti``, the the keyword argument ``reactionss`` was not
recognized. Seeing this message, it is clear that the problem is that
*reactions* is misspelled.

Cantera Errors
--------------

Now let's consider the other class of errors---ones that Cantera, not Python,
detects. Continuing the example above, suppose that the misspelling is
corrected, and the input file processed again. Again an error message results,
but this time it is from Cantera::

    cantera.error:
    Procedure: installSpecies
    Error: species C contains undeclared element C

The problem is that the phase definition specifies that all species are to be
imported from dataset gri30, but only the elements H and O are declared. The
gri30 datset contains species composed of the elements H, O, C, N, and Ar. If
the definition is modified to declare these additional elements::

    ideal_gas(name = 'gas',
              elements = 'H O C N Ar',
              species = 'gri30: all',
              reactions = 'gri30: all')

it may be imported successfully.

Errors of this type do not have to be fatal, as long as you tell Cantera how you
want to handle them. You can, for example, instruct Cantera to quietly skip
importing any species that contain undeclared elements, instead of flagging them
as errors. You can also specify that reactions containing undeclared species
(also usually an error) should be skipped. This allows you to very easily
extract a portion of a large reaction mechanism, as described in :ref:`sec-phase-options`.

.. _sec-ck-format-conversion:

Converting CK-format files
==========================

Many existing reaction mechanism files are in "CK format," by which we mean
the input file format developed for use with the Chemkin-II software package
as specified in the report describing the Chemkin software [SAND89]_.

Cantera comes with a converter utility program ``ck2cti`` (or ``ck2cti.py``)
that converts CK format into Cantera format. This program should be run from
the command line first to convert any CK files you plan to use into Cantera
format (CTI format).

Usage::

    ck2cti [--input=<filename>]
           [--thermo=<filename>]
           [--transport=<filename>]
           [--surface=<filename>]
           [--id=<phase-id>]
           [--output=<filename>]
           [--permissive]
           [-d | --debug]

Each of the terms in square brackets is an option that can be passed on the
command line to ``ck2cti``. ``--input`` is the chemistry input file, containing
a list of all the element names that are used, a list of all the species names,
and a list of all the reactions to be considered between the species. This file
can also optionally contain thermodynamic information for the species. If the
``--input`` file does not contain the thermodynamic data, a separate file
containing this information must be specified to the `--thermo`` option. Finally,
the ``--input`` file can also optionally contain transport information for the
species. If it does not, and the user wishes to use a part of Cantera that relies
on some transport properties, the ``--transport`` option must be used to specify
the file containing all the transport data for the species.

For the case of a surface mechanism, the gas phase input file should be
specified as ``--input`` and the surface phase input file should be specified as
``--surface``.

Example::

    ck2cti --input=chem.inp --thermo=therm.dat --transport=tran.dat

If the output file name is not given, an output file with the same name as the
input file, with the extension changed to '.cti'.

If the ck2cti script is not on your path but the Cantera Python module is,
ck2cti can also be used by running::

    python -m cantera.ck2cti --input=chem.inp --thermo=therm.dat --transport=tran.dat

An input file containing only species definitions (which can be referenced from
phase definitions in other input files) can be created by specifying only a
thermo file.

Many existing CK format files cause errors in ``ck2cti`` when they are
processed. Some of these errors may be avoided by specifying the
``--permissive`` option. This option allows certain recoverable parsing errors
(e.g. duplicate transport or thermodynamic data) to be ignored. Other errors
may be caused by incorrect formatting of lines in one or more of the input files.

Debugging common errors in CK files
-----------------------------------

When ``ck2cti`` encounters an error, it attempts to print the surrounding
information to help you to locate the error. Many of the most common errors
are due to an inconsistency of the input files from their standard, as defined
in the report for Chemkin referenced above. These errors include:

  * Each section of the input files must be started with a keyword representing that
    section and ending with the keyword ``END``. Keywords that may begin a section
    include:

    - ``ELEMENTS`` or ``ELEM``
    - ``SPECIES`` or ``SPEC``
    - ``THERMO`` or ``THERMO ALL``
    - ``REACTIONS`` or ``REAC``
    - ``TRANSPORT``

  * The thermodynamic data is read in a fixed format. This means that each
    column of the input has a particular meaning. *Many common errors are
    generated because information is missing or in the wrong column. Check
    thoroughly for extraneous or missing spaces.* The format for each
    thermodynamic entry should be as follows::

        N2                      N 2                 G200.000   6000.000  1000.00       1
         2.95258000E+00 1.39690000E-03-4.92632000E-07 7.86010000E-11-4.60755000E-15    2
        -9.23949000E+02 5.87189000E+00 3.53101000E+00-1.23661000E-04-5.02999000E-07    3
         2.43531000E-09-1.40881000E-12-1.04698000E+03 2.96747000E+00                   4

    The following table is adapted from the Chemkin manual [SAND89]_ to describe the
    column positioning of each required part of the entry. Empty columns should be
    filled with spaces.

    +---------+-------------------------------------+--------+
    |Line No. | Contents                            | Column |
    +=========+=====================================+========+
    | 1       | Species Name                        | 1--18  |
    +---------+-------------------------------------+--------+
    | 1       | Date (Optional)                     | 19--24 |
    +---------+-------------------------------------+--------+
    | 1       | Atomic Symbols and formula          | 25--44 |
    +---------+-------------------------------------+--------+
    | 1       | Phase of species (S, L, G)          | 45     |
    +---------+-------------------------------------+--------+
    | 1       | Low temperature                     | 46--55 |
    +---------+-------------------------------------+--------+
    | 1       | High temperature                    | 56--65 |
    +---------+-------------------------------------+--------+
    | 1       | Common temperature                  | 66--73 |
    +---------+-------------------------------------+--------+
    | 1       | Additional Atomic Symbols           | 74--78 |
    +---------+-------------------------------------+--------+
    | 1       | The integer ``1``                   | 80     |
    +---------+-------------------------------------+--------+
    | 2       | Coefficients :math:`a_1`            | 1--75  |
    |         | to :math:`a_5` for the upper        |        |
    |         | temperature interval                |        |
    +---------+-------------------------------------+--------+
    | 2       | The integer ``2``                   | 80     |
    +---------+-------------------------------------+--------+
    | 3       | Coefficients :math:`a_6,\ a_7`      | 1--75  |
    |         | for the upper temperature interval, |        |
    |         | and :math:`a_1,\ a_2,\ a_3` for     |        |
    |         | the lower temperature interval      |        |
    +---------+-------------------------------------+--------+
    | 3       | The integer ``3``                   | 80     |
    +---------+-------------------------------------+--------+
    | 4       | Coefficients :math:`a_4` through    | 1--60  |
    |         | :math:`a_7` for the lower           |        |
    |         | temperature interval                |        |
    +---------+-------------------------------------+--------+
    | 4       | The integer ``4``                   | 80     |
    +---------+-------------------------------------+--------+

    The first 18 columns are reserved for the species name. The name assigned
    to the species in the thermodynamic data must be the same as the species
    name defined in the ``SPECIES`` section. If the species name is shorter
    than 18 characters, the rest of the characters should be filled by spaces.
    The next six columns (columns 19--24) are typically used to write a date;
    they are not used further. The next 20 columns (25--44) are used to
    specify the elemental composition of the species. In column 45, the phase
    of the species (``S``, ``L``, or ``G`` for solid, liquid, or gas
    respectively) should be specified. The next 28 columns are reserved for
    the temperatures that delimit the ranges of the polynomials specified on
    the next several lines. The first two temperatures have a width of 10
    columns each (46--55 and 56--65), and represent the lowest temperature and
    highest temperature for which the polynomials are valid. The last
    temperature has a width of 8 columns (66--73) and is the "common"
    temperature, where the switch from low to high occurs. The next 5 columns
    (74--78) are reserved for atomic symbols and are usually left blank for
    the default behavior. Column 79 is blank and finally, the row is ended in
    column 80 with the integer ``1``.

    The next three lines of the thermodynamic entry have a similar format.
    They contain the coefficients of the polynomial described in
    :ref:`sec-thermo-models` for the NASA 7-coefficient polynomial formulation.
    The second row of the thermo entry (the first after the information row)
    contains the first five coefficients that apply the the temperature range
    between the midpoint and the upper limit. 15 columns are alloted for each
    coefficient (for a total of 75 columns), with no spaces between them.
    Although the entry above shows spaces between positive coefficients, it is
    to be noted that this is done only for formatting consistency with other
    lines that contain negative numbers. After the coefficients, four spaces
    in columns 76--79 are followed by the integer ``2`` in column 80. On the
    next line, the last two coefficients for the upper temperature range and
    the first three coefficients for the lower temperature range are
    specified. Once again, this takes up the first 75 columns, columns 76--79
    are blank, and the integer ``3`` is in column 80. Finally, on the last
    line of a particular entry, the last four coefficients of the lower
    temperature range are specified in columns 1--60, 19 blank spaces are
    present, and the integer ``4`` is in column 80. The 19 blank spaces in the
    last line are part of the standard. However, since the original Chemkin
    interpreter ignored those spaces, researchers began using that space to
    store additional information that was not necessary for the input file.
    Although these numbers create an error in ``ck2cti`` if present, they are
    harmless and can be ignored by using the ``--permissive`` option.

  * It may be the case that scientific formatted numbers are missing the ``E``.
    In this case, numbers often show up as ``1.1+01``, when they should be
    ``1.1E+01``. You can fix this with a simple Regular Expression find and
    replace::

        Find: (\d+\.\d+)([+-]\d+)
        Replace: \1E\2

  * The transport data file also has a specified format, as described in
    [SAND98]_, although the format is not as strict as for the thermodynamic
    entries. In particular, the first 15 columns of a line are reserved for
    the species name. *One common source of errors is a species that is present
    in the transport data file, but not in the thermodynamic data or in
    the species list; or a species that is present in the species list but
    not the transport data file.* The rest of the columns on a given line have
    no particular format, but must be present in the following order:

    +------------------+------------------------------------------------------+
    | Parameter Number | Parameter Name                                       |
    +==================+======================================================+
    | 1                | An integer with value 0, 1, or 2 indicating          |
    |                  | monatomic, linear, or non-linear molecular geometry. |
    +------------------+------------------------------------------------------+
    | 2                | The Lennard-Jones potential well depth               |
    |                  | :math:`\varepsilon/k_B` in Kelvin                    |
    +------------------+------------------------------------------------------+
    | 3                | The Lennard-Jones collision diameter :math:`\sigma`  |
    |                  | in Angstrom                                          |
    +------------------+------------------------------------------------------+
    | 4                | The dipole moment :math:`\mu` in Debye               |
    +------------------+------------------------------------------------------+
    | 5                | The polarizability :math:`\alpha` in Angstrom        |
    +------------------+------------------------------------------------------+
    | 6                | The rotational relaxation collision number           |
    |                  | :math:`Z_{rot}` at 298 K                             |
    +------------------+------------------------------------------------------+

    Another common error is if all 6 of these numbers are not present for every
    species.

.. [SAND89] See R. J. Kee, F. M. Rupley, and J. A. Miller, Sandia National
   Laboratories Report SAND89-8009 (1989).
   http://www.osti.gov/scitech/biblio/5681118

.. [SAND98] See R. J. Kee, G. Dixon-Lewis, J. Warnatz, M. E. Coltrin, J. A. Miller,
   H. K. Moffat, Sandia National Laboratories Report SAND86-8246B (1998).
