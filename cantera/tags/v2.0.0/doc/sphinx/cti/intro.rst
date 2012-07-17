************
Introduction
************

Virtually every Cantera simulation involves one or more phases of
matter. Depending on the calculation being performed, it may be necessary to
evaluate thermodynamic properties, transport properties, and/or homogeneous
reaction rates for the phase(s) present. In problems with multiple phases, the
properties of the interfaces between phases, and the heterogeneous reaction
rates at these interfaces, may also be required.

Before the properties can be evaluated, each phase must be defined, meaning that
the models to use to compute its properties and reaction rates must be
specified, along with any parameters the models require. For example, a solid
phase might be defined as being incompressible, with a specified density and
composition. A gaseous phase for a combustion simulation might be defined as an
ideal gas consisting of a mixture of many species that react with one another
via a specified set of reactions.

For phases containing multiple species and reactions, a large amount of data is
required to define the phase, since the contribution of each species to the
thermodynamic and transport properties must be specified, and rate information
must be given for each reaction. While this could be done directly in an
application program, a better approach is put the phase and interface
definitions in a text file that can be read by the application, so that a given
phase model can be re-used for other simulations.

This guide describes how to write such files to define phases and interfaces for
use in Cantera simulations. Section :ref:`sec-input-files` contains a summary of
some basic rules for writing input files, a discussion of how they are
processed, and of how errors are handled. In Section :ref:`sec-phases`, we will
go over how to define phases and interfaces, including how to import species and
reactions from external files. Then in :ref:`sec-species` and
:ref:`sec-reactions`, we'll look in depth at how to specify the component parts
of phase and interface models---the elements, species, and reactions.

.. In Section ##REF##, we'll put it all together, and present some complete,
   realistic example problems, showing the input file containing the definitions
   of all phases and interfaces, the application code to use the input file to
   solve a problem, and the resulting output.
