
*************************
A Very Simple C++ Program
*************************

A short C++ program that uses Cantera is shown below. This program reads in a
specification of a gas mixture from an input file, and then builds a new object
representing the mixture. It then sets the thermodynamic state and composition
of the gas mixture, and prints out a summary of its properties.

.. literalinclude:: demo1a.cpp
   :language: c++

This program produces the output below::

          temperature             500  K
             pressure          202650  Pa
              density        0.361163  kg/m^3
     mean mol. weight         7.40903  amu

                             1 kg            1 kmol
                          -----------      ------------
             enthalpy    -2.47725e+06       -1.835e+07     J
      internal energy    -3.03836e+06       -2.251e+07     J
              entropy         20700.1        1.534e+05     J/K
       Gibbs function    -1.28273e+07       -9.504e+07     J
    heat capacity c_p         3919.29        2.904e+04     J/K
    heat capacity c_v         2797.09        2.072e+04     J/K

                              X                 Y          Chem. Pot. / RT
                        -------------     ------------     ------------
                   H2            0.8         0.217667         -15.6441
                    H              0                0
                    O              0                0
                   O2              0                0
                   OH              0                0
                  H2O            0.1         0.243153         -82.9531
                  HO2              0                0
                 H2O2              0                0
                   AR            0.1          0.53918         -20.5027

As C++ programs go, this one is *very* short. It is the Cantera equivalent of
the "Hello, World" program most programming textbooks begin with. But it
illustrates some important points in writing Cantera C++ programs.

Catching :ct:`CanteraError` exceptions
======================================

The entire body of the program is put inside a function that is invoked within
a ``try`` block in the main program. In this way, exceptions thrown in the
function or in any procedure it calls may be caught. In this program, a
``catch`` block is defined for exceptions of type :ct:`CanteraError`. Cantera
throws exceptions of this type, so it is always a good idea to catch them. In
the ``catch`` block, function :ct:`showErrors` may be called to print the error
message associated with the exception.

The ``report`` function
=======================

The :ct:`report` function generates a nicely-formatted report of the properties of
a phase, including its composition in both mole (X) and mass (Y) units. For
each species present, the non-dimensional chemical potential is also printed.
This is handy particularly when doing equilibrium calculations. This function
is very useful to see at a glance the state of some phase.
