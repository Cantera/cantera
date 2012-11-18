
*******************
Language Interfaces
*******************

Although most of Cantera is written in C++, interfaces are provided to
allow users to work with Cantera from several different languages or
environments, including Fortran 90/95, Python, and MATLAB. Which
language should you choose? The basic rule of thumb is this: use
Python or MATLAB if possible; use C++ or Fortran if necessary.

Python
======

Python is a free scripting language that is designed to be easy to use. If you
are familiar with any other programming language, you can probably learn Python
in a couple of hours. It is also an elegant language, and provides a
user-friendly introduction to the concepts of object-oriented programming.
Python is great for solving problems quickly, and Cantera provides example
Python scripts to do calculations ranging from simple evaluation of
thermodynamic or transport properties, on up to chemical equilibrium in
multiphase mixtures, 1D laminar flames, reactor networks, and more.  If your
problem can be solved by using Cantera from Python, you'll almost certainly
solve it faster with Python than by writing programs in Fortran or C++. 

See http://www.python.org

Matlab
======

The comments above for Python apply to MATLAB too, except hat Python is free and
MATLAB isn't. If you have MATLAB already and are familiar with it, this is a
good choice for an environment from which to run Cantera. It is probably the
most popular Cantera application environment. http://www.mathworks.com.

C++
===

If you find that you need full access to the internals of Cantera, or want to
extend and customize Cantera, then C++ is the language for you. Most of Cantera
is itself written in C++, and so C++ application programs have more direct
access to Cantera's core functionality than do programs written in other
languages, which access Cantera through a library of C-like functions. From C++,
you can implement new equations of state, new models for transport properties,
and many other things that simply can't be done through the other language
interfaces. If you are doing substantial code development with Cantera, rather
than simply using it to solve a few problems, then you will probably want to use
it from C++.

Fortran
=======

Cantera provides an interface to Fortran 90/95, and can even be used from
Fortran 77 programs. Use this if you have existing Fortran code you want to port
to Cantera.
