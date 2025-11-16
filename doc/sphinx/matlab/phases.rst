.. default-domain:: mat
.. currentmodule:: Base.+ct

Objects Representing Phases
===========================

.. caution::
    The MATLAB toolbox is an experimental part of the Cantera API and may be changed
    without notice. It includes breaking changes from the legacy MATLAB API. While
    almost all features of the legacy MATLAB API are implemented, the toolbox does not
    include all functionality available for the C++ and Python interfaces.

.. contents::
   :local:

Composite Phase Objects
-----------------------
These classes are composite representations of a substance which has thermodynamic,
chemical kinetic, and (optionally) transport properties.

Solution
^^^^^^^^
.. autoclass:: Solution(src[, name][, transport_model])

Interface
^^^^^^^^^
.. autoclass:: Interface(src, name[, adj...])

.. _sec-matlab-purefluid:

Pure Fluid Phases
-----------------
The following convenience functions can be used to create `PureFluid` objects with the
indicated equation of state:

.. autofunction:: CarbonDioxide()
.. autofunction:: Heptane()
.. autofunction:: HFC134a()
.. autofunction:: Hydrogen()
.. autofunction:: Methane()
.. autofunction:: Nitrogen()
.. autofunction:: Oxygen()
.. autofunction:: Water()
