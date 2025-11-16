.. default-domain:: mat
.. currentmodule:: Base.+ct

Thermodynamic Properties
========================
These classes are used to describe the thermodynamic state of a system.

.. caution::
    The MATLAB toolbox is an experimental part of the Cantera API and may be changed
    without notice. It includes breaking changes from the legacy MATLAB API. While
    almost all features of the legacy MATLAB API are implemented, the toolbox does not
    include all functionality available for the C++ and Python interfaces.

Phases
------

ThermoPhase
^^^^^^^^^^^
.. autoclass:: ThermoPhase

Mixture
-------
.. autoclass:: Mixture(phases)
