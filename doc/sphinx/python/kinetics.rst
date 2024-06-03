.. py:currentmodule:: cantera

Chemical Kinetics
=================

.. contents::
   :local:

Kinetics Managers
-----------------

Kinetics
^^^^^^^^
.. autoclass:: Kinetics()

InterfaceKinetics
^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceKinetics()

Reactions
---------

This class contains the definition of a single reaction, independent of a specific
`Kinetics` object. Reaction rate evaluation is handled by `ReactionRate` objects.

Reaction
^^^^^^^^
.. autoclass:: Reaction(reactants=None, products=None, rate=None, *, equation=None, init=True, third_body=None)
   :no-undoc-members:


Reaction Rates
--------------

ReactionRate
^^^^^^^^^^^^
.. autoclass:: ReactionRate()

ArrheniusRateBase
^^^^^^^^^^^^^^^^^
.. autoclass:: ArrheniusRateBase(input_data)
   :no-undoc-members:

ArrheniusRate
^^^^^^^^^^^^^
.. autoclass:: ArrheniusRate(A, b, Ea)
   :no-undoc-members:

BlowersMaselRate
^^^^^^^^^^^^^^^^
.. autoclass:: BlowersMaselRate(A, b, Ea, w)
   :no-undoc-members:

TwoTempPlasmaRate
^^^^^^^^^^^^^^^^^
.. autoclass:: TwoTempPlasmaRate(A, b, Ea_gas, Ea_electron)
   :no-undoc-members:

ElectronCollisionPlasmaRate
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ElectronCollisionPlasmaRate(energy_levels=None, cross_sections=None, input_data=None)

FalloffRate
^^^^^^^^^^^
.. autoclass:: FalloffRate()
   :no-undoc-members:

LindemannRate
^^^^^^^^^^^^^
.. autoclass:: LindemannRate(low, high, falloff_coeffs)
   :no-undoc-members:

TroeRate
^^^^^^^^
.. autoclass:: TroeRate(low, high, falloff_coeffs)
   :no-undoc-members:

SriRate
^^^^^^^
.. autoclass:: SriRate(low, high, falloff_coeffs)
   :no-undoc-members:

TsangRate
^^^^^^^^^
.. autoclass:: TsangRate(low, high, falloff_coeffs)
   :no-undoc-members:

PlogRate
^^^^^^^^
.. autoclass:: PlogRate(rates)
   :no-undoc-members:

ChebyshevRate
^^^^^^^^^^^^^
.. autoclass:: ChebyshevRate(temperature_range, pressure_range, data)
   :no-undoc-members:

CustomRate
^^^^^^^^^^
.. autoclass:: CustomRate(k)
   :no-undoc-members:

ExtensibleRate
^^^^^^^^^^^^^^
.. autoclass:: ExtensibleRate()
   :no-undoc-members:

InterfaceRateBase
^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceRateBase()
   :no-undoc-members:

InterfaceArrheniusRate
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceArrheniusRate(A, b, Ea)
   :no-undoc-members:

InterfaceBlowersMaselRate
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceBlowersMaselRate(A, b, Ea0, w)
   :no-undoc-members:

StickRateBase
^^^^^^^^^^^^^
.. autoclass:: StickRateBase()
   :no-undoc-members:

StickingArrheniusRate
^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: StickingArrheniusRate(A, b, Ea)
   :no-undoc-members:

StickingBlowersMaselRate
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: StickingBlowersMaselRate(A, b, Ea0, w)
   :no-undoc-members:

Auxiliary Reaction Data
-----------------------

ExtensibleRateData
^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleRateData()
   :no-undoc-members:

ThirdBody
^^^^^^^^^
.. autoclass:: ThirdBody(collider="M", *, efficiencies=None, default_efficiency=None)
   :no-undoc-members:

Arrhenius
^^^^^^^^^
.. autoclass:: Arrhenius(A, b, E)

Reaction Path Analysis
----------------------

ReactionPathDiagram
^^^^^^^^^^^^^^^^^^^
.. autoclass:: ReactionPathDiagram(kin: Kinetics, element: str)
