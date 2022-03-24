.. py:currentmodule:: cantera

Chemical Kinetics
=================

.. contents::
   :local:

Kinetics Managers
-----------------

Kinetics
^^^^^^^^
.. autoclass:: Kinetics

InterfaceKinetics
^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceKinetics

Reactions
---------

These classes contain the definition of a single reaction, independent of a specific
`Kinetics` object. For legacy objects (CTI/XML input), each class integrates associated
rate expressions, whereas for the new, YAML-based implementation, reaction rate
evaluation is handled by dedicated `ReactionRate` objects.

Reaction
^^^^^^^^
.. autoclass:: Reaction
   :no-undoc-members:

ElementaryReaction
^^^^^^^^^^^^^^^^^^
.. autoclass:: ElementaryReaction
   :no-undoc-members:

ThreeBodyReaction
^^^^^^^^^^^^^^^^^
.. autoclass:: ThreeBodyReaction
   :no-undoc-members:

FalloffReaction
^^^^^^^^^^^^^^^
.. autoclass:: FalloffReaction
   :no-undoc-members:

ChemicallyActivatedReaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ChemicallyActivatedReaction
   :no-undoc-members:

PlogReaction
^^^^^^^^^^^^
.. autoclass:: PlogReaction
   :no-undoc-members:

ChebyshevReaction
^^^^^^^^^^^^^^^^^
.. autoclass:: ChebyshevReaction
   :no-undoc-members:

InterfaceReaction
^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceReaction
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

InterfaceRateBase
^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceRateBase
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
.. autoclass:: StickRateBase
   :no-undoc-members:

StickingArrheniusRate
^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: StickingArrheniusRate(A, b, Ea)
   :no-undoc-members:

StickingBlowersMaselRate
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: StickingBlowersMaselRate(A, b, Ea0, w)
   :no-undoc-members:

Auxiliary Reaction Data (legacy only)
-------------------------------------

Arrhenius
^^^^^^^^^
.. autoclass:: Arrhenius(A, b, E)

Falloff
^^^^^^^
.. autoclass:: Falloff(params=(), init=True)
   :no-undoc-members:

TroeFalloff
^^^^^^^^^^^
.. autoclass:: TroeFalloff(params=(), init=True)
   :no-undoc-members:

SriFalloff
^^^^^^^^^^
.. autoclass:: SriFalloff(params=(), init=True)
   :no-undoc-members:

Reaction Path Analysis
----------------------

ReactionPathDiagram
^^^^^^^^^^^^^^^^^^^
.. autoclass:: ReactionPathDiagram(Kinetics kin, str element)
