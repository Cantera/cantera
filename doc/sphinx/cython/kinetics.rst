.. py:currentmodule:: cantera

Chemical Kinetics
=================

.. contents::
   :local:

Kinetics Managers
-----------------

Kinetics
^^^^^^^^
.. autoclass:: Kinetics(infile='', phaseid='', phases=())

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
.. autoclass:: Reaction(reactants='', products='')
   :no-undoc-members:

ElementaryReaction
^^^^^^^^^^^^^^^^^^
.. autoclass:: ElementaryReaction(reactants='', products='')
   :no-undoc-members:

ThreeBodyReaction
^^^^^^^^^^^^^^^^^
.. autoclass:: ThreeBodyReaction(reactants='', products='')
   :no-undoc-members:

FalloffReaction
^^^^^^^^^^^^^^^
.. autoclass:: FalloffReaction(reactants='', products='')
   :no-undoc-members:

ChemicallyActivatedReaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ChemicallyActivatedReaction(reactants='', products='')
   :no-undoc-members:

PlogReaction
^^^^^^^^^^^^
.. autoclass:: PlogReaction(reactants='', products='')
   :no-undoc-members:

ChebyshevReaction
^^^^^^^^^^^^^^^^^
.. autoclass:: ChebyshevReaction(reactants='', products='')
   :no-undoc-members:

BlowersMaselReaction
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: BlowersMaselReaction(reactants='', products='')
   :no-undoc-members:

InterfaceReaction
^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceReaction(reactants='', products='')
   :no-undoc-members:

BlowersMaselInterfaceReaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: BlowersMaselInterfaceReaction(reactants='', products='')
   :no-undoc-members:

Reaction Rates
---------------------------------------

ReactionRate
^^^^^^^^^^^^
.. autoclass:: ReactionRate()

ArrheniusRate
^^^^^^^^^^^^^
.. autoclass:: ArrheniusRate(A, b, Ea)
   :no-undoc-members:

BlowersMaselRate
^^^^^^^^^^^^^^^^
.. autoclass:: BlowersMaselRate(A, b, Ea, w)
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

Auxilliary Reaction Data (legacy only)
--------------------------------------

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

BlowersMasel
^^^^^^^^^^^^
.. autoclass:: BlowersMasel(A, b, E0, w)

Reaction Path Analysis
----------------------

ReactionPathDiagram
^^^^^^^^^^^^^^^^^^^
.. autoclass:: ReactionPathDiagram(Kinetics kin, str element)
