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

These classes contain the definition of a single reaction and its associated
rate expression, independent of a specific `Kinetics` object.

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

InterfaceReaction
^^^^^^^^^^^^^^^^^
.. autoclass:: InterfaceReaction(reactants='', products='')
   :no-undoc-members:

Auxilliary Reaction Data
------------------------

Arrhenius
^^^^^^^^^
.. autoclass:: Arrhenius(A, b, E)

Falloff
^^^^^^^
.. autoclass:: Falloff(coeffs=(), init=True)
   :no-undoc-members:

TroeFalloff
^^^^^^^^^^^
.. autoclass:: TroeFalloff(coeffs=(), init=True)
   :no-undoc-members:

SriFalloff
^^^^^^^^^^
.. autoclass:: SriFalloff(coeffs=(), init=True)
   :no-undoc-members:

Reaction Path Analysis
----------------------

ReactionPathDiagram
^^^^^^^^^^^^^^^^^^^
.. autoclass:: ReactionPathDiagram(Kinetics kin, str element)
