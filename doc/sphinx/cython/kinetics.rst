.. py:currentmodule:: cantera

Chemical Kinetics
=================

.. autoclass:: Kinetics(infile='', phaseid='', phases=())

Reactions
---------

These classes contain the definition of a single reaction and its associated
rate expression, indepenent of a specific `Kinetics` object.

.. autoclass:: Reaction(reactants='', products='')
   :no-undoc-members:

.. autoclass:: ElementaryReaction(reactants='', products='')
   :no-undoc-members:

.. autoclass:: ThreeBodyReaction(reactants='', products='')
   :no-undoc-members:

.. autoclass:: FalloffReaction(reactants='', products='')
   :no-undoc-members:

.. autoclass:: ChemicallyActivatedReaction(reactants='', products='')
   :no-undoc-members:

.. autoclass:: PlogReaction(reactants='', products='')
   :no-undoc-members:

.. autoclass:: ChebyshevReaction(reactants='', products='')
   :no-undoc-members:

.. autoclass:: InterfaceReaction(reactants='', products='')
   :no-undoc-members:

Auxilliary Reaction Data
------------------------

.. autoclass:: Arrhenius(A, b, E)

.. autoclass:: Falloff(coeffs=(), init=True)
   :no-undoc-members:

.. autoclass:: TroeFalloff(coeffs=(), init=True)
   :no-undoc-members:

.. autoclass:: SriFalloff(coeffs=(), init=True)
   :no-undoc-members:

Reaction Path Analysis
----------------------

.. autoclass:: ReactionPathDiagram(Kinetics kin, str element)
