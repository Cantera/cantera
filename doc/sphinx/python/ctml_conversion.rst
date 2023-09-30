***********************
CTML to YAML conversion
***********************

.. seealso::
   For documentation of the command line interface, see the :ref:`ctml2yaml
   <sec-ctml2yaml>` section. For a tutorial, refer to the `Converting CTI and XML input
   files to YAML <https://cantera.org/tutorials/legacy2yaml.html>`_ pages.


Module-level documentation
==========================

.. py:module:: cantera.ctml2yaml
.. py:currentmodule:: cantera.ctml2yaml

The script ``ctml2yaml.py`` will convert files from the legacy CTML format to YAML input
format. The documentation below describes the classes and functions in the script. Each
function/method is annotated with the Python types that the function accepts.

Most users will access the functionality of this module via the command line with the
``ctml2yaml`` entry-point script. For programmatic access, the `main` and/or `convert`
functions should be used. `main` should be used when command line arguments must be
processed, while `convert` takes an input filename or a string containing the CTML file
to be converted, and optionally the name of the output file.

Module-level functions
----------------------

.. autofunction:: float2string
.. autofunction:: represent_float
.. autofunction:: get_float_or_quantity
.. autofunction:: split_species_value_string
.. autofunction:: clean_node_text
.. autofunction:: create_species_from_data_node
.. autofunction:: create_reactions_from_data_node
.. autofunction:: create_phases_from_data_node
.. autofunction:: convert
.. autofunction:: main

Conversion classes
------------------

.. autoclass:: Phase
   :no-undoc-members:
.. autoclass:: Species
   :no-undoc-members:
.. autoclass:: SpeciesThermo
   :no-undoc-members:
.. autoclass:: SpeciesTransport
   :no-undoc-members:
.. autoclass:: Reaction
   :no-undoc-members:

Exceptions
----------

.. autoexception:: MissingXMLNode
.. autoexception:: MissingXMLAttribute
.. autoexception:: MissingNodeText
