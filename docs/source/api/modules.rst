API Reference: Modules
======================

This section documents the analysis modules that power scVIZ.

Explore Dataset Module
-----------------------

.. automodule:: modules.explore_dataset
   :members:
   :undoc-members:
   :show-inheritance:

Visualize Dataset Module
-------------------------

.. automodule:: modules.visualize_dataset
   :members:
   :undoc-members:
   :show-inheritance:

Differential Expression Module
-------------------------------

.. automodule:: modules.differential_expression
   :members:
   :undoc-members:
   :show-inheritance:

Module Functions
----------------

All analysis modules implement a ``render()`` function that serves as the entry point:

.. code-block:: python

   def render(adata: AnnData) -> None:
       """Render the module UI.
       
       Parameters
       ----------
       adata : AnnData
           The loaded single-cell dataset
       """

This function is called by the main application when a module is selected.
