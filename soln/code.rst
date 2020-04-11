
Code Samples
~~~~~~~~~~~~

The following listings contain reference implementations of the modular
compartmental carbon cycle model I developed to solve this exercise. You're
encouraged to refine them and fine-tune them to your own needs.

Compartmental Carbon Cycle Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: carbon_model.py
   :language: python
   :caption: carbon_model.py - An implementation of the compartmental carbon model using a flexible interface for time integration.
   :linenos:

Modular Integrator Abstract Class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: integrator.py
  :language: python
  :caption: integrator.py - A modular time-integration scheme for ODEs
  :linenos:
