Preface
~~~~~~~

This is a reference solution for the **12.806 Computational Modeling
Exercise**, including versions of required and suggested figures for
reference and all code used to produce the solution. This solution was
implemented in Python, leaning on several widely-used and actively
developed and maintained scientific and numerical libraries.

For convenience, we'll use two very powerful libraries from the scientific
Python stack:

-  `NumPy <http://www.numpy.org>`_ - a numerical framework for Python
   centered around an :math:`n`-dimensional array object; includes many
   useful mathematical functions.
-  `pandas <http://pandas.pydata.org/>`_ - an analysis package based on
   data tables; allows you to quickly organize tabular data and run
   common statistics and analysis on them.

.. .. code-block:: python
..
..     import numpy as np
..     import pandas as pd

Additionally, we'll need to generate some plots, so we'll use
standard Python visualization libraries:

-  `matplotlib <http://matplotlib.org>`_ - Matlab-like visualization
   library
-  `seaborn <http://stanford.edu/~mwaskom/software/seaborn/>`_ -
   extension to matplotlib which generates quick statistical plots when
   data is packaged into pandas data structures; additionally includes
   aesthetic tweaks which greatly improve the matplotlib basics

.. .. code-block:: python
..
..     %matplotlib inline
..     import matplotlib.pyplot as plt
..     import seaborn as sns
..     # Set some plot aesthetics
..     plt.style.use(['seaborn-ticks', 'seaborn-talk'])

Finally, to read the Excel spreadsheet data source in the final part of
this notebook, you'll need the
`xlrd <https://xlrd.readthedocs.io/en/latest/>`_ package, which you can
install by invoking:

.. code-block:: shell

    $ pip install xlrd

All of the snippets shown here are culled from a companion `Jupyter Notebook <http://jupyter.org/>`_ which you can modify and run on your own personal machine.
