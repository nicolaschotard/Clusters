.. image:: https://readthedocs.org/projects/clusters/badge/?version=latest
   :target: http://clusters.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


.. warning::

   Package under development

____


.. include:: docs/source/overview

.. include:: docs/source/installation
    

Dependencies
------------

``Clusters`` has for now (too) many dependencies:

- The LSST DM `stack <https://developer.lsst.io/build-ci/lsstsw.html>`_)
- Python 2.7 and libraries
  
  - numpy
  - matplotlib
  - seaborn
  - astropy / astroquery
  - healpy
  - and probably other packages that will be replaced/listed at some point
    
- `LEPHARE <http://cesam.lam.fr/lephare/lephare.html>`_


Usage
-----

``Clusters`` consists of several command-line executables that you have
to run in the right order.

Get the input data and dump them in a pickle file (will change soon)::

  clusters_data config.yaml output.pkl

Correct the data for Milky Way extinction (output.pkl is the output of the previous step)::

  clusters_extinction.py config.yaml output.pkl --plot


Get the photometric redshift using LEPHARE::

  clusters_zphot.py config.yaml input_extcorr.pkl --plot

**Next parts will come soon**

Exctract background galaxies from the whole sample: remove the cluster
galaxies (red sequence) and other foreground galaxies using the
photometric redshifts::

  clusters_getbackground config.yaml input.pkl output.pkl

Etc.

With any command, you can run with ``-h`` or ``--help`` to see all the
optional arguments, e.g., ``clusters_data.py -h``.

Configuration file
------------------

All the scripts will take the same input YAML file. Keys are listed
below and are case-sensitive. Additional keys are simply ignored. You
can find examples of these comfiguration files in the
`config <https://github.com/nicolaschotard/Clusters/blob/master/configs>`_)
directory, or clicking `here <https://github.com/nicolaschotard/Clusters/blob/master/configs/MACSJ2243.3-0935.yaml>`_)
for MACSJ2243.3-0935.

+--------------------+--------+-------------------------------------------------------+
| Parameter          | Type   | Description [units]                                   |
+====================+========+=======================================================+
| ``"cluster"``      | string | Name of the cluster                                   |
+--------------------+--------+-------------------------------------------------------+
| ``"ra"``           | float  | RA coordinate of the cluster **[deg]**                |
+--------------------+--------+-------------------------------------------------------+
| ``"dec"``          | float  | DEC coordinate of the cluster **[deg]**               |
+--------------------+--------+-------------------------------------------------------+
| ``"redshift"``     | float  | Redshift the cluster                                  |
+--------------------+--------+-------------------------------------------------------+
| ``"filters"``      | string | Filter list to study, e.g., 'ugriz' (Megacam filters) |
+--------------------+--------+-------------------------------------------------------+
| ``"butler"``       | string | Absolute path to the intput data (butler)             |
+--------------------+--------+-------------------------------------------------------+
| ``"patches"``      | list   | List of patches to study                              |
+--------------------+--------+-------------------------------------------------------+
