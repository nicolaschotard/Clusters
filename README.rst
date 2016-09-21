.. image:: https://readthedocs.org/projects/clusters/badge/?version=latest
   :target: http://clusters.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://landscape.io/github/nicolaschotard/Clusters/master/landscape.svg?style=flat
   :target: https://landscape.io/github/nicolaschotard/Clusters/master
   :alt: Code Health

.. image:: https://travis-ci.org/nicolaschotard/Clusters.svg?branch=master
   :target: https://travis-ci.org/nicolaschotard/Clusters
   :alt: Travis CI build status (Linux)

.. image:: https://codecov.io/gh/nicolaschotard/Clusters/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/nicolaschotard/Clusters

____

**WARNING**: Package under development

____

.. inclusion-marker-do-not-remove

Clusters
--------

Python package wrapping up the ongoing cluster analysis of the french
LSST/DESC group. For more info, see the two following github
repositories:

- A collection of `notebooks <https://github.com/lsst-france/LSST_notebooks>`_ for LSST
- The `ReprocessingTaskForce <https://github.com/DarkEnergyScienceCollaboration/ReprocessingTaskForce>`_ repository

Installation
------------

To install::

  git clone https://github.com/nicolaschotard/Clusters.git
  pip install Clusters/

To install in a local directory ``mypath``, use::

  pip install --prefix='mypath' Clusters/

and do not forget to add it to your PYTHONPATH.

To upgrade to a new version (after a ``git pull`` or a local modification), use::

  pip install --upgrade (--prefix='mypath') Clusters/

To install a release version (no release version available yet)::

  pip install http://github.com/nicolaschotard/Cluster/archive/v0.1.tar.gz

Also works with the master::

  pip install (--upgrade) https://github.com/nicolaschotard/Clusters/archive/master.zip

In the future, release versions will be listed at this `location
<http://github.com/nicolaschotard/Clusters/releases>`_.
    

Dependencies
------------

``Clusters`` has for now the following dependencies:

- The LSST DM `stack <https://developer.lsst.io/build-ci/lsstsw.html>`_
- Python 2.7 and libraries listed in the `requirements <requirements.txt>`_ file
- `LEPHARE <http://cesam.lam.fr/lephare/lephare.html>`_


Usage
-----

``Clusters`` consists of several command-line executables that you
have to run in the right order.

- Get the input data and dump them in an hdf5 file containing astropy
  tables (see the `data format section
  <http://clusters.readthedocs.io/en/latest/data.html>`_ of the
  documentation for detail)::

    clusters_data config.yaml (--output output.hdf5)

- Correct the data for Milky Way extinction (input.hdf5 is the output of the previous step)::

    clusters_extinction.py config.yaml input.hdf5 --plot

- Get the photometric redshift using LEPHARE::

    clusters_zphot.py config.yaml input.hdf5 (--zpara zphot.para) --plot

**Next parts will come soon**

- Exctract background galaxies from the whole sample: remove the
  cluster galaxies (red sequence) and other foreground galaxies using
  the photometric redshifts::

    clusters_getbackground config.yaml input.hdf5 output.hdf5

Etc.

With any command, you can run with ``-h`` or ``--help`` to see all the
optional arguments, e.g., ``clusters_data.py -h``.

Configuration file
------------------

All the scripts will take the same input YAML file. Keys are listed
below and are case-sensitive. Additional keys are simply ignored. You
can find examples of these comfiguration files in the
`config <https://github.com/nicolaschotard/Clusters/blob/master/configs>`_
directory, or clicking `here <https://github.com/nicolaschotard/Clusters/blob/master/configs/MACSJ2243.3-0935.yaml>`_
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
